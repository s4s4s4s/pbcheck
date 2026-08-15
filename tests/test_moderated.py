"""The extraction of the moderated arm out of the probe must be numerically a no-op.

Amendment 2 Change 1 selects moderated eBayes on the strength of the 146-cell grid in
``pilot/testsel/summary.json``. That evidence is only evidence for the *shipped* arm if the shipped
code is the code that produced it. ``src/pbcheck/methods/moderated.py`` was lifted verbatim out of
``scripts/pb_calibration_probe.py``; these tests pin the lift by running a fixed-seed synthetic
pseudobulk matrix through **both** and asserting bit-level agreement of the p-values.

The probe is deliberately left untouched by that change (it is the diagnostic instrument the
selection was measured with), so it stays importable as an independent reference here.
"""

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from pbcheck.methods.moderated import (
    ebayes_from_pdata,
    fit_f_dist,
    moderated_t,
    trigamma_inverse,
    wls_two_group,
)

REPO = Path(__file__).resolve().parents[1]


def _probe():
    """Import ``scripts/pb_calibration_probe.py`` as a module (it is not on the import path)."""
    path = REPO / "scripts" / "pb_calibration_probe.py"
    if not path.exists():  # pragma: no cover
        pytest.skip("probe script not present")
    spec = importlib.util.spec_from_file_location("_pb_calibration_probe", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _pdata(seed=17, n_donors=8, n_genes=120):
    """A fixed-seed donor x gene pseudobulk count matrix with donor-level overdispersion.

    Not an oracle run: this is a self-contained matrix so the test stays in the seconds range and
    depends on nothing but numpy's generator.
    """
    rng = np.random.default_rng(seed)
    base = rng.gamma(shape=2.0, scale=200.0, size=n_genes)
    donor_eff = rng.lognormal(mean=0.0, sigma=0.5, size=(n_donors, n_genes))
    lam = base[None, :] * donor_eff
    X = rng.poisson(lam).astype(float)
    obs = pd.DataFrame(
        {"donor": [f"d{i}" for i in range(n_donors)],
         "condition": ["disease"] * (n_donors // 2) + ["ctrl"] * (n_donors - n_donors // 2)},
        index=[f"d{i}" for i in range(n_donors)],
    )
    var = pd.DataFrame(index=[f"g{j}" for j in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


def test_extraction_is_numerically_identical_to_the_probe():
    """The shipped arm must reproduce the probe's ebayes p-values exactly, not approximately."""
    probe = _probe()
    pdata = _pdata()
    universe = sorted(str(g) for g in pdata.var_names)
    cond = pdata.obs["condition"].to_numpy()

    ref = probe._ebayes_result(pdata[:, universe].copy(), cond, "disease", "ctrl", "ebayes")
    got = ebayes_from_pdata(pdata, universe=universe)

    assert list(got.table.index) == list(ref.table.index)
    for col in ("pval", "log2fc"):
        a = ref.table[col].to_numpy(dtype=float)
        b = got.table[col].to_numpy(dtype=float)
        assert np.allclose(a, b, rtol=0, atol=0, equal_nan=True), f"{col} diverged from the probe"

    # The realised prior — the fitType analog persisted per run (Amendment 2 Change 4 / spec C5) —
    # must match too, or the grid's prior_df_d0 / shrinkage columns do not describe this arm.
    for k in ("prior_df_d0", "shrinkage_factor_d0_over_d0_plus_d", "residual_df_d",
              "s0_squared", "prior_is_complete_pooling", "s0_squared_is_trended"):
        assert got.moderation[k] == ref._diag[k], k


def test_extraction_is_identical_for_the_trended_variant_too():
    """``trend=True`` is the variant Amendment 2 measured and rejected; it must still be faithful."""
    probe = _probe()
    pdata = _pdata(seed=23)
    universe = sorted(str(g) for g in pdata.var_names)
    cond = pdata.obs["condition"].to_numpy()

    ref = probe._ebayes_result(pdata[:, universe].copy(), cond, "disease", "ctrl", "ebayes_trend")
    got = ebayes_from_pdata(pdata, universe=universe, trend=True)

    assert np.allclose(ref.table["pval"].to_numpy(dtype=float),
                       got.table["pval"].to_numpy(dtype=float), rtol=0, atol=0, equal_nan=True)
    assert got.moderation["s0_squared_is_trended"] is True


def test_d0_zero_reduces_to_the_ordinary_pooled_t():
    """Smyth's anchor: no prior weight => the equal-variance two-sample Student t, exactly."""
    from scipy.stats import ttest_ind

    rng = np.random.default_rng(3)
    Y = rng.normal(size=(8, 200))
    x = np.array([1.0] * 4 + [0.0] * 4)
    beta, s2, d, uvar = wls_two_group(Y, x)
    _, p_mod, _, df_total = moderated_t(beta, s2, d, uvar, 0.0, 0.0)
    p_ref = ttest_ind(Y[x == 1], Y[x == 0], axis=0, equal_var=True).pvalue

    assert df_total == d == 6
    assert np.allclose(p_mod, p_ref)
    # and the unweighted design reduces to the textbook contrast / unscaled variance
    assert np.allclose(beta, Y[x == 1].mean(axis=0) - Y[x == 0].mean(axis=0))
    assert np.allclose(uvar, 1 / 4 + 1 / 4)


def test_trigamma_inverse_inverts_trigamma():
    from scipy.special import polygamma

    x = np.array([1e-8, 1e-3, 0.1, 1.0, 5.0, 1e3, 1e9])
    y = trigamma_inverse(x)
    assert np.max(np.abs(polygamma(1, y) / x - 1.0)) < 1e-6


def test_prior_is_complete_pooling_reported_not_clipped_on_homoscedastic_data():
    """Homoscedastic genes give d0 = +inf, which must surface as a flag, never as a silent number."""
    rng = np.random.default_rng(5)
    s2 = rng.chisquare(6, size=500) / 6.0  # exactly the null model with no between-gene variance
    d0, s02, info = fit_f_dist(s2, 6)
    assert np.isinf(d0) or d0 > 50, f"expected near-complete pooling, got d0={d0}"
    assert info["n_genes_used_in_prior"] == 500


def test_nonfinite_variance_yields_nan_pvalue_not_a_prior_derived_one():
    """A gene the arm cannot test must be dropped by BH, not handed a fabricated p-value."""
    beta = np.array([1.0, 1.0])
    s2 = np.array([np.nan, 0.5])
    uvar = np.array([0.5, 0.5])
    for d0 in (np.inf, 3.0):
        _, p, _, _ = moderated_t(beta, s2, 6, uvar, d0, 0.4)
        assert np.isnan(p[0]), f"d0={d0}: untestable gene got p={p[0]}"
        assert np.isfinite(p[1])
