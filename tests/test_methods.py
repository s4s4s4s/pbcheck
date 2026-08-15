"""The two arms' pinned settings and the moderated arm's calibration.

Both properties here have already been wrong once in this repo's history, which is why they are
tests rather than comments:

* Spec §2 pins the naive call **verbatim**, including ``tie_correct=True`` and ``pts=True`` — two
  arguments scanpy defaults to False. They were missing until R0 (``c13a21e``), and their absence
  is invisible in the output: the run succeeds and quietly reports a differently-calibrated
  Wilcoxon. Asserted here by capturing the actual call, not by reading the source.
* The pseudobulk arm's calibration is the binding validity gate (decision rule item 1 / §8(a)) and
  is what Amendment 1 found DESeq2-Wald failing. The moderated arm selected by Amendment 2 is
  checked here on a fast fixed-seed null, as a regression guard — the real evidence is the 146-cell
  grid, not these few hundred replicates.
"""

import warnings

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from pbcheck import gate_config as gc
from pbcheck import mtc
from pbcheck.methods.moderated import ebayes_from_pdata
from pbcheck.methods.naive import naive_de

warnings.filterwarnings("ignore")


def test_naive_arm_passes_the_settings_the_spec_pins_verbatim(monkeypatch, small_null_adata):
    """Spec §2: ``tie_correct=True`` and ``pts=True``, captured from the real call."""
    seen = {}
    real = sc.tl.rank_genes_groups

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(sc.tl, "rank_genes_groups", spy)

    naive_de(small_null_adata)

    assert seen.get("method") == "wilcoxon"
    # tie_correct: scRNA-seq counts are heavily tied (mostly zeros); without it the Wilcoxon normal
    # approximation is distorted and the distortion lands in lambda_naive as a tie artifact rather
    # than as the clustering miscalibration the study is measuring (spec §2, B5).
    assert seen.get("tie_correct") is True
    # pts: the per-group expressed fractions the min-expression sensitivity check needs (§6).
    assert seen.get("pts") is True


def test_naive_arm_omits_tie_correct_for_the_t_test_variant(monkeypatch, small_null_adata):
    """``tie_correct`` is defined only for the rank test; passing it to the t-test variant errors."""
    seen = {}
    real = sc.tl.rank_genes_groups
    monkeypatch.setattr(sc.tl, "rank_genes_groups",
                        lambda *a, **k: (seen.update(k), real(*a, **k))[1])

    naive_de(small_null_adata, method="t-test_overestim_var")

    assert "tie_correct" not in seen
    assert seen.get("pts") is True     # the robustness variant still carries pts


def test_naive_arm_returns_both_expressed_fractions(small_null_adata):
    """``pts=True`` must reach the result table on BOTH sides, not just the call.

    On scanpy 1.12.2 ``rank_genes_groups_df`` returns ``pct_nz_group`` but NOT ``pct_nz_reference``
    when an explicit ``reference`` is given, so the dataframe helper alone yields one side. The
    min-expression-floor sensitivity check (spec §6/B5) needs both, and both are recoverable from
    ``uns['rank_genes_groups']['pts']``, which ``naive_de`` reads directly.
    """
    res = naive_de(small_null_adata)
    assert {"pct_group", "pct_reference"} <= set(res.table.columns)
    for col in ("pct_group", "pct_reference"):
        frac = res.table[col]
        assert frac.notna().all(), f"{col} has holes after reindexing onto the ranked gene order"
        assert ((frac >= 0) & (frac <= 1)).all(), f"{col} is not a fraction"


def test_naive_arm_expressed_fractions_match_the_matrix(small_null_adata):
    """The two fractions must be the per-condition nonzero rates, not the same column twice.

    Guards the reindex: ``uns['pts']`` is in ``var_names`` order while the table is in ranked order,
    so a missing reindex would silently mis-assign fractions to genes.
    """
    res = naive_de(small_null_adata)
    a = small_null_adata
    cond = a.obs["condition"].to_numpy()
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    for level, col in (("disease", "pct_group"), ("ctrl", "pct_reference")):
        expected = pd.Series((X[cond == level] > 0).mean(axis=0), index=list(a.var_names))
        got = res.table[col]
        np.testing.assert_allclose(got.to_numpy(), expected.reindex(got.index).to_numpy(), atol=1e-12)


# ---------------------------------------------------------------------------
# The moderated arm's calibration on a fast fixed-seed null.
# ---------------------------------------------------------------------------

class _FakePdata:
    """Minimal donor x gene stand-in: the arm only ever touches ``.X`` and ``.var_names``."""

    def __init__(self, X, genes):
        self.X = X
        self.var_names = list(genes)
        self.n_vars = X.shape[1]

    def __getitem__(self, key):
        _, cols = key
        idx = [self.var_names.index(g) for g in cols]
        return _FakePdata(self.X[:, idx], [self.var_names[i] for i in idx])

    def copy(self):
        return _FakePdata(self.X.copy(), list(self.var_names))

    @property
    def obs_names(self):
        return [f"d{i}" for i in range(self.X.shape[0])]


def _null_donor_matrix(rng, n_donors=8, n_genes=300, donor_sigma=0.5):
    """Donor pseudobulk counts under the COMPLETE null: no condition effect anywhere.

    A donor-level log-normal random effect on an NB-ish mean, drawn directly at the donor level —
    no cells, no aggregation — so the test stays in the milliseconds and depends on nothing but
    numpy. Truth is exactly 0 DE genes by construction.
    """
    base = rng.gamma(shape=2.0, scale=150.0, size=n_genes)
    eff = rng.lognormal(mean=0.0, sigma=donor_sigma, size=(n_donors, n_genes))
    return rng.poisson(base[None, :] * eff).astype(float)


@pytest.mark.parametrize("donor_sigma", [0.0, 0.5])
def test_moderated_arm_is_calibrated_on_a_fixed_seed_null(donor_sigma):
    """Under the complete null: FP rate <= alpha and lambda in the pre-registered band.

    The criteria are the spec's own (decision rule item 1 / §8(a)) read off ``gate_config``, so if
    a threshold ever moves without an amendment this test moves with it and the drift is at least
    recorded in one place. ``donor_sigma = 0`` is the easy case; 0.5 is the gate's own operating
    point, and the one at which Amendment 1 found DESeq2-Wald at FP 0.35-0.50.
    """
    rng = np.random.default_rng(31)
    n_rep, n_donors, n_genes = 120, 8, 300
    genes = [f"g{j}" for j in range(n_genes)]
    cond = np.array(["disease"] * (n_donors // 2) + ["ctrl"] * (n_donors // 2), dtype=object)

    rejected, lambdas = 0, []
    for _ in range(n_rep):
        pdata = _FakePdata(_null_donor_matrix(rng, n_donors, n_genes, donor_sigma), genes)
        res = ebayes_from_pdata(pdata, condition_values=cond, universe=genes)
        res = mtc.bh_over_universe(res, genes, alpha=gc.ALPHA)
        rejected += int(res.n_significant(fdr=gc.ALPHA) >= 1)
        lambdas.append(_lambda(res.table["pval"].to_numpy()))

    fp_rate = rejected / n_rep
    lam = float(np.median(lambdas))
    lo, hi = gc.LAMBDA_BAND

    # Under the COMPLETE null, V = R, so FDR = P(R >= 1) and BH bounds it by alpha (spec B1).
    # The bound is on the true rate; over 120 replicates the binomial SE at alpha is ~0.02, so the
    # observed rate is allowed a couple of SE of headroom rather than being pinned at exactly alpha.
    assert fp_rate <= gc.ALPHA + 0.05, f"sigma={donor_sigma}: FP rate {fp_rate} against alpha {gc.ALPHA}"
    assert lo <= lam <= hi, f"sigma={donor_sigma}: lambda {lam} outside {gc.LAMBDA_BAND}"


def _lambda(pvals):
    from pbcheck.metrics import genomic_inflation
    return genomic_inflation(pvals)


def test_moderated_arm_detects_a_planted_effect():
    """The reverse control: an arm that reports nothing on real signal is not calibrated, it is dead.

    Spec §7's reverse control and §8(c)'s purpose in miniature — a conservative-to-the-point-of-
    useless arm would pass the null test above and break every inflation ratio as a denominator.

    Uses the repo's own ``positive_oracle`` rather than a hand-planted effect: injecting a fold
    change into a raw count matrix by hand also inflates those donors' library sizes, and the CPM
    transform then pushes every *other* gene down, manufacturing false positives that belong to the
    test's construction rather than to the arm. The oracle is the generative convention this study
    is calibrated against; a small, strong-effect instance of it runs in about a second.
    """
    from oracles import positive_oracle
    from pbcheck.gene_universe import frozen_universe
    from pbcheck.methods.pseudobulk import build_pseudobulk

    o = positive_oracle(n_de=40, log2fc=2.0, seed=41, n_genes=400, n_donors_per_group=6,
                        n_cells_per_donor=80, donor_sigma=0.2, dispersion=0.2)
    pdata = build_pseudobulk(o.adata)
    universe = frozen_universe(pdata, min_size=100)
    res = mtc.bh_over_universe(
        ebayes_from_pdata(pdata, universe=universe), universe, alpha=gc.ALPHA)

    hits = res.significant(fdr=gc.ALPHA)
    truth = set(o.de_genes) & set(universe)
    sensitivity = len(hits & truth) / max(len(truth), 1)
    fdp = len(hits - set(o.de_genes)) / max(len(hits), 1)

    assert sensitivity > 0.5, f"the arm must find a large shared effect; sensitivity {sensitivity}"
    assert fdp < 0.25, f"and must not be swamped by false positives; FDP {fdp}"
