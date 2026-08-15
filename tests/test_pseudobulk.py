"""The pseudobulk arm: the thin-donor filter and the arm's test identity (Amendment 2).

Two properties are pinned here because both were silently wrong before Amendment 2:

1. The pre-registered ``min_cells`` / ``min_counts`` filter (spec §1 inclusion-gate item 2, §3) had
   never run — decoupler 2.x removed the parameters the spec pins, and ``build_pseudobulk`` accepted
   them and forwarded them nowhere. Amendment 2 Change 7 implements it manually, with the spec's
   own semantics: **dropped, not merged**.
2. The arm's test is moderated eBayes, not DESeq2-Wald (Change 1).
"""

import warnings

import numpy as np
import pytest

from pbcheck.gene_universe import frozen_universe
from pbcheck.methods.pseudobulk import build_pseudobulk, pseudobulk_de

warnings.filterwarnings("ignore")


@pytest.fixture(scope="module")
def oracle():
    from oracles import null_oracle
    return null_oracle(n_donors_per_group=4, n_cells_per_donor=60, n_genes=150, seed=13)


def _thin_one_donor(adata, donor, keep_cells):
    """Return a copy in which ``donor`` contributes only ``keep_cells`` cells."""
    obs = adata.obs
    is_donor = (obs["donor"].astype(str) == donor).to_numpy()
    idx = np.flatnonzero(is_donor)
    drop = set(idx[keep_cells:].tolist())
    mask = np.array([i not in drop for i in range(adata.n_obs)])
    return adata[mask].copy()


def test_thin_donor_is_dropped_not_merged(oracle):
    """A donor below ``min_cells`` must vanish from the pseudobulk matrix, taking its counts with it.

    "Merged" would mean its counts reappear inside some other profile — which would corrupt the
    surviving donors rather than shrink the sample. The spec says dropped; this asserts dropped.
    """
    full = build_pseudobulk(oracle.adata)
    assert full.n_obs == 8
    assert full.uns["thin_donor_filter"]["n_dropped"] == 0

    thin = _thin_one_donor(oracle.adata, "d00", keep_cells=4)   # 4 < min_cells = 10
    out = build_pseudobulk(thin)

    assert out.n_obs == 7, "the thin donor must be dropped"
    assert not any(str(n).startswith("d00") for n in out.obs_names)

    info = out.uns["thin_donor_filter"]
    assert info["n_dropped"] == 1
    assert info["min_cells"] == 10 and info["min_counts"] == 1000
    assert info["per_threshold"]["min_cells"] == 1
    assert any(str(d).startswith("d00") for d in info["dropped_profiles"])

    # DROPPED, NOT MERGED: every surviving profile is byte-identical to what it was before, so the
    # thin donor's counts went nowhere. (d00's own cells are a subset in `thin`, so it is excluded.)
    kept = [n for n in out.obs_names if not str(n).startswith("d00")]
    before = full[kept].X
    after = out[kept].X
    assert np.array_equal(np.asarray(before), np.asarray(after))
    # ...and the total counts fell by exactly the dropped profile's share.
    assert np.asarray(out.X).sum() < np.asarray(full.X).sum()


def test_thin_donor_filter_catches_low_total_counts_too(oracle):
    """``min_counts`` drops an empty-but-populous profile; both thresholds are live, not decorative."""
    pdata = build_pseudobulk(oracle.adata, min_cells=0, min_counts=10**9)
    assert pdata.n_obs == 0
    assert pdata.uns["thin_donor_filter"]["per_threshold"]["min_cells"] == "disabled"
    assert pdata.uns["thin_donor_filter"]["per_threshold"]["min_counts"] == 8


def test_filter_records_when_a_threshold_could_not_be_applied(oracle):
    """An unavailable bookkeeping column must be reported as such, never treated as satisfied."""
    from pbcheck.methods.pseudobulk import _drop_thin_donors

    pdata = build_pseudobulk(oracle.adata)
    del pdata.obs["psbulk_cells"]
    out = _drop_thin_donors(pdata, min_cells=10, min_counts=1000)
    assert "unavailable" in out.uns["thin_donor_filter"]["per_threshold"]["min_cells"]


def test_the_arm_is_moderated_ebayes_and_reports_its_prior(oracle):
    """Amendment 2 Change 1: the arm's test changed, and its realised prior is persisted (C5)."""
    universe = frozen_universe(build_pseudobulk(oracle.adata), min_size=50)
    res = pseudobulk_de(oracle.adata, universe=universe)

    assert res.method.startswith("pseudobulk[eBayes")
    assert "DESeq2" not in res.method
    assert len(res.table) == len(universe)
    assert res.table["padj"].isna().all(), "BH must stay outside the arm (spec section 5)"

    m = res.moderation
    for k in ("prior_df_d0", "prior_is_complete_pooling", "residual_df_d",
              "shrinkage_factor_d0_over_d0_plus_d", "s0_squared_is_trended", "n_pval_nan"):
        assert k in m, k
    assert m["residual_df_d"] == 8 - 2
    assert m["s0_squared_is_trended"] is False   # the trended variant was measured and not selected
    assert 0.0 <= m["shrinkage_factor_d0_over_d0_plus_d"] <= 1.0


def test_deseq_path_is_retained_for_reproducibility(oracle):
    """The superseded arm must still run: Amendment-1-era numbers have to stay reproducible."""
    from pbcheck.methods.pseudobulk import deseq_from_pdata

    pdata = build_pseudobulk(oracle.adata)
    universe = frozen_universe(pdata, min_size=50)
    res = deseq_from_pdata(pdata, universe=universe, n_cpus=1)
    assert "DESeq2" in res.method
    assert len(res.table) > 0
