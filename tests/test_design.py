"""The design auditor reads only .obs and must flag the ingredients of pseudoreplication."""

import warnings

warnings.filterwarnings("ignore")

from oracles import null_oracle  # noqa: E402
from pbcheck.design import audit_design  # noqa: E402


def test_healthy_design_no_flags():
    rep = audit_design(null_oracle(n_donors_per_group=4, n_cells_per_donor=150, n_genes=300, seed=1).adata)
    assert rep.min_donors_per_group == 4
    assert rep.donor_nests_in_condition
    assert rep.usable_for_pseudobulk
    assert rep.flags == []


def test_too_few_donors_is_flagged():
    rep = audit_design(null_oracle(n_donors_per_group=2, n_cells_per_donor=150, n_genes=300, seed=1).adata)
    assert rep.min_donors_per_group == 2
    assert any("donor" in f and "smallest group" in f for f in rep.flags)


def test_group_imbalance_flagged():
    import numpy as np
    o = null_oracle(n_donors_per_group=4, n_cells_per_donor=200, n_genes=300, seed=2)
    a = o.adata
    # Drop 90% of ctrl cells to force a >5x imbalance.
    ctrl = a.obs["condition"].astype(str) == "ctrl"
    keep = ~ctrl | (np.arange(a.n_obs) % 10 == 0)
    a = a[keep].copy()
    rep = audit_design(a)
    assert rep.imbalance_ratio >= 5
    assert any("imbalance" in f for f in rep.flags)


def test_batch_confounding_detected():
    o = null_oracle(n_donors_per_group=4, n_cells_per_donor=100, n_genes=200, seed=3)
    a = o.adata
    # Make a batch column perfectly aligned with condition (worst case).
    a.obs["batch"] = a.obs["condition"].astype(str).map({"ctrl": "b1", "disease": "b2"})
    rep = audit_design(a, batch_cols=["batch"])
    assert rep.batch_confounded["batch"] >= 0.8
    assert any("confounded" in f for f in rep.flags)
