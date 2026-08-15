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
    assert rep.batch_separates_condition["batch"]
    assert any("PERFECTLY separates" in f for f in rep.flags)
    # A perfectly separated design cannot be attributed to pseudoreplication at all.
    assert not rep.usable_for_pseudobulk


def test_confound_is_measured_per_donor_not_per_cell():
    """A donor-level confound must not be hidden by cell-count weighting.

    This is the regression test for the auditor committing the very fallacy it audits: with the
    association measured over cells, a batch that contains only disease donors scored ~0.05 and went
    unflagged, because the donors carrying it happened to be small. Measured over donors it is 0.707.
    """
    import numpy as np

    o = null_oracle(n_donors_per_group=4, n_cells_per_donor=200, n_genes=200, seed=5)
    a = o.adata
    donors = list(dict.fromkeys(a.obs["donor"].astype(str)))
    cond = {d: a.obs.loc[a.obs["donor"].astype(str) == d, "condition"].astype(str).iloc[0] for d in donors}
    disease = [d for d in donors if cond[d] == "disease"]

    # One batch level containing ONLY disease donors — a genuine donor-level confound.
    confounded = set(disease[:2])
    a.obs["batch"] = ["bX" if d in confounded else "bY" for d in a.obs["donor"].astype(str)]

    # Starve exactly those donors of cells so a cell-weighted estimate would miss the confound.
    donor_arr = a.obs["donor"].astype(str).to_numpy()
    rank = np.zeros(a.n_obs, dtype=int)
    for d in donors:
        m = donor_arr == d
        rank[m] = np.arange(m.sum())
    keep = ~np.isin(donor_arr, list(confounded)) | (rank < 5)
    a = a[keep].copy()

    rep = audit_design(a, batch_cols=["batch"])
    donor_v = rep.batch_confounded["batch"]

    # What the old, cell-weighted implementation would have reported on this same design.
    from pbcheck.design import _cramers_v

    cell_v = _cramers_v(a.obs["batch"].astype(str), a.obs["condition"].astype(str))

    assert donor_v >= 0.5, f"donor-level confound scored only {donor_v}"
    assert donor_v > 3 * cell_v, (
        f"donor-level V={donor_v:.3f} is not clearly above the cell-weighted V={cell_v:.3f}; "
        "the association is still being driven by cell counts"
    )


def test_min_donors_threshold_reaches_the_verdict():
    """Regression: usable_for_pseudobulk hardcoded >= 3 regardless of the caller's min_donors.

    audit_design(..., min_donors=5) changed the warning text (via the flags list) but not the
    verdict, because usable_for_pseudobulk compared min_donors_per_group against a literal 3
    instead of self.min_donors. A 4-donors/group design must pass at the default/lenient threshold
    but fail once the caller asks for a stricter min_donors=5.
    """
    o = null_oracle(n_donors_per_group=4, n_cells_per_donor=100, n_genes=200, seed=9).adata
    lenient = audit_design(o, min_donors=3)
    strict = audit_design(o, min_donors=5)
    assert lenient.usable_for_pseudobulk
    assert not strict.usable_for_pseudobulk
    assert any("smallest group" in f for f in strict.flags)


def test_three_donor_rule_matches_the_spec_inclusion_gate():
    """Spec §1 item 1 requires >= 3 donors per group; 2 is not usable."""
    two = audit_design(null_oracle(n_donors_per_group=2, n_cells_per_donor=100, n_genes=200, seed=7).adata)
    three = audit_design(null_oracle(n_donors_per_group=3, n_cells_per_donor=100, n_genes=200, seed=7).adata)
    assert not two.usable_for_pseudobulk
    assert three.usable_for_pseudobulk
