"""The synthetic oracles are pbcheck's correctness spec.

A pseudoreplication auditor is only trustworthy if, on data where the truth is known:
  * it sees the naive test explode on a null that has donor structure, while pseudobulk stays clean;
  * that explosion DISAPPEARS when the donor structure is removed (the falsification control), proving
    the gap is pseudoreplication and not merely the naive test's larger sample size.

Both arms are scored through the same frozen, label-agnostic universe + identical BH (spec §3/§5).
"""

import warnings

import pytest

warnings.filterwarnings("ignore")

from oracles import null_oracle, no_donor_effect_oracle, positive_oracle  # noqa: E402
from pbcheck.methods import naive_de, pseudobulk_de  # noqa: E402
from pbcheck.methods.pseudobulk import build_pseudobulk  # noqa: E402
from pbcheck.gene_universe import frozen_universe  # noqa: E402
from pbcheck import mtc  # noqa: E402

# Every test in this module runs both arms end-to-end (naive + pseudobulk DE, BH, on real
# oracle-generated data) rather than exercising one function in isolation — the slowest tests in
# the suite, and the ones `pytest -m "not slow"` is for.
pytestmark = pytest.mark.slow

COMMON = dict(n_genes=800, n_donors_per_group=4, n_cells_per_donor=200, dispersion=0.2, seed=11)


def _run(oracle):
    pdata = build_pseudobulk(oracle.adata)
    uni = frozen_universe(pdata)
    pb = mtc.bh_over_universe(pseudobulk_de(oracle.adata, universe=uni), uni)
    nv = mtc.bh_over_universe(naive_de(oracle.adata, genes=uni), uni)
    return nv, pb, uni


def test_null_with_donor_effect_inflates_naive_only():
    nv, pb, _ = _run(null_oracle(donor_sigma=0.5, **COMMON))
    assert pb.n_significant() <= 5, f"pseudobulk should call ~0 on a null, got {pb.n_significant()}"
    assert nv.n_significant() >= 100, f"naive should inflate massively, got {nv.n_significant()}"
    assert nv.n_significant() > 20 * max(pb.n_significant(), 1)


def test_falsification_control_no_donor_effect():
    nv, pb, _ = _run(no_donor_effect_oracle(**COMMON))
    assert pb.n_significant() <= 5
    assert nv.n_significant() <= 20, (
        f"with no donor effect the naive test must stay clean, got {nv.n_significant()} "
        "(if this fails, the measured inflation could be a power artifact, not pseudoreplication)"
    )


def test_control_vs_donor_effect_gap():
    nv_donor, _, _ = _run(null_oracle(donor_sigma=0.5, **COMMON))
    nv_ctrl, _, _ = _run(no_donor_effect_oracle(**COMMON))
    assert nv_donor.n_significant() > 10 * max(nv_ctrl.n_significant(), 1)


def test_pseudobulk_recovers_true_signal_specifically():
    oracle = positive_oracle(n_de=80, log2fc=1.5, donor_sigma=0.4, **COMMON)
    _, pb, _ = _run(oracle)
    truth = set(oracle.de_genes)
    hits = pb.significant()
    false_pos = hits - truth
    assert len(false_pos) <= max(5, int(0.1 * len(hits) + 1)), f"too many pseudobulk false positives: {len(false_pos)}"
    assert len(hits & truth) >= 1
