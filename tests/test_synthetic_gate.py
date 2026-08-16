"""The gate's choice of which permutation floor is the headline — pinned, because it was wrong.

``scripts/synthetic_gate.py`` computed the headline naive floor as ``perm_floor`` over
``run_null``'s ``naive_ndeg``, an array that held paired-BH counts below index ``n_paired`` and the
naive arm's own BH above it. With ``gate_config.N_PERM == N_PERM_PB == 200`` the two halves
coincide and the committed artifact is unaffected; at spec §4's pre-registered sweep counts
(``n_perm`` 1000, ``n_perm_pb`` >= 200) 80 % of that array comes from the other convention and the
headline floor silently becomes a median over both.

These tests fail if a future caller re-mixes them: one checks the selection directly at counts
where the two series disagree, one greps the tree for the retired key, and the slow one drives the
real ``run_null`` at ``n_perm != n_perm_pb`` — a regime nothing in the suite exercised before.
"""

from __future__ import annotations

import ast
import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

from pbcheck import metrics

REPO = Path(__file__).resolve().parent.parent
GATE = REPO / "scripts" / "synthetic_gate.py"


def _load_gate():
    warnings.filterwarnings("ignore")
    spec = importlib.util.spec_from_file_location("_synthetic_gate", GATE)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _fake_null_res(paired_counts, solo_counts):
    return {
        "naive_ndeg_paired": metrics.NDegSeries(np.asarray(paired_counts), metrics.BH_PAIRED),
        "naive_ndeg_solo": metrics.NDegSeries(np.asarray(solo_counts), metrics.BH_SOLO),
    }


def test_naive_floors_reads_each_series_under_its_own_convention():
    """The headline floor is the paired series, and neither floor is a median over both."""
    gate = _load_gate()
    # Deliberately disjoint: paired sits at 100, the solo-only tail at 900. A median over the
    # concatenation would be 900 -- neither of the two legitimate answers.
    paired = [100, 100, 100]
    solo = [100, 100, 100] + [900] * 12
    floors = gate.naive_floors(_fake_null_res(paired, solo), n_genes=1000)

    assert floors["paired"]["bh_mode"] == metrics.BH_PAIRED
    assert floors["paired"]["n_perm"] == 3
    assert floors["paired"]["median_count"] == 100

    assert floors["solo"]["bh_mode"] == metrics.BH_SOLO
    assert floors["solo"]["n_perm"] == 15
    assert floors["solo"]["median_count"] == 900

    mixed = metrics.perm_floor(np.asarray(solo), n_genes=1000)["median_count"]
    assert floors["paired"]["median_count"] != mixed


def test_no_module_reads_the_retired_mixed_ndeg_key():
    """``naive_ndeg`` is gone, not deprecated: reintroducing it reintroduces the defect.

    Matched on the AST rather than on text, so ``naive_ndeg_paired`` and ``naive_ndeg_solo`` -- the
    two replacements, which contain the retired name as a substring -- do not trip it.
    """
    offenders = []
    for path in sorted((REPO / "src").rglob("*.py")) + sorted((REPO / "scripts").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if (isinstance(node, ast.Subscript)
                    and isinstance(node.slice, ast.Constant)
                    and node.slice.value == "naive_ndeg"):
                offenders.append(f"{path.relative_to(REPO)}:{node.lineno}")
    assert not offenders, (
        "the mixed-BH 'naive_ndeg' key is read again at "
        + ", ".join(offenders)
        + " -- use naive_ndeg_paired (cross-arm comparable) or naive_ndeg_solo (naive arm's own)"
    )


@pytest.mark.slow
def test_run_null_separates_the_conventions_at_unequal_permutation_counts():
    """The regime the defect needed and nothing exercised: ``n_perm`` != ``n_perm_pb``.

    Beyond ``n_paired`` there is no pseudobulk fit and therefore no common tested set, so no paired
    #DEG exists for those labelings at all. The paired series is sized to the paired permutations
    instead of being padded out with counts from the other BH.
    """
    warnings.filterwarnings("ignore")
    from oracles import null_oracle

    from pbcheck.gene_universe import frozen_universe
    from pbcheck.methods.pseudobulk import build_pseudobulk
    from pbcheck.permutation import run_null

    n_perm, n_perm_pb = 9, 4
    o = null_oracle(n_donors_per_group=4, n_cells_per_donor=60, n_genes=200, seed=13)
    uni = frozen_universe(build_pseudobulk(o.adata), min_size=50)
    r = run_null(o.adata, uni, n_perm=n_perm, n_perm_pb=n_perm_pb, fdr=0.05)

    assert "naive_ndeg" not in r, "the mixed-convention array is back"

    paired, solo, pb = r["naive_ndeg_paired"], r["naive_ndeg_solo"], r["pb_ndeg"]
    assert (paired.bh_mode, len(paired)) == (metrics.BH_PAIRED, n_perm_pb)
    assert (solo.bh_mode, len(solo)) == (metrics.BH_SOLO, n_perm)
    assert (pb.bh_mode, len(pb)) == (metrics.BH_PAIRED, n_perm_pb)
    assert r["n_perm_paired"] == n_perm_pb and r["n_perm_naive"] == n_perm

    # The Monte-Carlo block's floor is the paired one; its solo companion is the other series.
    mc = r["monte_carlo"]
    assert mc["naive_floor_median"] == float(np.median(paired.counts))
    assert mc["naive_floor_median_solo"] == float(np.median(solo.counts))

    # And the gate's headline floor is over exactly the paired permutations, not all of them.
    gate = _load_gate()
    floors = gate.naive_floors(r, len(uni))
    assert floors["paired"]["n_perm"] == n_perm_pb
    assert floors["solo"]["n_perm"] == n_perm
    assert floors["paired"]["median_count"] == mc["naive_floor_median"]
