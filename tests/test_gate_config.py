"""The pre-registered thresholds must match the frozen spec, and must not drift between copies.

Amendment 1 exists because a gate script quietly substituted an easier oracle and dropped two
binding criteria, and nothing noticed. These tests are the mechanical part of the fix: the values
are read back out of ``docs/PHASE0_SPEC.md`` itself, and the probe's independent copy is pinned
equal to the package's.
"""

import importlib.util
import re
import sys
from pathlib import Path

import pytest

from pbcheck import gate_config as gc

REPO = Path(__file__).resolve().parents[1]
SPEC = (REPO / "docs" / "PHASE0_SPEC.md").read_text(encoding="utf-8")
AMENDMENTS = (REPO / "docs" / "AMENDMENTS.md").read_text(encoding="utf-8")


def test_pre_registered_values_appear_in_the_frozen_spec():
    """Each pre-registered number must be findable in the spec text that pins it.

    A weak check by design — it cannot verify the *meaning* — but it is exactly strong enough to
    catch the failure that produced Amendment 1: a threshold silently replaced by an easier one.
    """
    assert "alpha = 0.05" in SPEC
    assert "[0.9, 1.1]" in SPEC
    assert re.search(r"sensitivity\s*(>=|≥)\s*0\.6", SPEC)
    assert re.search(r"log2FC\s*=\s*1\.0", SPEC) and re.search(r"K\s*=\s*200", SPEC)
    assert "default 200 genes" in SPEC                    # C5 / inclusion gate item 5
    assert "min_cells=10" in SPEC and "min_counts=1000" in SPEC
    assert re.search(r"lambda_naive\s*(>=|≥)\s*2\.0", SPEC)
    assert re.search(r"lambda_naive\s*(>=|≥)\s*3", SPEC)

    assert gc.ALPHA == 0.05
    assert gc.LAMBDA_BAND == (0.9, 1.1)
    assert gc.POWER_TARGET == 0.60
    assert (gc.ORACLE_LOG2FC, gc.ORACLE_K) == (1.0, 200)
    assert gc.MIN_UNIVERSE_SIZE == 200
    assert (gc.MIN_CELLS, gc.MIN_COUNTS) == (10, 1000)
    assert gc.LAMBDA_NAIVE_GO == 2.0 and gc.LAMBDA_NAIVE_STRONG_GO == 3.0


def test_amendment_3_values_are_pinned_by_the_amendment_that_introduces_them():
    """The evaluation regimes and the operating envelope are frozen by Amendment 3, not by the spec.

    ``PHASE0_SPEC.md`` never pinned a ``sigma_donor`` at all — that hole is exactly what Amendment 3
    closes — so these values cannot be checked against the spec text the way the §8 thresholds are.
    They are checked against the amendment that introduces them, which is the same guarantee at the
    same strength: a silent edit to either number breaks a test.
    """
    assert "## Amendment 3" in AMENDMENTS
    assert gc.CALIBRATION_EVAL_SIGMA == 0.5
    assert gc.POWER_EVAL_SIGMA == 0.35
    assert gc.ORACLE_SIM["donor_sigma"] == gc.CALIBRATION_EVAL_SIGMA

    # The envelope is Amendment 1's power frontier, adopted verbatim by Amendment 3 Change 1.
    frontier = [(r["sigma_donor"], r["min_donors_per_group"]) for r in gc.OPERATING_ENVELOPE]
    assert frontier == [(0.2, 4), (0.35, 8), (0.5, 13), (0.7, 23)]
    assert "n* = 4 / 8 / 13 / 23" in AMENDMENTS      # Amendment 1's own wording, still there
    # More donor variance can only ever demand more donors; a non-monotone envelope is a typo.
    assert frontier == sorted(frontier)

    # The structural invariant that makes the re-scoping honest: the point at which the power
    # criterion binds must itself lie INSIDE the declared envelope. If the gate's oracle had fewer
    # donors than the envelope demands at POWER_EVAL_SIGMA, the gate would be certifying a regime
    # its own envelope excludes.
    row = next(r for r in gc.OPERATING_ENVELOPE if r["sigma_donor"] == gc.POWER_EVAL_SIGMA)
    assert gc.ORACLE_SIM["n_donors_per_group"] >= row["min_donors_per_group"]


def test_the_two_permutation_counts_cannot_drift_apart():
    """``run_null`` pairs ``min(n_perm, n_perm_pb)``, so a larger ``N_PERM_PB`` alone does nothing.

    Amendment 3 Change 2 raised both to 200 for exactly this reason. Pinning them equal (rather than
    pinning the value, which is deliberately NOT pre-registered) guards the trap without freezing a
    convenience threshold.
    """
    assert gc.N_PERM == gc.N_PERM_PB


def test_the_probe_keeps_the_same_pre_registered_values():
    """The probe holds an independent copy on purpose; "independent" must not mean "divergent"."""
    path = REPO / "scripts" / "pb_calibration_probe.py"
    if not path.exists():  # pragma: no cover
        pytest.skip("probe script not present")
    spec = importlib.util.spec_from_file_location("_pb_probe_consts", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)

    assert mod.ALPHA == gc.ALPHA
    assert tuple(mod.LAMBDA_BAND) == gc.LAMBDA_BAND
    assert mod.POWER_TARGET == gc.POWER_TARGET
    assert mod.ORACLE_LOG2FC == gc.ORACLE_LOG2FC
    assert mod.ORACLE_K == gc.ORACLE_K


def test_the_two_blocks_stay_separate_in_the_manifest():
    """A convenience threshold must never be reportable as a pre-registered one."""
    m = gc.manifest()
    assert set(m["pre_registered"]) & set(m["instrument_sanity"]) == set()
    assert "NOT pre-registered" in m["instrument_sanity_source"]
    assert "AMENDMENTS" in m["pre_registered_source"]
    # The instrument-sanity naive threshold is looser than the spec's own; they must not be
    # confused, so assert the ordering that makes them distinguishable.
    assert gc.INSTRUMENT_LAMBDA_NAIVE_MIN < gc.LAMBDA_NAIVE_GO


def test_pre_registered_block_is_immutable():
    with pytest.raises(TypeError):
        gc.PRE_REGISTERED["alpha"] = 0.10
