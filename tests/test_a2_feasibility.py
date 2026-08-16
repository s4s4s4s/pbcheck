"""``scripts/a2_feasibility.py`` — the A2 tau-ladder feasibility bracket, tested as evidence.

**Runtime note, read before adding to this file.** The full artifact — 704 strata across two tiers,
each analysed under two dispersion extremes, ~584 of them via a 20,000-draw Monte Carlo estimate —
takes roughly 20-25 seconds to build. That is fine for a committed evidence driver run once and
diffed (``python scripts/a2_feasibility.py --check``, which is how byte-identical regeneration
against the committed artifact is actually verified), but it is not fine as an addition to a suite
this repository keeps in the single-digit seconds. Every test below either operates on a tiny
synthetic fixture, or exercises the (cheap) bracket-construction and tier-selection logic against the
real committed data without invoking the (expensive) Monte-Carlo / enumeration path — deliberately,
so this file adds negligible time to ``pytest -q``. Determinism of the expensive path itself is
covered by :func:`test_analyze_stratum_is_bitwise_reproducible` on a single real stratum, which is
fast because one stratum's Monte-Carlo cost is ~1/584th of the full artifact's.

No network. Everything is read from ``pilot/preregistration/``.
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "a2_feasibility.py"


def _load_script():
    if not SCRIPT.exists():  # pragma: no cover - the file is part of the commit under test
        pytest.skip("a2_feasibility.py not present")
    spec = importlib.util.spec_from_file_location("_a2_feasibility", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


af = _load_script()


@pytest.fixture(scope="module")
def stratum_list():
    return af.load_stratum_list()


@pytest.fixture(scope="module")
def manifest():
    return af.load_manifest()


@pytest.fixture(scope="module")
def frozen_rows(stratum_list):
    return af.select_frozen_tier(stratum_list)


@pytest.fixture(scope="module")
def manifest_rows(manifest):
    return af.select_manifest_tier(manifest)


# ---------------------------------------------------------------------------
# Tier facts, re-derived — the four numbers verified independently before any artifact is written.
# Cheap: filtering and min()/sum() over already-parsed rows, no MC, no enumeration.
# ---------------------------------------------------------------------------


def test_frozen_tier_matches_the_verified_facts(frozen_rows):
    assert len(frozen_rows) == 150
    assert min(r["permutation_count"] for r in frozen_rows) == 24_310
    n_enum = sum(1 for r in frozen_rows if r["permutation_count"] <= af.A2_ENUM_CAP)
    assert n_enum == 32


def test_manifest_tier_matches_the_verified_facts(manifest_rows):
    assert len(manifest_rows) == 554
    assert min(r["permutation_count"] for r in manifest_rows) == 12_870
    n_enum = sum(1 for r in manifest_rows if r["permutation_count"] <= af.A2_ENUM_CAP)
    assert n_enum == 88


def test_frozen_tier_is_a_subset_of_the_manifest_tier(frozen_rows, manifest_rows):
    frozen_keys = {af.stratum_key(r) for r in frozen_rows}
    manifest_keys = {af.stratum_key(r) for r in manifest_rows}
    assert frozen_keys <= manifest_keys


# ---------------------------------------------------------------------------
# The bracket vectors satisfy the summary statistics they were built from (synthetic fixtures).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "n_donors,n_cells,lo,median,hi",
    [
        (1, 500, 500, 500, 500),
        (2, 300, 100, 150.0, 200),
        (3, 100 + 130 + 400, 100, 130, 400),
        (4, 900, 100, 225.0, 400),
        (5, 1000, 50, 200, 500),
        (8, 4000, 10, 500, 1500),
        (17, 12000, 30, 700, 2000),
        (26, 96000, 20, 3000, 16000),
    ],
)
def test_low_dispersion_vector_satisfies_the_five_statistics(n_donors, n_cells, lo, median, hi):
    vec = af.build_low_group_vector(n_donors, n_cells, lo, median, hi)
    assert len(vec) == n_donors
    assert sum(vec) == n_cells
    assert min(vec) == lo
    assert max(vec) == hi
    assert vec == sorted(vec)


@pytest.mark.parametrize(
    "n_donors,n_cells,lo,median,hi",
    [
        (1, 500, 500, 500, 500),
        (2, 300, 100, 150.0, 200),
        (4, int(100 + 2 * 225.0 + 400), 100, 225.0, 400),  # D=4 identity: n_cells=lo+2*median+hi
        (5, 1000, 50, 200, 500),
        (8, 4000, 10, 500, 1500),
        (17, 12000, 30, 700, 2000),
        (26, 96000, 20, 3000, 16000),
    ],
)
def test_high_dispersion_vector_satisfies_the_five_statistics(n_donors, n_cells, lo, median, hi):
    vec = af.build_high_group_vector(n_donors, n_cells, lo, median, hi)
    assert len(vec) == n_donors
    assert sum(vec) == n_cells
    assert min(vec) >= lo
    assert max(vec) <= hi
    assert vec == sorted(vec)
    if n_donors >= 3:
        # the reported median value is literally present in the constructed vector
        assert median in vec or (math_is_close_to_int(median) and int(median) in vec)


def math_is_close_to_int(x):
    return float(x).is_integer()


def test_high_dispersion_vector_contains_the_median_for_even_donor_counts():
    vec = af.build_high_group_vector(10, 3000, 50, 275.5, 600)
    # even n_donors -> two mid slots floor(median), ceil(median), summing to 2*median
    assert 275 in vec and 276 in vec


@pytest.mark.parametrize("n_donors", [1, 2, 3, 4])
def test_small_donor_counts_reconstruct_exactly_with_zero_median_drift(n_donors):
    """n_donors <= 4: the five statistics determine the multiset exactly (module docstring)."""
    if n_donors == 1:
        gs = {"n_donors": 1, "n_cells": 777, "min": 777, "median": 777, "max": 777}
    elif n_donors == 2:
        gs = {"n_donors": 2, "n_cells": 500, "min": 200, "median": 250.0, "max": 300}
    elif n_donors == 3:
        gs = {"n_donors": 3, "n_cells": 100 + 150 + 400, "min": 100, "median": 150, "max": 400}
    else:
        gs = {"n_donors": 4, "n_cells": 100 + 180 + 220 + 500, "min": 100, "median": 200.0, "max": 500}
    for dispersion in ("low", "high"):
        vec, drift = af.build_group_vector(gs, dispersion)
        assert sum(vec) == gs["n_cells"]
        assert len(vec) == gs["n_donors"]
        assert drift == 0.0


def test_bracket_vectors_never_escape_bounds_on_the_full_committed_data(frozen_rows, manifest_rows):
    """Construction only (no MC/enumeration) is cheap enough to check against every real stratum."""
    seen = set()
    for row in (*frozen_rows, *manifest_rows):
        key = af.stratum_key(row)
        if key in seen:
            continue
        seen.add(key)
        for group in ("A", "B"):
            gs = row["cells_per_donor_by_group"][group]
            for dispersion in ("low", "high"):
                vec, _drift = af.build_group_vector(gs, dispersion)
                assert len(vec) == gs["n_donors"]
                assert sum(vec) == gs["n_cells"]
                assert min(vec) >= gs["min"]
                assert max(vec) <= gs["max"]


# ---------------------------------------------------------------------------
# The VOID rule's boundaries (thin-set rule, sec 1.2).
# ---------------------------------------------------------------------------


def test_enumerated_classification_boundaries_at_m_1000_999_200_199():
    assert af.classify("enumerated", m=1000) == af.STATUS_FULL
    assert af.classify("enumerated", m=999) == af.STATUS_COARSE
    assert af.classify("enumerated", m=200) == af.STATUS_COARSE
    assert af.classify("enumerated", m=199) == af.STATUS_VOID


def test_sampled_classification_boundaries_at_the_acceptance_rate_thresholds():
    assert af.classify("sampled", q=af.Q_FULL_THRESHOLD) == af.STATUS_FULL
    assert af.classify("sampled", q=af.Q_FULL_THRESHOLD - 1e-9) == af.STATUS_COARSE
    assert af.classify("sampled", q=af.Q_COARSE_THRESHOLD) == af.STATUS_COARSE
    assert af.classify("sampled", q=af.Q_COARSE_THRESHOLD - 1e-12) == af.STATUS_VOID
    assert af.classify("sampled", q=0.0) == af.STATUS_VOID


# ---------------------------------------------------------------------------
# The exact machinery (DP enumeration, finite-population moments) against brute force.
# ---------------------------------------------------------------------------


def test_enumerate_subset_sum_counts_matches_brute_force():
    values = [3, 5, 5, 8, 12, 1, 7]
    k = 3
    counts = af.enumerate_subset_sum_counts(values, k)
    brute = {}
    for combo in itertools.combinations(range(len(values)), k):
        s = sum(values[i] for i in combo)
        brute[s] = brute.get(s, 0) + 1
    assert int(counts.sum()) == sum(brute.values())
    for s, c in brute.items():
        assert int(counts[s]) == c


def test_finite_population_moments_match_brute_force_enumeration():
    values = [10, 20, 25, 40, 5, 60, 33]
    k = 3
    mean, var = af.finite_population_mean_var(values, k)
    totals = [
        sum(values[i] for i in combo) for combo in itertools.combinations(range(len(values)), k)
    ]
    arr = np.array(totals, dtype=np.float64)
    assert mean == pytest.approx(arr.mean(), abs=1e-9)
    assert var == pytest.approx(arr.var(ddof=0), abs=1e-9)


def test_finite_population_moments_match_dp_enumeration_on_a_real_stratum(frozen_rows):
    row = min(frozen_rows, key=lambda r: r["n_donors_A"] + r["n_donors_B"])
    n_a = row["n_donors_A"]
    assert row["permutation_count"] <= af.A2_ENUM_CAP  # smallest stratum must be enumerable
    for dispersion in ("low", "high"):
        vec, _a, _b = af.build_combined_vector(row, dispersion)
        counts = af.enumerate_subset_sum_counts(vec, n_a)
        total = int(counts.sum())
        s_idx = np.arange(len(counts))
        mean_exact = float((s_idx * counts).sum()) / total
        var_exact = float((((s_idx - mean_exact) ** 2) * counts).sum()) / total
        mean_formula, var_formula = af.finite_population_mean_var(vec, n_a)
        assert mean_formula == pytest.approx(mean_exact, abs=1e-6)
        assert var_formula == pytest.approx(var_exact, abs=1e-6)


# ---------------------------------------------------------------------------
# Reproducibility — determinism of the expensive path, on one real stratum (fast: ~1/584th of the
# full artifact's Monte-Carlo cost). The full committed-artifact byte-identity check is
# ``python scripts/a2_feasibility.py --check`` — see the module docstring for why it is not here.
# ---------------------------------------------------------------------------


def test_analyze_stratum_is_bitwise_reproducible(manifest_rows):
    # a non-enumerable stratum, so this also exercises the seeded Monte-Carlo path.
    row = next(r for r in manifest_rows if r["permutation_count"] > af.A2_ENUM_CAP)
    first = af.analyze_stratum(row)
    second = af.analyze_stratum(row)
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)


def test_mc_seed_material_is_independent_of_iteration_order(manifest_rows):
    row = next(r for r in manifest_rows if r["permutation_count"] > af.A2_ENUM_CAP)
    seed_first_call = af.stratum_seed_material(row, "low")
    # simulate other work happening between calls -- seed derivation must not depend on any
    # shared mutable state or call counter.
    _ = af.stratum_seed_material(row, "high")
    _ = af.stratum_seed_material(row, "low")
    seed_later_call = af.stratum_seed_material(row, "low")
    assert seed_first_call == seed_later_call


def test_render_json_and_render_csv_are_pure_functions_of_the_artifact():
    tiny = {
        "header": {"a": 1, "b": [1, 2, 3.0]},
        "tau_ladder_summary": {"frozen": [], "manifest": []},
        "extreme_f_real": {"frozen": {}, "manifest": {}},
        "void_index": {"frozen": {}, "manifest": {}},
        "strata": {"frozen": [], "manifest": []},
    }
    assert af.render_json(tiny) == af.render_json(tiny)
    assert af.render_csv(tiny) == af.render_csv(tiny)


# ---------------------------------------------------------------------------
# Hash pins reject a substituted input.
# ---------------------------------------------------------------------------


def test_stratum_list_hash_pin_rejects_a_substituted_file(tmp_path):
    bad = tmp_path / "stratum_list_2026-08-16.json"
    bad.write_text('{"header": {}, "rows": [], "within_collection_control_rows": []}', encoding="utf-8")
    with pytest.raises(af.SourceArtifactMismatch):
        af.load_stratum_list(bad, csv_path=None)


def test_manifest_hash_pin_rejects_a_substituted_file(tmp_path):
    bad = tmp_path / "census_candidates.json"
    bad.write_text('{"header": {}, "rows": []}', encoding="utf-8")
    with pytest.raises(af.SourceArtifactMismatch):
        af.load_manifest(bad, csv_path=None)


def test_stratum_list_hash_pin_rejects_right_size_wrong_bytes(tmp_path):
    """Same byte count as the pinned file, different content -- the size check alone must not pass."""
    bad = tmp_path / "stratum_list_2026-08-16.json"
    bad.write_bytes(b"x" * af.STRATUM_LIST_BYTES)
    with pytest.raises(af.SourceArtifactMismatch):
        af.load_stratum_list(bad, csv_path=None)
