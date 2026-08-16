"""Unit tests for the inflation metrics — fast, no DE runs."""

import numpy as np
import pytest

from pbcheck import metrics


def test_lambda_uniform_is_one():
    rng = np.random.default_rng(0)
    lam = metrics.genomic_inflation(rng.uniform(size=20000))
    assert abs(lam - 1.0) < 0.05, f"lambda on uniform p-values should be ~1, got {lam}"


def test_lambda_inflated_when_p_small():
    lam = metrics.genomic_inflation(np.full(1000, 1e-6))
    assert lam > 10, f"lambda on tiny p-values should be >> 1, got {lam}"


def test_lambda_ignores_nan():
    p = np.array([np.nan, 0.5, 0.5, 0.5, np.nan])
    lam = metrics.genomic_inflation(p)
    assert np.isfinite(lam)


def test_lambda_over_permutations_summary():
    rng = np.random.default_rng(1)
    mat = rng.uniform(size=(10, 2000))
    out = metrics.lambda_over_permutations(mat)
    assert out["n_perm"] == 10
    assert abs(out["lambda"] - 1.0) < 0.1


def test_lambda_over_permutations_empty_and_normal_paths_share_keys():
    """Regression: the empty/all-NaN path returned only 3 of the 4 keys the normal path returns.

    A caller that expects "per_perm" to always be present (as scripts/pb_calibration_probe.py's
    ``lam["lambda"]`` / ``lam.get("per_perm", [])`` mix does) could KeyError or silently branch on
    shape depending on which path ran. Both paths must expose the same key set.
    """
    rng = np.random.default_rng(2)
    normal = metrics.lambda_over_permutations(rng.uniform(size=(5, 100)))
    empty = metrics.lambda_over_permutations(np.empty((0, 100)))
    assert set(normal.keys()) == set(empty.keys())
    assert empty["n_perm"] == 0
    assert np.isnan(empty["lambda"])
    assert len(empty["per_perm"]) == 0


def test_empirical_perm_pvalues_are_uniform_under_exchangeability():
    """Spec B5's construction, pinned WITH its limitation (Amendment 2 Change 3).

    Under the sharp null the real labeling is exchangeable with every permuted one, so the rank of
    the real statistic is uniform BY CONSTRUCTION — for any test, however miscalibrated. That is
    why this is a check on the permutation machinery and NOT the calibration criterion, and it is
    pinned here so nobody later promotes it to one on the strength of it "looking calibrated".
    """
    rng = np.random.default_rng(4)
    G, n_perm = 4000, 200
    # A wildly anti-conservative arm: every p-value is tiny, real and permuted alike.
    null = rng.uniform(size=(n_perm, G)) * 1e-6
    real = rng.uniform(size=G) * 1e-6

    assert metrics.genomic_inflation(real) > 10          # the arm IS grossly inflated...
    p_emp = metrics.empirical_perm_pvalues(real, null)
    assert abs(metrics.genomic_inflation(p_emp) - 1.0) < 0.1   # ...and B5's lambda cannot see it

    assert p_emp.min() >= 1.0 / (1.0 + n_perm)           # resolution floor from the +1 correction
    assert p_emp.max() <= 1.0


def test_empirical_perm_pvalues_edges():
    real = np.array([0.001, 0.9, np.nan])
    null = np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]])
    p = metrics.empirical_perm_pvalues(real, null)
    assert p[0] == 1 / 4          # nothing in the null is as extreme as the real value
    assert p[1] == 4 / 4          # everything is
    assert np.isnan(p[2])         # an untested gene stays untested
    assert np.all(np.isnan(metrics.empirical_perm_pvalues(real, np.empty((0, 3)))))
    try:
        metrics.empirical_perm_pvalues(real, np.empty((5, 7)))
    except ValueError as e:
        assert "does not align" in str(e)
    else:
        raise AssertionError("a misaligned null matrix must raise, not broadcast")


def test_perm_floor_fraction():
    f = metrics.perm_floor(np.array([100, 120, 110, 90]), n_genes=1000)
    assert f["median_count"] == 105
    assert abs(f["median_frac"] - 0.105) < 1e-9
    # A bare array carries no convention, and the floor says so rather than guessing one.
    assert f["bh_mode"] is None


def test_ndeg_series_holds_one_bh_convention_and_cannot_be_extended():
    s = metrics.NDegSeries(np.array([3, 1, 2]), metrics.BH_PAIRED)
    assert len(s) == 3 and s.bh_mode == metrics.BH_PAIRED
    assert s.counts.dtype == np.int64
    # Read-only: the defect being prevented is exactly "append the other convention to this array".
    with pytest.raises(ValueError):
        s.counts[0] = 99
    with pytest.raises(ValueError):
        metrics.NDegSeries(np.array([1, 2]), "whatever")
    with pytest.raises(ValueError):
        metrics.NDegSeries(np.array([[1, 2], [3, 4]]), metrics.BH_SOLO)


def test_perm_floor_carries_the_bh_convention_and_refuses_a_contradiction():
    """The floor must say which BH produced it, and must not be relabelled at the call site.

    ``scripts/synthetic_gate.py`` used to take the headline naive floor over an array that was
    paired-BH below ``n_paired`` and solo-BH above it. Both halves are legitimate statistics; their
    median together is not one. Carrying the label into the floor dict is what makes the mistake
    visible in the run artifact instead of invisible in an intermediate array.
    """
    paired = metrics.NDegSeries(np.array([10, 12, 14]), metrics.BH_PAIRED)
    solo = metrics.NDegSeries(np.array([10, 12, 14, 100, 100]), metrics.BH_SOLO)

    f_paired = metrics.perm_floor(paired, n_genes=1000)
    f_solo = metrics.perm_floor(solo, n_genes=1000)
    assert f_paired["bh_mode"] == metrics.BH_PAIRED and f_paired["n_perm"] == 3
    assert f_solo["bh_mode"] == metrics.BH_SOLO and f_solo["n_perm"] == 5
    assert f_paired["median_count"] == 12 and f_solo["median_count"] == 14

    # Pooling the two would give yet a third number, which is the bug.
    pooled = metrics.perm_floor(np.concatenate([paired.counts, solo.counts[3:]]), n_genes=1000)
    assert pooled["median_count"] != f_paired["median_count"]

    with pytest.raises(ValueError, match="contradicts the series"):
        metrics.perm_floor(paired, n_genes=1000, bh_mode=metrics.BH_SOLO)
    with pytest.raises(ValueError, match="bh_mode must be one of"):
        metrics.perm_floor(np.array([1, 2]), n_genes=10, bh_mode="pared")
    # Restating the series' own label is allowed -- it is an assertion, not a contradiction.
    assert metrics.perm_floor(paired, 1000, bh_mode=metrics.BH_PAIRED)["bh_mode"] == "paired"


def test_signal_above_floor_near_one():
    assert abs(metrics.signal_above_floor(300, 295.0) - 300 / 295) < 1e-6


def test_jaccard():
    assert metrics.jaccard({"a", "b"}, {"b", "c"}) == 1 / 3
    assert np.isnan(metrics.jaccard(set(), set()))
