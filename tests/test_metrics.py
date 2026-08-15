"""Unit tests for the inflation metrics — fast, no DE runs."""

import numpy as np

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


def test_perm_floor_fraction():
    f = metrics.perm_floor(np.array([100, 120, 110, 90]), n_genes=1000)
    assert f["median_count"] == 105
    assert abs(f["median_frac"] - 0.105) < 1e-9


def test_signal_above_floor_near_one():
    assert abs(metrics.signal_above_floor(300, 295.0) - 300 / 295) < 1e-6


def test_jaccard():
    assert metrics.jaccard({"a", "b"}, {"b", "c"}) == 1 / 3
    assert np.isnan(metrics.jaccard(set(), set()))
