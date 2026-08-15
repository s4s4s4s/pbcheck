"""Inflation metrics (spec §6).

The primary, size-robust severity scalar is the **genomic-inflation factor** lambda, borrowed from
GWAS. Given a set of p-values that ought to be Uniform(0,1) under the null, convert each to a 1-df
chi-square deviate and compare the observed median to the theoretical one:

    lambda = median( chi2.isf(p, df=1) ) / chi2.ppf(0.5, df=1)

A calibrated test gives lambda ~ 1; a test whose p-values are systematically too small (as the naive
per-cell test is, because it treats correlated cells as independent) gives lambda >> 1. We compute it
on the p-values each arm produces **under the donor-permutation null**, where the correct answer is no
association — so lambda measures each arm's calibration, not real biology. It is used as a
binary/ordinal flag (calibrated vs inflated), not as a cross-dataset magnitude (spec B3).

The within-stratum companion is the **permutation false-positive floor**: the naive #DEG under the
permutation null (a sharp-null upper bound on the pseudoreplication false positives, spec A1/B2).
"""

from __future__ import annotations

import numpy as np
from scipy.stats import chi2

_CHI2_MEDIAN_1DF = float(chi2.ppf(0.5, df=1))  # 0.4549364...


def genomic_inflation(pvals: np.ndarray) -> float:
    """Genomic-inflation factor lambda from a vector of (null) p-values. ~1 = calibrated."""
    p = np.asarray(pvals, dtype=float)
    p = p[~np.isnan(p)]
    if p.size == 0:
        return float("nan")
    # Guard the tails so isf stays finite.
    p = np.clip(p, np.nextafter(0, 1), 1.0)
    stats = chi2.isf(p, df=1)
    return float(np.median(stats) / _CHI2_MEDIAN_1DF)


def lambda_over_permutations(pval_matrix: np.ndarray) -> dict:
    """Lambda summarized across permutation-null realizations. **The binding calibration measure.**

    ``pval_matrix`` is (n_perm x G): each row is one arm's own p-values under one permuted
    labeling. Returns the median lambda across permutations plus its spread.

    This is the quantity decision rule item 1 / §8(a) binds on. Note what it is: the arm's
    *nominal* p-values evaluated where the correct answer is "no association". For the naive
    per-cell arm the nominal p-value is (asymptotically) the **cell**-permutation p-value, so its
    lambda against the **donor**-permutation null measures exactly how much wider the true
    donor-level null is than the cell-level null the naive test assumes — i.e. the clustering
    miscalibration.

    CITATION CORRECTED (Amendment 2 Change 3). This docstring previously cited "(spec B5)" for the
    median-across-permutations step. That was an overclaim: taking the median mitigates p-value
    discreteness in a single realization, but it is not B5's construction, which is
    :func:`empirical_perm_pvalues`. B5's tie/discreteness concern is met here by ``tie_correct=True``
    in the naive arm (pinned by spec §2) together with this median.
    """
    per_perm = np.array([genomic_inflation(row) for row in pval_matrix], dtype=float)
    per_perm = per_perm[~np.isnan(per_perm)]
    if per_perm.size == 0:
        # Same key set as the non-empty path (regression: this branch used to return only 3 of
        # the 4 keys, so callers indexing "per_perm" directly would KeyError on an empty/all-NaN
        # permutation set). per_perm is an empty array, not a dropped key, matching how
        # scripts/pb_calibration_probe.py already consumes it via lam.get("per_perm", []).
        return {
            "lambda": float("nan"),
            "lambda_iqr": float("nan"),
            "n_perm": 0,
            "per_perm": per_perm,
        }
    return {
        "lambda": float(np.median(per_perm)),
        "lambda_iqr": float(np.subtract(*np.percentile(per_perm, [75, 25]))),
        "n_perm": int(per_perm.size),
        "per_perm": per_perm,
    }


def empirical_perm_pvalues(real_pvals: np.ndarray, null_pval_matrix: np.ndarray) -> np.ndarray:
    """Spec B5's construction: rank the observed statistic within the donor-permutation null.

    Per gene, ``p_emp = (1 + #{permutations at least as extreme as the real labeling}) /
    (1 + n_perm)``; "at least as extreme" is read off the arm's own p-values (a smaller p is a more
    extreme statistic), so this works for any arm without re-deriving its statistic. The ``+1``
    is the standard bias correction — it keeps the estimator valid and bounds the minimum attainable
    p-value at ``1/(1 + n_perm)``, which is also its resolution limit.

    WHAT THIS DOES AND DOES NOT MEASURE — read before using it (Amendment 2 Change 3).
    Under the donor-permutation **sharp null** the real labeling is exchangeable with every permuted
    one, so this rank is uniform **by construction**, for any test, however miscalibrated. Measured
    on the null oracle while Amendment 2 was written: lambda from this construction = 0.93 on data
    where the arm's lambda under the permutation null was 26.08 — and it returned the identical 0.93
    when fed a held-out *permuted* labeling instead of the real one. It therefore carries **zero
    information about the arm's calibration** and must never be used as the calibration criterion,
    despite B5's wording; :func:`lambda_over_permutations` is that criterion.

    What it is good for, and what it is reported as: a validity check on the permutation machinery
    itself. It must sit near 1; a departure means the permutation set is not exchangeable with the
    real labeling (an unbalanced or leaked permutation construction), which would invalidate every
    other number taken from the same null.
    """
    real = np.asarray(real_pvals, dtype=float)
    null = np.asarray(null_pval_matrix, dtype=float)
    if null.ndim != 2 or null.shape[1] != real.shape[0]:
        raise ValueError(f"null matrix {null.shape} does not align with {real.shape[0]} genes")
    n_perm = null.shape[0]
    if n_perm == 0:
        return np.full(real.shape[0], np.nan)
    with np.errstate(invalid="ignore"):
        at_least_as_extreme = np.nansum(null <= real[None, :], axis=0)
    p = (1.0 + at_least_as_extreme) / (1.0 + n_perm)
    return np.where(np.isnan(real), np.nan, p)


def perm_floor(ndeg_per_perm: np.ndarray, n_genes: int) -> dict:
    """Permutation false-positive floor = #DEG under the null (sharp-null upper bound, spec §6).

    Returned as an absolute count and as a fraction of the frozen universe, with spread.
    """
    x = np.asarray(ndeg_per_perm, dtype=float)
    return {
        "median_count": float(np.median(x)),
        "iqr_count": float(np.subtract(*np.percentile(x, [75, 25]))),
        "median_frac": float(np.median(x) / n_genes) if n_genes else float("nan"),
        "max_count": float(np.max(x)),
        "n_perm": int(x.size),
    }


def signal_above_floor(real_ndeg: int, floor_median_count: float) -> float:
    """naive #DEG(real) / naive #DEG(perm-null median). ~1 => calls are essentially the floor.

    Decision-relevant ONLY in high-donor strata (>= 8v8); leak-contaminated at 3v3 (spec A1/§6).
    """
    denom = max(floor_median_count, 1.0)
    return float(real_ndeg / denom)


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)
