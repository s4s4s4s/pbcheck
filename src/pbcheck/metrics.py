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
    """Lambda summarized across permutation-null realizations.

    ``pval_matrix`` is (n_perm x G): each row is one arm's p-values under one permuted labeling.
    Returns the median lambda across permutations plus its spread (robust to Wilcoxon p-value
    discreteness in any single realization, spec B5).
    """
    per_perm = np.array([genomic_inflation(row) for row in pval_matrix], dtype=float)
    per_perm = per_perm[~np.isnan(per_perm)]
    if per_perm.size == 0:
        return {"lambda": float("nan"), "lambda_iqr": float("nan"), "n_perm": 0}
    return {
        "lambda": float(np.median(per_perm)),
        "lambda_iqr": float(np.subtract(*np.percentile(per_perm, [75, 25]))),
        "n_perm": int(per_perm.size),
        "per_perm": per_perm,
    }


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
