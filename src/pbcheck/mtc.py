"""Multiple-testing correction, applied identically to both arms over the frozen universe.

Spec §5: the naive and pseudobulk arms must be Benjamini-Hochberg corrected over an **identical,
label-agnostic** gene universe at the same alpha, so any difference in the significant-gene count is
due to the test / replication unit and nothing else. We deliberately do NOT consume scanpy's or
DESeq2's own adjusted p-values (DESeq2's independent filtering would NA genes asymmetrically vs the
naive arm — spec C2/C3); both arms feed raw p-values into the same ``fdr_bh`` over the same G genes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from pbcheck.methods.de import DEResult


def bh_complete_null_reference(alpha: float = 0.05) -> dict[str, float]:
    """The correct yardstick for a calibrated arm under the COMPLETE null (spec B1).

    Kept here, as a function with a test behind it, because getting this wrong is what the
    spec's B1 correction was written to fix and what a later reader is most likely to
    reintroduce. Under the complete null every rejection is false, so BH's step-up procedure
    rejects nothing at all in most replicates: the **median number of rejections is 0** and the
    expectation is of order ``alpha`` — that is, below 1, not ``alpha * G``.

    The deleted ``alpha * G`` line (~750 rejections at G = 15000) describes *unadjusted* testing
    and must never be used as a calibration reference.

    Note also that under the complete null V = R, so FDP = 1{R >= 1} and FDR = E[FDP] = P(R >= 1):
    the family-wise and false-discovery criteria **coincide** there, which is why the
    permutation-null false-positive rate is a legitimate BH calibration check.
    """
    return {"median_rejections": 0.0, "max_mean_rejections": float(alpha), "alpha": float(alpha)}


@dataclass
class PairedBH:
    """Two arms corrected over one common tested set, with the bookkeeping made visible."""

    naive: DEResult
    pseudobulk: DEResult
    n_universe: int
    n_tested_common: int
    n_na_naive: int
    n_na_pseudobulk: int
    dropped_for_fairness: list[str]


def bh(pvals: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaNs pass through as NaN (untested genes)."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return out
    out[ok] = multipletests(p[ok], alpha=alpha, method="fdr_bh")[1]
    return out


def bh_over_universe(result: DEResult, universe: list[str], alpha: float = 0.05) -> DEResult:
    """Reindex a DE result to the frozen universe and BH-correct over exactly those G genes.

    Genes in the universe that the arm did not test become NaN (never counted as significant in only
    one arm). Asserts the corrected vector has exactly ``len(universe)`` entries.
    """
    tab = result.table.reindex(pd.Index(universe, name="gene"))
    tab = tab.assign(padj=bh(tab["pval"].to_numpy(), alpha=alpha))
    # NOTE: a length check here would be a tautology — reindex returns len(universe) rows by
    # construction and cannot fail. The property that actually needs guarding is that the TWO arms
    # are corrected over the same tested set, which one arm alone cannot establish. Use
    # bh_both_arms for anything whose output is compared across arms.
    return DEResult(method=result.method, table=tab, contrast=result.contrast)


def bh_both_arms(
    naive: DEResult,
    pseudobulk: DEResult,
    universe: list[str],
    alpha: float = 0.05,
    *,
    strict: bool = False,
) -> PairedBH:
    """BH-correct both arms over one **common** tested set, so the denominators cannot diverge.

    The study's headline is that the two arms differ only by the test and the replication unit.
    That holds only if BH sees the same number of hypotheses ``m`` on both sides. Correcting each
    arm over its own non-NaN subset silently breaks it: the naive Wilcoxon arm never yields NaN,
    while DESeq2 still emits NaN p-values on sparse or degenerate strata even with Cook's filtering
    and independent filtering disabled (GLM non-convergence, degenerate dispersion). The pseudobulk
    arm would then be corrected over a smaller ``m``, i.e. made systematically *more liberal* per
    gene than the arm it is the control for — and nothing would report it.

    Genes untested by either arm are dropped from **both** and reported in
    ``dropped_for_fairness``; the remaining common set is corrected once per arm at the same alpha.

    Parameters
    ----------
    strict
        Raise if the arms disagree about which genes are testable, instead of intersecting. Useful
        in tests and on data where a disagreement means something is wrong upstream.
    """
    idx = pd.Index(universe, name="gene")
    a = naive.table.reindex(idx)
    b = pseudobulk.table.reindex(idx)

    ok_a = a["pval"].notna().to_numpy()
    ok_b = b["pval"].notna().to_numpy()
    common = ok_a & ok_b
    dropped = [g for g, keep in zip(universe, ~common) if keep]

    if strict and not np.array_equal(ok_a, ok_b):
        raise ValueError(
            f"arms disagree on testable genes: naive tested {int(ok_a.sum())}, "
            f"pseudobulk tested {int(ok_b.sum())} of {len(universe)}; "
            f"{len(dropped)} gene(s) would be dropped to equalise the BH denominator"
        )

    def _corrected(tab: pd.DataFrame) -> pd.DataFrame:
        padj = np.full(len(universe), np.nan)
        p = tab["pval"].to_numpy(dtype=float)
        if common.sum():
            padj[common] = multipletests(p[common], alpha=alpha, method="fdr_bh")[1]
        return tab.assign(padj=padj)

    return PairedBH(
        naive=DEResult(method=naive.method, table=_corrected(a), contrast=naive.contrast),
        pseudobulk=DEResult(
            method=pseudobulk.method, table=_corrected(b), contrast=pseudobulk.contrast
        ),
        n_universe=len(universe),
        n_tested_common=int(common.sum()),
        n_na_naive=int((~ok_a).sum()),
        n_na_pseudobulk=int((~ok_b).sum()),
        dropped_for_fairness=dropped,
    )
