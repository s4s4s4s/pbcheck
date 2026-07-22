"""Multiple-testing correction, applied identically to both arms over the frozen universe.

Spec §5: the naive and pseudobulk arms must be Benjamini-Hochberg corrected over an **identical,
label-agnostic** gene universe at the same alpha, so any difference in the significant-gene count is
due to the test / replication unit and nothing else. We deliberately do NOT consume scanpy's or
DESeq2's own adjusted p-values (DESeq2's independent filtering would NA genes asymmetrically vs the
naive arm — spec C2/C3); both arms feed raw p-values into the same ``fdr_bh`` over the same G genes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from pbcheck.methods.de import DEResult


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
    assert len(tab) == len(universe), "BH vector length must equal the frozen universe size"
    return DEResult(method=result.method, table=tab, contrast=result.contrast)
