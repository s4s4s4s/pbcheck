"""The naive per-cell DE test — what pbcheck audits.

This treats individual cells as independent replicates (the pseudoreplication error): standard
``scanpy`` normalization + log1p, then ``rank_genes_groups`` across all cells of the two conditions.
It is deliberately the *default, careless* workflow, reproduced faithfully so the inflation we measure
is the inflation a real analyst would incur.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import scanpy as sc

from pbcheck.methods.de import DEResult


def naive_de(
    adata,
    *,
    condition_col: str = "condition",
    test_level: str = "disease",
    ref_level: str = "ctrl",
    genes=None,
    method: str = "wilcoxon",
    target_sum: float = 1e4,
) -> DEResult:
    """Run a per-cell DE test across ``condition``, cells as replicates.

    Parameters
    ----------
    adata
        Cells x genes with **raw counts** in ``.X`` and ``condition_col`` in ``.obs``.
    genes
        Optional gene universe to restrict to (for a fair comparison against pseudobulk, pass the
        genes that survived pseudobulk expression filtering).
    """
    a = adata if genes is None else adata[:, list(genes)]
    a = a.copy()
    sc.pp.normalize_total(a, target_sum=target_sum)
    sc.pp.log1p(a)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            a, groupby=condition_col, groups=[test_level], reference=ref_level, method=method
        )
    df = sc.get.rank_genes_groups_df(a, group=test_level)
    table = pd.DataFrame(
        {
            "pval": df["pvals"].to_numpy(),  # raw; BH deferred to mtc so both arms match (spec §5)
            "padj": np.nan,
            "log2fc": df["logfoldchanges"].to_numpy(),
        },
        index=pd.Index(df["names"], name="gene"),
    )
    return DEResult(
        method=f"naive[{method}]",
        table=table,
        contrast=(condition_col, test_level, ref_level),
    )
