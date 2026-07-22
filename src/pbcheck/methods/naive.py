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

    # The spec (§2) pins this call verbatim, including two arguments scanpy defaults to False.
    # tie_correct matters because scRNA-seq counts are heavily tied (mostly zeros), and without it
    # the Wilcoxon normal approximation is distorted — which would land in lambda_naive as a tie
    # artifact rather than as the clustering miscalibration we are trying to measure. It is only
    # defined for the rank test, so it is not passed to the t-test robustness variant.
    # pts adds the per-group expressed fractions, which the min-expression sensitivity check needs.
    kw = {"pts": True}
    if method == "wilcoxon":
        kw["tie_correct"] = True

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            a, groupby=condition_col, groups=[test_level], reference=ref_level, method=method, **kw
        )
    df = sc.get.rank_genes_groups_df(a, group=test_level)
    cols = {
        "pval": df["pvals"].to_numpy(),  # raw; BH deferred to mtc so both arms match (spec §5)
        "padj": np.nan,
        "log2fc": df["logfoldchanges"].to_numpy(),
    }
    for src, dst in (("pct_nz_group", "pct_group"), ("pct_nz_reference", "pct_reference")):
        if src in df.columns:
            cols[dst] = df[src].to_numpy()
    table = pd.DataFrame(cols, index=pd.Index(df["names"], name="gene"))
    return DEResult(
        method=f"naive[{method}]",
        table=table,
        contrast=(condition_col, test_level, ref_level),
    )
