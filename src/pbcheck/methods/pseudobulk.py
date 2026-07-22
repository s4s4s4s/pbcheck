"""The correct pseudobulk DE test — the reference pbcheck compares the naive test against.

Counts are aggregated to one profile per donor (per cell type) with decoupler, then tested across
donors with PyDESeq2. The unit of replication is the donor, which is what removes the
pseudoreplication inflation.

Two spec corrections are baked in:
  * **poscounts size factors** (C1): median-of-ratios silently drops any gene with a zero in any
    sample and degenerates on sparse per-cell-type pseudobulk, so we use ``poscounts``.
  * **cooks_filter=False** (C2): DESeq2's Cook's-distance / independent filtering would set some genes'
    p-values to NA asymmetrically, while the naive Wilcoxon arm filters nothing — that would make the
    two arms' gene sets differ. We pull raw p-values with no such filter, then BH-correct both arms
    identically over the frozen universe in :mod:`pbcheck.mtc`.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import decoupler as dc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from pbcheck.methods.de import DEResult


def build_pseudobulk(
    adata,
    *,
    donor_col: str = "donor",
    celltype_col: str | None = "cell_type",
    condition_col: str = "condition",
    mode: str = "sum",
    min_cells: int = 10,
    min_counts: int = 1000,
):
    """Aggregate raw counts to one pseudobulk profile per donor (per cell type).

    ``mode='sum'`` (not mean) because DESeq2 models library size and expects summed integer counts.
    """
    pdata = dc.pp.pseudobulk(
        adata, sample_col=donor_col, groups_col=celltype_col, mode=mode, raw=False
    )
    if condition_col not in pdata.obs.columns:
        raise ValueError(
            f"'{condition_col}' not carried into pseudobulk .obs; it must be constant within donor."
        )
    return pdata


def deseq_from_pdata(
    pdata,
    *,
    condition_col: str = "condition",
    test_level: str = "disease",
    ref_level: str = "ctrl",
    universe: list[str] | None = None,
    condition_values=None,
    fdr: float = 0.05,
    n_cpus: int = 4,
) -> DEResult:
    """Run PyDESeq2 on already-aggregated donor pseudobulk profiles.

    Split out from :func:`pseudobulk_de` so a permutation harness can aggregate donors **once** and
    only relabel donor->condition per permutation (donor profiles are invariant under relabeling).

    ``condition_values`` optionally overrides the per-donor condition vector (used by permutations);
    it must align with ``pdata.obs_names``.
    """
    if universe is not None:
        keep = [g for g in universe if g in set(pdata.var_names)]
        pdata = pdata[:, keep].copy()
    if pdata.n_vars == 0:
        raise ValueError("No genes to test after universe restriction.")

    counts = pd.DataFrame(
        np.asarray(pdata.X).round().astype(int),
        index=pdata.obs_names,
        columns=pdata.var_names,
    )
    cond = (pd.Series(condition_values, index=pdata.obs_names) if condition_values is not None
            else pdata.obs[condition_col])
    metadata = pd.DataFrame({condition_col: cond.astype(str).to_numpy()}, index=pdata.obs_names)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design=f"~{condition_col}",
            ref_level=[condition_col, ref_level],
            size_factors_fit_type="poscounts",  # C1
            refit_cooks=True,
            quiet=True,
            n_cpus=n_cpus,
        )
        dds.deseq2()
        stats = DeseqStats(
            dds,
            contrast=[condition_col, test_level, ref_level],
            alpha=fdr,
            cooks_filter=False,  # C2: no asymmetric NA vs the naive arm
            independent_filter=False,
            quiet=True,
            n_cpus=n_cpus,
        )
        stats.summary()

    res = stats.results_df
    table = pd.DataFrame(
        {
            "pval": res["pvalue"].to_numpy(),
            "padj": np.nan,  # deferred to mtc.bh_over_universe (spec §5)
            "log2fc": res["log2FoldChange"].to_numpy(),
        },
        index=pd.Index(res.index, name="gene"),
    )
    return DEResult(
        method="pseudobulk[DESeq2:poscounts]",
        table=table,
        contrast=(condition_col, test_level, ref_level),
    )


def pseudobulk_de(
    adata,
    *,
    donor_col: str = "donor",
    celltype_col: str | None = "cell_type",
    condition_col: str = "condition",
    test_level: str = "disease",
    ref_level: str = "ctrl",
    universe: list[str] | None = None,
    min_count: int = 10,
    min_total_count: int = 15,
    fdr: float = 0.05,
    n_cpus: int = 4,
) -> DEResult:
    """Correct DE: pseudobulk per donor -> PyDESeq2 ``~condition`` across donors.

    If ``universe`` is given, the test is restricted to exactly those genes (the frozen label-agnostic
    universe, spec §3); otherwise a self-contained expression filter is applied for standalone use.
    Returns raw p-values + log2FC; call :func:`pbcheck.mtc.bh_over_universe` to get significance.
    """
    pdata = build_pseudobulk(
        adata, donor_col=donor_col, celltype_col=celltype_col, condition_col=condition_col
    )
    if universe is None:
        dc.pp.filter_by_expr(
            pdata, group=condition_col, min_count=min_count, min_total_count=min_total_count
        )
    return deseq_from_pdata(
        pdata, condition_col=condition_col, test_level=test_level, ref_level=ref_level,
        universe=universe, fdr=fdr, n_cpus=n_cpus,
    )
