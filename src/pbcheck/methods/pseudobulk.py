"""The correct pseudobulk DE test — the reference pbcheck compares the naive test against.

Counts are aggregated to one profile per donor (per cell type) with decoupler, then tested across
donors. The unit of replication is the donor, which is what removes the pseudoreplication inflation.

**The arm's test is moderated eBayes** (:mod:`pbcheck.methods.moderated`) as of Amendment 2
(2026-08-15) Change 1, which closes the replacement-test question Amendment 1 left open.
DESeq2-Wald — which Amendment 1 showed fails the binding validity gate of decision rule item 1 /
§8(a) — is retired from the arm; :func:`deseq_from_pdata` is retained, superseded, so that every
Amendment-1-era number stays reproducible from this repository.

Spec corrections still live here:
  * **thin-donor filtering** (§1 inclusion-gate item 2, §3): donor profiles below ``min_cells``
    cells or ``min_counts`` summed counts are **dropped, not merged**. Implemented manually per
    Amendment 2 Change 7 — decoupler 2.x removed the parameters the spec pins, so this filter had
    silently never run.
  * **cooks_filter=False** (C2) on the retained DESeq2 path. C2's underlying requirement — no
    asymmetric NA between the arms — now binds for *any* source of NA and is carried by
    :func:`pbcheck.mtc.bh_both_arms` (Amendment 2 Change 2), which is strictly stronger.
  * **poscounts size factors** (C1) on the retained DESeq2 path. C1 is DESeq2-specific and lapses
    with it; the moderated arm normalises by library size directly through log2(CPM + 1).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import decoupler as dc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from pbcheck.methods.de import DEResult
from pbcheck.methods.moderated import ebayes_from_pdata

#: Pre-registered thin-donor thresholds (spec §1 inclusion-gate item 2 and the §3 aggregation call).
#: Changing them is a protocol change and requires an AMENDMENTS.md entry.
MIN_CELLS = 10
MIN_COUNTS = 1000


def build_pseudobulk(
    adata,
    *,
    donor_col: str = "donor",
    celltype_col: str | None = "cell_type",
    condition_col: str = "condition",
    mode: str = "sum",
    min_cells: int = MIN_CELLS,
    min_counts: int = MIN_COUNTS,
):
    """Aggregate raw counts to one pseudobulk profile per donor (per cell type), then thin-filter.

    ``mode='sum'`` (not mean) because the downstream models library size and expect summed integer
    counts.

    **Thin-donor filter (Amendment 2 Change 7).** Spec §1 inclusion-gate item 2 pins "≥ 10 cells per
    donor in that cell_type (donors below threshold are **dropped, not merged**)" and §3 pins
    ``dc.pp.pseudobulk(..., min_cells=10, min_counts=1000)``. decoupler 2.x removed both parameters,
    so from the port onwards the pre-registered filter silently never ran; R0 documented that hole
    and deferred the fix to an amendment, which is Change 7. The filter is applied here on
    decoupler's own bookkeeping columns (``psbulk_cells`` / ``psbulk_counts``) after aggregation and
    before the universe is frozen, exactly where the spec places it. A profile below either
    threshold is **removed**; nothing is merged into another donor and nothing is back-filled.

    The surviving profiles carry ``pdata.uns['thin_donor_filter']`` recording the thresholds, the
    counts dropped for each reason and the dropped profile names, so a stratum that loses donors is
    visible in the artifact rather than silently smaller. The spec's separate "≥ 3 pseudosamples per
    group post-aggregation, else SKIP" rule (§3) applies to what survives and is enforced by the
    caller.

    Passing ``min_cells=0, min_counts=0`` disables the filter; it is not a way to soften the
    protocol, only to inspect what the filter removed.
    """
    pdata = dc.pp.pseudobulk(
        adata, sample_col=donor_col, groups_col=celltype_col, mode=mode, raw=False
    )
    if condition_col not in pdata.obs.columns:
        raise ValueError(
            f"'{condition_col}' not carried into pseudobulk .obs; it must be constant within donor."
        )
    return _drop_thin_donors(pdata, min_cells=min_cells, min_counts=min_counts)


def _drop_thin_donors(pdata, *, min_cells: int, min_counts: int):
    """Drop (never merge) donor profiles below the pre-registered thin-sample thresholds.

    ``psbulk_cells`` / ``psbulk_counts`` are decoupler's own per-profile bookkeeping. Where they are
    absent (a hand-built ``pdata`` in a test, say) the corresponding threshold cannot be applied and
    is recorded as not applied, rather than being silently treated as satisfied — the failure mode
    this whole change exists to remove.
    """
    n_before = pdata.n_obs
    keep = np.ones(n_before, dtype=bool)
    applied = {}

    for col, thresh, label in (("psbulk_cells", min_cells, "min_cells"),
                               ("psbulk_counts", min_counts, "min_counts")):
        if thresh <= 0:
            applied[label] = "disabled"
            continue
        if col not in pdata.obs.columns:
            applied[label] = "unavailable: decoupler did not emit " + col
            continue
        vals = pd.to_numeric(pdata.obs[col], errors="coerce").to_numpy(dtype=float)
        fails = ~(vals >= thresh)  # NaN fails: an unmeasurable profile is not a passing one
        applied[label] = int(fails.sum())
        keep &= ~fails

    dropped = [str(n) for n, k in zip(pdata.obs_names, keep) if not k]
    out = pdata[keep].copy() if not keep.all() else pdata
    out.uns["thin_donor_filter"] = {
        "min_cells": int(min_cells),
        "min_counts": int(min_counts),
        "n_profiles_before": int(n_before),
        "n_profiles_after": int(out.n_obs),
        "n_dropped": int(n_before - out.n_obs),
        "dropped_profiles": dropped,
        "per_threshold": applied,
        "semantics": "dropped, not merged (spec section 1 item 2; Amendment 2 Change 7)",
    }
    return out


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

    **SUPERSEDED by Amendment 2 (2026-08-15) Change 1 — kept for reproducibility of
    Amendment-1-era results.** This is no longer the pseudobulk arm's test: Amendment 1 showed
    DESeq2-Wald fails the binding validity gate (λ = 1.21–1.25 outside [0.9, 1.1]; permutation-null
    FP 0.35–0.50 against α = 0.05; power 0.47–0.57 against ≥ 0.60), and the 146-cell selection grid
    reproduced the failure independently (FP = 0.26 and 0.56 at 8v8, binomial P = 1.1e-17 and
    5.8e-68). The arm's test is now :func:`pbcheck.methods.moderated.ebayes_from_pdata`.

    It is retained, not deleted, because every number in Amendment 1 and ``PILOT_FINDINGS.md`` was
    produced by this function: retiring a method by removing the code that produced the retired
    results would make the record unfalsifiable. New callers must use the moderated arm.

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
    min_cells: int = MIN_CELLS,
    min_counts: int = MIN_COUNTS,
    trend: bool = False,
) -> DEResult:
    """Correct DE: pseudobulk per donor -> **moderated eBayes** ``~condition`` across donors.

    The arm's test as of Amendment 2 Change 1 (Smyth 2004 moderated t on log2(CPM + 1) donor
    profiles; see :mod:`pbcheck.methods.moderated` for the estimator and the selection evidence).
    Thin donor profiles are dropped before the test per Amendment 2 Change 7.

    If ``universe`` is given, the test is restricted to exactly those genes (the frozen
    label-agnostic universe, spec §3); otherwise a self-contained expression filter is applied for
    standalone use. Returns raw p-values + log2FC with the realised prior on ``.moderation``; call
    :func:`pbcheck.mtc.bh_both_arms` (not the per-arm path) for anything compared across arms.
    """
    pdata = build_pseudobulk(
        adata, donor_col=donor_col, celltype_col=celltype_col, condition_col=condition_col,
        min_cells=min_cells, min_counts=min_counts,
    )
    if universe is None:
        dc.pp.filter_by_expr(
            pdata, group=condition_col, min_count=min_count, min_total_count=min_total_count
        )
    return ebayes_from_pdata(
        pdata, condition_col=condition_col, test_level=test_level, ref_level=ref_level,
        universe=universe, trend=trend,
    )
