"""The single, label-agnostic gene universe shared by both DE arms and every permutation.

Correction A3 in ``docs/PHASE0_SPEC.md``: if the tested gene set were built with the real condition
labels, it would shift under permutation and the naive-vs-pseudobulk comparison would not be on an
identical footing. So the universe is filtered on expression **pooled across all donors, ignoring
condition** — an edgeR ``filterByExpr``-style rule with no grouping. This frozen set is used by the
naive arm, the pseudobulk arm, and every permutation-null realization alike (spec §3, §5).
"""

from __future__ import annotations

import numpy as np


def frozen_universe(
    pdata,
    *,
    min_total_count: int = 15,
    min_prop: float = 0.5,
    min_size: int | None = None,
) -> list[str]:
    """Build the label-agnostic common universe from donor-pseudobulk profiles.

    Parameters
    ----------
    pdata
        Donor-pseudobulk ``AnnData`` (one row per donor per cell type), raw summed counts in ``.X``.
    min_total_count
        A gene must have at least this many total counts across all donors.
    min_prop
        A gene must be detected (count > 0) in at least this fraction of donors.
    min_size
        If given, **enforced**: raises :class:`UniverseTooSmall` when fewer genes survive (spec C5,
        inclusion gate item 5). Default ``None`` leaves the gate off, because the right place to
        SKIP a degenerate stratum is the stratum driver, which can record the skip in the manifest.
        It defaults to off rather than to 200 so that the gate is never applied silently — but a
        caller that passes a number gets a real check, not a decorative argument.

    Returns
    -------
    Sorted list of gene names forming the frozen universe.

    Raises
    ------
    UniverseTooSmall
        Only when ``min_size`` is given and the surviving set is smaller. DESeq2's dispersion trend
        degenerates on tiny gene sets, so such a stratum must be skipped rather than contributing a
        meaningless lambda and floor to the aggregate.
    """
    X = np.asarray(pdata.X)
    n_donors = X.shape[0]
    total = X.sum(axis=0)
    detected = (X > 0).sum(axis=0)
    keep = (total >= min_total_count) & (detected >= np.ceil(min_prop * n_donors))
    genes = sorted(str(g) for g, k in zip(pdata.var_names, keep) if k)
    if min_size is not None:
        require_min_size(genes, min_size)
    return genes


class UniverseTooSmall(ValueError):
    """Raised when the frozen universe is below the minimum viable gene count (spec C5)."""


def require_min_size(genes: list[str], min_size: int = 200) -> list[str]:
    if len(genes) < min_size:
        raise UniverseTooSmall(
            f"frozen universe has {len(genes)} genes (< {min_size}); stratum too sparse — SKIP."
        )
    return genes
