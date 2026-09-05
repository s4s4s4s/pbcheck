"""pbcheck - an auditor of pseudoreplication in single-cell RNA-seq differential expression.

v0.1.0 ships the single-stratum audit (``pbcheck.audit``), its Markdown/HTML report
(``pbcheck.render``) and the ``pbcheck`` CLI, all outside the Phase 0 protocol. The Phase 0
real-data harness (``controls``, ``decision``, the spec section 9 ``report``) is specified but
not built; no ``risk_score`` exists and none is promised (see README's "What exists, and what
does not"). No ``py.typed`` marker: annotation coverage across the package is uneven and a
deliberate call, not an oversight.
"""

from pbcheck.design import DesignReport, audit_design
from pbcheck.gene_universe import UniverseTooSmall, frozen_universe
from pbcheck.methods import ebayes_from_pdata, naive_de, pseudobulk_de
from pbcheck.metrics import (
    BH_PAIRED,
    BH_SOLO,
    NDegSeries,
    empirical_perm_pvalues,
    genomic_inflation,
    jaccard,
    lambda_over_permutations,
    perm_floor,
    signal_above_floor,
)
from pbcheck.mtc import PairedBH, bh_both_arms, bh_over_universe
from pbcheck.permutation import build_perms, labels_for, run_null

__version__ = "0.1.0"

__all__ = [
    "__version__",
    # design.py — the metadata-only auditor
    "audit_design",
    "DesignReport",
    # gene_universe.py — the frozen, label-agnostic universe both arms share
    "frozen_universe",
    "UniverseTooSmall",
    # permutation.py — the donor-permutation null
    "run_null",
    "build_perms",
    "labels_for",
    # metrics.py — the inflation metrics (spec section 6)
    "genomic_inflation",
    "lambda_over_permutations",
    "empirical_perm_pvalues",
    "perm_floor",
    "signal_above_floor",
    "jaccard",
    # ...and the labelled #DEG series, so the two BH conventions cannot be pooled by accident
    "NDegSeries",
    "BH_PAIRED",
    "BH_SOLO",
    # mtc.py — the shared BH correction (spec section 5)
    "bh_both_arms",
    "PairedBH",
    "bh_over_universe",
    # methods — the naive arm and the current (moderated eBayes) pseudobulk arm
    "naive_de",
    "ebayes_from_pdata",
    "pseudobulk_de",
]
