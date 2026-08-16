"""pbcheck — an auditor of pseudoreplication in single-cell RNA-seq differential expression.

Phase 0 (pilot): the public surface is curated rather than exhaustive while the measurement
engine is validated on synthetic oracles — everything below is exercised by the test suite. The
rest of the real-data harness (``controls``, ``decision``, ``report``), the ``risk_score``, the
HTML report and the no-raw-counts mode are specified but not yet built (see
``docs/PHASE0_SPEC.md`` and README's "What exists, and what does not"). ``pbcheck.census_select``
(obs only) and ``pbcheck.io_counts`` (the counts gate) are built and are imported directly rather
than re-exported here: they are drivers for the real-data sweep, and between them they emit
*candidates* — not an admission, and not the §1 pre-registration of the stratum list. No
``py.typed`` marker: annotation coverage across the package is uneven and a deliberate call, not
an oversight.
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

__version__ = "0.0.1.dev0"

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
