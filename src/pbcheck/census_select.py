"""Real-data stratum selection over the CELLxGENE Census — the ``.obs``-only first module of the
Phase 0 harness (spec §1; §9 item 1).

What it does: opens the Census **at the pinned version and only that version**, reads the ``obs``
frame the spec pins, forms the candidate strata — ``(dataset_id × cell_type)``, one binary
``<disease> vs normal`` contrast per disease term, **never** a pooled "any disease" arm (§1, Unit of
analysis) — applies the parts of the §1 inclusion gate that are decidable from metadata, runs the
donor-level confound pre-screen (§1, §7.1, D3/D4) and writes an auditable candidate manifest.

What it reads: ``obs`` and nothing else. Counts (``X``) are ``io_counts.py``'s job (§9 item 2), which
is also where the two gate items that need counts get settled — see the pending fields below. A
module that reads no counts can screen the whole Census cheaply, which is the point.

**Nothing here pre-registers anything.** §1 requires the stratum list to be pre-registered before any
metric is computed; this module *proposes* candidates. Freezing a list is a separate, deliberate act
(a commit of the manifest, plus the §1 "span the outcome space" curation of 8–12 datasets, which is a
judgement about biology and not a filter). Every row carries ``admitted_to_sweep=False``.

**What this does NOT settle** (each one is a real hole, named rather than papered over):

* **Inclusion-gate items 4 and 5** — raw integer counts and frozen-universe size ≥ C5 — cannot be
  evaluated from ``obs``. They are emitted as the explicit ``pending`` columns ``integer_check`` and
  ``frozen_universe_size``, to be filled by ``io_counts.py`` / :mod:`pbcheck.gene_universe`. No value
  is guessed, and a row with pending columns is not a passing row.
* **``sigma_donor`` and envelope membership** (Amendment 3 Change 1). Stratum inclusion in the real
  sweep now additionally requires a per-stratum ``sigma_donor`` estimate and membership in
  :data:`pbcheck.gate_config.OPERATING_ENVELOPE` at that stratum's donors-per-group. Amendment 3
  supplies a *mechanism* (the moderated arm's ``s0^2``, times ln 2) and explicitly **not** an anchor:
  that quantity is an *upper bound* on ``donor_sigma``, the conversion is unvalidated, and "the
  ``sigma_donor`` anchoring demanded by Amendment 1 therefore remains OPEN". So
  ``sigma_donor_estimate`` and ``envelope_min_donors_per_group`` are pending here.
  ``envelope_max_sigma_supported`` is the one envelope number this module *can* compute honestly —
  the largest ``sigma_donor`` whose donor requirement this stratum's donor count already meets, i.e.
  a statement of the form "inside the envelope **if** its unmeasured σ is at most this" — and it is
  never a pass: no stratum may be admitted "on the strength of this entry alone".
* **Cell-type ontology depth** (D5). ``obs`` exposes the ontology *term id*, not its depth in the CL
  graph; computing depth needs the ontology itself, which this module deliberately does not fetch.
  The term id is recorded so ``controls.annotation_depth()`` (§9 item 9) can fill the depth later.
* **Pooling, where the Census exposes no library/pool identifier.** D3 requires pooled designs to be
  flagged and the donor-pseudobulk-is-calibrated claim restricted to non-pooled data. Where no
  library/pool/multiplex column exists in ``obs`` the answer is not "not pooled" — it is
  :data:`POOLING_UNRESOLVED`, which is D3's "where pooling cannot be resolved" state and carries its
  "donor-pseudobulk is a lower bound" caveat.
* **Nothing about carrying a confounded covariate into the design.** §1 would admit one "only if
  C4's df rule allows", but Amendment 2 Change 1 retired DESeq2 from the arm and recorded that "C1
  and C4 are DESeq2-specific and lapse with it"; the shipped moderated arm fits ``~ 1 + x`` and has
  no covariate slot. The manifest therefore carries ``residual_df`` (§1's own column, and the
  moderated arm's own ``d = n − 2``) and deliberately carries **no** covariate-affordability column:
  publishing one would advertise a mechanism this package does not currently have. A partially
  confounded stratum is tagged in ``covariate_candidates`` and is neutralised only by the
  permutation null (§7.1), which is what the tag means and all it means.
* **Whether any stratum is analyzable at all.** This module reports; it does not conclude. D4's
  excluded fraction and its reasons are aggregated into the manifest header for exactly that reason.

Donor-level confounding is measured with :mod:`pbcheck.design`'s helpers, imported rather than
reimplemented: the unit of observation for a condition-vs-covariate association is the **donor**, and
that decision, its rationale and its regression test already live there (see
``tests/test_design.py::test_confound_is_measured_per_donor_not_per_cell``). A second copy of that
logic here is exactly how the two would drift apart.
"""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _package_version
from pathlib import Path
from types import MappingProxyType

import pandas as pd

from pbcheck import gate_config as gc
from pbcheck.design import _cramers_v, _donor_level_table, _perfectly_separates

# ---------------------------------------------------------------------------
# Frozen by spec §1. Changing any of these is a protocol change (docs/AMENDMENTS.md).
# ---------------------------------------------------------------------------

#: The Census version, pinned by §1: ``open_soma(census_version="2025-01-30")`` — and the spec's own
#: next four words are "never ``latest``". A Census "version" that resolves differently on two days
#: is not a version, and every downstream number would silently belong to a different dataset.
CENSUS_VERSION = "2025-01-30"

#: The two aliases the Census resolves *at call time*. They are named here so they can be REJECTED,
#: which is the only place either string may legitimately appear in this package.
MUTABLE_CENSUS_ALIASES = ("latest", "stable")

#: obs ``value_filter`` verbatim from §1.
VALUE_FILTER = "is_primary_data == True and organism == 'Homo sapiens' and disease != 'na'"

CENSUS_ORGANISM = "homo_sapiens"

#: Group B of every contrast (§1: "Group A = one disease term; Group B = ``normal``").
REFERENCE_LEVEL = "normal"

#: The obs columns §1 materializes.
OBS_COLUMNS = (
    "dataset_id", "donor_id", "disease", "cell_type", "assay", "suspension_type",
    "tissue_general", "self_reported_ethnicity", "sex", "development_stage",
)

#: Without these four there is no stratum and no contrast, so their absence is fatal rather than
#: recorded. The rest of :data:`OBS_COLUMNS` is reported as missing and the run continues.
REQUIRED_OBS_COLUMNS = ("dataset_id", "donor_id", "disease", "cell_type")

#: Requested when present, not required: the ontology term id (D5, depth filled later) and the
#: per-cell count summaries the Census computes (used for the §1 sequencing-depth bin and
#: ``median_counts_per_cell_by_group``).
OBS_COLUMNS_OPTIONAL = ("cell_type_ontology_term_id", "raw_sum", "nnz", "soma_joinid")

#: Candidate names for the "library/pool/batch/multiplex identifier (D3)" §1 asks for *if the Census
#: exposes one*. As of the pinned version its obs schema exposes none of these, so the honest outcome
#: on real data is :data:`POOLING_UNRESOLVED` rather than "not pooled"; the list and the pattern
#: exist so that a dataset- or version-specific column is picked up rather than missed.
POOL_COLUMN_CANDIDATES = (
    "library_id", "library", "pool_id", "pool", "multiplex_id", "batch", "batch_id",
    "sample_id", "lane", "hashing_id",
)
_POOL_COLUMN_PATTERN = re.compile(r"(library|pool|multiplex|lane|hashtag|hashing|batch)", re.I)
_NEVER_POOL_COLUMNS = frozenset(OBS_COLUMNS) | frozenset(OBS_COLUMNS_OPTIONAL)

#: Per-cell columns the Census may expose that stand in for sequencing depth; the first present one
#: drives the §1 pre-screen's depth bin. Absent → the depth bin is pending, never invented.
DEPTH_COLUMN_CANDIDATES = ("raw_sum", "nnz", "raw_mean_nnz")

#: The only one of those that is a **count**. ``median_counts_per_cell_by_group`` is computed from
#: this column and from no other: ``nnz`` is a count of detected genes, not of molecules, and
#: writing it into a column named "counts" would be a fabricated number wearing the right label.
COUNTS_COLUMN = "raw_sum"

DEPTH_BIN_COL = "sequencing_depth_bin"
N_DEPTH_BINS = 3

#: The derived two-level column every downstream module keys on (§9 item 2's ``attach_obs``).
CONDITION_COL = "condition"

#: §1 inclusion gate item 1: "≥ 3 distinct ``donor_id`` in EACH group after filtering". Matches
#: :func:`pbcheck.design.audit_design`'s default ``min_donors``; the thin-donor threshold of item 2
#: is taken from :data:`pbcheck.gate_config.MIN_CELLS` rather than restated here.
MIN_DONORS_PER_GROUP = 3

#: §1 pre-screen: "perfectly separated by assay/suspension/pool (V ≈ 1) → EXCLUDE". ``≈`` is read as
#: "separates, or V at least this" — separation is checked structurally (a perfectly separating but
#: unbalanced table can score below 1, see :func:`pbcheck.design._perfectly_separates`), and the V
#: cut is the belt to that braces.
PERFECT_CONFOUND_V = 0.999

#: Tag threshold for "nearly confounded", NOT an exclusion rule and NOT pre-registered: §1 pins only
#: the V ≈ 1 exclusion. Same value as :mod:`pbcheck.design`'s flag threshold, on purpose.
PARTIAL_CONFOUND_V = 0.8

#: Covariates whose perfect separation from condition EXCLUDES the stratum. §1 names exactly three —
#: "assay/suspension/pool" — even though the pre-screen measures five. ``tissue_general`` and the
#: sequencing-depth bin are therefore screened and TAGGED but never excluded on: this is the reading
#: of §1 as written, recorded here because the alternative (excluding on all five) would silently
#: drop strata the frozen spec keeps. Depth in particular is addressed by §7.2 through
#: depth-matching and the permutation null, not by exclusion.
EXCLUDING_COVARIATES = ("assay", "suspension_type")

#: The sentinel written into a manifest column this module cannot fill. It is a string, not ``None``,
#: so that "not computed here" can never be read as "computed, and empty".
PENDING = "pending"

#: Third state for ``pooled_flag``, distinct from both True and False (D3).
POOLING_UNRESOLVED = "unresolved"

#: Every pending column, with the module that owes it. An artifact that says which of its holes are
#: known is auditable; one that writes zeros is not.
PENDING_FIELDS = MappingProxyType({
    "cell_type_ontology_depth":
        "controls.annotation_depth() (§9 item 9, D5) — obs carries the CL term id, not its depth",
    "median_counts_per_cell_by_group":
        "io_counts.py — unless the Census obs exposes its per-cell count column (raw_sum), in "
        "which case it is computed here from obs; nnz counts genes, not molecules, and is never "
        "substituted for it",
    "integer_check":
        "io_counts.assert_integer_counts() — inclusion-gate item 4; Census raw is not guaranteed "
        "integer (§10 risk 6) and integrality is a property of X, not of obs",
    "frozen_universe_size":
        "gene_universe.frozen_universe() with MIN_UNIVERSE_SIZE — inclusion-gate item 5 / C5",
    "fitType":
        "the pseudobulk arm — C5. Under Amendment 2 the arm is moderated eBayes, whose recorded "
        "analog of DESeq2's fitType is the realised prior (prior df / shrinkage)",
    "sigma_donor_estimate":
        "OPEN — Amendment 3 Change 1 supplies a mechanism (sqrt(s0^2) * ln 2 from the moderated "
        "fit) and states it is an UPPER BOUND, not an estimate; the anchor is not delivered",
    "envelope_min_donors_per_group":
        "requires sigma_donor_estimate; the envelope demand is a function of that unmeasured knob",
    "envelope_membership":
        "requires sigma_donor_estimate — 'no stratum may be admitted to the real sweep on the "
        "strength of this entry alone' (Amendment 3 Change 1)",
})

#: One row per stratum-contrast; the §1 manifest columns verbatim, plus the Amendment 3 envelope
#: columns and the bookkeeping that makes an exclusion auditable (D4).
MANIFEST_FIELDS = (
    "dataset_id",
    "cell_type",
    "cell_type_ontology_term_id",
    "cell_type_ontology_depth",
    "disease",
    "reference",
    "n_donors_A",
    "n_donors_B",
    "n_cells",
    "cells_per_donor_by_group",
    "median_counts_per_cell_by_group",
    "permutation_count",
    "confound_flags",
    "confound_cramers_v",
    "pooled_flag",
    "residual_df",
    "integer_check",
    "frozen_universe_size",
    "fitType",
    "sigma_donor_estimate",
    "envelope_min_donors_per_group",
    "envelope_max_sigma_supported",
    "envelope_membership",
    "gate_status",
    "gate_failures",
    "dropped_thin_donors",
    "dropped_missing_donor_cells",
    "confound_exclusion_reasons",
    "covariate_candidates",
    "admitted_to_sweep",
    "admission_blockers",
    "census_version",
)

#: Packages whose exact versions §1 requires in the run manifest.
MANIFEST_PACKAGES = (
    "cellxgene-census", "anndata", "scanpy", "decoupler", "pydeseq2",
    "statsmodels", "numpy", "scipy", "pandas",
)

_NON_DISEASE_TERMS = frozenset({"na"})


class CensusNotInstalled(ImportError):
    """``cellxgene-census`` is an optional extra and is not installed."""


@dataclass(frozen=True)
class Contrast:
    """One binary ``disease vs normal`` comparison inside one ``(dataset_id × cell_type)`` stratum.

    §1: "If a dataset has > 2 disease terms, run one binary ``disease-vs-normal`` contrast per term;
    never pool into 'any disease'." Each term therefore gets its own :class:`Contrast`, and the two
    share nothing but the stratum they sit in.
    """

    dataset_id: str
    cell_type: str
    disease: str
    reference: str = REFERENCE_LEVEL

    @property
    def key(self) -> str:
        return f"{self.dataset_id}|{self.cell_type}|{self.disease}_vs_{self.reference}"


@dataclass
class GateResult:
    """Outcome of the obs-decidable half of the §1 inclusion gate for one contrast."""

    contrast: Contrast
    obs: pd.DataFrame  # surviving cells, with the derived CONDITION_COL; thin donors already dropped
    n_donors: dict[str, int]
    n_cells: dict[str, int]
    cells_per_donor: dict[str, dict[str, float]]
    median_counts_per_cell: dict[str, float] | str
    dropped_thin_donors: dict[str, int] = field(default_factory=dict)  # donor -> cells it had
    dropped_missing_donor_cells: int = 0  # cells whose donor_id is null — §1 item 3, "present"
    failures: list[str] = field(default_factory=list)
    pending: dict[str, str] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        """Passed *the checks obs can decide*. Items 4 and 5 are pending, not passed — see
        :data:`PENDING_FIELDS` and the module docstring."""
        return not self.failures


@dataclass
class ConfoundReport:
    """Donor-level Cramér's V + perfect separation + pooling for one gated stratum (§1, D3, D4).

    ``cramers_v`` and ``separates`` are measured over the covariate's **donor-sharing levels only**
    (those held by ≥ 2 donors) — see :func:`confound_prescreen`. ``None`` / ``False`` therefore mean
    "no association measurable among the levels that carry information beyond ``donor_id``", which
    is not the same as "measured and zero"; ``flags`` records how many levels were set aside.
    """

    cramers_v: dict[str, float]
    separates: dict[str, bool]
    flags: list[str]
    exclude: bool
    exclusion_reasons: list[str]
    covariate_candidates: list[str]
    pooled: bool | str
    pool_columns: list[str]
    pool_evidence: dict[str, dict[str, int]]
    pending: dict[str, str] = field(default_factory=dict)


def open_census(*, census_version: str = CENSUS_VERSION):
    """Open the Census at the pinned version (§1). Lazy-imports the ``[census]`` extra.

    The version parameter exists to be *checked*, not to be varied: anything other than
    :data:`CENSUS_VERSION` raises, and the two mutable aliases raise with their own message, because
    a run whose Census version is decided at call time cannot be reproduced. Pointing this package at
    a different Census release is a protocol change and goes through ``docs/AMENDMENTS.md`` first.
    """
    if census_version in MUTABLE_CENSUS_ALIASES:
        raise ValueError(
            f"census_version={census_version!r} is a MUTABLE alias that resolves differently over "
            f"time; spec §1 pins {CENSUS_VERSION!r} and forbids it explicitly."
        )
    if census_version != CENSUS_VERSION:
        raise ValueError(
            f"census_version={census_version!r} is not the pinned version {CENSUS_VERSION!r} "
            "(spec §1). Changing it is a protocol change and requires a docs/AMENDMENTS.md entry "
            "written and committed first."
        )
    try:
        import cellxgene_census
    except ImportError as exc:  # pragma: no cover - exercised via a stubbed sys.modules in tests
        raise CensusNotInstalled(
            "cellxgene-census is not installed. It is an optional extra on purpose — the whole "
            'synthetic engine validation runs without it. Install with: pip install -e ".[census]"'
        ) from exc
    return cellxgene_census.open_soma(census_version=census_version)


def _schema_names(soma_obs) -> list[str]:
    schema = getattr(soma_obs, "schema", None)
    names = getattr(schema, "names", None)
    return [str(n) for n in names] if names else []


def detect_pool_columns(names) -> list[str]:
    """Library/pool/batch/multiplex identifiers among ``names`` (D3), by name.

    Explicit candidates first, then a pattern, because the Census schema is not ours to control and
    a column we fail to notice becomes a pooled design silently scored as non-pooled — which is
    precisely the failure D3 exists to prevent. Nothing in :data:`OBS_COLUMNS` can be picked up here.
    """
    names = [str(n) for n in names]
    hits = [n for n in names if n in POOL_COLUMN_CANDIDATES and n not in _NEVER_POOL_COLUMNS]
    hits += [
        n for n in names
        if n not in hits and n not in _NEVER_POOL_COLUMNS and _POOL_COLUMN_PATTERN.search(n)
    ]
    return hits


def _depth_column(columns) -> str | None:
    for c in DEPTH_COLUMN_CANDIDATES:
        if c in set(columns):
            return c
    return None


def query_obs(
    census,
    *,
    organism: str = CENSUS_ORGANISM,
    value_filter: str = VALUE_FILTER,
    columns=None,
    extra_columns=(),
) -> pd.DataFrame:
    """Read the §1 obs slice into a DataFrame. Reads obs only — never ``X`` (that is ``io_counts``).

    Columns §1 names but the Census does not expose are recorded in ``frame.attrs['missing_columns']``
    and left out of the request rather than faked; the four columns without which there is no stratum
    (:data:`REQUIRED_OBS_COLUMNS`) raise instead. Any library/pool identifier the schema does expose
    is added to the request automatically (D3).

    ``frame.attrs`` carries the provenance the manifest header needs: census version, the filter
    string actually used, requested/missing columns, detected pool columns and the depth column.
    """
    soma_obs = census["census_data"][organism].obs
    available = _schema_names(soma_obs)
    wanted = list(dict.fromkeys([*(columns or OBS_COLUMNS), *OBS_COLUMNS_OPTIONAL, *extra_columns]))

    if available:
        missing_required = [c for c in REQUIRED_OBS_COLUMNS if c not in available]
        if missing_required:
            raise ValueError(
                f"the Census obs schema at version {CENSUS_VERSION} lacks {missing_required}; "
                "stratum formation (dataset_id x cell_type, donor, disease) is impossible without it"
            )
        pool_columns = detect_pool_columns(available)
        requested = [c for c in [*wanted, *pool_columns] if c in available]
        missing = [c for c in wanted if c not in available]
    else:  # a handle that does not expose a schema: ask for what the spec names and find out.
        pool_columns = []
        requested = wanted
        missing = []

    frame = soma_obs.read(column_names=requested, value_filter=value_filter).concat().to_pandas()
    frame.attrs.update({
        "census_version": CENSUS_VERSION,
        "organism": organism,
        "value_filter": value_filter,
        "requested_columns": list(requested),
        "missing_columns": [c for c in missing if c in OBS_COLUMNS],
        "missing_optional_columns": [c for c in missing if c not in OBS_COLUMNS],
        "pool_columns": [c for c in pool_columns if c in frame.columns],
        "depth_column": _depth_column(frame.columns),
    })
    return frame


def _require_columns(obs: pd.DataFrame, cols) -> None:
    missing = [c for c in cols if c not in obs.columns]
    if missing:
        raise ValueError(f"obs is missing required column(s) {missing}")


def build_contrasts(
    obs: pd.DataFrame,
    *,
    dataset_col: str = "dataset_id",
    celltype_col: str = "cell_type",
    disease_col: str = "disease",
    reference: str = REFERENCE_LEVEL,
) -> list[Contrast]:
    """Candidate stratum-contrasts: ``(dataset_id × cell_type)`` × one disease term vs ``normal``.

    A stratum with no ``normal`` cells yields nothing — there is no within-stratum control to compare
    against, and borrowing controls from another dataset would confound the contrast with the dataset
    itself (D2). A stratum with three disease terms yields three separate contrasts, never one pooled
    "any disease" arm (§1).

    Cells with a null ``dataset_id``, ``cell_type`` or ``disease`` define no stratum and no contrast,
    so they take part in none (under pandas 3.0 a null survives ``.astype(str)`` rather than becoming
    the string "nan", which would otherwise produce a "nan" stratum).
    """
    _require_columns(obs, [dataset_col, celltype_col, disease_col])
    keys = (
        obs[[dataset_col, celltype_col, disease_col]]
        .dropna()
        .astype(str)
        .drop_duplicates()
    )
    contrasts: list[Contrast] = []
    for (dataset_id, cell_type), sub in keys.groupby([dataset_col, celltype_col], sort=True):
        terms = sorted(set(sub[disease_col]))
        if reference not in terms:
            continue
        for term in terms:
            if term == reference or term in _NON_DISEASE_TERMS:
                continue
            contrasts.append(Contrast(str(dataset_id), str(cell_type), term, reference))
    return contrasts


def _per_group_stats(sub: pd.DataFrame, donor_col: str, counts_col: str | None):
    per_donor = sub.groupby(donor_col, observed=True).size()
    stats = {
        "n_donors": int(per_donor.size),
        "n_cells": int(sub.shape[0]),
        "median": float(per_donor.median()) if per_donor.size else 0.0,
        "min": int(per_donor.min()) if per_donor.size else 0,
        "max": int(per_donor.max()) if per_donor.size else 0,
    }
    counts = float(sub[counts_col].median()) if counts_col and sub.shape[0] else None
    return stats, counts


def apply_inclusion_gate(
    obs: pd.DataFrame,
    contrast: Contrast,
    *,
    dataset_col: str = "dataset_id",
    celltype_col: str = "cell_type",
    disease_col: str = "disease",
    donor_col: str = "donor_id",
    min_cells_per_donor: int = gc.MIN_CELLS,
    min_donors_per_group: int = MIN_DONORS_PER_GROUP,
) -> GateResult:
    """The obs-decidable part of the §1 inclusion gate, for one contrast.

    Applied in the spec's own order, after one prerequisite: **cells whose ``donor_id`` is null are
    removed first** and counted in ``dropped_missing_donor_cells``. Item 3 begins with "``donor_id``
    present", and a cell without a donor has no replication unit — it can be neither thin-filtered
    nor permuted. It must be dropped explicitly rather than left to a null-skipping aggregation.

    1. **item 2 first** — donors with fewer than ``min_cells_per_donor`` cells *in this cell type*
       are **dropped, not merged**, and the dropped ids and their cell counts are kept in the result
       so a stratum that lost donors is visible rather than merely smaller;
    2. **item 1** — ≥ ``min_donors_per_group`` distinct ``donor_id`` in EACH group *after* that drop;
    3. **item 3** — ``donor_id`` present, non-constant, and nested within condition. A donor with
       cells under both conditions breaks nesting (this is not a paired design and the permutation
       unit would be ill-defined); a group carrying exactly one donor is the "donor collinear with
       condition, n = 1" case §1 names explicitly. That last check is strictly weaker than item 1 and
       is reported separately anyway, so the manifest says *which* rule a stratum died on.

    Cells are counted **inside the contrast**: a donor's cells belonging to some third disease term
    in the same stratum are not tested by this contrast and so do not help it clear item 2.

    Items 4 (integer counts) and 5 (frozen-universe size) are **not** evaluated — obs cannot decide
    them — and are returned in ``pending``. ``passed`` therefore means "passed everything obs can
    check", never "included".
    """
    _require_columns(obs, [dataset_col, celltype_col, disease_col])
    levels = [contrast.disease, contrast.reference]
    mask = (
        (obs[dataset_col].astype(str) == contrast.dataset_id)
        & (obs[celltype_col].astype(str) == contrast.cell_type)
        & (obs[disease_col].astype(str).isin(levels))
    )
    sub = obs.loc[mask].copy()
    sub[CONDITION_COL] = sub[disease_col].astype(str)
    counts_col = COUNTS_COLUMN if COUNTS_COLUMN in sub.columns else None

    pending = {
        "integer_check": PENDING_FIELDS["integer_check"],
        "frozen_universe_size": PENDING_FIELDS["frozen_universe_size"],
    }
    empty = {"A": 0, "B": 0}

    if donor_col not in sub.columns:
        return GateResult(
            contrast=contrast, obs=sub.iloc[:0], n_donors=dict(empty), n_cells=dict(empty),
            cells_per_donor={}, median_counts_per_cell=PENDING,
            failures=["donor_id_missing"], pending=pending,
        )

    # §1 item 3 opens with "donor_id present", and a cell whose donor is null has no replication
    # unit: it cannot be a donor, cannot be thin-filtered, and cannot be permuted. It must be
    # removed explicitly. Under pandas 3.0 (``future.infer_string``) ``.astype(str)`` does NOT
    # render null as the string "nan" — the null survives, and then ``value_counts`` / ``nunique``
    # skip it while ``isin(thin)`` is False for it, so such cells used to pass the gate unseen and
    # travel on into the stratum. The count is kept, like the thin-donor drop, so a stratum that
    # lost cells is visible rather than merely smaller.
    missing_donor = sub[donor_col].isna()
    dropped_missing = int(missing_donor.sum())
    if dropped_missing:
        sub = sub.loc[~missing_donor].copy()

    donors = sub[donor_col].astype(str)
    sub[donor_col] = donors
    per_donor = donors.value_counts()
    thin = per_donor[per_donor < min_cells_per_donor]
    dropped = {str(d): int(n) for d, n in thin.items()}
    sub = sub.loc[~sub[donor_col].isin(thin.index)].copy()

    group_of = {contrast.disease: "A", contrast.reference: "B"}
    n_donors: dict[str, int] = {}
    n_cells: dict[str, int] = {}
    cells_per_donor: dict[str, dict[str, float]] = {}
    median_counts: dict[str, float] = {}
    for level, tag in group_of.items():
        part = sub.loc[sub[CONDITION_COL] == level]
        stats, counts = _per_group_stats(part, donor_col, counts_col)
        n_donors[tag] = stats["n_donors"]
        n_cells[tag] = stats["n_cells"]
        cells_per_donor[tag] = stats
        if counts is not None:
            median_counts[tag] = counts

    failures: list[str] = []
    donor_condition = sub[[donor_col, CONDITION_COL]].astype(str).drop_duplicates()
    n_distinct_donors = int(donor_condition[donor_col].nunique())

    if n_distinct_donors == 0 and dropped_missing and not dropped:
        failures.append(
            f"donor_id is null on all {dropped_missing} cell(s) of the stratum "
            "(§1 inclusion gate item 3 requires donor_id present)"
        )
    elif n_distinct_donors == 0:
        failures.append(f"no donor clears the >= {min_cells_per_donor}-cells-per-donor rule")
    elif n_distinct_donors < 2:
        failures.append("donor_id is constant (1 donor in the whole stratum)")

    spanning = donor_condition.groupby(donor_col)[CONDITION_COL].nunique()
    if n_distinct_donors and bool((spanning > 1).any()):
        failures.append(
            "donor_id is not nested within condition: "
            f"{int((spanning > 1).sum())} donor(s) carry cells under both conditions"
        )

    for tag in ("A", "B"):
        if n_donors[tag] == 1:
            failures.append(f"group {tag} has exactly 1 donor — donor collinear with condition (n=1)")
    for tag in ("A", "B"):
        if n_donors[tag] < min_donors_per_group:
            failures.append(
                f"group {tag} has {n_donors[tag]} donor(s) after the thin-donor drop "
                f"(< {min_donors_per_group}; §1 inclusion gate item 1)"
            )

    return GateResult(
        contrast=contrast,
        obs=sub,
        n_donors=n_donors,
        n_cells=n_cells,
        cells_per_donor=cells_per_donor,
        median_counts_per_cell=median_counts if median_counts else PENDING,
        dropped_thin_donors=dropped,
        dropped_missing_donor_cells=dropped_missing,
        failures=failures,
        pending=pending,
    )


def _depth_bins(sub: pd.DataFrame, donor_col: str, depth_col: str) -> pd.Series:
    """Donor-level sequencing-depth tertiles, as a per-cell column.

    The bin is a property of the DONOR (its median per-cell depth), not of the cell, because the
    association it feeds is measured per donor — see :func:`pbcheck.design._donor_level_table`.
    """
    donor_depth = sub.groupby(donor_col, observed=True)[depth_col].median()
    if donor_depth.nunique() < 2:
        codes = pd.Series(0, index=donor_depth.index)
    else:
        n_bins = min(N_DEPTH_BINS, int(donor_depth.nunique()))
        codes = pd.Series(pd.qcut(donor_depth, q=n_bins, labels=False, duplicates="drop"),
                          index=donor_depth.index).fillna(0)
    labels = codes.map(lambda c: f"depth_bin_{int(c) + 1}")
    return sub[donor_col].map(labels)


def confound_prescreen(
    obs: pd.DataFrame,
    *,
    condition_col: str = CONDITION_COL,
    donor_col: str = "donor_id",
    covariate_cols=("assay", "suspension_type", "tissue_general"),
    pool_cols=None,
    depth_col=None,
    partial_v: float = PARTIAL_CONFOUND_V,
) -> ConfoundReport:
    """§1's confound pre-screen, at the DONOR level: Cramér's V + perfect separation + pooling flag.

    Screened: ``assay``, ``suspension_type``, ``tissue_general``, the sequencing-depth bin and any
    library/pool id — the five §1 lists. Perfect separation (or V ≥ :data:`PERFECT_CONFOUND_V`) on
    **assay, suspension_type or a pool id** EXCLUDES the stratum: with condition determined by the
    covariate, an inflation reading cannot be attributed to pseudoreplication rather than to the
    covariate, so the stratum is not merely noisy but uninterpretable. Separation on
    ``tissue_general`` or on the depth bin is TAGGED and not excluded — §1's exclusion clause names
    three covariates, not five, and depth is handled by §7.2's depth-matching plus the permutation
    null instead (see :data:`EXCLUDING_COVARIATES`).

    Partial confounding is retained and tagged in ``covariate_candidates``, and that is a **tag, not
    a plan**: §1 would carry such a covariate into the DESeq2 design "only if C4's df rule allows",
    but Amendment 2 Change 1 retired DESeq2 from the arm and states that "C1 and C4 are
    DESeq2-specific and lapse with it". The shipped moderated arm fits ``~ 1 + x`` and has no
    covariate slot at all, so a partially confounded stratum is neutralised **only** by the
    permutation null (§7.1: "dropped and neutralized by the permutation null, which breaks
    donor↔condition"), not modelled. The tag records what was seen; nothing downstream acts on it.

    Association and separation are measured **only over levels that group two or more donors**. A
    level held by a single donor is a relabeling of ``donor_id``, which §1 item 3 already requires
    to be nested within condition, so "separation" by it is the design working as specified rather
    than a confound; the count set aside is reported in ``flags`` (see the inline note below).

    Pooling (D3): donors of the stratum sharing a pool/library id → ``pooled=True``. With no pool
    column exposed at all the value is :data:`POOLING_UNRESOLVED`, never ``False``.

    The association is measured once per donor throughout. Cell-weighting it would let a
    many-celled donor dominate — the cells-as-replicates fallacy this package exists to detect.
    """
    _require_columns(obs, [condition_col])
    work = obs.copy()
    pending: dict[str, str] = {}

    if donor_col not in work.columns:
        raise ValueError(f"obs is missing required column '{donor_col}'")
    if bool(work[donor_col].isna().any()):
        # Never silently: a null donor would be dropped by every groupby here and would quietly
        # change the denominator of an association measured per donor.
        raise ValueError(
            f"'{donor_col}' is null on {int(work[donor_col].isna().sum())} cell(s); run "
            "apply_inclusion_gate() first — it removes them and records the count (§1 item 3)"
        )
    work[donor_col] = work[donor_col].astype(str)

    if pool_cols is None:
        pool_cols = [c for c in (obs.attrs.get("pool_columns") or []) if c in work.columns]
        if not pool_cols:
            pool_cols = detect_pool_columns(work.columns)
    pool_cols = [c for c in pool_cols if c in work.columns]

    depth_col = depth_col or obs.attrs.get("depth_column") or _depth_column(work.columns)
    screened = [c for c in covariate_cols if c in work.columns]
    if depth_col and work.shape[0]:
        work[DEPTH_BIN_COL] = _depth_bins(work, donor_col, depth_col)
        screened.append(DEPTH_BIN_COL)
    else:
        pending[DEPTH_BIN_COL] = (
            "no per-cell count column (raw_sum/nnz) in this obs frame; the sequencing-depth bin "
            "of the §1 pre-screen is deferred to io_counts.py, which sees X"
        )
    screened += [c for c in pool_cols if c not in screened]
    if not pool_cols:
        pending["pooled_flag"] = (
            "no library/pool/multiplex identifier exposed in obs — D3's 'where pooling cannot be "
            "resolved' state: donor-pseudobulk is a LOWER BOUND on the replication unit here"
        )

    cramers: dict[str, float] = {}
    separates: dict[str, bool] = {}
    flags: list[str] = []
    exclusion_reasons: list[str] = []
    covariate_candidates: list[str] = []

    excluding = set(EXCLUDING_COVARIATES) | set(pool_cols)
    for col in screened:
        table = _donor_level_table(work, donor_col, [col, condition_col])

        # A level held by exactly one donor — a per-donor library id, a tertile bin over 3 donors —
        # is a relabeling of that donor, and donor IS nested within condition by design (§1 item 3
        # requires it). Such a level "separates" condition trivially, so it is set aside BEFORE the
        # association is measured rather than after: doing it after (testing the full table and
        # only skipping when EVERY level is donor-unique) let two donors who happen to share one
        # library and one condition, among 38 donors on unique libraries, score V = 1.0 and exclude
        # the stratum. Real Census library ids are near-donor-unique, so that would have thrown out
        # usable strata systematically. Only levels that group >= 2 donors carry information beyond
        # donor_id, and they are what the screen is computed on.
        donors_per_level = table.groupby(col)[donor_col].nunique()
        shared_levels = donors_per_level[donors_per_level > 1].index
        n_singleton = int((donors_per_level == 1).sum())
        informative = table.loc[table[col].isin(shared_levels)]

        if informative.empty:
            cramers[col] = None
            separates[col] = False
            flags.append(f"{col}:donor-unique levels — no information beyond donor_id, not screened")
            continue
        if n_singleton:
            flags.append(
                f"{col}:{n_singleton} donor-unique level(s) set aside before screening; "
                f"{len(shared_levels)} level(s) group >= 2 donors and carry the test"
            )

        v = _cramers_v(informative[col], informative[condition_col])
        sep = _perfectly_separates(informative[col], informative[condition_col])
        cramers[col] = None if math.isnan(v) else round(float(v), 3)
        separates[col] = bool(sep)

        confounded = sep or (not math.isnan(v) and v >= PERFECT_CONFOUND_V)
        if confounded:
            flags.append(f"{col}:perfect_separation")
            if col in excluding:
                exclusion_reasons.append(
                    f"condition is perfectly separated by '{col}' at the donor level — inflation "
                    "here is uninterpretable, not merely confounded (§1 pre-screen)"
                )
        elif not math.isnan(v) and v >= partial_v:
            flags.append(f"{col}:near_confound_v={cramers[col]}")
            covariate_candidates.append(col)
        elif not math.isnan(v) and v > 0:
            flags.append(f"{col}:partial_confound_v={cramers[col]}")
            covariate_candidates.append(col)

    pool_evidence: dict[str, dict[str, int]] = {}
    for col in pool_cols:
        table = work[[donor_col, col]].astype(str).drop_duplicates()
        shared = table.groupby(col)[donor_col].nunique()
        shared = shared[shared > 1]
        if shared.size:
            pool_evidence[col] = {str(k): int(v) for k, v in shared.items()}
    if not pool_cols:
        pooled: bool | str = POOLING_UNRESOLVED
    else:
        pooled = bool(pool_evidence)
    if pooled is True:
        flags.append("pooled:donors share a pool/library id (D3)")
    elif pooled == POOLING_UNRESOLVED:
        flags.append("pooled:unresolved — no pool/library id in obs (D3 lower bound)")

    return ConfoundReport(
        cramers_v=cramers,
        separates=separates,
        flags=flags,
        exclude=bool(exclusion_reasons),
        exclusion_reasons=exclusion_reasons,
        covariate_candidates=covariate_candidates,
        pooled=pooled,
        pool_columns=list(pool_cols),
        pool_evidence=pool_evidence,
        pending=pending,
    )


def envelope_max_sigma_supported(donors_per_group: int) -> float | None:
    """Largest ``sigma_donor`` in the declared envelope this donor count already satisfies.

    Reads as: *"inside the operating envelope only if this stratum's (unmeasured) ``sigma_donor`` is
    at most X"*. ``None`` means even the envelope's easiest row (σ = 0.2) demands more donors than
    the stratum has, i.e. it is outside the envelope at every σ the envelope states.

    This is **not** an admission test — ``sigma_donor`` itself is unanchored (Amendment 3's own
    closing item) and "no stratum may be admitted to the real sweep on the strength of this entry
    alone". It is the envelope's necessary condition, computed from a number obs does know.
    """
    supported = [
        float(row["sigma_donor"]) for row in gc.OPERATING_ENVELOPE
        if donors_per_group >= int(row["min_donors_per_group"])
    ]
    return max(supported) if supported else None


def package_versions() -> dict[str, str | None]:
    """Exact installed versions of the packages §1 requires in the manifest.

    ``None`` means *not installed in the environment that produced this manifest* — which is the
    honest reading for ``cellxgene-census`` when a manifest is built from a cached obs frame.
    Queried through the metadata, never by importing, so building a manifest costs no imports.
    """
    versions: dict[str, str | None] = {}
    for name in MANIFEST_PACKAGES:
        try:
            versions[name] = _package_version(name)
        except PackageNotFoundError:
            versions[name] = None
    try:
        versions["pbcheck"] = _package_version("pbcheck")
    except PackageNotFoundError:  # pragma: no cover - editable installs always resolve
        versions["pbcheck"] = None
    return versions


def _manifest_row(contrast: Contrast, gate: GateResult, confound: ConfoundReport | None) -> dict:
    n_a, n_b = gate.n_donors.get("A", 0), gate.n_donors.get("B", 0)
    total_donors = n_a + n_b
    term_ids = (
        sorted(set(gate.obs["cell_type_ontology_term_id"].astype(str)))
        if "cell_type_ontology_term_id" in gate.obs.columns and gate.obs.shape[0]
        else []
    )
    if not gate.passed:
        status = "excluded_inclusion_gate"
    elif confound is not None and confound.exclude:
        status = "excluded_confound"
    else:
        status = "candidate"

    blockers = ["integer_check", "frozen_universe_size", "sigma_donor_estimate",
                "envelope_membership"]
    if status != "candidate":
        blockers.insert(0, status)

    return {
        "dataset_id": contrast.dataset_id,
        "cell_type": contrast.cell_type,
        "cell_type_ontology_term_id": term_ids[0] if len(term_ids) == 1 else (term_ids or None),
        "cell_type_ontology_depth": PENDING,
        "disease": contrast.disease,
        "reference": contrast.reference,
        "n_donors_A": n_a,
        "n_donors_B": n_b,
        "n_cells": int(gate.n_cells.get("A", 0) + gate.n_cells.get("B", 0)),
        "cells_per_donor_by_group": gate.cells_per_donor,
        "median_counts_per_cell_by_group": gate.median_counts_per_cell,
        "permutation_count": math.comb(total_donors, n_a) if total_donors else 0,
        "confound_flags": confound.flags if confound else ["not_screened:failed_inclusion_gate"],
        "confound_cramers_v": confound.cramers_v if confound else {},
        "pooled_flag": confound.pooled if confound else PENDING,
        "residual_df": total_donors - 2,
        "integer_check": PENDING,
        "frozen_universe_size": PENDING,
        "fitType": PENDING,
        "sigma_donor_estimate": PENDING,
        "envelope_min_donors_per_group": PENDING,
        "envelope_max_sigma_supported": envelope_max_sigma_supported(min(n_a, n_b)),
        "envelope_membership": PENDING,
        "gate_status": status,
        "gate_failures": list(gate.failures),
        "dropped_thin_donors": dict(gate.dropped_thin_donors),
        "dropped_missing_donor_cells": int(gate.dropped_missing_donor_cells),
        "confound_exclusion_reasons": confound.exclusion_reasons if confound else [],
        "covariate_candidates": confound.covariate_candidates if confound else [],
        "admitted_to_sweep": False,
        "admission_blockers": blockers,
        "census_version": CENSUS_VERSION,
    }


def screen_strata(
    obs: pd.DataFrame,
    *,
    dataset_col: str = "dataset_id",
    celltype_col: str = "cell_type",
    disease_col: str = "disease",
    donor_col: str = "donor_id",
    reference: str = REFERENCE_LEVEL,
    min_cells_per_donor: int = gc.MIN_CELLS,
    min_donors_per_group: int = MIN_DONORS_PER_GROUP,
) -> list[dict]:
    """Every candidate contrast in ``obs``, gated and screened — one manifest row each.

    Strata that fail the inclusion gate are kept in the output with their failure reasons rather
    than dropped: D4 asks for "the count and characteristics of strata excluded", and a list of
    survivors cannot answer that. The confound pre-screen is skipped for them (a Cramér's V over a
    degenerate design is a number with no meaning), which the row records as ``not_screened``.
    """
    rows: list[dict] = []
    for contrast in build_contrasts(
        obs, dataset_col=dataset_col, celltype_col=celltype_col,
        disease_col=disease_col, reference=reference,
    ):
        gate = apply_inclusion_gate(
            obs, contrast, dataset_col=dataset_col, celltype_col=celltype_col,
            disease_col=disease_col, donor_col=donor_col,
            min_cells_per_donor=min_cells_per_donor, min_donors_per_group=min_donors_per_group,
        )
        gate.obs.attrs.update(obs.attrs)
        confound = confound_prescreen(gate.obs, donor_col=donor_col) if gate.passed else None
        rows.append(_manifest_row(contrast, gate, confound))
    return rows


def _counts(rows) -> dict:
    by_status: dict[str, int] = {}
    gate_reasons: dict[str, int] = {}
    confound_reasons: dict[str, int] = {}
    for row in rows:
        by_status[row["gate_status"]] = by_status.get(row["gate_status"], 0) + 1
        for reason in row["gate_failures"]:
            gate_reasons[reason] = gate_reasons.get(reason, 0) + 1
        for reason in row["confound_exclusion_reasons"]:
            confound_reasons[reason] = confound_reasons.get(reason, 0) + 1
    n = len(rows)
    n_gate = by_status.get("excluded_inclusion_gate", 0)
    n_confound = by_status.get("excluded_confound", 0)
    return {
        "n_contrasts": n,
        "by_gate_status": by_status,
        "excluded_fraction": round((n_gate + n_confound) / n, 4) if n else 0.0,
        "excluded_confound_fraction": round(n_confound / n, 4) if n else 0.0,
        "inclusion_gate_failure_reasons": gate_reasons,
        "confound_exclusion_reasons": confound_reasons,
        "n_pooled": sum(1 for r in rows if r["pooled_flag"] is True),
        "n_pooling_unresolved": sum(1 for r in rows if r["pooled_flag"] == POOLING_UNRESOLVED),
    }


def _csv_value(value):
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=isinstance(value, dict))
    if value is None:
        return ""
    return value


def emit_manifest(
    rows,
    *,
    out_dir=None,
    stem: str = "census_strata",
    obs: pd.DataFrame | None = None,
    notes=None,
) -> dict:
    """The §1 manifest: one row per stratum-contrast, provenance in the header. JSON + CSV.

    The header carries what §1 demands be persisted — the exact census version and the exact package
    versions — plus the pre-registered config block (:func:`pbcheck.gate_config.manifest`), the
    declared operating envelope, the pending-field map and D4's excluded fraction with its reasons.
    The grid summaries under ``pilot/testsel/`` are a bare row list because their provenance lives in
    per-cell JSONs; here there is no second file to carry it, so the JSON is
    ``{"header": ..., "rows": [...]}`` and the CSV holds the rows alone (with dict/list cells
    JSON-encoded). ``census_version`` is repeated per row so the CSV is self-describing.

    ``out_dir=None`` writes nothing and just returns the manifest.

    **This artifact is a candidate list.** Committing it is not the §1 pre-registration, and no row
    of it is an admission to the sweep — see ``admitted_to_sweep`` (always ``False``) and the module
    docstring.
    """
    attrs = dict(obs.attrs) if obs is not None else {}
    header = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "census_version": CENSUS_VERSION,
        "census_version_source":
            "spec §1 pins open_soma(census_version='2025-01-30'); the mutable aliases "
            f"{list(MUTABLE_CENSUS_ALIASES)} are rejected by open_census()",
        "value_filter": attrs.get("value_filter", VALUE_FILTER),
        "organism": attrs.get("organism", CENSUS_ORGANISM),
        "obs_columns_spec": list(OBS_COLUMNS),
        "obs_columns_requested": attrs.get("requested_columns"),
        "obs_columns_missing": attrs.get("missing_columns"),
        "pool_columns_detected": attrs.get("pool_columns"),
        "depth_column": attrs.get("depth_column"),
        "package_versions": package_versions(),
        "config": gc.manifest(),
        "operating_envelope": [dict(row) for row in gc.OPERATING_ENVELOPE],
        "operating_envelope_source":
            "docs/AMENDMENTS.md Amendment 3 Change 1. SYNTHETIC: sigma_donor is an unanchored "
            "simulator knob, so envelope membership is PENDING for every row here and no row may "
            "be admitted to the sweep on the strength of the envelope entry alone.",
        "pending_fields": dict(PENDING_FIELDS),
        "counts": _counts(rows),
        "scope":
            "CANDIDATES ONLY. This is not the §1 pre-registration of the stratum list, and no row "
            "is an admission: inclusion-gate items 4 and 5 are pending (io_counts.py / "
            "gene_universe), sigma_donor and envelope membership are pending (Amendment 3), and "
            "the '8-12 datasets spanning the outcome space' curation is a separate judgement.",
        "thresholds_applied": {
            "min_cells_per_donor": gc.MIN_CELLS,
            "min_donors_per_group": MIN_DONORS_PER_GROUP,
            "perfect_confound_v": PERFECT_CONFOUND_V,
            "partial_confound_v_tag": PARTIAL_CONFOUND_V,
            "excluding_covariates": list(EXCLUDING_COVARIATES) + ["<any pool/library id>"],
        },
    }
    if notes:
        header["notes"] = notes
    manifest = {"header": header, "rows": [dict(row) for row in rows]}

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / f"{stem}.json", "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=1, ensure_ascii=False)
        with open(out_dir / f"{stem}.csv", "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(MANIFEST_FIELDS))
            writer.writeheader()
            for row in manifest["rows"]:
                writer.writerow({k: _csv_value(row.get(k)) for k in MANIFEST_FIELDS})
    return manifest
