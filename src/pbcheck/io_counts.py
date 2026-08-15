"""Counts I/O for one real-data stratum — the ``X``-reading second module of the Phase 0 harness
(spec §9 item 2), the counterpart to :mod:`pbcheck.census_select`.

What it does: takes the cells × genes count matrix of ONE ``(dataset_id × cell_type)`` stratum and
settles the two §1 inclusion-gate items that ``obs`` cannot decide — item 4 (raw integer counts
confirmed at load) and, through :mod:`pbcheck.gene_universe`, item 5 (frozen-universe size ≥ C5) —
then fills the manifest columns ``census_select`` was obliged to write as ``pending`` because it
reads ``obs`` and nothing else.

**The core reads an AnnData and never the network.** Every check below works on a matrix obtained
however the caller likes — a Census fetch, a cached ``.h5ad``, a simulator. The Census access is one
thin wrapper (:func:`load_stratum`) with a lazy import and the pinned
:data:`~pbcheck.census_select.CENSUS_VERSION` imported from ``census_select`` rather than restated,
so there is exactly one place in this package where the Census version lives. That split is
deliberate: the integrality gate is the part that must be testable, and a gate whose tests need a
network is a gate that stops being run.

**A non-integer stratum is DROPPED, never repaired.** §1 item 4 says "non-integer strata dropped",
§10 risk 6 says "Census raw not integer for some datasets … assert integrality at load, drop
offenders". Rounding, rescaling or ``expm1``-ing a matrix back toward integers would manufacture
counts that were never measured, and the resulting stratum would enter the sweep wearing a
provenance it does not have. So the failure path here writes a reason into the manifest and stops.
Note what is *not* a failure: **integral values stored in a float dtype**. Census ``raw`` is
routinely ``float32`` holding 0.0, 1.0, 7.0; those are integer counts and pass. What fails is a
genuinely fractional value — the fingerprint of ``normalize_total`` / ``log1p`` — together with
non-finite and negative values, neither of which any count matrix may contain.

**The obs snapshot and the load are two separate acts, and this module keeps them separate.**
``census_select`` counted donors and cells from an ``obs`` query; the load happens later, possibly
against a different Census build, possibly after a different filter. Where the two disagree —
fewer cells, a donor absent, per-cell depth that does not match the ``raw_sum`` the snapshot used —
:func:`reconcile_with_obs_snapshot` records a **discrepancy flag**. It never quietly adopts the
newer number, because a manifest row that silently re-derives itself cannot be audited against the
artifact that was committed.

**What this module does NOT settle** (each a real hole, named rather than papered over):

* **``fitType``** stays ``pending``. It belongs to the pseudobulk arm: Amendment 2 Change 1 retired
  DESeq2 and Change 4 replaced ``fitType`` with the moderated arm's realised prior (``d0``, the
  shrinkage factor, the residual df), which exists only once the arm has been fitted. This module
  aggregates donor profiles to *size* the universe and reads nothing off the fit.
* **``cell_type_ontology_depth``** stays ``pending`` (D5). Depth is a property of the CL graph, not
  of ``X``; it is ``controls.annotation_depth()``'s (§9 item 9).
* **Admission to the sweep.** Filling items 4 and 5 removes two blockers and admits nothing:
  Amendment 3 Change 1 additionally requires a per-stratum ``sigma_donor`` estimate and membership
  in the declared operating envelope, and its own closing section records that the ``sigma_donor``
  anchor "remains OPEN". ``admitted_to_sweep`` is therefore ``False`` on every row this module
  returns, exactly as it was on every row ``census_select`` produced.
* **Cell-count stratification of the permutation set (A2).** The per-cell depth vector computed here
  feeds §1's confound pre-screen and is what §7.2's depth-matching will downsample against, but A2
  itself is deferred to Phase 1 (Amendment 2 Change 6) and nothing here changes that.
* **§3's "≥ 3 pseudosamples per group post-aggregation" decision.** The number is measured and
  reported (a donor can clear §1 item 2's ≥ 10 cells and still fail the thin-donor filter's
  ``min_counts``), and the row is marked, but the arm's own SKIP is the arm's to take.
* **Pooling.** ``X`` carries no library/pool identifier, so D3's ``pooled_flag`` — including its
  ``unresolved`` third state — is left exactly as ``census_select`` found it in ``obs``.

Thresholds are imported, never restated: :data:`pbcheck.gate_config.MIN_CELLS` /
:data:`~pbcheck.gate_config.MIN_COUNTS` (the pre-registered thin-donor pair, §1 item 2 / §3,
implemented as "dropped, not merged" by Amendment 2 Change 7 and applied here through the arm's own
:func:`pbcheck.methods.pseudobulk.build_pseudobulk`) and
:data:`pbcheck.gate_config.MIN_UNIVERSE_SIZE` (C5). A second copy of any of them here is how the
gate and the arm would come to disagree about what was filtered.
"""

from __future__ import annotations

import copy
import csv
import inspect
import json
from contextlib import ExitStack
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from types import MappingProxyType

import numpy as np
from scipy import sparse as sp

from pbcheck import census_select as cs
from pbcheck import gate_config as gc
from pbcheck.census_select import CENSUS_VERSION, CensusNotInstalled
from pbcheck.gene_universe import UniverseTooSmall, frozen_universe, require_min_size

# ---------------------------------------------------------------------------
# Frozen by spec §1/§3/C5. Everything numeric here is imported from gate_config.
# ---------------------------------------------------------------------------

#: §1: "X is loaded from the Census **raw** layer". Any other layer is already normalised, which is
#: precisely what inclusion-gate item 4 exists to refuse — so the layer name is pinned, not defaulted
#: by the client library.
RAW_LAYER = "raw"

#: The Census measurement holding gene expression.
MEASUREMENT_NAME = "RNA"

#: Value written into the manifest's ``integer_check`` column when item 4 is satisfied. A failure is
#: written as ``"failed: <reason>"`` so that the manifest carries *why*, not just *that*.
INTEGER_CHECK_PASSED = "passed"
INTEGER_CHECK_FAILED_PREFIX = "failed: "

#: Distinct from :data:`pbcheck.census_select.PENDING`, and the distinction is load-bearing:
#: ``pending`` means "owed by a module that has not run yet", ``not_computed`` means "deliberately
#: not measured, because the stratum had already failed an earlier gate and the number would have no
#: meaning". Writing ``pending`` for the second case would advertise work that nobody owes.
NOT_COMPUTED = "not_computed"

#: ``gate_status`` values this module can set, each naming the clause it comes from.
STATUS_DROPPED_NON_INTEGER = "dropped_non_integer_counts"        # §1 item 4 ("dropped")
STATUS_SKIPPED_UNIVERSE_TOO_SMALL = "skipped_universe_too_small"  # §1 item 5 / C5 ("SKIP")
STATUS_SKIPPED_THIN_PSEUDOBULK = "skipped_thin_pseudobulk"        # §3 ("else SKIP; never fudge")

#: Precedence when more than one trips: the earliest gate in the spec's own order wins the status,
#: and every reason is recorded in ``gate_failures`` regardless.
IO_COUNTS_STATUSES = (
    STATUS_DROPPED_NON_INTEGER,
    STATUS_SKIPPED_UNIVERSE_TOO_SMALL,
    STATUS_SKIPPED_THIN_PSEUDOBULK,
)

#: §3: "Require ≥ 3 pseudosamples per group post-aggregation (else SKIP; never fudge)."
MIN_PROFILES_PER_GROUP = 3

#: §3's label-agnostic universe filter (A3), edgeR ``filterByExpr`` semantics grouped on nothing.
#: Same defaults as :func:`pbcheck.gene_universe.frozen_universe`; named here so the manifest can
#: record what was actually applied.
UNIVERSE_MIN_TOTAL_COUNT = 15
UNIVERSE_MIN_PROP = 0.5

#: Rows scanned at a time when ``X`` is dense. Bounds the temporary the integrality test allocates;
#: for sparse ``X`` no chunking is needed at all (see :func:`check_integer_counts`).
DENSE_CHUNK_ROWS = 8192

#: How many offending values a failure carries into the manifest. Enough to recognise a log1p matrix
#: on sight, few enough that the artifact stays readable.
MAX_EXAMPLES = 5

#: Relative tolerance for the cross-check between the X-derived per-cell depth and the ``raw_sum``
#: medians the obs snapshot recorded. Not exact equality: the Census stores ``raw_sum`` as a float
#: and a float32 row sum over ~60k genes carries accumulation error of that order. A disagreement
#: beyond it means the load is not the matrix the snapshot described, which is a discrepancy, not a
#: rounding artifact.
COUNTS_AGREEMENT_RTOL = 1e-3

#: The columns this module adds to a ``census_select`` manifest row. They are a superset, not a
#: replacement: :data:`MANIFEST_FIELDS` below is what a CSV writer must use, because
#: ``census_select.emit_manifest`` writes its own (shorter) field list and would drop these silently.
IO_COUNTS_FIELDS = (
    "integer_check_detail",
    "frozen_universe_detail",
    "load_vs_obs_snapshot",
    "counts_provenance",
)

#: The full row schema after this module has run: §1's manifest columns, then the counts columns.
MANIFEST_FIELDS = tuple(cs.MANIFEST_FIELDS) + IO_COUNTS_FIELDS

#: What each added column is for, mirroring :data:`pbcheck.census_select.PENDING_FIELDS`' contract
#: that an artifact should say what its own columns mean.
IO_COUNTS_FIELD_NOTES = MappingProxyType({
    "integer_check_detail":
        "the §1 item 4 / §10 risk 6 evidence: dtype, values scanned, and the counts of "
        "non-integer / negative / non-finite entries with examples. A drop is auditable only if "
        "the reason travels with it",
    "frozen_universe_detail":
        "the §1 item 5 / C5 evidence: measured universe size against MIN_UNIVERSE_SIZE, the "
        "surviving pseudobulk profiles per group (§3's >= 3 rule) and the thin-donor filter's own "
        "bookkeeping (dropped, not merged — Amendment 2 Change 7)",
    "load_vs_obs_snapshot":
        "donors and cells as the obs snapshot counted them vs as the load delivered them. A "
        "disagreement is FLAGGED, never silently adopted: the two are separate acts",
    "counts_provenance":
        "where X came from (census version + layer, or the caller), when, and which derived "
        "quantities were taken off it",
})


class NonIntegerCounts(ValueError):
    """``X`` is not a raw integer count matrix (§1 inclusion-gate item 4; §10 risk 6)."""


# ---------------------------------------------------------------------------
# Reading X without densifying it.
# ---------------------------------------------------------------------------


def _as_matrix(x):
    """The count matrix of ``x``, which may be an ``AnnData`` or the matrix itself."""
    if hasattr(x, "X") and hasattr(x, "obs"):
        return x.X
    return x


def _iter_value_blocks(X, *, chunk_rows: int = DENSE_CHUNK_ROWS):
    """Yield the values of ``X`` as 1-D arrays, without ever materializing it dense.

    For a sparse matrix that is exactly ``.data`` — the stored entries — and nothing else. The
    implicit zeros need no test: they are exactly 0, which is integral, non-negative and finite. On
    a real Census stratum this is the whole point, since scanning the dense form of a
    50k cells × 60k genes matrix costs 12 GB to learn something its ~2% stored entries already say.

    A dense array is walked in row blocks so the temporary the caller allocates per block is bounded
    by ``chunk_rows``, not by the matrix.

    Anything else — a backed/lazy ``X``, a dask array — raises rather than being coerced: a silent
    ``np.asarray`` on a lazy handle is exactly the densification this function exists to avoid.
    """
    if sp.issparse(X):
        data = getattr(X, "data", None)
        if not isinstance(data, np.ndarray) or data.dtype == object:  # lil/dok store lists
            data = X.tocsr().data
        yield np.asarray(data)
        return
    if isinstance(X, np.ndarray):
        arr = np.asarray(X)
        if arr.ndim <= 1:
            yield arr.ravel()
            return
        for start in range(0, arr.shape[0], max(int(chunk_rows), 1)):
            yield np.asarray(arr[start:start + chunk_rows]).ravel()
        return
    raise TypeError(
        f"unsupported X of type {type(X).__name__}: expected a numpy array or a scipy.sparse "
        "matrix. A backed or lazy X must be materialized by the caller (e.g. adata.to_memory()) — "
        "this module will not densify one implicitly."
    )


@dataclass(frozen=True)
class IntegerCheck:
    """The §1 item 4 verdict on one stratum's ``X``, with the evidence attached.

    ``passed`` is the gate. ``reason`` is what goes in the manifest next to a failure, and the
    counts below are what makes a drop auditable months later without re-fetching the matrix.
    """

    passed: bool
    reason: str | None
    dtype: str
    layer: str
    n_cells: int
    n_genes: int
    n_values_checked: int
    n_noninteger: int
    n_negative: int
    n_nonfinite: int
    max_abs_value: float | None
    examples: tuple[float, ...] = ()
    notes: tuple[str, ...] = ()

    @property
    def manifest_value(self) -> str:
        """The string written into the manifest's ``integer_check`` column."""
        if self.passed:
            return INTEGER_CHECK_PASSED
        return f"{INTEGER_CHECK_FAILED_PREFIX}{self.reason}"

    def as_dict(self) -> dict:
        return {
            "passed": bool(self.passed),
            "reason": self.reason,
            "dtype": self.dtype,
            "layer": self.layer,
            "n_cells": int(self.n_cells),
            "n_genes": int(self.n_genes),
            "n_values_checked": int(self.n_values_checked),
            "n_noninteger": int(self.n_noninteger),
            "n_negative": int(self.n_negative),
            "n_nonfinite": int(self.n_nonfinite),
            "max_abs_value": self.max_abs_value,
            "examples": list(self.examples),
            "notes": list(self.notes),
            "rule": "spec §1 inclusion gate item 4; §10 risk 6 — non-integer strata are DROPPED, "
                    "never rounded. Integral values in a float dtype PASS.",
        }


def check_integer_counts(
    x,
    *,
    layer: str = RAW_LAYER,
    max_examples: int = MAX_EXAMPLES,
    chunk_rows: int = DENSE_CHUNK_ROWS,
) -> IntegerCheck:
    """Confirm ``X`` is a raw integer count matrix (§1 item 4). Reports; never raises on a failure.

    Three things disqualify a matrix, and they are reported separately because they mean different
    things:

    1. **fractional values** — the fingerprint of ``normalize_total`` / ``log1p`` / a CPM matrix.
       This is §10 risk 6 ("Census raw not integer for some datasets") and the stratum is dropped.
    2. **non-finite values** (NaN/inf) — no count is NaN; something upstream already computed on it.
    3. **negative values** — a count cannot be negative; the matrix is a residual or a z-score.

    **A float dtype is not a failure.** Census ``raw`` is commonly ``float32`` carrying 0.0, 1.0,
    7.0; those are integer counts and pass. Integrality is tested on the values, not on the dtype —
    testing the dtype instead would drop most of the Census for a storage decision.

    On a sparse matrix only the **stored** entries are scanned (see :func:`_iter_value_blocks`); the
    implicit zeros are integral by construction. A matrix whose stored set is empty (all-zero X)
    therefore passes integrality and carries a note saying so — it is not this gate's job to call a
    stratum degenerate, and the universe check (item 5 / C5) will refuse it on its own terms.

    An ``AnnData`` may be passed directly; its ``.X`` is used.
    """
    X = _as_matrix(x)
    shape = tuple(int(v) for v in (getattr(X, "shape", ()) or ()))
    n_cells, n_genes = (shape + (0, 0))[:2]
    dtype = str(getattr(X, "dtype", "unknown"))

    n_checked = n_noninteger = n_negative = n_nonfinite = 0
    max_abs: float | None = None
    examples: list[float] = []

    for block in _iter_value_blocks(X, chunk_rows=chunk_rows):
        if block.size == 0:
            continue
        n_checked += int(block.size)
        if block.dtype.kind in "iub":  # an integer/bool dtype cannot hold a fractional value
            finite = block
        else:
            is_finite = np.isfinite(block)
            n_nonfinite += int((~is_finite).sum())
            finite = block[is_finite]
            fractional = finite != np.floor(finite)
            n_frac = int(fractional.sum())
            n_noninteger += n_frac
            if n_frac and len(examples) < max_examples:
                room = max_examples - len(examples)
                examples.extend(float(v) for v in finite[fractional][:room])
        n_negative += int((finite < 0).sum())
        if finite.size:
            block_max = float(np.max(np.abs(finite)))
            max_abs = block_max if max_abs is None else max(max_abs, block_max)

    notes: list[str] = []
    reasons: list[str] = []
    if n_cells == 0 or n_genes == 0:
        reasons.append(f"X is empty (shape {n_cells} x {n_genes}) — nothing was loaded to check")
    if n_nonfinite:
        reasons.append(f"non-finite values in X ({n_nonfinite} of {n_checked} stored)")
    if n_noninteger:
        shown = ", ".join(repr(round(v, 6)) for v in examples)
        reasons.append(
            f"non-integer values in X ({n_noninteger} of {n_checked} stored, e.g. {shown}) — this "
            "is a normalised and/or log-transformed matrix, not the Census raw layer"
        )
    if n_negative:
        reasons.append(f"negative values in X ({n_negative} of {n_checked} stored)")
    if n_checked == 0 and n_cells and n_genes:
        notes.append(
            "X stores no values at all (every entry is a structural zero); integrality holds "
            "trivially and the frozen-universe gate (item 5 / C5) is what will refuse this stratum"
        )

    return IntegerCheck(
        passed=not reasons,
        reason="; ".join(reasons) if reasons else None,
        dtype=dtype,
        layer=layer,
        n_cells=n_cells,
        n_genes=n_genes,
        n_values_checked=n_checked,
        n_noninteger=n_noninteger,
        n_negative=n_negative,
        n_nonfinite=n_nonfinite,
        max_abs_value=max_abs,
        examples=tuple(examples),
        notes=tuple(notes),
    )


def assert_integer_counts(x, **kwargs) -> IntegerCheck:
    """:func:`check_integer_counts`, raising :class:`NonIntegerCounts` on a failure (§9 item 2).

    The **stratum-level** response to a failed check is a DROP recorded in the manifest — see
    :func:`update_manifest_row` — not an exception: a sweep that aborts on the first non-integer
    dataset cannot report D4's excluded fraction. This raising form exists for callers that want to
    fail fast (a one-stratum script, a test) and for the assertion §9 item 2 names by that word.
    """
    check = check_integer_counts(x, **kwargs)
    if not check.passed:
        raise NonIntegerCounts(
            f"{check.reason}. Spec §1 inclusion-gate item 4 / §10 risk 6: such a stratum is "
            "DROPPED. Do not round, rescale or expm1 it back into shape — the counts would then "
            "be manufactured rather than measured."
        )
    return check


def counts_per_cell(x, *, chunk_rows: int = DENSE_CHUNK_ROWS) -> np.ndarray:
    """Per-cell total counts (row sums of ``X``) as a float array — the sequencing-depth vector.

    This is the ``raw_sum`` the Census exposes in ``obs`` when it exposes one, recomputed from the
    matrix that was actually loaded. It feeds §1's sequencing-depth bin in the confound pre-screen,
    §7.2's depth-matching, and ``median_counts_per_cell_by_group`` in the §1 manifest. Sparse input
    is summed sparsely.
    """
    X = _as_matrix(x)
    if sp.issparse(X):
        return np.asarray(X.sum(axis=1)).ravel().astype(float)
    if isinstance(X, np.ndarray):
        arr = np.asarray(X)
        if arr.ndim <= 1:
            return arr.astype(float).ravel()
        out = np.empty(arr.shape[0], dtype=float)
        step = max(int(chunk_rows), 1)
        for start in range(0, arr.shape[0], step):
            stop = start + step
            out[start:stop] = np.asarray(arr[start:stop]).sum(axis=1, dtype=float)
        return out
    raise TypeError(
        f"unsupported X of type {type(X).__name__}: expected a numpy array or a scipy.sparse matrix"
    )


def median_counts_per_cell_by_group(
    adata,
    *,
    condition_col: str = cs.CONDITION_COL,
    test_level: str | None = None,
    ref_level: str | None = None,
    counts=None,
) -> dict[str, float]:
    """§1's ``median_counts_per_cell_by_group``, computed from ``X`` rather than from ``obs``.

    Keyed ``{"A": ..., "B": ...}`` when ``test_level`` / ``ref_level`` are given, matching
    ``census_select``'s own convention (A = the disease term, B = ``normal``); keyed by the raw
    condition values otherwise. ``census_select`` fills this column from ``obs`` when the Census
    exposes ``raw_sum`` and leaves it ``pending`` when it does not — this is the ``pending`` case.
    """
    obs = adata.obs
    if condition_col not in obs.columns:
        raise ValueError(f"'{condition_col}' is not in .obs; call attach_obs() first")
    depth = counts_per_cell(adata) if counts is None else np.asarray(counts, dtype=float)
    if depth.shape[0] != obs.shape[0]:
        raise ValueError(
            f"counts vector has {depth.shape[0]} entries for {obs.shape[0]} cells"
        )
    tags = {}
    if test_level is not None:
        tags[str(test_level)] = "A"
    if ref_level is not None:
        tags[str(ref_level)] = "B"

    values = obs[condition_col].astype(str).to_numpy()
    out: dict[str, float] = {}
    for level in dict.fromkeys(values):
        mask = values == level
        if not mask.any():
            continue
        out[tags.get(level, level)] = float(np.median(depth[mask]))
    return out


# ---------------------------------------------------------------------------
# attach_obs (§9 item 2) — the derived columns the rest of the harness keys on.
# ---------------------------------------------------------------------------


def attach_obs(
    adata,
    contrast=None,
    *,
    donor_col: str = "donor_id",
    disease_col: str = "disease",
    condition_col: str = cs.CONDITION_COL,
    counts_col: str = cs.COUNTS_COLUMN,
    add_counts_column: bool = True,
    drop_missing_donor_cells: bool = True,
    inplace: bool = True,
):
    """Derive the columns the harness keys on: ``condition``, and the per-cell depth from ``X``.

    §9 item 2 lists ``attach_obs() (donor_id/condition/cell_type/assay/pool)``. ``donor_id``,
    ``cell_type``, ``assay`` and any pool identifier come from the Census obs and are left exactly as
    loaded; what has to be *derived* is the two-level ``condition`` column (from ``disease``, the
    same derivation ``census_select.apply_inclusion_gate`` makes) and — because ``obs`` may expose no
    count column at all — the per-cell ``raw_sum``, taken off ``X``.

    ``contrast`` (a :class:`pbcheck.census_select.Contrast`) is checked, not applied: if the frame
    carries a disease term outside ``{contrast.disease, contrast.reference}`` this raises, because a
    third level silently present in a two-group test is not something to fix by subsetting here —
    the loader's own filter was supposed to have done it, and quietly repairing it would hide that
    it did not.

    Cells whose ``donor_id`` is null have no replication unit (§1 item 3 opens with "donor_id
    present"), so ``drop_missing_donor_cells`` removes them — the same act
    ``apply_inclusion_gate`` performs on the obs side — and records the count in
    ``uns['pbcheck_obs']['n_missing_donor_cells']`` rather than leaving it to be discovered as a
    changed denominator downstream.

    ``inplace=True`` (the default) writes the derived columns into ``adata.obs`` and returns the same
    object: on a real stratum ``X`` is gigabytes and copying it to add two columns is not a neutral
    default. Pass ``inplace=False`` for a copy. Either way the returned object is a **new** one
    whenever cells had to be removed, so always use the return value.
    """
    if donor_col not in adata.obs.columns:
        raise ValueError(
            f"'{donor_col}' is not in .obs — a cell without a donor has no replication unit "
            "(§1 inclusion gate item 3) and can be neither thin-filtered nor permuted"
        )

    out = adata if inplace else adata.copy()
    n_missing = int(out.obs[donor_col].isna().sum())
    if n_missing and drop_missing_donor_cells:
        out = out[~out.obs[donor_col].isna().to_numpy()].copy()

    obs = out.obs
    if condition_col not in obs.columns:
        if disease_col not in obs.columns:
            raise ValueError(
                f"neither '{condition_col}' nor '{disease_col}' is in .obs; there is no contrast "
                "to derive"
            )
        obs[condition_col] = obs[disease_col].astype(str)
    else:
        obs[condition_col] = obs[condition_col].astype(str)

    if contrast is not None:
        allowed = {str(contrast.disease), str(contrast.reference)}
        seen = set(obs[condition_col].astype(str))
        extra = sorted(seen - allowed)
        if extra:
            raise ValueError(
                f"condition carries level(s) {extra} outside the contrast {sorted(allowed)}; a "
                "binary disease-vs-normal test (§1) must not be handed a third level. Fix the load "
                "filter — this function will not subset silently."
            )

    if add_counts_column and counts_col not in obs.columns:
        obs[counts_col] = counts_per_cell(out)

    uns = getattr(out, "uns", None)
    if uns is not None:
        uns["pbcheck_obs"] = {
            "condition_col": condition_col,
            "donor_col": donor_col,
            "counts_column": counts_col if add_counts_column else None,
            "counts_column_source": "row sums of the loaded X (not the Census obs raw_sum)",
            "n_missing_donor_cells": n_missing,
            "missing_donor_cells_dropped": bool(n_missing and drop_missing_donor_cells),
        }
    return out


# ---------------------------------------------------------------------------
# Inclusion-gate item 5 / C5 — the frozen universe, sized on real counts.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class UniverseCheck:
    """§1 item 5 / C5 for one stratum, plus what the thin-donor filter did on the way there."""

    size: int
    min_size: int
    skip_statuses: tuple[str, ...]
    reasons: tuple[str, ...]
    n_profiles: int
    profiles_by_group: dict[str, int]
    thin_donor_filter: dict
    genes: tuple[str, ...] = field(default=(), repr=False)

    @property
    def passed(self) -> bool:
        return not self.skip_statuses

    def as_dict(self) -> dict:
        return {
            "size": int(self.size),
            "min_size": int(self.min_size),
            "passed": bool(self.passed),
            "skip_statuses": list(self.skip_statuses),
            "reasons": list(self.reasons),
            "n_profiles": int(self.n_profiles),
            "profiles_by_group": dict(self.profiles_by_group),
            "min_profiles_per_group": MIN_PROFILES_PER_GROUP,
            "filter": {
                "min_total_count": UNIVERSE_MIN_TOTAL_COUNT,
                "min_prop": UNIVERSE_MIN_PROP,
                "label_agnostic": True,
            },
            "thin_donor_filter": dict(self.thin_donor_filter),
            "rule": "spec §1 inclusion gate item 5 / C5 (min universe size) and §3 (>= 3 "
                    "pseudosamples per group post-aggregation); thin donor profiles are dropped, "
                    "not merged (Amendment 2 Change 7)",
        }


def frozen_universe_check(
    adata,
    *,
    donor_col: str = "donor_id",
    celltype_col: str | None = "cell_type",
    condition_col: str = cs.CONDITION_COL,
    test_level: str | None = None,
    ref_level: str | None = None,
    min_size: int = gc.MIN_UNIVERSE_SIZE,
    min_cells: int = gc.MIN_CELLS,
    min_counts: int = gc.MIN_COUNTS,
    min_total_count: int = UNIVERSE_MIN_TOTAL_COUNT,
    min_prop: float = UNIVERSE_MIN_PROP,
    min_profiles_per_group: int = MIN_PROFILES_PER_GROUP,
) -> UniverseCheck:
    """Size the label-agnostic frozen universe for this stratum and apply C5 (§1 item 5).

    The universe is built the way the arm will build it, by calling the arm's own code rather than a
    look-alike: :func:`pbcheck.methods.pseudobulk.build_pseudobulk` (``mode='sum'`` per §3, then the
    pre-registered thin-donor filter — dropped, never merged, Amendment 2 Change 7) followed by
    :func:`pbcheck.gene_universe.frozen_universe` (pooled across donors, ignoring condition — A3).
    A second implementation here would be a second answer to "how big is this stratum's universe",
    and the manifest would record the wrong one.

    C5 is applied through :func:`pbcheck.gene_universe.require_min_size`, but the size is measured
    **first and kept either way**: a manifest row saying "skipped, universe too small" without the
    number is not auditable, and the number is what tells a reader whether the stratum missed by 5
    genes or by 195.

    §3's "≥ 3 pseudosamples per group post-aggregation" is measured here too and reported as its own
    skip status, because it can only be seen after the thin-donor filter has run on real counts: a
    donor that cleared §1 item 2's ≥ 10 cells can still fall below ``min_counts`` = 1000.

    ``decoupler`` and ``pydeseq2`` are imported lazily, inside this function, so that importing
    :mod:`pbcheck.io_counts` — and therefore the integrality gate, which needs neither — stays cheap.
    """
    from pbcheck.methods.pseudobulk import build_pseudobulk  # heavy: decoupler + pydeseq2

    groups_col = celltype_col if celltype_col and celltype_col in adata.obs.columns else None
    pdata = build_pseudobulk(
        adata,
        donor_col=donor_col,
        celltype_col=groups_col,
        condition_col=condition_col,
        min_cells=min_cells,
        min_counts=min_counts,
    )
    genes = frozen_universe(pdata, min_total_count=min_total_count, min_prop=min_prop)

    tags = {}
    if test_level is not None:
        tags[str(test_level)] = "A"
    if ref_level is not None:
        tags[str(ref_level)] = "B"
    profiles_by_group: dict[str, int] = {}
    if condition_col in pdata.obs.columns:
        for level, n in pdata.obs[condition_col].astype(str).value_counts().items():
            profiles_by_group[tags.get(str(level), str(level))] = int(n)

    skip_statuses: list[str] = []
    reasons: list[str] = []
    try:
        require_min_size(genes, min_size)
    except UniverseTooSmall as exc:
        skip_statuses.append(STATUS_SKIPPED_UNIVERSE_TOO_SMALL)
        reasons.append(f"{exc} (§1 inclusion gate item 5 / C5)")

    thin = [
        f"group {tag} has {n} pseudobulk profile(s) after the thin-donor filter "
        f"(< {min_profiles_per_group}; §3)"
        for tag, n in sorted(profiles_by_group.items())
        if n < min_profiles_per_group
    ]
    expected_tags = {"A", "B"} if tags else set()
    thin += [
        f"group {tag} has 0 pseudobulk profile(s) after the thin-donor filter "
        f"(< {min_profiles_per_group}; §3)"
        for tag in sorted(expected_tags - set(profiles_by_group))
    ]
    if thin:
        skip_statuses.append(STATUS_SKIPPED_THIN_PSEUDOBULK)
        reasons.extend(thin)

    return UniverseCheck(
        size=len(genes),
        min_size=int(min_size),
        skip_statuses=tuple(skip_statuses),
        reasons=tuple(reasons),
        n_profiles=int(pdata.n_obs),
        profiles_by_group=profiles_by_group,
        thin_donor_filter=dict(pdata.uns.get("thin_donor_filter", {})),
        genes=tuple(genes),
    )


# ---------------------------------------------------------------------------
# The §1 pre-screen's sequencing-depth bin, which needs X when obs carries no counts.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DepthScreen:
    """The sequencing-depth half of §1's confound pre-screen, computed from ``X``."""

    counts_column: str
    bin_column: str
    cramers_v: float | None
    separates: bool
    flags: tuple[str, ...]
    source: str


def depth_screen(
    adata,
    *,
    donor_col: str = "donor_id",
    condition_col: str = cs.CONDITION_COL,
    counts=None,
) -> DepthScreen:
    """Screen condition against the sequencing-depth bin, with the depth taken off ``X``.

    §1's pre-screen lists "sequencing-depth bin" among the five covariates it measures Cramér's V
    against. ``census_select`` can only supply it when the Census obs exposes ``raw_sum``/``nnz``,
    and marks it pending otherwise — this is that pending case, answered from the loaded matrix.

    The screen itself is :func:`pbcheck.census_select.confound_prescreen`, called with **no** other
    covariates and **no** pool columns so that it computes the depth bin and nothing else. Calling
    the real screen rather than re-deriving the association keeps one implementation of the
    donor-level Cramér's V, of perfect separation, and — most easily got wrong — of the rule that
    sets aside covariate levels held by a single donor. Its ``pooled`` verdict is discarded here:
    pooling is an ``obs`` question that ``census_select`` already answered, and ``X`` has nothing to
    add to it.

    Depth **never excludes a stratum** (:data:`pbcheck.census_select.EXCLUDING_COVARIATES` names
    assay and suspension_type only); §7.2 handles it by depth-matching and by the permutation null,
    which preserves each donor's depth profile and randomizes it across condition. A separation here
    is a flag to carry, not a verdict.
    """
    obs = adata.obs.copy()
    counts_col = cs.COUNTS_COLUMN
    obs[counts_col] = counts_per_cell(adata) if counts is None else np.asarray(counts, dtype=float)
    report = cs.confound_prescreen(
        obs,
        condition_col=condition_col,
        donor_col=donor_col,
        covariate_cols=(),
        pool_cols=[],
        depth_col=counts_col,
    )
    bin_col = cs.DEPTH_BIN_COL
    return DepthScreen(
        counts_column=counts_col,
        bin_column=bin_col,
        cramers_v=report.cramers_v.get(bin_col),
        separates=bool(report.separates.get(bin_col, False)),
        flags=tuple(f for f in report.flags if f.startswith(f"{bin_col}:")),
        source="row sums of the loaded X",
    )


# ---------------------------------------------------------------------------
# The obs snapshot vs the actual load — flagged, never reconciled by fiat.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LoadDiscrepancy:
    """How the loaded stratum differs from the ``obs`` snapshot ``census_select`` measured."""

    agrees: bool
    flags: tuple[str, ...]
    expected: dict
    loaded: dict
    donors_missing: dict
    donors_unexpected: dict
    donors_dropped_by_gate: dict
    thin_donors_no_longer_thin: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "agrees": bool(self.agrees),
            "flags": list(self.flags),
            "obs_snapshot": dict(self.expected),
            "loaded": dict(self.loaded),
            "donors_missing_from_load": {k: list(v) for k, v in self.donors_missing.items()},
            "donors_absent_from_snapshot": {k: list(v) for k, v in self.donors_unexpected.items()},
            "donors_dropped_by_the_obs_gate": {
                k: list(v) for k, v in self.donors_dropped_by_gate.items()
            },
            "donors_dropped_as_thin_but_not_thin_in_the_load": {
                k: dict(v) for k, v in self.thin_donors_no_longer_thin.items()
            },
            "rule": "the obs snapshot and the load are separate acts (spec §1 vs §9 item 2): a "
                    "disagreement is FLAGGED here and never silently adopted into the row",
        }


def _as_mapping(value):
    """A dict from ``value``, which may already be one or may be its JSON form (a CSV round-trip)."""
    if isinstance(value, dict):
        return value
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
        except (TypeError, ValueError):
            return None
        return parsed if isinstance(parsed, dict) else None
    return None


def _snapshot_from_row(row) -> dict:
    per_group = _as_mapping(row.get("cells_per_donor_by_group")) or {}
    cells: dict[str, int] = {}
    for tag, stats in per_group.items():
        stats = _as_mapping(stats) or {}
        if "n_cells" in stats:
            cells[str(tag)] = int(stats["n_cells"])
    return {
        "n_cells": int(row["n_cells"]) if row.get("n_cells") is not None else None,
        "n_donors": {"A": int(row.get("n_donors_A", 0)), "B": int(row.get("n_donors_B", 0))},
        "n_cells_by_group": cells,
        "median_counts_per_cell_by_group": _as_mapping(row.get("median_counts_per_cell_by_group")),
        "dropped_thin_donors": {
            str(k): int(v) for k, v in (_as_mapping(row.get("dropped_thin_donors")) or {}).items()
        },
    }


def _snapshot_from_gate(gate) -> dict:
    medians = gate.median_counts_per_cell
    return {
        "n_cells": int(sum(gate.n_cells.values())),
        "n_donors": {k: int(v) for k, v in gate.n_donors.items()},
        "n_cells_by_group": {k: int(v) for k, v in gate.n_cells.items()},
        "median_counts_per_cell_by_group": medians if isinstance(medians, dict) else None,
        "dropped_thin_donors": {str(k): int(v) for k, v in gate.dropped_thin_donors.items()},
        "donors_by_group": {
            tag: sorted(set(gate.obs.loc[gate.obs[cs.CONDITION_COL] == level, "donor_id"]
                            .astype(str)))
            for tag, level in (("A", gate.contrast.disease), ("B", gate.contrast.reference))
            if "donor_id" in gate.obs.columns and cs.CONDITION_COL in gate.obs.columns
        },
    }


def reconcile_with_obs_snapshot(
    snapshot,
    adata,
    *,
    donor_col: str = "donor_id",
    condition_col: str = cs.CONDITION_COL,
    test_level: str | None = None,
    ref_level: str | None = None,
    counts=None,
    min_cells: int = gc.MIN_CELLS,
    compare_magnitudes: bool = True,
) -> LoadDiscrepancy:
    """Compare what the load delivered against what the ``obs`` snapshot said, and FLAG the gap.

    ``snapshot`` is either a ``census_select`` manifest row (a dict, possibly straight out of the
    committed CSV) or a :class:`pbcheck.census_select.GateResult`, which additionally carries the
    donor ids and so supports an id-level comparison rather than only a count-level one.

    Nothing here corrects anything. If the load has fewer cells, or a donor the snapshot counted is
    absent, the row keeps the snapshot's numbers and gains a flag: the manifest is an artifact of a
    query made at a moment in time, and a row that silently re-derives itself against a later load
    cannot be checked against the commit that published it. Which of the two is *right* is a
    question for whoever reads the flag — it may be that the Census build moved, that the value
    filter differs, or that the snapshot was taken before a gate step.

    Two expected differences are recognised and not flagged as surprises:

    * **donors the obs gate dropped as thin** (§1 item 2) reappear in a raw load, because a load
      filters on dataset/cell type/disease and not on donor. The comparison is therefore made
      against the load **net of the donors the snapshot itself recorded as dropped** — applying the
      snapshot's own recorded decision, not inventing a new one — with the raw totals kept beside it
      under ``n_cells_including_gate_dropped_donors`` and the donors listed under
      ``donors_dropped_by_the_obs_gate``;
    * a per-cell depth that matches the snapshot's own ``raw_sum`` medians to within
      :data:`COUNTS_AGREEMENT_RTOL`.

    **That first allowance is bounded, and the bound is checked.** Carrying the snapshot's drop
    across is only legitimate while the load agrees the donor is thin. A donor recorded as dropped
    at 4 cells that arrives with 500 is not the same donor-decision reapplied — it is the load and
    the snapshot disagreeing about the very fact the gate was applied to, and netting it out would
    turn the largest possible divergence into the quietest possible output. Such donors are compared
    against ``min_cells`` (:data:`pbcheck.gate_config.MIN_CELLS`, the threshold §1 item 2 pins and
    the one aggregation will apply again), listed under
    ``donors_dropped_as_thin_but_not_thin_in_the_load``, **flagged, and ``agrees`` is False.** The
    downstream science is safe either way — the thin-donor filter re-runs at aggregation (Amendment
    2 Change 7) and would now keep the donor — but auditability is not, and auditability is what
    this function is for.

    ``compare_magnitudes=False`` reports donors and cells and **nothing derived from the values of**
    ``X``: no loaded medians, no depth cross-check. :func:`update_manifest_row` passes it when the
    integrality gate failed. Two reasons, both concrete: a median row sum of a matrix this study has
    just refused is a number with no referent, and the depth cross-check would re-report that same
    refusal under a different name — a ``×0.5`` rescale would surface as "obs snapshot 1809, X row
    sums 904.5", a "depth discrepancy" that is really the non-integrality already recorded in
    ``integer_check``. Counts of donors and cells stay meaningful whatever the values are, so they
    are still compared.
    """
    if hasattr(snapshot, "contrast"):  # a GateResult
        expected = _snapshot_from_gate(snapshot)
        test_level = test_level or snapshot.contrast.disease
        ref_level = ref_level or snapshot.contrast.reference
    else:
        expected = _snapshot_from_row(snapshot)
        test_level = test_level or snapshot.get("disease")
        ref_level = ref_level or snapshot.get("reference") or cs.REFERENCE_LEVEL

    obs = adata.obs
    if condition_col not in obs.columns:
        raise ValueError(f"'{condition_col}' is not in .obs; call attach_obs() first")
    tags = {str(test_level): "A", str(ref_level): "B"}
    levels = obs[condition_col].astype(str)
    donors = obs[donor_col].astype(str) if donor_col in obs.columns else None

    thin_snapshot = expected.get("dropped_thin_donors") or {}
    thin = set(thin_snapshot)
    depth = None
    if compare_magnitudes:
        depth = counts_per_cell(adata) if counts is None else np.asarray(counts, dtype=float)
    keep = np.ones(obs.shape[0], dtype=bool)
    dropped_by_gate: dict[str, list[str]] = {}
    no_longer_thin: dict[str, dict[str, int]] = {}
    cells_in_load: dict[str, int] = {}
    if donors is not None and thin:
        in_thin = donors.isin(thin)
        keep = ~in_thin.to_numpy()
        cells_in_load = {str(d): int(n) for d, n in donors[in_thin].value_counts().items()}

    loaded_donors: dict[str, list[str]] = {}
    loaded_cells: dict[str, int] = {}
    loaded_medians: dict[str, float] = {}
    for level, tag in tags.items():
        in_group = (levels == level).to_numpy()
        mask = in_group & keep
        loaded_cells[tag] = int(mask.sum())
        if depth is not None and mask.any():
            loaded_medians[tag] = float(np.median(depth[mask]))
        if donors is not None:
            loaded_donors[tag] = sorted(set(donors[mask]))
            gate_dropped = sorted(set(donors[in_group & ~keep]))
            if gate_dropped:
                dropped_by_gate[tag] = gate_dropped
            revived = {d: cells_in_load.get(d, 0) for d in gate_dropped
                       if cells_in_load.get(d, 0) >= min_cells}
            if revived:
                no_longer_thin[tag] = revived
    loaded = {
        "n_cells": int(keep.sum()),
        "n_donors": {tag: len(v) for tag, v in loaded_donors.items()},
        "n_cells_by_group": loaded_cells,
        "levels_present": sorted(set(levels)),
        "median_counts_per_cell_by_group": loaded_medians if compare_magnitudes else NOT_COMPUTED,
    }
    if not compare_magnitudes:
        loaded["magnitudes_note"] = (
            "X failed the integer-counts gate: nothing derived from its values is reported here "
            "and the depth cross-check was not run — it would re-report that same failure as a "
            "'depth discrepancy'. Donor and cell counts are unaffected and were compared."
        )
    if int(keep.sum()) != obs.shape[0]:
        loaded["n_cells_including_gate_dropped_donors"] = int(obs.shape[0])

    flags: list[str] = []
    missing: dict[str, list[str]] = {}
    unexpected: dict[str, list[str]] = {}

    if expected["n_cells"] is not None and expected["n_cells"] != loaded["n_cells"]:
        flags.append(
            f"n_cells:obs snapshot {expected['n_cells']}, load {loaded['n_cells']} "
            f"(delta {loaded['n_cells'] - expected['n_cells']:+d})"
        )
    for tag in ("A", "B"):
        want = expected["n_cells_by_group"].get(tag)
        got = loaded_cells.get(tag)
        if want is not None and got is not None and want != got:
            flags.append(
                f"n_cells_{tag}:obs snapshot {want}, load {got} (delta {got - want:+d})"
            )
        want_d = expected["n_donors"].get(tag)
        got_d = loaded["n_donors"].get(tag)
        if want_d is not None and got_d is not None and want_d != got_d:
            flags.append(f"n_donors_{tag}:obs snapshot {want_d}, load {got_d}")

    snapshot_donors = expected.get("donors_by_group")
    if snapshot_donors and donors is not None:
        for tag, want_ids in snapshot_donors.items():
            got_ids = set(loaded_donors.get(tag, []))
            gone = sorted(set(want_ids) - got_ids)
            extra = sorted(got_ids - set(want_ids))
            if gone:
                missing[tag] = gone
                flags.append(f"donors_{tag}:in the obs snapshot but not in the load {gone}")
            if extra:
                unexpected[tag] = extra
                flags.append(f"donors_{tag}:in the load but not in the obs snapshot {extra}")

    for tag, revived in sorted(no_longer_thin.items()):
        for donor_id, n_in_load in sorted(revived.items()):
            was = thin_snapshot.get(donor_id)
            had = f"{was} cell(s)" if was is not None else "too few cells"
            flags.append(
                f"donors_{tag}:'{donor_id}' was dropped as thin by the obs snapshot ({had}) but "
                f"carries {n_in_load} cell(s) in the load (>= {min_cells}) — the two acts disagree "
                "about the fact §1 item 2 was applied to, so the drop was not carried across as a "
                "silent agreement"
            )

    unknown = sorted(set(levels) - set(tags))
    if unknown:
        flags.append(
            f"condition:level(s) {unknown} outside the contrast "
            f"['{test_level}', '{ref_level}'] are present in the load"
        )

    want_med = expected.get("median_counts_per_cell_by_group")
    got_med = loaded.get("median_counts_per_cell_by_group")
    if compare_magnitudes and isinstance(want_med, dict) and isinstance(got_med, dict):
        for tag, want_v in want_med.items():
            got_v = got_med.get(tag)
            if got_v is None or want_v in (None, 0):
                continue
            if abs(float(got_v) - float(want_v)) > COUNTS_AGREEMENT_RTOL * abs(float(want_v)):
                flags.append(
                    f"median_counts_per_cell_{tag}:obs snapshot {float(want_v):.6g}, "
                    f"X row sums {float(got_v):.6g}"
                )

    return LoadDiscrepancy(
        agrees=not flags,
        flags=tuple(flags),
        expected=expected,
        loaded=loaded,
        donors_missing=missing,
        donors_unexpected=unexpected,
        donors_dropped_by_gate=dropped_by_gate,
        thin_donors_no_longer_thin=no_longer_thin,
    )


# ---------------------------------------------------------------------------
# Filling the manifest row.
# ---------------------------------------------------------------------------


def update_manifest_row(
    row,
    adata,
    *,
    donor_col: str = "donor_id",
    celltype_col: str | None = "cell_type",
    condition_col: str = cs.CONDITION_COL,
    min_size: int = gc.MIN_UNIVERSE_SIZE,
    min_cells: int = gc.MIN_CELLS,
    min_counts: int = gc.MIN_COUNTS,
    min_profiles_per_group: int = MIN_PROFILES_PER_GROUP,
    screen_depth: bool = True,
    source: str | None = None,
) -> dict:
    """One ``census_select`` manifest row + this stratum's counts → a **new** row. Never mutates.

    Fills what ``obs`` could not decide and this module can:

    * ``integer_check`` — ``"passed"`` or ``"failed: <reason>"`` (§1 item 4), with the evidence in
      ``integer_check_detail``;
    * ``frozen_universe_size`` — the measured size (§1 item 5 / C5), with the thin-donor filter's
      bookkeeping and the per-group profile counts in ``frozen_universe_detail``;
    * ``median_counts_per_cell_by_group`` — from ``X``, but **only when it was pending**. When
      ``census_select`` already filled it from the Census's own ``raw_sum``, the snapshot's value
      stays and the X-derived one is cross-checked against it (a disagreement becomes a flag);
    * the ``sequencing_depth_bin`` entry of ``confound_flags`` / ``confound_cramers_v`` — §1's
      pre-screen covariate that needs a count column;
    * ``gate_status`` / ``gate_failures`` / ``admission_blockers``, so a dropped or skipped stratum
      says which clause it died on and is still in the manifest for D4's excluded fraction.

    **Nothing derived from the magnitudes of a failed matrix is written, at any depth of the row.**
    If the integrality check fails, the universe, the depth bin and the counts medians are recorded
    as :data:`NOT_COMPUTED` rather than computed: the median row sum of a ``normalize_total`` matrix
    is a fact about the normalisation, not about sequencing depth, and putting it in a column named
    "counts" would be a fabricated number wearing the right label. The same applies **inside**
    ``load_vs_obs_snapshot``: the reconciliation still runs, because donor and cell counts stay
    meaningful whatever the values are, but it runs with ``compare_magnitudes=False`` so that no
    loaded median is reported and the depth cross-check cannot fire — otherwise a ``×0.5`` rescale
    would be reported twice, once truthfully as a failed integer check and once as a spurious
    "depth discrepancy" against the obs snapshot.

    **An exclusion already on the row wins.** A row that failed the obs-side inclusion gate or the
    confound pre-screen keeps that status; the counts findings are appended to ``gate_failures``
    rather than overwriting the reason it first failed.

    ``admitted_to_sweep`` stays ``False`` on every row: Amendment 3 Change 1 requires a per-stratum
    ``sigma_donor`` estimate and envelope membership on top of §1, and its own text records that the
    anchor "remains OPEN".
    """
    out = copy.deepcopy(dict(row))
    test_level = str(out.get("disease")) if out.get("disease") is not None else None
    ref_level = str(out.get("reference") or cs.REFERENCE_LEVEL)

    if condition_col not in adata.obs.columns:
        raise ValueError(
            f"'{condition_col}' is not in .obs — call attach_obs() first; this function fills a "
            "manifest row and deliberately derives no columns on the caller's AnnData"
        )

    # The integrality verdict comes FIRST and gates every magnitude below it, including the ones
    # nested inside the reconciliation report — a row that says "not_computed" at the top level
    # while carrying medians of the same refused matrix two dict levels down is not honest, it is
    # merely tidy.
    integer = check_integer_counts(adata)
    out["integer_check"] = integer.manifest_value
    out["integer_check_detail"] = integer.as_dict()
    counts = counts_per_cell(adata) if integer.passed else None

    discrepancy = reconcile_with_obs_snapshot(
        out, adata, donor_col=donor_col, condition_col=condition_col,
        test_level=test_level, ref_level=ref_level, counts=counts,
        min_cells=min_cells, compare_magnitudes=integer.passed,
    )
    out["load_vs_obs_snapshot"] = discrepancy.as_dict()

    new_failures: list[str] = []
    statuses: list[str] = []
    provenance = {
        "source": source or _default_source(adata),
        "layer": RAW_LAYER,
        "census_version": CENSUS_VERSION,
        "measured_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_cells_loaded": int(adata.n_obs),
        "n_genes_loaded": int(adata.n_vars),
        "counts_per_cell_source": "row sums of the loaded X",
        "thresholds": {
            "min_universe_size": int(min_size),
            "min_cells": int(min_cells),
            "min_counts": int(min_counts),
            "min_profiles_per_group": int(min_profiles_per_group),
        },
    }

    if not integer.passed:
        statuses.append(STATUS_DROPPED_NON_INTEGER)
        new_failures.append(
            f"raw integer counts not confirmed at load: {integer.reason} "
            "(§1 inclusion gate item 4 / §10 risk 6 — DROPPED, never rounded)"
        )
        out["frozen_universe_size"] = NOT_COMPUTED
        out["frozen_universe_detail"] = {
            "status": NOT_COMPUTED,
            "reason": "the stratum was dropped at the integer-counts gate; a universe frozen over "
                      "non-integer X would describe a matrix this study may not use",
        }
        if out.get("median_counts_per_cell_by_group") == cs.PENDING:
            out["median_counts_per_cell_by_group"] = NOT_COMPUTED
        provenance["derived_from_x"] = ["integer_check", "load_vs_obs_snapshot"]
    else:
        universe = frozen_universe_check(
            adata, donor_col=donor_col, celltype_col=celltype_col, condition_col=condition_col,
            test_level=test_level, ref_level=ref_level, min_size=min_size, min_cells=min_cells,
            min_counts=min_counts, min_profiles_per_group=min_profiles_per_group,
        )
        out["frozen_universe_size"] = int(universe.size)
        out["frozen_universe_detail"] = universe.as_dict()
        statuses.extend(universe.skip_statuses)
        new_failures.extend(universe.reasons)

        medians = median_counts_per_cell_by_group(
            adata, condition_col=condition_col, test_level=test_level, ref_level=ref_level,
            counts=counts,
        )
        derived = ["integer_check", "frozen_universe_size", "load_vs_obs_snapshot"]
        if out.get("median_counts_per_cell_by_group") == cs.PENDING:
            out["median_counts_per_cell_by_group"] = medians
            derived.append("median_counts_per_cell_by_group")
        else:
            provenance["median_counts_per_cell_by_group_from_X"] = medians
            provenance["median_counts_per_cell_by_group_kept"] = (
                "the obs snapshot's own raw_sum medians; the X-derived value is recorded beside it "
                "and any disagreement is flagged in load_vs_obs_snapshot"
            )
        if screen_depth:
            _merge_depth_screen(out, adata, donor_col=donor_col, condition_col=condition_col,
                                counts=counts, derived=derived)
        provenance["derived_from_x"] = derived

    out["counts_provenance"] = provenance
    if new_failures:
        out["gate_failures"] = list(out.get("gate_failures") or []) + new_failures

    # An exclusion the obs side already recorded happened first and keeps the status; the counts
    # findings are still carried, as reasons in gate_failures and as blockers below.
    tripped = [s for s in IO_COUNTS_STATUSES if s in statuses]
    if tripped and out.get("gate_status") in (None, "candidate"):
        out["gate_status"] = tripped[0]

    settled = {"integer_check"} if integer.passed else set()
    if isinstance(out.get("frozen_universe_size"), int):
        settled.add("frozen_universe_size")
    blockers = [
        b for b in list(out.get("admission_blockers") or [])
        if b not in settled and b not in IO_COUNTS_STATUSES
    ]
    for status in reversed(tripped):
        blockers.insert(0, status)
    out["admission_blockers"] = blockers
    out["admitted_to_sweep"] = False  # Amendment 3: sigma_donor and the envelope are still OPEN
    return out


def _default_source(adata) -> str:
    load = getattr(adata, "uns", {}).get("pbcheck_load") if hasattr(adata, "uns") else None
    if isinstance(load, dict) and load.get("source"):
        return str(load["source"])
    return "caller-supplied AnnData (not loaded by pbcheck.io_counts.load_stratum)"


def _merge_depth_screen(out, adata, *, donor_col, condition_col, counts, derived) -> None:
    """Add §1's depth-bin covariate to an existing confound pre-screen, without redoing the rest.

    If ``census_select`` already screened the depth bin (the Census exposed ``raw_sum`` in obs) the
    entry is left alone: re-deriving it from ``X`` and overwriting would replace a screened value
    with a differently-sourced one for no stated reason.
    """
    cramers = out.get("confound_cramers_v")
    if not isinstance(cramers, dict):
        cramers = _as_mapping(cramers) or {}
    if cs.DEPTH_BIN_COL in cramers:
        return
    screen = depth_screen(adata, donor_col=donor_col, condition_col=condition_col, counts=counts)
    cramers[cs.DEPTH_BIN_COL] = screen.cramers_v
    out["confound_cramers_v"] = cramers
    flags = list(out.get("confound_flags") or [])
    flags.extend(screen.flags)
    if screen.separates:
        flags.append(f"{cs.DEPTH_BIN_COL}:perfect_separation")
    out["confound_flags"] = flags
    derived.append(cs.DEPTH_BIN_COL)


# ---------------------------------------------------------------------------
# Emitting the extended manifest.
# ---------------------------------------------------------------------------


def emit_manifest(
    rows,
    *,
    out_dir=None,
    stem: str = "census_strata",
    obs=None,
    notes=None,
    fields=MANIFEST_FIELDS,
) -> dict:
    """:func:`pbcheck.census_select.emit_manifest`, widened to the columns this module adds.

    The header is built by ``census_select`` — one implementation of the provenance block, not two —
    and then gains a ``counts_gate`` section, because ``census_select``'s own ``counts["excluded_
    fraction"]`` is an obs-only accounting that knows nothing of a stratum dropped for non-integer
    counts. D4 asks for the excluded fraction *and its characteristics*; a reader who saw only the
    obs-side number would understate both.

    This exists because ``census_select.emit_manifest`` writes exactly its own ``MANIFEST_FIELDS``
    to the CSV, so the four columns added here would be present in the JSON and **silently absent**
    from the CSV — the kind of quiet gap this repository keeps finding in its own history.
    """
    manifest = cs.emit_manifest(rows, out_dir=None, obs=obs, notes=notes)
    manifest["header"]["counts_gate"] = _counts_gate(manifest["rows"])
    manifest["header"]["io_counts_fields"] = dict(IO_COUNTS_FIELD_NOTES)

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / f"{stem}.json", "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=1, ensure_ascii=False)
        with open(out_dir / f"{stem}.csv", "w", newline="", encoding="utf-8") as fh:
            # census_select's own cell encoder, imported rather than re-written: two ways of
            # writing a dict cell would give two CSV dialects for one manifest.
            writer = csv.DictWriter(fh, fieldnames=list(fields))
            writer.writeheader()
            for row in manifest["rows"]:
                writer.writerow({k: cs._csv_value(row.get(k)) for k in fields})
    return manifest


def _counts_gate(rows) -> dict:
    by_status: dict[str, int] = {}
    for row in rows:
        status = row.get("gate_status")
        if status in IO_COUNTS_STATUSES:
            by_status[status] = by_status.get(status, 0) + 1
    n = len(rows)
    n_dropped = sum(by_status.values())
    return {
        "n_rows": n,
        "by_counts_gate_status": by_status,
        "n_integer_check_pending": sum(
            1 for r in rows if r.get("integer_check") == cs.PENDING
        ),
        "n_integer_check_failed": sum(
            1 for r in rows
            if isinstance(r.get("integer_check"), str)
            and r["integer_check"].startswith(INTEGER_CHECK_FAILED_PREFIX)
        ),
        "n_load_discrepancies": sum(
            1 for r in rows
            if isinstance(r.get("load_vs_obs_snapshot"), dict)
            and not r["load_vs_obs_snapshot"].get("agrees", True)
        ),
        "counts_excluded_fraction": round(n_dropped / n, 4) if n else 0.0,
        "note":
            "counts-gate exclusions (§1 items 4 and 5, §3) are NOT included in "
            "header['counts']['excluded_fraction'], which census_select computes from obs alone.",
    }


# ---------------------------------------------------------------------------
# The Census fetch wrapper — the only part of this module that touches a network.
# ---------------------------------------------------------------------------


def _import_census():
    """Lazy-import ``cellxgene_census``, mirroring :func:`pbcheck.census_select.open_census`."""
    try:
        import cellxgene_census
    except ImportError as exc:  # pragma: no cover - exercised via a stubbed sys.modules in tests
        raise CensusNotInstalled(
            "cellxgene-census is not installed. It is an optional extra on purpose — every check "
            "in this module runs on a caller-supplied AnnData. Install with: "
            'pip install -e ".[census]"'
        ) from exc
    return cellxgene_census


def _filter_literal(value) -> str:
    """Quote ``value`` for a SOMA ``value_filter``, refusing anything that would need escaping.

    Cell-type labels carry commas, parentheses and hyphens, all of which quote cleanly; a quote or a
    backslash does not, and a mis-escaped filter does not fail — it silently returns the wrong
    cells. So this raises instead of guessing.
    """
    text = str(value)
    if "'" in text or '"' in text or "\\" in text:
        raise ValueError(
            f"cannot build a Census value_filter for {text!r}: it contains a quote or backslash. "
            "Query it by soma_joinid instead — a silently mis-escaped filter returns the wrong "
            "cells rather than an error."
        )
    return f"'{text}'"


def stratum_value_filter(
    dataset_id: str,
    cell_type: str,
    *,
    diseases=None,
    value_filter: str = cs.VALUE_FILTER,
) -> str:
    """§1's obs ``value_filter``, narrowed to one ``(dataset_id × cell_type)`` stratum.

    The §1 filter is kept verbatim as the first conjunct — it is what defines the queried population
    (primary data, Homo sapiens, disease known) — and the stratum keys are appended. The disease
    terms, when given, restrict the load to the two levels of one binary contrast (§1: one
    ``disease-vs-normal`` contrast per term, never a pooled "any disease" arm).
    """
    parts = [value_filter] if value_filter else []
    parts.append(f"dataset_id == {_filter_literal(dataset_id)}")
    parts.append(f"cell_type == {_filter_literal(cell_type)}")
    if diseases:
        terms = ", ".join(_filter_literal(d) for d in diseases)
        parts.append(f"disease in [{terms}]")
    return " and ".join(parts)


def load_stratum(
    dataset_id: str,
    cell_type: str,
    *,
    diseases=None,
    contrast=None,
    census=None,
    census_version: str = CENSUS_VERSION,
    organism: str = cs.CENSUS_ORGANISM,
    value_filter: str = cs.VALUE_FILTER,
    obs_columns=None,
    x_name: str = RAW_LAYER,
    measurement_name: str = MEASUREMENT_NAME,
    attach: bool = True,
) -> object:
    """Fetch one ``(dataset_id × cell_type)`` stratum from the Census, raw layer (§9 item 2).

    A thin wrapper and nothing more: it builds §1's value filter narrowed to the stratum, asks
    ``cellxgene_census.get_anndata`` for the **raw** ``X`` (§1: "X is loaded from the Census raw
    layer"), derives the ``condition`` and per-cell-depth columns via :func:`attach_obs`, and stamps
    the provenance into ``uns['pbcheck_load']``. Every check in this module then runs on the result
    — and equally on an AnnData from anywhere else, which is why none of them live in here.

    The Census version is :data:`pbcheck.census_select.CENSUS_VERSION`, imported rather than
    restated, and is validated by :func:`pbcheck.census_select.open_census`, which refuses the
    mutable aliases. Pass an already-open ``census`` to reuse a handle across strata; one opened
    here is closed here.

    **The integrality check is run and recorded, not enforced.** The result lands in
    ``uns['pbcheck_integer_check']`` and no exception is raised, because §1 item 4's response to a
    non-integer stratum is a DROP recorded in the manifest (:func:`update_manifest_row`), and a
    loader that raised would abort the sweep instead of reporting D4's excluded fraction.
    """
    if contrast is not None:
        if str(dataset_id) != str(contrast.dataset_id) or str(cell_type) != str(contrast.cell_type):
            raise ValueError(
                f"contrast {contrast.key} is not the ({dataset_id}, {cell_type}) stratum being "
                "loaded; one of the two is wrong and guessing which would load the wrong cells"
            )
        if diseases is None:
            diseases = (contrast.disease, contrast.reference)

    obs_value_filter = stratum_value_filter(
        dataset_id, cell_type, diseases=diseases, value_filter=value_filter
    )
    wanted = list(obs_columns) if obs_columns is not None else list(
        dict.fromkeys([*cs.OBS_COLUMNS, *cs.OBS_COLUMNS_OPTIONAL])
    )

    cellxgene_census = _import_census()
    with ExitStack() as stack:
        if census is None:
            census = cs.open_census(census_version=census_version)
            if hasattr(census, "__enter__"):
                stack.enter_context(census)
        kwargs = _get_anndata_kwargs(
            cellxgene_census.get_anndata,
            obs_value_filter=obs_value_filter,
            obs_columns=wanted,
            x_name=x_name,
            measurement_name=measurement_name,
        )
        adata = cellxgene_census.get_anndata(census, organism, **kwargs)

    if hasattr(adata, "uns"):
        adata.uns["pbcheck_load"] = {
            "source": f"cellxgene-census {census_version}, {measurement_name}/{x_name} layer",
            "census_version": census_version,
            "measurement_name": measurement_name,
            "x_name": x_name,
            "organism": organism,
            "obs_value_filter": obs_value_filter,
            "obs_columns_requested": wanted,
            "dataset_id": dataset_id,
            "cell_type": cell_type,
            "diseases": list(diseases) if diseases else None,
            "loaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
    if attach:
        adata = attach_obs(adata, contrast)
    if hasattr(adata, "uns"):
        adata.uns["pbcheck_integer_check"] = check_integer_counts(adata, layer=x_name).as_dict()
    return adata


def load_contrast(contrast, **kwargs):
    """:func:`load_stratum` for a :class:`pbcheck.census_select.Contrast` — both of its levels."""
    return load_stratum(
        contrast.dataset_id, contrast.cell_type, contrast=contrast, **kwargs
    )


def _get_anndata_kwargs(get_anndata, *, obs_value_filter, obs_columns, x_name, measurement_name):
    """Build ``get_anndata``'s keyword arguments for the installed client version.

    ``cellxgene-census`` renamed the obs column selector: newer releases take
    ``obs_column_names=[...]`` while older ones take ``column_names={"obs": [...]}``. Both are
    accepted by some versions and only one by others, and passing the wrong one is not a silent
    degradation — it either raises or returns every obs column of a multi-million-cell slice. The
    installed signature is therefore read rather than guessed, which is the same discipline
    ``docs/ENV_NOTES.md`` records for decoupler 2.x and PyDESeq2.
    """
    kwargs = {
        "measurement_name": measurement_name,
        "X_name": x_name,
        "obs_value_filter": obs_value_filter,
    }
    try:
        params = set(inspect.signature(get_anndata).parameters)
    except (TypeError, ValueError):  # pragma: no cover - builtins have no signature
        params = set()
    if "obs_column_names" in params:
        kwargs["obs_column_names"] = list(obs_columns)
    elif "column_names" in params:
        kwargs["column_names"] = {"obs": list(obs_columns)}
    return kwargs
