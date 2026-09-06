"""Audit orchestration: one stratum of one file, from ``.obs`` to a ``pbcheck-audit/1`` payload.

This is the only module of the product that calls the measurement engine. It prepares the stratum
from the caller's ``AnnData`` (never mutating it), runs the design audit, resolves which matrix
holds raw counts, freezes the gene universe, runs the two DE arms and the donor-permutation null,
and assembles the payload that :mod:`pbcheck.audit_schema` validates and :mod:`pbcheck.render`
renders.

What it is not: a Phase 0 measurement. The engine it drives was calibrated on synthetic oracles
inside the operating envelope of Amendment 3; this module runs it on a user's file, outside the
pre-registered protocol, and the payload says so in every report it feeds (the caveat block).

Every threshold this module owns is a product value, not a pre-registered one, and is defined at
the top of this file with that fact stated next to it. The pre-registered constants it reads
(``gate_config.ALPHA``, ``gate_config.LAMBDA_BAND``, ``gate_config.MIN_UNIVERSE_SIZE``,
``gate_config.MIN_CELLS``, ``gate_config.MIN_COUNTS``, ``gate_config.OPERATING_ENVELOPE``) are
copied into ``settings.protocol_constants`` and ``provenance`` so a reader can tell the two kinds
apart. The permutation counts are deliberately not taken from ``gate_config``: the gate script's
permutation budget answers a different question than a user's audit does.
"""

from __future__ import annotations

import platform
import re
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from math import comb
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse as sp

from pbcheck import (
    __version__,
    audit_schema,
    gate_config,
    gene_universe,
    io_counts,
    metrics,
    mtc,
    product_constants,
)
from pbcheck.design import audit_design
from pbcheck.methods.naive import naive_de
from pbcheck.methods.naive_engine import NaiveRelabelEngine
from pbcheck.methods.pseudobulk import build_pseudobulk, pseudobulk_de
from pbcheck.permutation import build_perms, run_null
from pbcheck.render.sections import STATUS_REASON_WORDS
from pbcheck.render.text import (
    LAMBDA_CLASS_WORDS,
    caveat_text,
    envelope_rows,
    readout_caveat_text,
    sentence_text,
)

# ---------------------------------------------------------------------------
# Product constants. None of these is pre-registered: they are this tool's own defaults and may be
# changed by an ordinary engineering change, unlike anything in pbcheck.gate_config.
# ---------------------------------------------------------------------------

#: Requested donor permutations for the naive arm. PRODUCT VALUE, not a protocol constant: it is
#: this module's own budget and never the gate script's pre-registered permutation count, which
#: answers a different question and is not read anywhere here (``pbcheck.gate_config`` is out of
#: bounds for any runtime decision of the product). Set by the runtime rule of the release plan
#: (``scripts/measure_audit_runtime.py``): the reference shape (``null_oracle(seed=1, n_genes=8000,
#: n_donors_per_group=8, n_cells_per_donor=625)``) at (1000, 200) measured 228.6 s wall, under the
#: 300 s cap, on a Windows 11 laptop (Intel Core, 20 logical CPUs) on 2026-09-05.
PRODUCT_N_PERM = 1000

#: Requested donor permutations for the pseudobulk arm; the paired series is
#: ``min(n_perm, n_perm_pb)`` permutations long. PRODUCT VALUE, not a protocol constant. Set by the
#: same 2026-09-05 measurement as :data:`PRODUCT_N_PERM`.
PRODUCT_N_PERM_PB = 200

#: Minimum donors per condition, before aggregation, below which the DE arms are not run at all.
#: PRODUCT VALUE: it is passed to :func:`pbcheck.design.audit_design` as its ``min_donors``, so the
#: design report's own verdict and this module's gate are judged against one number rather than
#: two that could drift apart. Numerically the spec's inclusion threshold, kept here because the
#: product decides for itself what it refuses to run on.
PRODUCT_MIN_DONORS_PER_GROUP = 3

#: Minimum surviving donor pseudobulk profiles per condition after the thin-donor filter, below
#: which the pseudobulk arm is dropped. PRODUCT VALUE: it is numerically the spec's
#: "at least 3 pseudosamples per group post-aggregation" rule, which the engine documents as the
#: caller's job and enforces nowhere on the product path, so this module enforces it itself.
PRODUCT_MIN_PROFILES_PER_GROUP = 3

#: Fraction of donors a gene must be detected in for the non-integer fallback universe.
#: PRODUCT VALUE, numerically equal to the pseudobulk universe's detection rule but computed on
#: cells of the matrix as found rather than on donor pseudobulk profiles.
FALLBACK_UNIVERSE_MIN_PROP = 0.5

#: Minimum size of the fallback universe before the run is refused. PRODUCT VALUE, numerically
#: equal to ``gate_config.MIN_UNIVERSE_SIZE`` but not that constant: the protocol's gate is about
#: the frozen pseudobulk universe, which does not exist on this path.
FALLBACK_UNIVERSE_MIN_SIZE = 200

#: Total counts a gene needs across donor pseudobulk profiles to enter the frozen universe.
#: PRODUCT VALUE, passed explicitly to :func:`pbcheck.gene_universe.frozen_universe` so the payload
#: never reports a filter parameter the call did not actually use. Bound to the engine's own name
#: rather than re-typed, so the product and the engine cannot drift apart unnoticed.
PRODUCT_UNIVERSE_MIN_TOTAL_COUNT = io_counts.UNIVERSE_MIN_TOTAL_COUNT

#: Fraction of donor pseudobulk profiles a gene must be detected in. PRODUCT VALUE, same reason.
PRODUCT_UNIVERSE_MIN_PROP = io_counts.UNIVERSE_MIN_PROP

#: Largest fraction of the file's cells that may be dropped for missing donor/condition/cell-type
#: values before the input is refused outright. PRODUCT VALUE: past this point the audited stratum
#: is not the file the user thinks they handed over, so a silent partial run would mislead.
MAX_MISSING_OBS_FRACTION = 0.5

#: Achieved naive permutation count below which the null is called coarse (note N9).
#: PRODUCT VALUE, a display threshold only: nothing is gated on it.
COARSE_NULL_THRESHOLD = 100

#: The clause the floor read-out line carries when the achieved permutation count is below
#: :data:`COARSE_NULL_THRESHOLD`. It is prose, and its home is ``pbcheck.render.text`` next to the
#: sentence it is interpolated into; it lives here until the renderer's remaining fixed sentences
#: are moved there, because that module cannot import this one (this one imports it).
COARSE_NULL_NOTE = ", coarse because few distinct donor splits exist"

#: Donors per group below which floors are called coarse and cross-file comparison is warned
#: against (note N4). Defined in :mod:`pbcheck.product_constants` with its origin, and re-exported
#: here under its established name: the report prose interpolates the same constant, and
#: ``pbcheck.render.text`` cannot import this module (this module imports it).
FEW_DONORS_THRESHOLD = product_constants.FEW_DONORS_THRESHOLD

#: The constant cell-type column added when the caller selects no cell type, so that the engine's
#: ``celltype_col`` argument is never ``None`` and the pseudobulk aggregation has one group.
STRATUM_COL = "_pbcheck_stratum"

#: Value of :data:`STRATUM_COL`, and the cell-type name of a pooled stratum.
STRATUM_VALUE = "all_cells"

#: Columns whose name looks like a cell-type annotation (note N8).
CELLTYPE_LIKE_PATTERN = re.compile(r"cell.?type|celltype|annotation", re.IGNORECASE)

#: The universe rule rendered verbatim in the report next to ``universe.builder`` on the ordinary
#: path. Built from the constants above so the prose and the call cannot drift apart.
PSEUDOBULK_UNIVERSE_RULE = (
    "A gene enters the frozen universe when its total count across the surviving donor pseudobulk "
    f"profiles is at least {PRODUCT_UNIVERSE_MIN_TOTAL_COUNT} and it is detected (count greater "
    f"than 0) in at least ceil({PRODUCT_UNIVERSE_MIN_PROP} * n_profiles) of them. The rule ignores "
    "the condition labels, so the tested gene set cannot shift when the labels are permuted "
    "(pbcheck.gene_universe.frozen_universe, applied after the thin-donor filter)."
)

#: The universe rule used when no matrix passes the raw-count check, rendered verbatim.
FALLBACK_UNIVERSE_RULE = (
    "The count matrix did not pass the raw-count check, so the pseudobulk universe does not exist "
    "for this run. A gene is kept when it has a value greater than 0 in at least one cell of at "
    f"least ceil({FALLBACK_UNIVERSE_MIN_PROP} * n_donors) donors, computed on the matrix as found. "
    "There is no total-count rule: a threshold on summed counts has no meaning on a normalised "
    "matrix. The rule ignores the condition labels, as the pseudobulk one does."
)

#: Packages whose installed version is recorded in ``provenance.packages``.
_PROVENANCE_PACKAGES = ("numpy", "scipy", "pandas", "anndata", "scanpy", "statsmodels",
                        "decoupler", "pydeseq2")

#: The three ``readout.lambda_*_class`` values, taken from the keys of the renderer's phrase table
#: so the payload and the prose cannot drift apart. The names describe a position against the
#: donor-pseudobulk arm's band and award no property to the arm they describe.
_LAMBDA_IN_BAND, _LAMBDA_ABOVE_BAND, _LAMBDA_BELOW_BAND = LAMBDA_CLASS_WORDS


class AuditInputError(ValueError):
    """The caller's file or settings cannot be audited as asked.

    Raised for a missing ``.obs`` column, a condition or cell-type level that is not in the data,
    a ``--counts-layer`` that the file does not have, duplicated ``var_names``, an empty stratum,
    or more than :data:`MAX_MISSING_OBS_FRACTION` of the cells missing a required ``.obs`` value.
    The message names the column and lists the values that are available, so the caller can fix
    the command without opening the file in a notebook.
    """


@dataclass(frozen=True)
class AuditSettings:
    """Everything the audit needs beyond the data itself.

    ``alpha`` defaults to the pre-registered ``gate_config.ALPHA``; the permutation counts default
    to this module's own product constants and are never read from ``gate_config``.
    """

    donor_col: str
    condition_col: str
    test_level: str
    ref_level: str
    celltype_col: str | None = None
    celltype_value: str | None = None
    batch_cols: tuple[str, ...] = ()
    counts_layer: str | None = None
    n_perm: int = PRODUCT_N_PERM
    n_perm_pb: int = PRODUCT_N_PERM_PB
    seed: int = 0
    alpha: float = gate_config.ALPHA
    design_only: bool = False
    top_n: int = 25


# ---------------------------------------------------------------------------
# Small helpers.
# ---------------------------------------------------------------------------


@contextmanager
def _stage(timings: dict, name: str):
    """Add one stage's seconds to ``timings[name]``; recorded even when the stage raises.

    Additive rather than assigning: a stage that runs in two pieces (the permutation null and the
    summarisation of its output, the two universe builders) reports the work it actually did, so
    the per-stage table sums to the run's wall clock instead of losing whichever piece came last.
    """
    started = time.perf_counter()
    try:
        yield
    finally:
        timings[name] = float(timings.get(name, 0.0) + (time.perf_counter() - started))


def _available(values) -> str:
    """The distinct values of an ``.obs`` column, for an error message."""
    seen = pd.Index(pd.Series(values).dropna().astype(str)).unique().tolist()
    shown = ", ".join(repr(v) for v in seen[:20])
    return shown + (", ..." if len(seen) > 20 else "")


def _lambda_class(value: float | None) -> str | None:
    """The report's word for a genomic-inflation factor, or ``None`` if the arm did not run.

    The three class names are the keys of ``pbcheck.render.text.LAMBDA_CLASS_WORDS``. The band is
    ``gate_config.LAMBDA_BAND``, the donor-pseudobulk arm's band, with the inclusive bounds of
    ``band()`` in ``scripts/synthetic_gate.py``; where the naive arm's class is rendered, the
    sentence says that the band is shown to describe the number, not as that arm's own criterion.
    """
    if value is None or not np.isfinite(value):
        return None
    low, high = gate_config.LAMBDA_BAND
    if value < low:
        return _LAMBDA_BELOW_BAND
    if value > high:
        return _LAMBDA_ABOVE_BAND
    return _LAMBDA_IN_BAND


def _floor_mc_se(counts: np.ndarray) -> float:
    """Monte-Carlo standard error of a permutation floor.

    The same estimator ``run_null`` uses for ``monte_carlo["naive_floor_mc_se"]`` (the standard
    deviation of the per-permutation counts over the square root of their number), reproduced here
    because the engine computes it for the paired series only and the headline floor of this report
    is the solo one. ``tests/test_audit.py`` pins the two against each other on the paired series.
    """
    x = np.asarray(counts, dtype=float)
    if x.size < 2:
        return float("nan")
    return float(np.std(x, ddof=1) / np.sqrt(x.size))


def _floor_block(series, n_genes: int) -> dict:
    """``metrics.perm_floor`` plus the Monte-Carlo SE of the same series."""
    block = metrics.perm_floor(series, n_genes)
    counts = series.counts if isinstance(series, metrics.NDegSeries) else series
    block["mc_se"] = _floor_mc_se(counts)
    return block


def _lambda_block(pval_matrix: np.ndarray, real_pvals: np.ndarray) -> dict:
    """Lambda over permutations plus the B5 empirical-permutation-p machinery check.

    B5's construction is a validity check on the permutation machinery, never a calibration
    measure: see ``metrics.empirical_perm_pvalues``, which says so at length and measured it.
    """
    lam = metrics.lambda_over_permutations(pval_matrix)
    b5 = metrics.empirical_perm_pvalues(real_pvals, pval_matrix)
    return {
        "lambda": float(lam["lambda"]),
        "lambda_iqr": float(lam["lambda_iqr"]),
        "n_perm": int(lam["n_perm"]),
        "b5_lambda_empirical": float(metrics.genomic_inflation(b5)),
    }


def _top_genes(table: pd.DataFrame, top_n: int, columns: tuple[str, ...]) -> list[dict]:
    """The ``top_n`` genes of a DE table by raw p-value, as plain JSON-able records."""
    tab = table[table["pval"].notna()].sort_values("pval", kind="mergesort").head(top_n)
    out = []
    for gene, row in tab.iterrows():
        record = {"gene": str(gene)}
        for column in columns:
            value = row[column] if column in tab.columns else float("nan")
            record[column] = float(value)
        out.append(record)
    return out


def _donor_condition_map(obs: pd.DataFrame, donor_col: str, condition_col: str) -> pd.Series:
    """Donor to condition, one row per donor.

    The same derivation ``permutation._true_map`` makes, reproduced from public operations rather
    than imported, so this module does not depend on a private engine name. ``run_null`` builds its
    donor list and true test set from exactly this series, and ``naive_null`` must build the same
    ones for its permutations to be the engine's permutations.
    """
    return (
        obs[[donor_col, condition_col]].astype(str).drop_duplicates()
        .set_index(donor_col)[condition_col]
    )


def _n_distinct_splits(n_donors: int, n_test: int) -> int:
    """Balanced donor splits available to ``build_perms``: all of them minus the ones it excludes.

    ``build_perms`` computes ``comb(n_donors, n_test) - (2 if complement != true_test_set else 1)``:
    it skips the true test set and its exact complement, and counts them once when they are the
    same set. They can only be the same set if both are empty, which a stratum that passed the
    donor gate never is (the two groups partition the donors and each has at least
    :data:`PRODUCT_MIN_DONORS_PER_GROUP` members), so the subtraction here is 2 unless one group is
    empty. Clamped at 0 for the degenerate shapes the gates above already refuse.
    """
    if n_donors <= 0 or n_test <= 0 or n_test >= n_donors:
        return 0
    return max(int(comb(n_donors, n_test)) - 2, 0)


def _package_versions() -> dict[str, str]:
    """Installed versions of the packages the numbers depend on."""
    versions = {}
    for name in _PROVENANCE_PACKAGES:
        try:
            versions[name] = importlib_metadata.version(name)
        except importlib_metadata.PackageNotFoundError:
            # Recorded rather than raised: pydeseq2 is only needed by the retired DESeq2 path, and
            # an artifact that says "not installed" is more useful than a crash at report time.
            versions[name] = "not installed"
    return versions


def _working_adata(matrix, obs: pd.DataFrame, var: pd.DataFrame):
    """The object every engine call receives: one matrix, the stratum's ``obs``, its source's ``var``."""
    return ad.AnnData(X=matrix, obs=obs.copy(), var=var.copy())


# ---------------------------------------------------------------------------
# Stratum preparation.
# ---------------------------------------------------------------------------


def prepare_stratum(adata, settings: AuditSettings) -> tuple[ad.AnnData, dict]:
    """Subset the caller's data to the audited stratum and describe what was dropped.

    Returns the working copy (the caller's object is never touched) and the payload's ``input``
    block with every drop counted. When no cell type is selected a constant :data:`STRATUM_COL`
    column is added, so the engine's ``celltype_col`` argument is never ``None``.
    """
    obs = adata.obs
    for column in (settings.donor_col, settings.condition_col):
        if column not in obs.columns:
            raise AuditInputError(
                f"obs has no column {column!r}; available columns: {', '.join(map(str, obs.columns))}"
            )
    if settings.celltype_col is not None and settings.celltype_col not in obs.columns:
        raise AuditInputError(
            f"obs has no cell-type column {settings.celltype_col!r}; available columns: "
            f"{', '.join(map(str, obs.columns))}"
        )
    missing_batch = [c for c in settings.batch_cols if c not in obs.columns]
    if missing_batch:
        raise AuditInputError(
            f"obs has no batch column(s) {', '.join(repr(c) for c in missing_batch)}; "
            f"available columns: {', '.join(map(str, obs.columns))}"
        )

    duplicated = pd.Index(adata.var_names)[pd.Index(adata.var_names).duplicated()].unique()
    if len(duplicated):
        shown = ", ".join(repr(str(g)) for g in duplicated[:5])
        raise AuditInputError(
            f"var_names has {len(duplicated)} duplicated name(s) (e.g. {shown}); the gene universe "
            "and both arms index by gene name, so call adata.var_names_make_unique() first"
        )

    n_cells_loaded = int(adata.n_obs)
    missing_donor = obs[settings.donor_col].isna().to_numpy()
    missing_condition = obs[settings.condition_col].isna().to_numpy()
    if settings.celltype_col is not None:
        missing_celltype = obs[settings.celltype_col].isna().to_numpy()
    else:
        missing_celltype = np.zeros(n_cells_loaded, dtype=bool)

    any_missing = missing_donor | missing_condition | missing_celltype
    if n_cells_loaded and any_missing.mean() > MAX_MISSING_OBS_FRACTION:
        raise AuditInputError(
            f"{int(any_missing.sum())} of {n_cells_loaded} cells are missing a value in "
            f"{settings.donor_col!r}, {settings.condition_col!r} or the cell-type column, which is "
            f"more than {MAX_MISSING_OBS_FRACTION:.0%} of the file; the audited stratum would not "
            "be the data you asked about"
        )

    condition = obs[settings.condition_col].astype(str).to_numpy()
    if settings.test_level == settings.ref_level:
        # Refused here rather than deep inside the DE arms: both levels exist, the donor gate reads
        # the one group twice and passes, and the run would only fail after paying for itself.
        raise AuditInputError(
            f"condition column {settings.condition_col!r} was given the same value "
            f"{settings.test_level!r} as both the test and the reference level, so there is one "
            f"group, not a contrast; available values: {_available(obs[settings.condition_col])}"
        )
    levels = {settings.test_level, settings.ref_level}
    for level in (settings.test_level, settings.ref_level):
        if not bool((condition[~missing_condition] == level).any()):
            raise AuditInputError(
                f"condition column {settings.condition_col!r} has no value {level!r}; "
                f"available values: {_available(obs[settings.condition_col])}"
            )
    in_levels = np.isin(condition, list(levels)) & ~missing_condition

    if settings.celltype_col is not None:
        if settings.celltype_value is None:
            raise AuditInputError(
                f"a cell-type column {settings.celltype_col!r} was given without a cell-type "
                f"value; available values: {_available(obs[settings.celltype_col])}"
            )
        celltype = obs[settings.celltype_col].astype(str).to_numpy()
        if not bool((celltype[~missing_celltype] == settings.celltype_value).any()):
            raise AuditInputError(
                f"cell-type column {settings.celltype_col!r} has no value "
                f"{settings.celltype_value!r}; available values: "
                f"{_available(obs[settings.celltype_col])}"
            )
        in_celltype = (celltype == settings.celltype_value) & ~missing_celltype
    else:
        if settings.celltype_value is not None:
            # Symmetric to the check above: a cell-type value with no column to read it from would
            # be recorded in the payload as a subset that was never taken.
            raise AuditInputError(
                f"a cell-type value {settings.celltype_value!r} was given without a cell-type "
                f"column; name the column with celltype_col, or drop the value to pool every cell "
                f"type; columns that look like cell-type columns in this file: "
                + (", ".join(str(c) for c in adata.obs.columns
                             if CELLTYPE_LIKE_PATTERN.search(str(c))) or "(none)")
            )
        in_celltype = np.ones(n_cells_loaded, dtype=bool)

    kept = ~any_missing & in_levels & in_celltype
    dropped_other_condition = int((~any_missing & ~in_levels).sum())
    dropped_other_celltype = int((~any_missing & in_levels & ~in_celltype).sum())

    if not kept.any():
        raise AuditInputError(
            f"no cells left after selecting {settings.condition_col!r} in "
            f"{sorted(levels)!r}"
            + (f" and {settings.celltype_col!r} == {settings.celltype_value!r}"
               if settings.celltype_col is not None else "")
        )

    work = adata[kept].copy()
    for column in work.obs.columns:
        if isinstance(work.obs[column].dtype, pd.CategoricalDtype):
            work.obs[column] = work.obs[column].cat.remove_unused_categories()
    if settings.celltype_col is None:
        work.obs[STRATUM_COL] = pd.Categorical([STRATUM_VALUE] * work.n_obs)

    celltype_like = [str(c) for c in adata.obs.columns if CELLTYPE_LIKE_PATTERN.search(str(c))]

    input_block = {
        "path": None,
        "n_cells_loaded": n_cells_loaded,
        "n_genes_loaded": int(adata.n_vars),
        "n_cells_audited": int(work.n_obs),
        "n_cells_dropped_other_condition": dropped_other_condition,
        "n_cells_dropped_other_celltype": dropped_other_celltype,
        "n_cells_dropped_missing_donor": int(missing_donor.sum()),
        "n_cells_dropped_missing_condition": int(missing_condition.sum()),
        "n_cells_dropped_missing_celltype": int(missing_celltype.sum()),
        "donor_col": settings.donor_col,
        "condition_col": settings.condition_col,
        "test_level": settings.test_level,
        "ref_level": settings.ref_level,
        "celltype_col": settings.celltype_col,
        "celltype_value": settings.celltype_value,
        "celltype_like_columns_found": celltype_like,
        "batch_cols": [str(c) for c in settings.batch_cols],
        "counts_source": None,
    }
    return work, input_block


# ---------------------------------------------------------------------------
# Counts resolution.
# ---------------------------------------------------------------------------


def _counts_candidates(work, settings: AuditSettings) -> list[tuple[str, object, pd.DataFrame]]:
    """(label, matrix, var) for every place raw counts may live, in the order they are tried."""
    if settings.counts_layer is not None:
        if settings.counts_layer not in work.layers:
            raise AuditInputError(
                f"no layer {settings.counts_layer!r} in this file; available layers: "
                + (", ".join(map(str, work.layers.keys())) or "(none)")
            )
        return [(f"layers:{settings.counts_layer}", work.layers[settings.counts_layer], work.var)]

    candidates: list[tuple[str, object, pd.DataFrame]] = [("X", work.X, work.var)]
    for name in ("counts", io_counts.RAW_LAYER):
        if name in work.layers and f"layers:{name}" not in [label for label, _, _ in candidates]:
            candidates.append((f"layers:{name}", work.layers[name], work.var))
    if work.raw is not None:
        raw_adata = work.raw.to_adata()
        candidates.append(("raw.X", raw_adata.X, raw_adata.var))
    return candidates


def resolve_counts(work, settings: AuditSettings):
    """Find the matrix that will be audited and build the working object from it.

    Returns ``(working AnnData or None, counts source label or None, IntegerCheck)``. The check is
    the first candidate's when nothing passes, which is the check of ``X`` unless the caller named
    a layer, in which case it is that layer's: those are the matrices the report then talks about.
    Nothing is ever rounded or rescaled to make a matrix pass.
    """
    candidates = _counts_candidates(work, settings)
    first_check = None
    for label, matrix, var in candidates:
        check = io_counts.check_integer_counts(matrix, layer=label)
        if first_check is None:
            first_check = check
        if check.passed:
            return _working_adata(matrix, work.obs, var), label, check
    return None, None, first_check


# ---------------------------------------------------------------------------
# The naive half of the permutation null.
# ---------------------------------------------------------------------------


def naive_null(work, universe: list[str], settings: AuditSettings) -> dict:
    """The donor-permutation null for the naive arm alone, from public engine names.

    ``permutation.run_null`` always runs both arms and has no flag to skip the pseudobulk one, so
    the two statuses that have no pseudobulk arm (a non-integer matrix, or too few profiles after
    the thin-donor filter) need this path. It is not a second implementation of the statistic: the
    donor list, the true test set, the permutation set and the per-permutation test are the same
    public calls ``run_null`` makes, in the same order, and
    ``test_naive_null_matches_run_null_bitwise`` pins the p-value matrix and the solo #DEG series
    against the engine's on data where both paths can run.

    The permutation set is built with ``n_perm=max(settings.n_perm, settings.n_perm_pb)`` and then
    truncated to ``settings.n_perm``, because that is what ``run_null`` does and the sampled set
    depends on the requested count. The pseudobulk arm does not run here, but the naive
    permutations stay the ones the engine would have drawn.
    """
    tmap = _donor_condition_map(work.obs, settings.donor_col, settings.condition_col)
    donors = list(tmap.index)
    true_test = set(tmap.index[tmap == settings.test_level])

    perms = build_perms(donors, true_test,
                        n_perm=max(settings.n_perm, settings.n_perm_pb), seed=settings.seed)
    perms = perms[:settings.n_perm]

    uni_index = pd.Index(universe, name="gene")
    engine = NaiveRelabelEngine.from_adata(work, donor_col=settings.donor_col, genes=universe)
    real = engine.test(true_test, condition_col=settings.condition_col,
                       test_level=settings.test_level, ref_level=settings.ref_level)

    naive_pvals = np.full((len(perms), len(universe)), np.nan)
    ndeg_solo = np.zeros(len(perms), dtype=np.int64)
    cells_per_donor = work.obs.groupby(settings.donor_col, observed=True).size()
    perm_test_cells = np.zeros(len(perms), dtype=float)
    for i, test_donors in enumerate(perms):
        result = engine.test(test_donors, condition_col=settings.condition_col,
                             test_level=settings.test_level, ref_level=settings.ref_level)
        naive_pvals[i] = result.table["pval"].reindex(uni_index).to_numpy()
        ndeg_solo[i] = mtc.bh_over_universe(
            result, universe, alpha=settings.alpha).n_significant(fdr=settings.alpha)
        perm_test_cells[i] = float(cells_per_donor[list(test_donors)].sum())

    real_test_cells = float(cells_per_donor[list(true_test)].sum())
    inside = bool(perm_test_cells.size
                  and perm_test_cells.min() <= real_test_cells <= perm_test_cells.max())
    percentile = (float((perm_test_cells <= real_test_cells).mean()) if perm_test_cells.size
                  else float("nan"))

    return {
        "naive_real": real,
        "naive_pvals": naive_pvals,
        "naive_pvals_real": real.table["pval"].reindex(uni_index).to_numpy(),
        "naive_ndeg_solo": metrics.NDegSeries(ndeg_solo, metrics.BH_SOLO),
        "n_perm_naive": len(perms),
        "n_donors": len(donors),
        "n_test_donors": len(true_test),
        "real_split_inside_perm_range": inside,
        "real_split_percentile_in_perms": percentile,
    }


# ---------------------------------------------------------------------------
# Universes.
# ---------------------------------------------------------------------------


def _fallback_universe(work, donor_col: str) -> list[str]:
    """The non-integer path's gene universe: detection in at least half the donors, label-agnostic.

    Defined in full by :data:`FALLBACK_UNIVERSE_RULE`, which the report renders verbatim.
    """
    donors = work.obs[donor_col].astype(str).to_numpy()
    unique_donors = np.unique(donors)
    detected = np.zeros(work.n_vars, dtype=np.int64)
    matrix = work.X
    for donor in unique_donors:
        block = matrix[donors == donor]
        if sp.issparse(block):
            present = np.asarray((block > 0).sum(axis=0)).ravel() > 0
        else:
            present = np.asarray(block > 0).any(axis=0)
        detected += present.astype(np.int64)
    needed = int(np.ceil(FALLBACK_UNIVERSE_MIN_PROP * unique_donors.size))
    return sorted(str(gene) for gene, keep in zip(work.var_names, detected >= needed) if keep)


# ---------------------------------------------------------------------------
# The audit itself.
# ---------------------------------------------------------------------------


def _design_block(report, work, donor_col: str) -> dict:
    """The payload's ``design`` block.

    ``cells_per_donor`` is the per-donor cell count the schema names, computed here rather than
    copied from ``DesignReport.cells_per_donor``, which is a min/median/max summary of a different
    shape.
    """
    per_donor = work.obs.groupby(donor_col, observed=True).size()
    return {
        "condition_col": report.condition_col,
        "donor_col": report.donor_col,
        "n_cells": int(report.n_cells),
        "groups": {str(k): int(v) for k, v in report.groups.items()},
        "donors_per_group": {str(k): int(v) for k, v in report.donors_per_group.items()},
        "cells_per_donor": {str(k): int(v) for k, v in per_donor.items()},
        "donor_nests_in_condition": bool(report.donor_nests_in_condition),
        "min_donors_per_group": int(report.min_donors_per_group),
        "imbalance_ratio": float(report.imbalance_ratio),
        "batch_confounded": {str(k): float(v) for k, v in report.batch_confounded.items()},
        "batch_separates_condition": {str(k): bool(v)
                                      for k, v in report.batch_separates_condition.items()},
        "flags": [str(f) for f in report.flags],
        "min_donors": int(report.min_donors),
        "usable_for_pseudobulk": bool(report.usable_for_pseudobulk),
    }


def _settings_block(settings: AuditSettings) -> dict:
    """The payload's ``settings`` block: the tool's own values and the pre-registered ones apart."""
    return {
        "tool": {
            "n_perm_requested": int(settings.n_perm),
            "n_perm_pb_requested": int(settings.n_perm_pb),
            "seed": int(settings.seed),
            "top_n": int(settings.top_n),
            "design_only": bool(settings.design_only),
            "naive_method": "wilcoxon",
            "naive_engine": "fast",
            "pseudobulk_method": "moderated_ebayes",
            "trend": False,
            "min_donors_per_group": PRODUCT_MIN_DONORS_PER_GROUP,
            "universe_min_total_count": PRODUCT_UNIVERSE_MIN_TOTAL_COUNT,
            "universe_min_prop": PRODUCT_UNIVERSE_MIN_PROP,
            "fallback_universe_min_prop": FALLBACK_UNIVERSE_MIN_PROP,
            "fallback_universe_min_size": FALLBACK_UNIVERSE_MIN_SIZE,
            "min_profiles_per_group_after_thin_filter": PRODUCT_MIN_PROFILES_PER_GROUP,
        },
        "protocol_constants": {
            "alpha": float(gate_config.ALPHA),
            "lambda_band": [float(v) for v in gate_config.LAMBDA_BAND],
            "min_universe_size": int(gate_config.MIN_UNIVERSE_SIZE),
            "min_cells": int(gate_config.MIN_CELLS),
            "min_counts": int(gate_config.MIN_COUNTS),
        },
    }


def _provenance_block() -> dict:
    """Where every frozen number in this payload came from, and what produced the rest."""
    manifest = gate_config.manifest()
    return {
        "pre_registered": dict(manifest["pre_registered"]),
        "pre_registered_source": str(manifest["pre_registered_source"]),
        "operating_envelope": [
            {
                "sigma_donor": float(row["sigma_donor"]),
                "min_donors_per_group": int(row["min_donors_per_group"]),
                "grid_support": str(row["grid_support"]),
            }
            for row in gate_config.OPERATING_ENVELOPE
        ],
        "platform": platform.platform(),
        "python": platform.python_version(),
        "packages": _package_versions(),
    }


def _caveats(payload: dict, settings: AuditSettings) -> list[dict]:
    """The caveat block: N1, N2, N3 and N7 always, the rest exactly under their own condition.

    Every text comes from ``pbcheck.render.text`` and every interpolated value from the payload
    built above, so no number in the prose is typed by hand here. The conditions are the payload's
    own flags where the schema defines one (``readout.few_donors`` for N4,
    ``readout.paired_floor_shown`` for R1's pseudobulk clause), so validation cannot find a note
    whose condition this same payload denies.
    """
    design = payload["design"]
    readout = payload["readout"]
    protocol_names = ", ".join(payload["settings"]["protocol_constants"])

    out = [
        {"id": "N1", "text": caveat_text("N1", version=payload["pbcheck_version"])},
        {"id": "N2", "text": caveat_text("N2", envelope_rows=envelope_rows())},
        {"id": "N3", "text": caveat_text("N3")},
    ]

    if readout["few_donors"]:
        out.append({"id": "N4", "text": caveat_text("N4")})

    separating = [col for col, sep in design["batch_separates_condition"].items() if sep]
    if separating:
        out.append({"id": "N5", "text": caveat_text("N5", cols=", ".join(separating))})

    if payload["status_reason"] == "non_integer_counts":
        reason = (payload["counts_check"] or {}).get("reason") or "no reason recorded"
        out.append({"id": "N6", "text": caveat_text("N6", reason=reason)})

    out.append({"id": "N7", "text": caveat_text("N7", protocol_constant_names=protocol_names)})

    if settings.celltype_col is None and payload["input"]["celltype_like_columns_found"]:
        out.append({"id": "N8", "text": caveat_text(
            "N8", col=payload["input"]["celltype_like_columns_found"][0])})

    achieved = readout["n_perm_naive_achieved"]
    if achieved is not None and achieved < COARSE_NULL_THRESHOLD:
        floor = readout["naive_floor_solo"] or {}
        out.append({"id": "N9", "text": caveat_text(
            "N9", n=achieved, requested=readout["n_perm_naive_requested"],
            se=f"{floor.get('mc_se', float('nan')):.3g}")})

    if payload["status"] == "complete":
        floor = readout["naive_floor_solo"]
        paired_floor_shown = bool(readout["paired_floor_shown"])
        pb_floor = payload["permutation_null"]["pseudobulk"]["floor"]
        out.append({"id": "R1", "text": readout_caveat_text(
            floor_solo=f"{floor['median_count']:.0f}",
            universe_size=payload["universe"]["size"],
            alpha=f"{settings.alpha:g}",
            floor_pct=f"{100.0 * floor['median_frac']:.1f}",
            real_solo=readout["naive_real_solo"],
            few_donors=bool(readout["few_donors"]),
            paired_floor_shown=paired_floor_shown,
            pb_real=(payload["real_label"]["pseudobulk"]["n_significant_paired"]
                     if paired_floor_shown else None),
            pb_floor=f"{pb_floor['median_count']:.0f}" if paired_floor_shown else None,
        )})
    return out


def format_scalar(value: object) -> str:
    """Render a payload scalar for a read-out sentence: ``None`` as ``"n/a"``, a bool as
    ``"yes"``/``"no"``, a float to four significant figures, everything else via ``str``.

    The report's number convention, applied here because the read-out lines are payload fields
    (``readout.sentences``) built by this module and printed verbatim by the renderer.
    ``tests/test_audit.py`` pins it against the renderer's own cell formatter, so a table cell and
    a sentence can never show the same number differently.
    """
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def _readout_sentences(payload: dict) -> list[str]:
    """The plain-language read-out lines, filled from the payload's own blocks.

    Written into ``readout.sentences`` and printed verbatim by ``pbcheck.render``: the numbers a
    reader sees in the read-out paragraph are the numbers this module put in the payload, not a
    second derivation made at render time. Empty when no permutation null ran (a ``design_only``
    status has nothing to read out).
    """
    null = payload["permutation_null"]
    if null is None:
        return []

    design = payload["design"]
    input_block = payload["input"]
    universe = payload["universe"]
    readout = payload["readout"]
    real = payload["real_label"]
    naive_null = null["naive"]
    floor = naive_null["floor_solo"]
    achieved = readout["n_perm_naive_achieved"]

    lines = [
        sentence_text(
            "floor_solo",
            median_count=format_scalar(floor["median_count"]),
            median_frac_pct=format_scalar(floor["median_frac"] * 100),
            universe_size=universe["size"],
            mc_se=format_scalar(floor["mc_se"]),
            n_perm_achieved=format_scalar(achieved),
            coarse_note=(COARSE_NULL_NOTE
                         if achieved is not None and achieved < COARSE_NULL_THRESHOLD else ""),
            real_solo=format_scalar(readout["naive_real_solo"]),
        ),
        sentence_text(
            "lambda_naive",
            **{"lambda": format_scalar(naive_null["lambda"])},
            iqr=format_scalar(naive_null["lambda_iqr"]),
            class_word=LAMBDA_CLASS_WORDS[readout["lambda_naive_class"]]
            if readout["lambda_naive_class"] is not None else format_scalar(None),
            band_lo=gate_config.LAMBDA_BAND[0],
            band_hi=gate_config.LAMBDA_BAND[1],
        ),
    ]

    pseudobulk_null = null["pseudobulk"]
    if pseudobulk_null is None:
        lines.append(sentence_text(
            "pseudobulk_not_run",
            reason_text=STATUS_REASON_WORDS.get(
                payload["status_reason"] or "", "no reason recorded"),
        ))
    else:
        pb_real = real["pseudobulk"] if real is not None else None
        pb_class = readout["lambda_pseudobulk_class"]
        lines.append(sentence_text(
            "pseudobulk",
            **{"lambda": format_scalar(pseudobulk_null["lambda"])},
            class_word=(LAMBDA_CLASS_WORDS[pb_class] if pb_class is not None
                        else format_scalar(None)),
            band_lo=gate_config.LAMBDA_BAND[0],
            band_hi=gate_config.LAMBDA_BAND[1],
            fp_rate=format_scalar(pseudobulk_null["fp_rate"]),
            se=format_scalar(pseudobulk_null["fp_rate_mc_se"]),
            median_count=format_scalar(pseudobulk_null["floor"]["median_count"]),
            real_paired=format_scalar(
                pb_real["n_significant_paired"] if pb_real is not None else None),
        ))

    donors_per_group = design["donors_per_group"]
    lines.append(sentence_text(
        "donors",
        n_test=format_scalar(donors_per_group.get(input_block["test_level"])),
        test_level=input_block["test_level"],
        n_ref=format_scalar(donors_per_group.get(input_block["ref_level"])),
        ref_level=input_block["ref_level"],
        n_distinct_splits=format_scalar(null["n_distinct_splits"]),
    ))

    profiles = universe["profiles_per_group_after_thin_filter"]
    if profiles is not None:
        protocol = payload["settings"]["protocol_constants"]
        lines.append(sentence_text(
            "profiles",
            min_cells=protocol["min_cells"],
            min_counts=protocol["min_counts"],
            p_test=format_scalar(profiles.get(input_block["test_level"])),
            p_ref=format_scalar(profiles.get(input_block["ref_level"])),
        ))
    return lines


def run_audit(adata, settings: AuditSettings) -> dict:
    """Audit one stratum of ``adata`` and return the validated ``pbcheck-audit/1`` payload.

    The caller's object is never mutated: everything runs on a copy made by
    :func:`prepare_stratum`. The status transitions are the ones the product contract fixes, in
    this order: ``--design-only``, a donor measured under both conditions, too few donors, a matrix
    that is not raw counts, a universe below the minimum, too few pseudobulk profiles after the
    thin-donor filter. Any engine exception other than
    :class:`pbcheck.gene_universe.UniverseTooSmall` propagates unchanged; nothing here turns a
    failed run into a plausible-looking number.
    """
    started = time.perf_counter()
    payload = audit_schema.empty_payload()
    timings: dict = {}
    payload["pbcheck_version"] = __version__
    payload["settings"] = _settings_block(settings)
    payload["provenance"] = _provenance_block()

    with _stage(timings, "prepare"):
        work, input_block = prepare_stratum(adata, settings)
    payload["input"] = input_block
    celltype_col = settings.celltype_col if settings.celltype_col is not None else STRATUM_COL

    with _stage(timings, "design"):
        report = audit_design(
            work,
            condition_col=settings.condition_col,
            donor_col=settings.donor_col,
            batch_cols=list(settings.batch_cols),
            min_donors=PRODUCT_MIN_DONORS_PER_GROUP,
        )
    payload["design"] = _design_block(report, work, settings.donor_col)
    donors_per_group = payload["design"]["donors_per_group"]
    min_donors_per_group = min(donors_per_group.get(settings.test_level, 0),
                               donors_per_group.get(settings.ref_level, 0))
    payload["readout"]["min_donors_per_group"] = min_donors_per_group
    payload["readout"]["few_donors_threshold"] = int(FEW_DONORS_THRESHOLD)
    payload["readout"]["few_donors"] = bool(min_donors_per_group < FEW_DONORS_THRESHOLD)
    payload["readout"]["n_perm_naive_requested"] = int(settings.n_perm)
    payload["readout"]["n_perm_pb_requested"] = int(settings.n_perm_pb)

    def _finish(status: str, reason: str | None) -> dict:
        payload["status"] = status
        payload["status_reason"] = reason
        payload["schema_version"] = audit_schema.AUDIT_SCHEMA_VERSION
        payload["generated_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        payload["runtime_by_stage_seconds"] = {
            key: timings.get(key) for key in audit_schema.RUNTIME_STAGE_KEYS
        }
        payload["runtime_seconds"] = float(time.perf_counter() - started)
        payload["caveats"] = _caveats(payload, settings)
        payload["readout"]["sentences"] = _readout_sentences(payload)
        audit_schema.validate(payload)
        return payload

    if settings.design_only:
        return _finish("design_only", "design_only_requested")
    if not report.donor_nests_in_condition:
        return _finish("design_only", "donor_spans_conditions")
    if min_donors_per_group < report.min_donors:
        return _finish("design_only", "too_few_donors")

    groups = payload["design"]["groups"]
    if len(groups) != 2:
        # Belt and braces for the arms: ``prepare_stratum`` keeps exactly the two levels, so this
        # can only fire if a future selection path lets a single-group stratum through. Refused
        # with the same error the input checks use rather than handed to the DE engines.
        raise AuditInputError(
            f"the audited stratum has {len(groups)} group(s) in {settings.condition_col!r} "
            f"({sorted(groups)!r}); a contrast needs exactly the two levels "
            f"{settings.test_level!r} and {settings.ref_level!r}"
        )

    with _stage(timings, "counts"):
        counts_work, counts_source, counts_check = resolve_counts(work, settings)
    payload["counts_check"] = counts_check.as_dict()
    payload["input"]["counts_source"] = counts_source

    if counts_source is None:
        # Non-integer matrix: the pseudobulk arm is dropped, never rounded, and the naive arm runs
        # on the matrix as found (the source that was checked), with its own normalisation on top.
        source_label, matrix, var = _counts_candidates(work, settings)[0]
        del source_label
        naive_work = _working_adata(matrix, work.obs, var)
        with _stage(timings, "universe"):
            universe = _fallback_universe(naive_work, settings.donor_col)
            payload["universe"].update({
                "size": len(universe),
                "min_size": FALLBACK_UNIVERSE_MIN_SIZE,
                "builder": "naive_detection_fallback",
                "builder_rule": FALLBACK_UNIVERSE_RULE,
            })
        if len(universe) < FALLBACK_UNIVERSE_MIN_SIZE:
            return _finish("design_only", "universe_too_small")
        _run_naive_only(payload, naive_work, universe, settings, timings)
        return _finish("naive_only", "non_integer_counts")

    with _stage(timings, "pseudobulk_build"):
        pdata = build_pseudobulk(counts_work, donor_col=settings.donor_col,
                                 celltype_col=celltype_col, condition_col=settings.condition_col)
        profiles = pdata.obs[settings.condition_col].astype(str).value_counts().to_dict()
        profiles_per_group = {str(k): int(v) for k, v in profiles.items()}

    with _stage(timings, "universe"):
        payload["universe"].update({
            "min_size": int(gate_config.MIN_UNIVERSE_SIZE),
            "builder": "pseudobulk_frozen",
            "builder_rule": PSEUDOBULK_UNIVERSE_RULE,
            "thin_donor_filter": dict(pdata.uns.get("thin_donor_filter", {})) or None,
            "profiles_per_group_after_thin_filter": profiles_per_group,
        })
        try:
            universe = gene_universe.frozen_universe(
                pdata,
                min_total_count=PRODUCT_UNIVERSE_MIN_TOTAL_COUNT,
                min_prop=PRODUCT_UNIVERSE_MIN_PROP,
                min_size=gate_config.MIN_UNIVERSE_SIZE,
            )
        except gene_universe.UniverseTooSmall:
            # The one engine exception this module answers instead of propagating: it is the
            # protocol's own SKIP verdict on a degenerate stratum, and the status carries it.
            universe = gene_universe.frozen_universe(
                pdata,
                min_total_count=PRODUCT_UNIVERSE_MIN_TOTAL_COUNT,
                min_prop=PRODUCT_UNIVERSE_MIN_PROP,
            )
            payload["universe"]["size"] = len(universe)
            return _finish("design_only", "universe_too_small")
        payload["universe"]["size"] = len(universe)

    min_profiles = min(profiles_per_group.get(settings.test_level, 0),
                       profiles_per_group.get(settings.ref_level, 0))
    payload["readout"]["min_profiles_per_group_after_thin_filter"] = min_profiles
    if min_profiles < PRODUCT_MIN_PROFILES_PER_GROUP:
        _run_naive_only(payload, counts_work, universe, settings, timings)
        return _finish("naive_only", "too_few_profiles_after_thin_filter")

    _run_both_arms(payload, counts_work, pdata, universe, settings, timings, celltype_col)
    return _finish("complete", None)


def _readout_from(payload: dict, universe: list[str]) -> None:
    """Fill the read-out from the blocks already in ``payload``. No new thresholds live here."""
    readout = payload["readout"]
    real = payload["real_label"]
    null = payload["permutation_null"]

    readout["lambda_naive_class"] = _lambda_class(null["naive"]["lambda"])
    floor_solo = null["naive"]["floor_solo"]
    readout["naive_floor_solo"] = {
        "median_count": float(floor_solo["median_count"]),
        "median_frac": float(floor_solo["median_frac"]),
        "mc_se": float(floor_solo["mc_se"]),
        "n_perm": int(floor_solo["n_perm"]),
    }
    naive_real_solo = int(real["naive"]["n_significant_solo"])
    readout["naive_real_solo"] = naive_real_solo
    readout["naive_real_over_floor_solo"] = metrics.signal_above_floor(
        naive_real_solo, floor_solo["median_count"])
    readout["paired_floor_shown"] = bool(real["paired_bh"]["n_na_pseudobulk"] == 0
                                         and null["naive"]["floor_paired"] is not None)
    readout["n_perm_naive_achieved"] = int(null["n_perm_naive_achieved"])
    readout["n_perm_pb_achieved"] = null["n_perm_pb_achieved"]
    readout["n_distinct_splits"] = int(null["n_distinct_splits"])

    pseudobulk = null["pseudobulk"]
    if pseudobulk is None:
        readout["lambda_pseudobulk_class"] = None
        readout["pseudobulk_real_over_floor"] = None
    else:
        readout["lambda_pseudobulk_class"] = _lambda_class(pseudobulk["lambda"])
        readout["pseudobulk_real_over_floor"] = metrics.signal_above_floor(
            int(real["pseudobulk"]["n_significant_paired"]), pseudobulk["floor"]["median_count"])
    del universe


def _run_naive_only(payload: dict, work, universe: list[str], settings: AuditSettings,
                    timings: dict) -> None:
    """The naive arm and its donor-permutation null, with no pseudobulk counterpart.

    Used by both statuses that have no pseudobulk arm. The real-label naive result is the one the
    permutation engine produced for the true labeling: it is the same statistic
    ``pbcheck.methods.naive_de`` computes, bit for bit (``tests/test_naive_engine.py``), so
    ranking the stratum a second time would only cost time.
    """
    with _stage(timings, "permutation_null"):
        null = naive_null(work, universe, settings)

    with _stage(timings, "real_label"):
        real_solo = mtc.bh_over_universe(null["naive_real"], universe, alpha=settings.alpha)
        n_tested = int(len(real_solo.genes_tested))
        payload["real_label"] = {
            "naive": {
                "n_significant_paired": None,
                "n_significant_solo": int(real_solo.n_significant(fdr=settings.alpha)),
                "n_tested": n_tested,
                "top": _top_genes(real_solo.table, settings.top_n,
                                  ("pval", "padj", "log2fc", "pct_group", "pct_reference")),
            },
            "pseudobulk": None,
            # No gene is testable in both arms when only one arm ran, and the paired bookkeeping
            # says exactly that rather than pretending the pseudobulk arm tested everything.
            "paired_bh": {
                "n_universe": len(universe),
                "n_tested_common": 0,
                "n_na_naive": len(universe) - n_tested,
                "n_na_pseudobulk": len(universe),
                "n_dropped_for_fairness": len(universe),
                "pseudobulk_na_free": False,
            },
            "consistent_with_permutation_path": None,
        }

    # Same reason as in :func:`_run_both_arms`: the lambda block runs the empirical permutation
    # p-values over the whole matrix, which is the null's work and belongs to the null's stage.
    with _stage(timings, "permutation_null"):
        payload["permutation_null"] = {
            "naive": {
                **_lambda_block(null["naive_pvals"], null["naive_pvals_real"]),
                "floor_solo": _floor_block(null["naive_ndeg_solo"], len(universe)),
                "floor_paired": None,
            },
            "pseudobulk": None,
            "monte_carlo": None,
            "real_split_inside_perm_range": null["real_split_inside_perm_range"],
            "real_split_percentile_in_perms": null["real_split_percentile_in_perms"],
            "n_donors": int(null["n_donors"]),
            "n_distinct_splits": _n_distinct_splits(null["n_donors"], null["n_test_donors"]),
            "n_perm_naive_achieved": int(null["n_perm_naive"]),
            "n_perm_pb_achieved": None,
            "n_perm_paired": None,
            "engine_path": "naive_null",
        }
        _readout_from(payload, universe)


def _run_both_arms(payload: dict, work, pdata, universe: list[str], settings: AuditSettings,
                   timings: dict, celltype_col: str) -> None:
    """Both arms on the real labels, then both arms under the donor-permutation null."""
    with _stage(timings, "real_label"):
        naive = naive_de(work, condition_col=settings.condition_col,
                         test_level=settings.test_level, ref_level=settings.ref_level,
                         genes=universe)
        naive_solo = mtc.bh_over_universe(naive, universe, alpha=settings.alpha)
        pseudobulk = pseudobulk_de(work, donor_col=settings.donor_col, celltype_col=celltype_col,
                                   condition_col=settings.condition_col,
                                   test_level=settings.test_level, ref_level=settings.ref_level,
                                   universe=universe)
        paired = mtc.bh_both_arms(naive, pseudobulk, universe, alpha=settings.alpha)
        payload["real_label"] = {
            "naive": {
                "n_significant_paired": int(paired.naive.n_significant(fdr=settings.alpha)),
                "n_significant_solo": int(naive_solo.n_significant(fdr=settings.alpha)),
                "n_tested": int(len(naive_solo.genes_tested)),
                "top": _top_genes(naive_solo.table, settings.top_n,
                                  ("pval", "padj", "log2fc", "pct_group", "pct_reference")),
            },
            "pseudobulk": {
                "n_significant_paired": int(paired.pseudobulk.n_significant(fdr=settings.alpha)),
                "n_tested": int(len(paired.pseudobulk.genes_tested)),
                "top": _top_genes(paired.pseudobulk.table, settings.top_n,
                                  ("pval", "padj", "log2fc")),
                "moderation": dict(getattr(pseudobulk, "moderation", {}) or {}),
            },
            "paired_bh": {
                "n_universe": int(paired.n_universe),
                "n_tested_common": int(paired.n_tested_common),
                "n_na_naive": int(paired.n_na_naive),
                "n_na_pseudobulk": int(paired.n_na_pseudobulk),
                "n_dropped_for_fairness": int(len(paired.dropped_for_fairness)),
                "pseudobulk_na_free": bool(paired.n_na_pseudobulk == 0),
            },
            "consistent_with_permutation_path": None,
        }

    with _stage(timings, "permutation_null"):
        result = run_null(
            work, universe,
            donor_col=settings.donor_col, condition_col=settings.condition_col,
            test_level=settings.test_level, ref_level=settings.ref_level,
            celltype_col=celltype_col,
            n_perm=settings.n_perm, n_perm_pb=settings.n_perm_pb, fdr=settings.alpha,
            seed=settings.seed, naive_engine="fast",
        )

    # Summarising the null is part of the null's cost: the empirical permutation p-values behind
    # every lambda are computed over the whole permutation matrix here, so the stage is reopened
    # instead of leaving that work outside the per-stage table.
    with _stage(timings, "permutation_null"):
        payload["real_label"]["consistent_with_permutation_path"] = bool(
            paired.n_tested_common == result["paired_bh_real"]["n_tested_common"]
        )

        monte_carlo = result["monte_carlo"]
        floor_paired = _floor_block(result["naive_ndeg_paired"], len(universe))
        floor_paired["mc_se"] = float(monte_carlo["naive_floor_mc_se"])
        pb_floor = _floor_block(result["pb_ndeg"], len(universe))
        pb_floor["mc_se"] = float(monte_carlo["pb_floor_mc_se"])

        tmap = _donor_condition_map(work.obs, settings.donor_col, settings.condition_col)
        n_test_donors = int((tmap == settings.test_level).sum())

        payload["permutation_null"] = {
            "naive": {
                **_lambda_block(result["naive_pvals"], result["naive_pvals_real"]),
                "floor_solo": _floor_block(result["naive_ndeg_solo"], len(universe)),
                "floor_paired": floor_paired,
            },
            "pseudobulk": {
                **_lambda_block(result["pb_pvals"], result["pb_pvals_real"]),
                "floor": pb_floor,
                "fp_rate": float(monte_carlo["pb_fp_rate"]),
                "fp_rate_mc_se": float(monte_carlo["pb_fp_rate_mc_se"]),
            },
            "monte_carlo": {k: (float(v) if isinstance(v, (int, float)) else v)
                            for k, v in monte_carlo.items()},
            "real_split_inside_perm_range": bool(result["real_split_inside_perm_range"]),
            "real_split_percentile_in_perms": float(result["real_split_percentile_in_perms"]),
            "n_donors": int(result["n_donors"]),
            "n_distinct_splits": _n_distinct_splits(int(result["n_donors"]), n_test_donors),
            "n_perm_naive_achieved": int(result["n_perm_naive"]),
            "n_perm_pb_achieved": int(result["n_perm_pb"]),
            "n_perm_paired": int(result["n_perm_paired"]),
            "engine_path": "run_null",
        }
        _readout_from(payload, universe)


def audit_h5ad(path: str | Path, settings: AuditSettings) -> dict:
    """Read an ``.h5ad`` file and audit it; the payload records the path and the load time."""
    file_path = Path(path)
    started = time.perf_counter()
    adata = ad.read_h5ad(file_path)
    load_seconds = float(time.perf_counter() - started)

    payload = run_audit(adata, settings)
    payload["input"]["path"] = str(file_path)
    payload["runtime_by_stage_seconds"]["load"] = load_seconds
    payload["runtime_seconds"] = float(payload["runtime_seconds"] + load_seconds)
    audit_schema.validate(payload)
    return payload
