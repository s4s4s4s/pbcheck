"""Candidate stratum manifest over the whole pinned Census — the streaming, two-pass driver.

Runs :mod:`pbcheck.census_select` against the real CELLxGENE Census and writes the candidate
manifest (``census_candidates.json`` / ``.csv``). **Nothing here is a method.** Every gate, every
threshold, every column and every exclusion reason comes from ``census_select``; this file only
decides what to read, in what order, and how to survive a multi-hour network job on a 16 GB runner.
No function in ``census_select`` is modified, wrapped or re-implemented, and no CLI option of this
driver can move a pre-registered number (§1 thresholds, the ``value_filter``, the Census pin) —
that is a protocol change and goes through ``docs/AMENDMENTS.md`` first.

**No row of the output is an admission.** ``admitted_to_sweep`` is ``False`` on every row,
inclusion-gate items 4 and 5 stay ``pending``, and the §1 pre-registration of the stratum list is a
separate, deliberate act. This is a list of *candidates*.

Why two passes and not one query
--------------------------------
The §1 ``value_filter`` leaves roughly 50-65 M cells of the pinned Census.
:func:`pbcheck.census_select.query_obs` materialises its slice (``.concat().to_pandas()``), and at
that row count the peak — the concat plus the pandas conversion — lands somewhere between 8 GB and
30-something GB depending on whether the Census hands ``donor_id`` back dictionary-encoded or as
Python objects. A public runner has 16 GB and cannot be *told* which it will be, and an OOM kill
(exit 137) five hours in tells you nothing. So the driver removes the question rather than betting
on the answer:

**Pass 1 — streaming, coarse.** Four columns (:data:`KEY_COLUMNS`, exactly
``census_select.REQUIRED_OBS_COLUMNS``) read off the SOMA reader **without** ``.concat()``, one
Arrow batch at a time, each batch folded immediately into a
``(dataset_id, cell_type, disease, donor_id) -> cells`` counter and then dropped. The accumulator
grows with the number of distinct keys (10^5-10^6, i.e. tens to hundreds of MB), not with the
number of cells. Addition is commutative and associative, so the totals do not depend on where the
Census happens to cut its batches — which is what makes the run reproducible at all: the pinned
Census pins the *data*, not the batch boundaries.

**Pass 2 — per dataset, exact.** Pass 1 says which datasets could possibly hold a stratum that
clears the gate; each of those is re-read with one scoped ``query_obs`` and handed to
``screen_strata`` unchanged. Everything the manifest reports — the gate, the confound pre-screen,
the pooling flag, D4's excluded fractions — is computed by ``census_select`` on real per-cell rows.

Why pass 2 re-reads the cells instead of reusing pass 1's counts
----------------------------------------------------------------
:func:`pbcheck.census_select.apply_inclusion_gate` counts **rows**
(``groupby(donor_id).size()`` / ``value_counts()``). Feeding it the pass-1 aggregate — one row per
donor — would tell it that every donor has exactly one cell, every donor would fall to the
thin-donor rule, and the entire Census would report as failing the gate. The aggregate decides
*which datasets to read*; only per-cell frames ever reach the gate.

What the coarse filter may and may not do
-----------------------------------------
:func:`surviving_datasets` implements one thing: the §1 item-2/item-1 pair (donors with at least
``gate_config.MIN_CELLS`` cells, at least ``census_select.MIN_DONORS_PER_GROUP`` of them in each
arm), on the same ``(dataset_id, cell_type, disease)`` split the gate uses. It is deliberately
*conservative* — it may pass a dataset the gate later rejects (item 3's nesting rule, for one, is
not checked here), and it must never do the reverse. ``tests/test_census_candidates.py`` cross-
checks that implication against the real gate on the same synthetic rows, because a filter that
silently drops a usable stratum would be invisible in the output.

Not resumable, on purpose: the Census is pinned, so a re-dispatch reproduces the run exactly, and a
resume file is one more thing that can be stale. Pass 1 gets one full retry (the SOMA reader is
stateful — a broken read cannot be resumed mid-stream, only restarted); each pass-2 query gets two,
with backoff, and a dataset that still fails is recorded in the manifest's ``notes`` rather than
quietly missing from it.
"""

from __future__ import annotations

import argparse
import signal
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from pbcheck import census_select as cs  # noqa: E402
from pbcheck import gate_config as gc  # noqa: E402

#: Written under ``pilot/results/``, which ``.gitignore`` excludes: in CI this leaves as an
#: artifact, and a candidate manifest is not committed as a side effect of a run.
DEFAULT_OUT_DIR = REPO / "pilot" / "results" / "census_candidates"

#: Not ``census_strata`` (``emit_manifest``'s default): the file name should say on sight that this
#: is the candidate list and not the pre-registered one.
MANIFEST_STEM = "census_candidates"

#: The four columns pass 1 streams — the stratum key, the donor, the disease term. Same set as
#: :data:`pbcheck.census_select.REQUIRED_OBS_COLUMNS` (pinned by a test), ordered as the counter
#: key rather than as the spec lists them.
KEY_COLUMNS = ("dataset_id", "cell_type", "disease", "donor_id")

#: Terms that name no disease and therefore form no contrast — mirrors ``build_contrasts``. The
#: §1 ``value_filter`` already excludes ``'na'``; this is the second copy that keeps the coarse
#: filter's arithmetic aligned with the gate's, and reads the real one when it can.
NON_DISEASE_TERMS = frozenset(getattr(cs, "_NON_DISEASE_TERMS", frozenset({"na"})))

#: A dry run reads five datasets. Not the five biggest and not the five smallest — see
#: :func:`spread_sample`: the point of a dry run is to calibrate the timing of the real one.
DRY_RUN_MAX_DATASETS = 5

#: Attempts per pass-2 query (one try plus two retries) and the waits between them. Pass 1 is
#: one retry of the WHOLE pass, because a stateful reader cannot resume a chunk.
PASS2_ATTEMPTS = 3
RETRY_BACKOFF_SECONDS = (2.0, 8.0)

#: Exceptions that mean "the code, the schema or the machine is wrong", not "the network hiccuped".
#: Retrying these burns the job's budget and buries the message that would explain the failure.
FATAL_EXCEPTIONS = (
    ValueError, TypeError, KeyError, IndexError, AttributeError,
    ImportError, MemoryError, NotImplementedError,
)

TOP_N_SUMMARY = 15

#: Read-buffer budget per column, in bytes, handed to the Census as ``py.init_buffer_bytes`` /
#: ``soma.init_buffer_bytes``.
#:
#: cellxgene-census's own ``DEFAULT_TILEDB_CONFIGURATION`` sets **1 GiB each**, which is a fine
#: default for a workstation and fatal here: the reader allocates that per column before the first
#: batch arrives, and pass 1 asks for four string columns (each of which needs data, offset and
#: validity buffers). The second live dry run died 22 seconds in with SIGTERM (exit 143 — the
#: hosted runner's memory-pressure kill, not the kernel OOM killer's SIGKILL) before a single
#: progress line, which is the shape of an allocation that never got to read anything.
#:
#: 128 MiB is the value cellxgene-census's own documentation uses in its example of overriding
#: these keys, and it is a budget, not a correctness knob: buffer size cannot change which cells
#: come back, only how many arrive per batch. Exposed as ``--reader-buffer-mb`` so the next
#: dispatch can go lower without a code change.
READER_BUFFER_BYTES = 128 * 1024**2

#: Progress lines for the first few batches, then every ``--progress-every``. Calibrating
#: cells/second should not require surviving twenty batches — the first one is the interesting one.
PROGRESS_FIRST_N = 5
PROGRESS_EVERY = 20


# ---------------------------------------------------------------------------
# Pass 1 — streaming aggregation. Pure functions below iter_obs_batches, so the arithmetic is
# testable on a machine that has no cellxgene-census (which is the development machine).
# ---------------------------------------------------------------------------

def reader_tiledb_config(buffer_bytes: int = READER_BUFFER_BYTES) -> dict:
    """The TileDB overrides pass 1 opens the Census with — a memory budget and nothing else.

    Both keys are set because both exist and cellxgene-census's own default sets both: ``py.*`` is
    the TileDB-Py path and ``soma.*`` the SOMA one. ``open_soma(tiledb_config=...)`` merges these
    over :data:`cellxgene_census.DEFAULT_TILEDB_CONFIGURATION` rather than replacing it, so the S3
    region, the anonymous credentials and the user-agent all survive untouched.

    Nothing here can change *which* cells are read: a buffer size decides how many rows arrive per
    batch, and pass 1's arithmetic is batch-boundary independent by construction.
    """
    buffer_bytes = int(buffer_bytes)
    return {"py.init_buffer_bytes": buffer_bytes, "soma.init_buffer_bytes": buffer_bytes}


def _rss_bytes():
    """Resident set size, from ``/proc`` on Linux. ``None`` anywhere else — no new dependency.

    The runner kills on RSS, so RSS is the number worth printing. It is read rather than estimated
    because the interesting allocations here (TileDB read buffers, Arrow buffers) live outside
    Python's own allocator and are invisible to anything that counts Python objects.
    """
    try:
        with open("/proc/self/status", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        return None
    return None


def _peak_rss_bytes():
    """High-water RSS for this process, or ``None`` on a platform without ``resource``."""
    try:
        import resource
    except ImportError:  # Windows, where this driver is developed but never run for real
        return None
    # Linux reports ru_maxrss in kibibytes (macOS in bytes); CI is Linux and that is the number
    # this line exists to report.
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def _fmt_bytes(n) -> str:
    if n is None:
        return "n/a"
    return f"{n / 1024**2:,.1f} MiB"


def _install_fatal_signal_logging() -> None:
    """Print RSS when the runner kills us, so the next post-mortem starts from a number.

    A hosted runner under memory pressure sends **SIGTERM** (exit 143), which is catchable — that
    is how the second dry run died, silently, 22 seconds in. (The kernel OOM killer's SIGKILL is
    not catchable and shows up as 137; if that is what happens next, the distinction is itself
    diagnostic.) The handler raises rather than exiting quietly, so the ``finally`` blocks close
    the Census handle and the exit status still says "killed by signal".
    """
    def handler(signum, _frame):
        print(f"[census-candidates] FATAL: killed by signal {signum} "
              f"({signal.Signals(signum).name}) — rss={_fmt_bytes(_rss_bytes())} "
              f"peak_rss={_fmt_bytes(_peak_rss_bytes())}", file=sys.stderr, flush=True)
        raise SystemExit(128 + signum)

    for sig in (signal.SIGTERM, signal.SIGINT):
        try:
            signal.signal(sig, handler)
        except (ValueError, OSError, AttributeError):  # not the main thread, or no such signal
            pass


def iter_obs_batches(
    census,
    *,
    organism: str = cs.CENSUS_ORGANISM,
    value_filter: str = cs.VALUE_FILTER,
    columns=KEY_COLUMNS,
):
    """Yield the §1 obs slice one **Arrow** batch at a time. **No** ``.concat()``, no pandas.

    ``soma_obs.read(...)`` returns a reader that is already an iterator of ``pyarrow.Table``
    chunks; ``query_obs`` calls ``.concat()`` on it, which is the full materialisation this pass
    exists to avoid. Iterating it directly costs nothing extra and never holds more than one chunk.

    The tables are handed on **as Arrow**. Calling ``.to_pandas()`` here was the driver's second
    memory bug: four string columns become Python ``object`` arrays, i.e. one heap-allocated
    ``str`` per cell per column — the most expensive representation available, for data that is
    about to be reduced to a few thousand group counts. :func:`aggregate_batches` groups in Arrow
    instead and only ever materialises the grouped keys.

    Only :data:`KEY_COLUMNS` are requested. The optional columns (``raw_sum`` and friends) are not
    needed to decide *which datasets to read* and would multiply the bytes crossing the network by
    the least useful factor available.
    """
    soma_obs = census["census_data"][organism].obs
    reader = soma_obs.read(column_names=list(columns), value_filter=value_filter)
    yield from reader


#: One message for both places that meet a batch they cannot read, so an unexpected reader type
#: says the same thing whichever of them notices it first.
_BATCH_TYPE_MESSAGE = (
    "cannot aggregate a batch of type {kind}: expected a pyarrow.Table (what the SOMA reader "
    "yields) or a pandas DataFrame"
)


def _batch_rows(batch) -> int:
    rows = getattr(batch, "num_rows", None)  # pyarrow.Table
    if rows is None:
        rows = getattr(batch, "shape", (None,))[0]  # pandas.DataFrame
    if rows is None:
        raise TypeError(_BATCH_TYPE_MESSAGE.format(kind=type(batch).__name__))
    return int(rows)


def _arrow_key_counts(table, key_columns):
    """``((key tuple), cells)`` for one Arrow batch, grouped in Arrow.

    Three things happen here and each of them is load-bearing:

    * **dictionary columns are decoded to their values.** The Census may hand back
      ``dictionary<values=string, indices=int32>``, and the indices are assigned *per batch* — code
      7 in one batch and code 7 in the next are unrelated. An accumulator keyed on codes would
      silently merge distinct donors, so the keys can only be the values. The decode goes through
      ``column.type.value_type`` rather than a hard-coded string type, so a ``large_string``
      dictionary stays a ``large_string``.
    * **nulls go before the grouping**, not after: ``Table.drop_null()`` drops a row with a null in
      *any* of the four columns, which is exactly what the gate does (a cell with no donor has no
      replication unit; a cell with no stratum key forms no stratum).
    * **only the grouped result becomes Python.** A batch of half a million cells collapses to at
      most a few thousand ``(dataset, cell type, disease, donor)`` keys before anything is turned
      into a ``str``.
    """
    keyed = table.select(list(key_columns))
    for i, name in enumerate(key_columns):
        column = keyed.column(name)
        value_type = getattr(column.type, "value_type", None)
        if value_type is not None:
            keyed = keyed.set_column(i, name, column.cast(value_type))

    keyed = keyed.drop_null()
    if not keyed.num_rows:
        return []

    grouped = keyed.group_by(list(key_columns)).aggregate([(key_columns[0], "count")])
    count_name = f"{key_columns[0]}_count"
    counts = (grouped.column(count_name) if count_name in grouped.schema.names
              else grouped.column(grouped.num_columns - 1)).to_pylist()
    keys = [grouped.column(name).to_pylist() for name in key_columns]
    return zip(zip(*keys), counts)


def _frame_key_counts(frame, key_columns):
    """The same, for a pandas frame — kept so the accumulator stays testable without pyarrow."""
    keyed = frame.loc[:, list(key_columns)].dropna()
    if not keyed.shape[0]:
        return []
    sizes = keyed.astype(str).groupby(list(key_columns), observed=True, sort=False).size()
    return sizes.items()


def _key_counts(batch, key_columns):
    if hasattr(batch, "group_by"):  # pyarrow.Table — the streaming path
        return _arrow_key_counts(batch, key_columns)
    if isinstance(batch, pd.DataFrame):
        return _frame_key_counts(batch, key_columns)
    raise TypeError(_BATCH_TYPE_MESSAGE.format(kind=type(batch).__name__))


def aggregate_batches(batches, *, key_columns=KEY_COLUMNS, progress_every: int = PROGRESS_EVERY,
                      label="pass1"):
    """Fold an iterable of obs batches into ``{(dataset, cell_type, disease, donor): cells}``.

    Returns ``(counts, n_cells_seen, n_batches)``. Pure apart from the progress print: it takes an
    iterable of batches — Arrow tables from the reader, or pandas frames from a test — so the whole
    of pass 1's arithmetic can be exercised with hand-built data and no network.

    Two properties this function must have, both tested:

    * **order- and chunking-independence.** Only ``+`` is used across batches, so re-cutting or
      re-ordering the stream cannot change the result. Batch boundaries are an implementation
      detail of the Census reader and are not pinned by the Census version; anything that depended
      on them would be irreproducible for a reason nobody would think to look for.
    * **nulls are dropped BEFORE grouping.** A cell with a null ``donor_id`` has no replication
      unit and is dropped by ``apply_inclusion_gate`` (§1 item 3); a null in any stratum column
      forms no stratum (``build_contrasts`` drops those too). Dropping them here keeps the coarse
      counter aligned with the gate — and under pandas 3.0 a null survives ``.astype(str)`` rather
      than becoming the string ``"nan"``, so a "nan" donor would otherwise be counted as a real
      one. ``n_cells_seen`` counts every row read, dropped ones included, so the discard is
      visible rather than silent.
    """
    counts: dict[tuple[str, str, str, str], int] = {}
    key_columns = list(key_columns)
    n_cells_seen = 0
    n_batches = 0
    t0 = time.time()

    for batch in batches:
        n_batches += 1
        rows = _batch_rows(batch)
        n_cells_seen += rows

        # The first batch is the one the last run never lived to report. Its rows, its bytes and
        # the RSS on either side of folding it are what turn "died at 22 s" into a cause.
        rss_before = _rss_bytes() if n_batches == 1 else None

        for key, n in _key_counts(batch, key_columns):
            counts[key] = counts.get(key, 0) + int(n)

        if not progress_every:  # 0 = quiet, for the tests that only want the arithmetic
            continue
        if n_batches == 1:
            print(f"[{label}] first batch: {rows:,} rows | "
                  f"{_fmt_bytes(getattr(batch, 'nbytes', None))} arrow | "
                  f"{len(counts):,} keys | rss {_fmt_bytes(rss_before)} -> "
                  f"{_fmt_bytes(_rss_bytes())} | {time.time() - t0:.1f}s", flush=True)
        elif n_batches <= PROGRESS_FIRST_N or n_batches % progress_every == 0:
            elapsed = max(time.time() - t0, 1e-9)
            print(
                f"[{label}] {n_batches} batches | {n_cells_seen:,} cells | {len(counts):,} keys | "
                f"{elapsed / 60:.1f} min | {n_cells_seen / elapsed:,.0f} cells/s | "
                f"rss {_fmt_bytes(_rss_bytes())}",
                flush=True,
            )
    return counts, n_cells_seen, n_batches


def surviving_datasets(
    counts,
    *,
    min_cells_per_donor: int = gc.MIN_CELLS,
    min_donors_per_group: int = cs.MIN_DONORS_PER_GROUP,
    reference: str = cs.REFERENCE_LEVEL,
):
    """Datasets that *could* hold a stratum clearing the §1 gate, plus each dataset's cell total.

    Returns ``(survivors, dataset_total_cells)``. A dataset survives when some
    ``(dataset_id, cell_type)`` has a disease term with at least ``min_donors_per_group`` donors of
    at least ``min_cells_per_donor`` cells **and** at least that many such ``normal`` donors in the
    same stratum — §1 inclusion-gate items 2 and 1, on the gate's own split.

    **Conservative by construction, and it has to be.** Item 3 (donor nested within condition) is
    not checked: a donor with cells under both conditions is counted as thin-or-fat in each, which
    can only *add* datasets here, and the gate rejects the contrast in pass 2 with its reason
    recorded. Nothing in this function can subtract a dataset the gate would have passed — for a
    contrast to pass, its donors are nested, so each donor's cells inside the contrast are exactly
    its cells in one ``(dataset, cell_type, disease)`` cell of this counter, which is what is
    thresholded here. That implication is the one thing the cross-check test pins.

    ``dataset_total_cells`` counts donor-labelled cells only (the null-donor rows never entered the
    counter). It sizes the pass-2 query, where a slight underestimate is immaterial next to the
    2 M-cell default cap; it is not reported as a dataset's cell count anywhere.
    """
    fat_donors: dict[tuple[str, str, str], int] = {}
    dataset_total_cells: dict[str, int] = {}

    for (dataset_id, cell_type, disease, _donor_id), n_cells in counts.items():
        n_cells = int(n_cells)
        dataset_total_cells[dataset_id] = dataset_total_cells.get(dataset_id, 0) + n_cells
        if n_cells >= min_cells_per_donor:
            stratum = (dataset_id, cell_type, disease)
            fat_donors[stratum] = fat_donors.get(stratum, 0) + 1

    survivors: set[str] = set()
    for (dataset_id, cell_type, disease), n_donors in fat_donors.items():
        if disease == reference or disease in NON_DISEASE_TERMS:
            continue
        if n_donors < min_donors_per_group:
            continue
        if fat_donors.get((dataset_id, cell_type, reference), 0) < min_donors_per_group:
            continue
        survivors.add(dataset_id)
    return survivors, dataset_total_cells


def dataset_cell_type_cells(counts, dataset_id) -> dict[str, int]:
    """``{cell_type: cells}`` for one dataset, from the pass-1 counter."""
    per_cell_type: dict[str, int] = {}
    for (ds_id, cell_type, _disease, _donor_id), n_cells in counts.items():
        if ds_id == dataset_id:
            per_cell_type[cell_type] = per_cell_type.get(cell_type, 0) + int(n_cells)
    return per_cell_type


def plan_queries(counts, dataset_id, *, total_cells: int, max_cells_per_query: int):
    """The pass-2 query scopes for one dataset: ``[None]``, or one cell type per query.

    A dataset bigger than ``max_cells_per_query`` is **split, not skipped**: pass 1 already knows
    its cell types, and a stratum is ``(dataset_id × cell_type)``, so cutting on ``cell_type``
    partitions the dataset's strata exactly — no contrast spans two cell types, and the union of
    the per-cell-type manifests is the manifest of the whole dataset (pinned by a test). That is
    also why the split stops there: a cut on ``disease`` would sever a contrast from its ``normal``
    arm, and a cut on ``donor_id`` would sever a group from its donors.

    A single cell type over the cap is queried anyway. The stratum is the atom of this analysis;
    dropping one because it is large would silently bias the candidate list toward small strata.
    """
    if max_cells_per_query <= 0 or total_cells <= max_cells_per_query:
        return [None]
    return sorted(dataset_cell_type_cells(counts, dataset_id)) or [None]


def spread_sample(ordered, n: int):
    """``n`` items spread evenly across ``ordered``, both ends included, order preserved.

    Used for ``--dry-run``: the point of a dry run is to calibrate the real one, and five datasets
    off either end of the size distribution calibrate nothing. Deterministic — with the Census
    pinned, the same dry run picks the same five datasets.
    """
    ordered = list(ordered)
    if n <= 0 or n >= len(ordered):
        return ordered
    if n == 1:
        return ordered[:1]
    last = len(ordered) - 1
    picks = sorted({round(i * last / (n - 1)) for i in range(n)})
    return [ordered[i] for i in picks]


# ---------------------------------------------------------------------------
# Pass 2 — one scoped query per dataset (or per cell type), screened by census_select as is.
# ---------------------------------------------------------------------------

def _escape(value) -> str:
    """``value`` as a **Python** string literal, quotes included, for a SOMA ``value_filter``.

    TileDB-SOMA parses a value filter with ``ast.parse(expression, mode="eval")`` — the expression
    is Python source, not SQL. The SQL escape (doubling the quote) is therefore not merely wrong
    here, it is *silently* wrong: ``'O''Brien'`` is Python's implicit concatenation of two adjacent
    literals and evaluates to ``"OBrien"``, a value that matches nothing. The query would succeed,
    return zero rows, and the stratum would be missing from the manifest with no error anywhere.
    ``repr`` emits a literal that Python's own parser reads back as the original string, which is
    exactly the contract the filter needs. Cell-type labels are free ontology text and may carry an
    apostrophe; dataset ids are UUIDs and need nothing, but go through the same door.

    :func:`_assert_scope` is the belt to this brace: if a scoped query ever comes back empty or
    carrying the wrong keys, the run says so loudly instead of writing a quietly short manifest.
    """
    return repr(str(value))


def scope_filter(dataset_id, cell_type=None, *, value_filter: str = cs.VALUE_FILTER) -> str:
    """The §1 filter, narrowed to one dataset (and optionally one cell type).

    The base filter is used verbatim and only ever has conjunctions appended: this driver narrows
    *what is read*, never *what qualifies*.
    """
    scope = f"{value_filter} and dataset_id == {_escape(dataset_id)}"
    if cell_type is not None:
        scope += f" and cell_type == {_escape(cell_type)}"
    return scope


def _assert_scope(obs, dataset_id, cell_type, *, expected_cells=None, scope="") -> None:
    """Fail loudly if a scoped read did not return what it was scoped to.

    Two failure modes this catches, both of which are otherwise silent: a filter clause that did
    not apply (the frame carries other datasets — the next stop would be an OOM, or worse, a
    manifest with rows the header does not describe), and a filter clause that matched nothing
    although pass 1 counted cells there (the escaping is wrong for this Census's parser).
    """
    for column, wanted in (("dataset_id", dataset_id), ("cell_type", cell_type)):
        if wanted is None or column not in obs.columns or not obs.shape[0]:
            continue
        seen = set(obs[column].astype(str).unique()) - {str(wanted)}
        if seen:
            raise ValueError(
                f"the value_filter did not scope the read: asked for {column} == {wanted!r}, got "
                f"{len(seen)} other value(s), e.g. {sorted(seen)[:3]}. Filter was: {scope}"
            )
    if expected_cells and not obs.shape[0]:
        raise ValueError(
            f"scoped read returned 0 rows, but pass 1 counted {expected_cells:,} cells for "
            f"dataset {dataset_id!r}"
            + (f" / cell_type {cell_type!r}" if cell_type is not None else "")
            + f". The filter matched nothing — check the literal escaping. Filter was: {scope}"
        )


def requery_dataset(census, dataset_id, *, cell_type=None, expected_cells=None,
                    value_filter: str = cs.VALUE_FILTER):
    """One scoped ``query_obs`` + ``screen_strata``, both unmodified. Returns
    ``(rows, attrs, n_cells)``.

    ``query_obs`` brings the optional columns along (``raw_sum`` for the depth bin and the median
    counts, the ontology term id, any pool/library column the schema exposes), so the gate and the
    confound pre-screen see exactly what they would see in a single-dataset run. The frame itself
    is not returned — it is the large object here, and it is dropped as soon as the rows are out;
    only its ``.attrs`` (provenance for the manifest header) and its row count survive.
    """
    scope = scope_filter(dataset_id, cell_type, value_filter=value_filter)
    obs = cs.query_obs(census, value_filter=scope)
    _assert_scope(obs, dataset_id, cell_type, expected_cells=expected_cells, scope=scope)
    rows = cs.screen_strata(obs)
    return rows, dict(obs.attrs), int(obs.shape[0])


def _header_carrier(attrs) -> pd.DataFrame:
    """An empty frame carrying pass-2 provenance, with the base §1 ``value_filter`` restored.

    :func:`pbcheck.census_select.emit_manifest` reads the header's ``value_filter`` off
    ``obs.attrs``, and every pass-2 frame's attrs carry the *scoped* filter (``... and dataset_id
    == '...'``) — true of that one query, false of the run. The header must state the filter the
    candidate set was drawn under, which is §1's; the scoping is recorded in ``notes`` instead.
    ``emit_manifest`` reads nothing but ``.attrs``, so an empty frame carries it exactly.
    """
    carrier = pd.DataFrame()
    carrier.attrs.update(dict(attrs))
    carrier.attrs["value_filter"] = cs.VALUE_FILTER
    return carrier


class _CensusHandle:
    """The open Census plus the one thing a long network job needs from it: reopening.

    A SOMA reader is stateful and a connection that dropped mid-read cannot be resumed, so every
    retry in this driver starts from a fresh handle rather than reusing one whose state after the
    error is unknown.
    """

    def __init__(self, *, tiledb_config=None):
        self.tiledb_config = tiledb_config
        self.handle = cs.open_census(tiledb_config=tiledb_config)

    def reopen(self) -> None:
        self.close()
        self.handle = cs.open_census(tiledb_config=self.tiledb_config)

    def close(self) -> None:
        closer = getattr(self.handle, "close", None)
        if closer is None:
            return
        try:
            closer()
        except Exception as exc:  # closing a dead handle must not mask the real failure
            print(f"[census-candidates] closing the census handle failed: {exc!r}",
                  file=sys.stderr, flush=True)


def _is_retryable(exc: BaseException) -> bool:
    return not isinstance(exc, FATAL_EXCEPTIONS)


def _retry_note(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {exc}"


def _requery_with_retries(handle, dataset_id, *, cell_type, expected_cells,
                          attempts: int = PASS2_ATTEMPTS):
    """:func:`requery_dataset` with backoff. Raises the last exception if every attempt fails."""
    for attempt in range(1, max(attempts, 1) + 1):
        try:
            return requery_dataset(handle.handle, dataset_id, cell_type=cell_type,
                                   expected_cells=expected_cells)
        except Exception as exc:
            if attempt >= attempts or not _is_retryable(exc):
                raise
            wait = RETRY_BACKOFF_SECONDS[min(attempt - 1, len(RETRY_BACKOFF_SECONDS) - 1)]
            scope = dataset_id if cell_type is None else f"{dataset_id} / {cell_type}"
            print(f"[pass2] {scope}: {_retry_note(exc)} — retry {attempt}/{attempts - 1} "
                  f"in {wait:.0f}s", file=sys.stderr, flush=True)
            time.sleep(wait)
            handle.reopen()
    raise AssertionError("unreachable: attempts must be >= 1")


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

def _run_pass1(handle, *, label="pass1"):
    """Pass 1 with one full retry. Returns ``(counts, n_cells_seen, n_batches, seconds)``.

    The retry restarts the *whole* pass: the reader is stateful, so a chunk that failed cannot be
    re-requested on its own, and half a stream is not half an answer. A second failure is raised —
    a driver that quietly writes a manifest from a truncated read would be worse than one that
    dies.
    """
    for attempt in (1, 2):
        t0 = time.time()
        try:
            counts, n_cells_seen, n_batches = aggregate_batches(
                iter_obs_batches(handle.handle), label=label
            )
            return counts, n_cells_seen, n_batches, time.time() - t0
        except Exception as exc:
            if attempt == 2 or not _is_retryable(exc):
                raise
            print(f"[{label}] FAILED after {(time.time() - t0) / 60:.1f} min: {_retry_note(exc)} "
                  "— restarting the pass from scratch (a stateful reader cannot resume a chunk); "
                  "the next failure is fatal", file=sys.stderr, flush=True)
            handle.reopen()
    raise AssertionError("unreachable")


def run(args) -> dict:
    """Both passes, then the manifest. Returns the manifest dict ``emit_manifest`` built."""
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[census-candidates] census_version={cs.CENSUS_VERSION} "
          f"organism={cs.CENSUS_ORGANISM}", flush=True)
    print(f"[census-candidates] value_filter (applied to obs): {cs.VALUE_FILTER}", flush=True)
    print(f"[census-candidates] value_filter (§1 as written): {cs.SPEC_VALUE_FILTER}", flush=True)
    print(f"[census-candidates] organism clause realised by: {cs.ORGANISM_REALIZED_BY}", flush=True)
    print(f"[census-candidates] thresholds: min_cells_per_donor={gc.MIN_CELLS} "
          f"min_donors_per_group={cs.MIN_DONORS_PER_GROUP} (frozen; no CLI option moves them)",
          flush=True)
    print(f"[census-candidates] out_dir={out_dir} dry_run={bool(args.dry_run)} "
          f"max_datasets={args.max_datasets} "
          f"max_cells_per_dataset_query={args.max_cells_per_dataset_query:,}", flush=True)

    tiledb_config = reader_tiledb_config(int(args.reader_buffer_mb) * 1024**2)
    print(f"[census-candidates] reader buffers: {args.reader_buffer_mb} MiB per column "
          f"({sorted(tiledb_config)}) — cellxgene-census defaults to 1 GiB each, which four "
          "string columns do not fit into a 16 GB runner", flush=True)
    print(f"[census-candidates] rss at start: {_fmt_bytes(_rss_bytes())}", flush=True)

    t_start = time.time()
    handle = _CensusHandle(tiledb_config=tiledb_config)
    try:
        counts, n_cells_seen, n_batches, pass1_seconds = _run_pass1(handle)
        if not n_cells_seen:
            raise RuntimeError(
                "pass 1 read 0 cells from the pinned Census — the filter, the organism key or the "
                "handle is wrong; refusing to write a manifest that would say 'no candidates'"
            )

        survivors, dataset_total_cells = surviving_datasets(counts)
        # Largest first: deterministic given the pin, and it puts the queries whose timing matters
        # most at the front of the log, where a stalled run is visible early.
        ordered = sorted(survivors, key=lambda d: (-dataset_total_cells.get(d, 0), d))
        cap = min([n for n in (DRY_RUN_MAX_DATASETS if args.dry_run else 0, args.max_datasets)
                   if n > 0] or [0])
        if args.dry_run:
            selected = spread_sample(ordered, cap)
            selection_rule = f"--dry-run: {len(selected)} datasets spread across the size range"
        elif cap:
            selected = ordered[:cap]
            selection_rule = f"--max-datasets {cap}: the {len(selected)} largest surviving datasets"
        else:
            selected = ordered
            selection_rule = "every surviving dataset"

        print(f"[pass1] done: {n_batches} batches | {n_cells_seen:,} cells | {len(counts):,} keys "
              f"| {len(dataset_total_cells):,} datasets | {pass1_seconds / 60:.1f} min | "
              f"{n_cells_seen / max(pass1_seconds, 1e-9):,.0f} cells/s", flush=True)
        print(f"[pass1] {len(survivors):,} dataset(s) could hold a stratum clearing the gate; "
              f"{selection_rule}", flush=True)

        rows: list[dict] = []
        failed_datasets: dict[str, str] = {}
        split_datasets: dict[str, int] = {}
        header_attrs: dict = {}
        n_queries = 0
        n_pass2_cells = 0
        t_pass2 = time.time()

        for i, dataset_id in enumerate(selected, 1):
            total_cells = dataset_total_cells.get(dataset_id, 0)
            scopes = plan_queries(counts, dataset_id, total_cells=total_cells,
                                  max_cells_per_query=args.max_cells_per_dataset_query)
            per_cell_type = {}
            if scopes != [None]:
                per_cell_type = dataset_cell_type_cells(counts, dataset_id)
                split_datasets[dataset_id] = len(scopes)
                print(f"[pass2] {i}/{len(selected)} {dataset_id}: {total_cells:,} cells > cap "
                      f"{args.max_cells_per_dataset_query:,} — splitting into {len(scopes)} "
                      "per-cell_type queries", flush=True)
                oversized = {c: n for c, n in per_cell_type.items()
                             if n > args.max_cells_per_dataset_query}
                for cell_type, n_cells in sorted(oversized.items(), key=lambda kv: -kv[1]):
                    print(f"[pass2] {dataset_id} / {cell_type}: {n_cells:,} cells is over the cap "
                          "on its own — querying anyway (a stratum is never split below "
                          "dataset x cell_type)", flush=True)

            t_ds = time.time()
            ds_rows: list[dict] = []
            ds_cells = 0
            for cell_type in scopes:
                expected = per_cell_type.get(cell_type, total_cells) if cell_type else total_cells
                try:
                    q_rows, attrs, n_cells = _requery_with_retries(
                        handle, dataset_id, cell_type=cell_type, expected_cells=expected
                    )
                except Exception as exc:
                    scope = dataset_id if cell_type is None else f"{dataset_id} / {cell_type}"
                    failed_datasets[dataset_id] = f"{scope}: {_retry_note(exc)}"
                    print(f"[pass2] {i}/{len(selected)} FAILED {scope}: {_retry_note(exc)} — the "
                          "whole dataset is dropped from the manifest (a partial dataset would "
                          "corrupt D4's excluded fractions) and recorded in notes.failed_datasets",
                          file=sys.stderr, flush=True)
                    ds_rows = []
                    break
                ds_rows.extend(q_rows)
                ds_cells += n_cells
                n_queries += 1
                header_attrs = header_attrs or attrs
            else:
                rows.extend(ds_rows)
                n_pass2_cells += ds_cells
                elapsed = max(time.time() - t_ds, 1e-9)
                print(f"[pass2] {i}/{len(selected)} {dataset_id} | {ds_cells:,} cells | "
                      f"{len(ds_rows)} contrast(s) | {elapsed:.1f}s | "
                      f"{ds_cells / elapsed:,.0f} cells/s", flush=True)

        pass2_seconds = time.time() - t_pass2
        notes = {
            "driver": "scripts/census_candidates.py",
            "strategy":
                "two passes over the pinned Census. Pass 1 streams the four stratum columns "
                "(no .concat()) and folds them into per-(dataset, cell_type, disease, donor) cell "
                "counts; pass 2 re-reads each surviving dataset with a scoped value_filter and "
                "runs census_select.screen_strata on the PER-CELL frame, unmodified. The pass-1 "
                "aggregate decides only WHICH datasets to read: the inclusion gate counts rows, "
                "so a donor-level aggregate would report every donor as holding one cell.",
            "coverage":
                "CANDIDATES FROM THE SELECTED DATASETS ONLY. The counts block below (D4's excluded "
                "fraction and its reasons) is computed over the contrasts of the datasets pass 2 "
                "actually read, not over the whole Census: datasets that cannot hold a "
                "gate-clearing stratum are never queried, so their strata appear in no denominator "
                "here.",
            "pass1": {
                "columns_streamed": list(KEY_COLUMNS),
                "n_batches": n_batches,
                "n_cells_seen": n_cells_seen,
                "n_keys": len(counts),
                "n_datasets_seen": len(dataset_total_cells),
                "n_datasets_surviving": len(survivors),
                "seconds": round(pass1_seconds, 1),
                "cells_per_second": round(n_cells_seen / max(pass1_seconds, 1e-9), 1),
                "reader_tiledb_config": dict(tiledb_config),
                "aggregation":
                    "Arrow-native: each batch is grouped with pyarrow (dictionary columns decoded "
                    "to their values first, since the codes are per-batch) and only the grouped "
                    "keys are ever turned into Python strings. Nothing calls .to_pandas() on a "
                    "batch and nothing calls .concat() on the reader.",
                "peak_rss_bytes": _peak_rss_bytes(),
                "coarse_filter":
                    f"a dataset survives when some (dataset_id, cell_type) has a disease term with "
                    f">= {cs.MIN_DONORS_PER_GROUP} donors of >= {gc.MIN_CELLS} cells and as many "
                    f"such '{cs.REFERENCE_LEVEL}' donors (§1 inclusion-gate items 2 and 1). "
                    "Conservative: item 3 (donor nested within condition) is left to pass 2, so "
                    "the filter may admit a dataset the gate then rejects, never the reverse.",
            },
            "pass2": {
                "selection_rule": selection_rule,
                "n_datasets_selected": len(selected),
                "n_datasets_failed": len(failed_datasets),
                "n_queries": n_queries,
                "n_cells_read": n_pass2_cells,
                "seconds": round(pass2_seconds, 1),
                "max_cells_per_dataset_query": int(args.max_cells_per_dataset_query),
                "datasets_split_per_cell_type": split_datasets,
                "value_filter_scoping":
                    "the header's value_filter is the executable §1 filter, verbatim (the "
                    "pre-registered text is beside it as value_filter_spec; §1's organism clause "
                    "is realised by the experiment key, see organism_realized_by). Each pass-2 "
                    "read appended "
                    "\"and dataset_id == '<id>'\" (and, for a split dataset, "
                    "\"and cell_type == '<label>'\") to narrow WHAT IS READ; nothing about what "
                    "qualifies was changed.",
                "split_caveat":
                    "in a per-cell_type split, a cell type present only on cells whose donor_id is "
                    "null is not queried (pass 1 drops null-donor rows before counting). Such a "
                    "stratum would have failed the gate on 'donor_id null on all cells'; it is "
                    "absent from this manifest rather than present as an excluded row."
                    if split_datasets else "no dataset needed the per-cell_type split",
            },
            "failed_datasets": failed_datasets,
            "dry_run": bool(args.dry_run),
            "total_seconds": round(time.time() - t_start, 1),
        }
        manifest = cs.emit_manifest(rows, out_dir=out_dir, stem=MANIFEST_STEM,
                                    obs=_header_carrier(header_attrs), notes=notes)
    finally:
        handle.close()

    print(f"[census-candidates] manifest -> {out_dir / (MANIFEST_STEM + '.json')} "
          f"(+ .csv), {len(manifest['rows'])} row(s), "
          f"{(time.time() - t_start) / 60:.1f} min total", flush=True)
    return manifest


def print_summary(manifest) -> None:
    """The run summary §1's D4 asks about, plus the largest candidates, on stdout."""
    header = manifest["header"]
    counts = header["counts"]
    notes = header.get("notes") or {}

    print("", flush=True)
    print("=== candidate manifest ===", flush=True)
    print(f"census_version         : {header['census_version']}", flush=True)
    print(f"value_filter           : {header['value_filter']}", flush=True)
    print(f"n_contrasts            : {counts['n_contrasts']}", flush=True)
    for status, n in sorted(counts["by_gate_status"].items()):
        print(f"  {status:<26}: {n}", flush=True)
    print(f"excluded_fraction      : {counts['excluded_fraction']}", flush=True)
    print(f"excluded_confound_frac : {counts['excluded_confound_fraction']}", flush=True)
    print(f"pooled / unresolved    : {counts['n_pooled']} / {counts['n_pooling_unresolved']}",
          flush=True)

    reasons = sorted(counts["inclusion_gate_failure_reasons"].items(), key=lambda kv: -kv[1])
    if reasons:
        print("top inclusion-gate failure reasons:", flush=True)
        for reason, n in reasons[:5]:
            print(f"  {n:>6}  {reason}", flush=True)

    failed = notes.get("failed_datasets") or {}
    if failed:
        print(f"FAILED datasets (absent from the rows above): {len(failed)}", flush=True)
        for dataset_id, why in sorted(failed.items()):
            print(f"  {dataset_id}: {why}", flush=True)

    candidates = [r for r in manifest["rows"] if r["gate_status"] == "candidate"]
    top = sorted(candidates, key=lambda r: -int(r["n_cells"]))[:TOP_N_SUMMARY]
    print(f"top {len(top)} candidate(s) by n_cells (of {len(candidates)}):", flush=True)
    for row in top:
        print(f"  {int(row['n_cells']):>9,}  {row['n_donors_A']}v{row['n_donors_B']} donors  "
              f"{row['dataset_id']}  {row['cell_type']} | {row['disease']} vs {row['reference']}",
              flush=True)
    print("admitted_to_sweep is False on every row: this is a candidate list, not the §1 "
          "pre-registration.", flush=True)


def build_arg_parser() -> argparse.ArgumentParser:
    """The four knobs, and deliberately no fifth.

    Nothing here can move a pre-registered number: the Census version, the ``value_filter`` and the
    §1 thresholds are read from :mod:`pbcheck.census_select` / :mod:`pbcheck.gate_config` and are
    not parameters of this driver. Changing one is a protocol change (``docs/AMENDMENTS.md``), not
    a command-line flag — pinned by ``tests/test_census_candidates.py``.
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--dry-run", action="store_true",
        help=f"read only {DRY_RUN_MAX_DATASETS} datasets, spread across the size range, to "
             "calibrate the timing of the full run. Pass 1 still streams the whole Census — that "
             "is the part whose cost is unknown.",
    )
    parser.add_argument(
        "--max-datasets", type=int, default=0,
        help="cap the number of datasets pass 2 reads (the largest first); 0 = no cap",
    )
    parser.add_argument(
        "--max-cells-per-dataset-query", type=int, default=2_000_000,
        help="a dataset with more cells than this is read one cell type at a time instead of in "
             "one query (never skipped); 0 = never split",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUT_DIR,
        help="where census_candidates.json/.csv are written (gitignored; the CI artifact)",
    )
    parser.add_argument(
        "--reader-buffer-mb", type=int, default=READER_BUFFER_BYTES // 1024**2,
        help="per-column read buffer handed to the Census as py./soma.init_buffer_bytes. A memory "
             "budget, not a data knob: it decides how many rows arrive per batch and nothing "
             "about which rows exist. cellxgene-census's own default is 1024 (1 GiB per column), "
             "which is what killed the second dry run on a 16 GB runner.",
    )
    return parser


def exit_code(manifest) -> int:
    """0 for a run that produced a candidate list, 1 for one that only looks like it did.

    Two failures are indistinguishable from success *in the artifact* and must not be
    indistinguishable in CI, because both write a well-formed manifest:

    * **nothing survived pass 1.** Zero datasets out of the whole Census clearing a filter that is
      only §1 items 2 and 1 is not an empty result, it is a regression — a moved threshold, a
      renamed stratum column, a :data:`pbcheck.census_select.REFERENCE_LEVEL` that no longer
      matches the Census's ``normal``. The manifest would be valid, nearly empty, and green.
    * **every selected dataset failed in pass 2.** A transport or scoping problem that hits all of
      them, reported as "0 candidates".

    One flaky dataset out of sixty is deliberately *not* one of these: that run still carries its
    result, the dataset is named in ``notes.failed_datasets``, and failing the job would discard
    five hours of work over a dropped connection.
    """
    notes = manifest["header"].get("notes") or {}
    pass1 = notes.get("pass1") or {}
    pass2 = notes.get("pass2") or {}
    failed = notes.get("failed_datasets") or {}
    selected = int(pass2.get("n_datasets_selected", 0))
    n_datasets_seen = int(pass1.get("n_datasets_seen", 0))
    n_cells_seen = int(pass1.get("n_cells_seen", 0))

    if not selected and (n_datasets_seen or n_cells_seen):
        print(
            "[census-candidates] 0 datasets survived the coarse filter across the whole Census "
            f"({n_datasets_seen:,} dataset(s) and {n_cells_seen:,} cells read in pass 1) — a "
            "regression, not an empty result. The §1 thresholds, the stratum column names and "
            f"REFERENCE_LEVEL ({cs.REFERENCE_LEVEL!r}) are where to look. The manifest is "
            "well-formed and nearly empty, which is exactly why this run exits non-zero.",
            file=sys.stderr, flush=True,
        )
        return 1
    if selected and len(failed) >= selected:
        print(f"[census-candidates] every one of the {selected} selected dataset(s) failed — "
              "this is a broken run, not an empty result", file=sys.stderr, flush=True)
        return 1
    return 0


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)
    _install_fatal_signal_logging()
    try:
        manifest = run(args)
        print_summary(manifest)
        return exit_code(manifest)
    finally:
        # Always, on every path including the fatal ones: the high-water mark is the number that
        # says whether the next run can afford a bigger buffer or needs a smaller one.
        print(f"[census-candidates] peak rss: {_fmt_bytes(_peak_rss_bytes())}",
              file=sys.stderr, flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
