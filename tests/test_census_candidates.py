"""The Census candidate driver, tested on a machine that has no Census.

``scripts/census_candidates.py`` is the only place in this repo that reads all ~60 M cells of the
pinned Census, and its two most dangerous failure modes are both **silent**:

* the coarse pass-1 filter drops a dataset that would have cleared the §1 inclusion gate — the
  stratum is simply absent from the manifest, with no row, no reason and no count;
* a scoped ``value_filter`` is escaped the way SQL would escape it, matches nothing, and the query
  returns zero rows without an error (TileDB-SOMA parses value filters with ``ast.parse``, so
  ``'O''Brien'`` is Python's implicit string concatenation — ``"OBrien"`` — not an escaped quote).

Neither would show up as a failure in CI, so both are pinned here. The conservativeness of the
coarse filter is not asserted against a hand-written expectation but **cross-checked against the
real** :func:`pbcheck.census_select.apply_inclusion_gate` on the same synthetic per-cell rows: for
every scenario, every dataset the real gate passes must survive the filter.

No network and no ``cellxgene-census`` anywhere in this file — including a test that the driver
itself does not import it at module level, which is what keeps this suite runnable at all.
"""

from __future__ import annotations

import ast
import importlib.util
import inspect
import json
import sys
import types
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from pbcheck import census_select as cs  # noqa: E402
from pbcheck import gate_config as gc  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DRIVER = REPO / "scripts" / "census_candidates.py"


def _load_driver():
    """Import the driver as a module (``scripts/`` is not on the import path)."""
    if not DRIVER.exists():  # pragma: no cover
        pytest.skip("census candidate driver not present")
    spec = importlib.util.spec_from_file_location("_census_candidates", DRIVER)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


cc = _load_driver()


# ---------------------------------------------------------------------------
# Synthetic obs frames — same shape as tests/test_census_select.py's, per CELL.
# ---------------------------------------------------------------------------

DEFAULT_CELL = {
    "assay": "10x 3' v3",
    "suspension_type": "cell",
    "tissue_general": "blood",
    "self_reported_ethnicity": "European",
    "sex": "female",
    "development_stage": "adult",
    "cell_type_ontology_term_id": "CL:0000084",
}


def obs_frame(donors, *, dataset_id="ds1", cell_type="T cell", **constants) -> pd.DataFrame:
    """A per-cell obs frame from a per-donor spec (``donor_id``, ``disease``, ``n_cells``)."""
    base = {**DEFAULT_CELL, **constants}
    rows = []
    for spec in donors:
        spec = dict(spec)
        n_cells = int(spec.pop("n_cells", 20))
        row = {**base, "dataset_id": dataset_id, "cell_type": cell_type, **spec}
        rows.extend(dict(row) for _ in range(n_cells))
    return pd.DataFrame(rows)


def healthy_donors(n_per_group=3, n_cells=20, disease="COVID-19", **extra):
    donors = [{"donor_id": f"A{i}", "disease": disease, "n_cells": n_cells, **extra}
              for i in range(n_per_group)]
    donors += [{"donor_id": f"B{i}", "disease": "normal", "n_cells": n_cells, **extra}
               for i in range(n_per_group)]
    return donors


def batches(frame, size):
    """Cut a frame into fixed-size chunks, the way a SOMA reader hands back Arrow batches."""
    return [frame.iloc[i:i + size].copy() for i in range(0, frame.shape[0], size)]


def counts_of(frame, size=17):
    counts, _seen, _n = cc.aggregate_batches(batches(frame, size), progress_every=0)
    return counts


# ---------------------------------------------------------------------------
# Pass 1: the aggregation must not depend on how the Census cut its batches.
# ---------------------------------------------------------------------------

def test_the_streamed_columns_are_the_stratum_columns_census_select_requires():
    assert set(cc.KEY_COLUMNS) == set(cs.REQUIRED_OBS_COLUMNS)
    assert inspect.signature(cc.iter_obs_batches).parameters["value_filter"].default \
        == cs.VALUE_FILTER
    assert inspect.signature(cc.iter_obs_batches).parameters["organism"].default \
        == cs.CENSUS_ORGANISM


def test_aggregation_does_not_depend_on_batch_boundaries_or_batch_order():
    """Batch boundaries are not pinned by the Census version; the data is. A result that moved
    with the chunking would be irreproducible for a reason nobody would think to look for."""
    frame = pd.concat([
        obs_frame(healthy_donors(n_per_group=4, n_cells=23), dataset_id="ds1"),
        obs_frame(healthy_donors(n_per_group=3, n_cells=11), dataset_id="ds1",
                  cell_type="B cell"),
        obs_frame(healthy_donors(n_per_group=5, n_cells=7), dataset_id="ds2"),
    ], ignore_index=True)

    fine, seen_fine, n_fine = cc.aggregate_batches(batches(frame, 7), progress_every=0)
    coarse, seen_coarse, n_coarse = cc.aggregate_batches(batches(frame, 1000), progress_every=0)

    shuffled = frame.sample(frac=1.0, random_state=0).reset_index(drop=True)
    chunks = batches(shuffled, 13)
    chunks.reverse()
    mixed, seen_mixed, _ = cc.aggregate_batches(chunks, progress_every=0)

    assert fine == coarse == mixed
    assert seen_fine == seen_coarse == seen_mixed == frame.shape[0]
    assert n_fine > n_coarse >= 1
    assert sum(fine.values()) == frame.shape[0]
    assert fine[("ds1", "T cell", "COVID-19", "A0")] == 23


def test_cells_without_a_donor_or_a_stratum_are_dropped_before_grouping():
    """§1 item 3 requires ``donor_id`` present and ``build_contrasts`` forms no stratum from a
    null — and under pandas 3.0 a null survives ``.astype(str)`` instead of becoming ``"nan"``,
    so grouping first would invent a donor called "nan" holding thousands of cells."""
    frame = pd.concat([
        obs_frame(healthy_donors()),
        obs_frame([{"donor_id": None, "disease": "COVID-19", "n_cells": 15}]),
        obs_frame([{"donor_id": "d9", "disease": None, "n_cells": 5}]),
    ], ignore_index=True)

    counts, seen, _ = cc.aggregate_batches(batches(frame, 9), progress_every=0)

    assert seen == frame.shape[0], "every row read is counted, including the dropped ones"
    assert sum(counts.values()) == frame.shape[0] - 20
    assert not any("nan" in key or "None" in key for key in counts)
    assert all(all(part is not None for part in key) for key in counts)


# ---------------------------------------------------------------------------
# The cross-check: the coarse filter must never subtract what the real gate would pass.
# ---------------------------------------------------------------------------

def _two_datasets():
    return pd.concat([
        obs_frame(healthy_donors(), dataset_id="ds1"),
        obs_frame(healthy_donors(n_cells=4), dataset_id="ds2"),
    ], ignore_index=True)


def _donor_in_both_arms():
    donors = healthy_donors()
    donors.append({"donor_id": "A0", "disease": "normal", "n_cells": 20})
    return obs_frame(donors)


def _two_disease_terms():
    donors = healthy_donors()
    donors += [{"donor_id": f"C{i}", "disease": "influenza", "n_cells": 20} for i in range(2)]
    return obs_frame(donors)


SCENARIOS = {
    # name: (frame, datasets the gate can pass, datasets the coarse filter must reject)
    "clean 3v3": (obs_frame(healthy_donors()), {"ds1"}, set()),
    "exactly at the thin-donor threshold": (
        obs_frame(healthy_donors(n_cells=gc.MIN_CELLS)), {"ds1"}, set()),
    "one cell below it": (obs_frame(healthy_donors(n_cells=gc.MIN_CELLS - 1)), set(), {"ds1"}),
    "two donors in group A": (
        obs_frame([{"donor_id": f"A{i}", "disease": "COVID-19"} for i in range(2)]
                  + [{"donor_id": f"B{i}", "disease": "normal"} for i in range(4)]),
        set(), {"ds1"}),
    "no normal arm": (
        obs_frame([{"donor_id": f"A{i}", "disease": "COVID-19"} for i in range(4)]),
        set(), {"ds1"}),
    "a donor with cells under both conditions": (_donor_in_both_arms(), set(), set()),
    "a second disease term that is too thin": (_two_disease_terms(), {"ds1"}, set()),
    "null donors alongside a usable design": (
        pd.concat([obs_frame(healthy_donors()),
                   obs_frame([{"donor_id": None, "disease": "normal", "n_cells": 40}])],
                  ignore_index=True), {"ds1"}, set()),
    "two datasets, one usable": (_two_datasets(), {"ds1"}, {"ds2"}),
}


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_the_coarse_filter_never_drops_a_dataset_the_real_gate_would_pass(name):
    """The one property pass 1 owes pass 2, checked against the real gate, not against a guess.

    ``surviving_datasets`` decides which datasets are read at all, so a dataset it drops is gone
    from the manifest without a row, a reason or a count — the failure D4 exists to make
    impossible. It may over-admit (item 3's nesting rule is deliberately left to pass 2, which is
    why "a donor with cells under both conditions" expects nothing in either direction here).
    """
    frame, must_survive, must_not_survive = SCENARIOS[name]
    survivors, dataset_cells = cc.surviving_datasets(counts_of(frame))

    passed_the_real_gate = {
        contrast.dataset_id for contrast in cs.build_contrasts(frame)
        if cs.apply_inclusion_gate(frame, contrast).passed
    }
    assert passed_the_real_gate <= survivors, (
        f"{name}: the coarse filter dropped {passed_the_real_gate - survivors}, which the §1 "
        "inclusion gate passes"
    )
    assert must_survive <= survivors
    assert not (must_not_survive & survivors), f"{name}: the filter is admitting the unusable"
    for dataset_id in set(frame["dataset_id"]):
        assert dataset_cells[dataset_id] > 0


def test_the_coarse_filter_reads_its_thresholds_from_the_frozen_config():
    """Not restated here, and not reachable from the command line: they are pre-registered."""
    params = inspect.signature(cc.surviving_datasets).parameters
    assert params["min_cells_per_donor"].default == gc.MIN_CELLS == 10
    assert params["min_donors_per_group"].default == cs.MIN_DONORS_PER_GROUP == 3
    assert params["reference"].default == cs.REFERENCE_LEVEL == "normal"


def test_no_command_line_option_can_move_a_pre_registered_number():
    """Changing a threshold, the value_filter or the Census pin is a protocol change
    (docs/AMENDMENTS.md), so the driver must not offer a flag that does it."""
    options = {opt for action in cc.build_arg_parser()._actions for opt in action.option_strings}
    assert options == {
        "-h", "--help", "--dry-run", "--max-datasets",
        "--max-cells-per-dataset-query", "--output-dir",
    }
    # ...and no call in the driver passes a census_version at all: open_census() defaults to the
    # §1 pin, and a driver that named the version could point the run at a different Census.
    tree = ast.parse(DRIVER.read_text(encoding="utf-8"))
    keywords = {kw.arg for node in ast.walk(tree) if isinstance(node, ast.Call)
                for kw in node.keywords}
    assert "census_version" not in keywords


def test_the_driver_feeds_the_gate_per_cell_rows_never_the_pass_one_aggregate():
    """Why pass 2 re-reads the cells at all: ``apply_inclusion_gate`` counts ROWS.

    Handing it the pass-1 aggregate (one row per donor) reports every donor as holding a single
    cell, so every donor falls to the thin-donor rule and the whole Census reads as failing the
    gate. This is the mistake the two-pass design exists to make impossible; it is pinned here so
    that "just use the counts we already have" is a red test rather than a plausible shortcut.
    """
    frame = obs_frame(healthy_donors())
    contrast = cs.build_contrasts(frame)[0]
    assert cs.apply_inclusion_gate(frame, contrast).passed

    donor_aggregate = pd.DataFrame([
        {"dataset_id": ds, "cell_type": ct, "disease": disease, "donor_id": donor,
         "n_cells": n_cells}
        for (ds, ct, disease, donor), n_cells in counts_of(frame).items()
    ])
    gate = cs.apply_inclusion_gate(donor_aggregate, contrast)
    assert not gate.passed
    assert any("cells-per-donor" in failure for failure in gate.failures)


# ---------------------------------------------------------------------------
# The scoped value_filter: escaping, and the belt-and-braces scope check.
# ---------------------------------------------------------------------------

def _literals(expression: str) -> set:
    """Every literal in a value filter, read with the parser TileDB-SOMA itself uses."""
    tree = ast.parse(expression, mode="eval")
    return {node.value for node in ast.walk(tree) if isinstance(node, ast.Constant)}


def test_a_scoped_filter_only_narrows_what_is_read():
    scope = cc.scope_filter("dataset-uuid", "T cell")
    assert scope.startswith(cs.VALUE_FILTER + " and ")
    assert scope.count(cs.VALUE_FILTER) == 1
    assert cc.scope_filter("dataset-uuid") == f"{cs.VALUE_FILTER} and dataset_id == 'dataset-uuid'"


def test_a_scoped_filter_survives_an_apostrophe_in_a_free_text_cell_type():
    """Cell-type labels are ontology free text. A label with an apostrophe must round-trip."""
    cell_type = "cell of Bowman's capsule"
    scope = cc.scope_filter("ds-1", cell_type)
    literals = _literals(scope)  # parses at all -> the expression is valid for the SOMA parser
    assert cell_type in literals and "ds-1" in literals


def test_the_sql_style_escape_would_have_been_silently_wrong():
    """Documents why ``_escape`` is ``repr`` and not the SQL doubling rule.

    TileDB-SOMA parses a value filter with ``ast.parse(expr, mode="eval")`` — Python source, not
    SQL. ``'O''Brien'`` is therefore implicit concatenation of two adjacent literals: it parses
    without error, means ``"OBrien"``, matches nothing, and the stratum vanishes from the manifest
    with nothing anywhere to say so.
    """
    assert _literals("cell_type == 'O''Brien'") == {"OBrien"}
    assert "O'Brien" in _literals(cc.scope_filter("ds-1", "O'Brien"))
    assert "\\" in cc._escape("back\\slash") and "back\\slash" in _literals(
        f"cell_type == {cc._escape('back\\slash')}"
    )


# ---------------------------------------------------------------------------
# Pass 2 against a SOMA stand-in that actually applies the filter (no network).
# ---------------------------------------------------------------------------

def _equality_clauses(expression: str):
    """The ``name == literal`` clauses of a value filter, via ``ast`` — as SOMA would read it."""
    tree = ast.parse(expression, mode="eval")
    clauses = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Compare) and len(node.ops) == 1
                and isinstance(node.ops[0], ast.Eq)
                and isinstance(node.left, ast.Name)
                and isinstance(node.comparators[0], ast.Constant)):
            clauses.append((node.left.id, node.comparators[0].value))
    return clauses


class _FakeBatch:
    def __init__(self, frame):
        self._frame = frame

    def to_pandas(self):
        return self._frame


class _FakeRead:
    """A SOMA reader: iterable in chunks (pass 1) *and* concat-able (pass 2, via query_obs)."""

    def __init__(self, frame, *, batch_size=16):
        self._frame = frame
        self._batch_size = batch_size

    def __iter__(self):
        return iter(_FakeBatch(chunk) for chunk in batches(self._frame, self._batch_size))

    def concat(self):
        return self

    def to_pandas(self):
        return self._frame


class _FakeSomaObs:
    """A SOMA obs stand-in. Unlike a mock, it *applies* the equality clauses of the filter, so a
    driver that mis-escaped a literal would come back empty here exactly as it would on S3."""

    def __init__(self, frame, *, honour_filter=True, fail_on=None):
        self._frame = frame
        self.honour_filter = honour_filter
        self.fail_on = fail_on  # substring of the value_filter whose read should blow up
        self.schema = types.SimpleNamespace(names=list(frame.columns))
        self.calls = []
        self.failures = 0

    def read(self, column_names=None, value_filter=None):
        self.calls.append({"column_names": list(column_names), "value_filter": value_filter})
        if self.fail_on and self.fail_on in (value_filter or ""):
            self.failures += 1
            raise RuntimeError(f"simulated transport failure for {self.fail_on!r}")
        frame = self._frame
        if self.honour_filter:
            for column, wanted in _equality_clauses(value_filter):
                if column in frame.columns:
                    frame = frame.loc[frame[column].astype(str) == str(wanted)]
        return _FakeRead(frame[[c for c in column_names if c in frame.columns]].copy())


class _FakeCensus:
    def __init__(self, soma_obs):
        self._soma_obs = soma_obs
        self.closed = 0

    def __getitem__(self, key):
        assert key == "census_data"
        return {cs.CENSUS_ORGANISM: types.SimpleNamespace(obs=self._soma_obs)}

    def close(self):
        self.closed += 1


def _two_dataset_frame():
    frame = pd.concat([
        obs_frame(healthy_donors(), dataset_id="ds1"),
        obs_frame(healthy_donors(), dataset_id="ds1", cell_type="B cell"),
        obs_frame(healthy_donors(), dataset_id="ds2"),
    ], ignore_index=True)
    frame["raw_sum"] = 2000.0
    return frame


def test_requery_dataset_scopes_the_read_and_screens_the_per_cell_frame():
    soma = _FakeSomaObs(_two_dataset_frame())
    rows, attrs, n_cells = cc.requery_dataset(_FakeCensus(soma), "ds1", expected_cells=240)

    assert soma.calls[0]["value_filter"].startswith(cs.VALUE_FILTER)
    assert n_cells == 240
    assert {row["dataset_id"] for row in rows} == {"ds1"}
    assert {row["cell_type"] for row in rows} == {"T cell", "B cell"}
    assert all(row["gate_status"] == "candidate" for row in rows)
    assert all(row["admitted_to_sweep"] is False for row in rows)
    # The attrs carry the SCOPED filter, which is why the manifest header cannot use them as is.
    assert attrs["value_filter"] != cs.VALUE_FILTER


def test_requery_dataset_can_be_scoped_to_one_cell_type():
    soma = _FakeSomaObs(_two_dataset_frame())
    rows, _attrs, n_cells = cc.requery_dataset(
        _FakeCensus(soma), "ds1", cell_type="B cell", expected_cells=120
    )
    assert n_cells == 120
    assert {(row["dataset_id"], row["cell_type"]) for row in rows} == {("ds1", "B cell")}


def test_requery_dataset_fails_loudly_when_the_filter_did_not_apply():
    soma = _FakeSomaObs(_two_dataset_frame(), honour_filter=False)
    with pytest.raises(ValueError, match="did not scope"):
        cc.requery_dataset(_FakeCensus(soma), "ds1")


def test_requery_dataset_fails_loudly_when_the_filter_matched_nothing():
    """The escaping tripwire: pass 1 counted cells there, so an empty read is a broken filter,
    not an empty dataset — and a quietly short manifest is the failure mode that hides it."""
    soma = _FakeSomaObs(_two_dataset_frame())
    with pytest.raises(ValueError, match="0 rows"):
        cc.requery_dataset(_FakeCensus(soma), "ds-absent", expected_cells=1234)


# ---------------------------------------------------------------------------
# The per-cell_type fallback for datasets over the query cap.
# ---------------------------------------------------------------------------

BIG_COUNTS = {
    ("ds1", "T cell", "COVID-19", "d1"): 1_500_000,
    ("ds1", "B cell", "normal", "d2"): 900_000,
    ("ds2", "T cell", "normal", "d3"): 10,
}


def test_a_dataset_over_the_cap_is_split_per_cell_type_never_skipped():
    assert cc.plan_queries(BIG_COUNTS, "ds1", total_cells=2_400_000,
                           max_cells_per_query=2_000_000) == ["B cell", "T cell"]
    assert cc.plan_queries(BIG_COUNTS, "ds1", total_cells=1_000_000,
                           max_cells_per_query=2_000_000) == [None]
    assert cc.plan_queries(BIG_COUNTS, "ds1", total_cells=2_400_000,
                           max_cells_per_query=0) == [None], "0 disables the split"
    # A single cell type over the cap is still queried: the stratum is the atom of the analysis.
    assert cc.plan_queries(BIG_COUNTS, "ds1", total_cells=2_400_000,
                           max_cells_per_query=1) == ["B cell", "T cell"]
    assert cc.dataset_cell_type_cells(BIG_COUNTS, "ds1") == {
        "T cell": 1_500_000, "B cell": 900_000}


def _row_keys(rows):
    return sorted(json.dumps(row, sort_keys=True, default=str) for row in rows)


def test_splitting_a_dataset_by_cell_type_does_not_change_a_single_manifest_row():
    """The property that makes the fallback safe rather than merely cheaper.

    A stratum is ``(dataset_id x cell_type)``, so a cut on ``cell_type`` partitions the dataset's
    strata exactly: no contrast spans two cell types, and the union of the per-cell-type manifests
    is the manifest of the whole dataset. (A cut on ``disease`` would sever a contrast from its
    ``normal`` arm; a cut on ``donor_id`` would sever a group from its donors. Neither is offered.)
    """
    frame = _two_dataset_frame()
    frame = frame.loc[frame["dataset_id"] == "ds1"].reset_index(drop=True)

    whole = cs.screen_strata(frame)
    split = []
    for cell_type in sorted(set(frame["cell_type"])):
        part = frame.loc[frame["cell_type"] == cell_type].reset_index(drop=True)
        split.extend(cs.screen_strata(part))

    assert len(whole) == 2
    assert _row_keys(whole) == _row_keys(split)


def _run_args(tmp_path, **overrides):
    args = dict(output_dir=tmp_path, dry_run=False, max_datasets=0,
                max_cells_per_dataset_query=2_000_000)
    args.update(overrides)
    return types.SimpleNamespace(**args)


def test_both_passes_end_to_end_produce_the_manifest_and_its_provenance(tmp_path, monkeypatch):
    """The whole driver against a SOMA stand-in: stream, aggregate, re-query, screen, emit."""
    frame = _two_dataset_frame()
    census = _FakeCensus(_FakeSomaObs(frame))
    monkeypatch.setattr(cs, "open_census", lambda **kwargs: census)

    manifest = cc.run(_run_args(tmp_path))
    cc.print_summary(manifest)  # must survive a real manifest

    assert (tmp_path / f"{cc.MANIFEST_STEM}.json").exists()
    assert (tmp_path / f"{cc.MANIFEST_STEM}.csv").exists()
    assert census.closed == 1

    header = manifest["header"]
    notes = header["notes"]
    assert header["value_filter"] == cs.VALUE_FILTER
    assert notes["failed_datasets"] == {}
    assert notes["pass1"]["n_cells_seen"] == frame.shape[0]
    assert notes["pass1"]["n_datasets_surviving"] == 2
    assert notes["pass2"]["n_queries"] == 2, "one query per dataset when nothing is over the cap"
    assert notes["pass2"]["datasets_split_per_cell_type"] == {}
    assert len(manifest["rows"]) == 3
    assert {row["dataset_id"] for row in manifest["rows"]} == {"ds1", "ds2"}
    assert all(row["admitted_to_sweep"] is False for row in manifest["rows"])


def test_the_split_path_end_to_end_produces_the_same_rows_as_the_single_query(
        tmp_path, monkeypatch):
    """A cap low enough to split every dataset must change the number of QUERIES, not the rows."""
    frame = _two_dataset_frame()
    monkeypatch.setattr(cs, "open_census", lambda **kwargs: _FakeCensus(_FakeSomaObs(frame)))

    whole = cc.run(_run_args(tmp_path / "whole"))
    split = cc.run(_run_args(tmp_path / "split", max_cells_per_dataset_query=1))

    assert split["header"]["notes"]["pass2"]["n_queries"] == 3
    assert split["header"]["notes"]["pass2"]["datasets_split_per_cell_type"] == {"ds1": 2, "ds2": 1}
    assert _row_keys(split["rows"]) == _row_keys(whole["rows"])


# ---------------------------------------------------------------------------
# The exit code: which broken runs must not be green. Both failures below write a well-formed
# manifest, so the artifact alone cannot tell them from a real result.
# ---------------------------------------------------------------------------

@pytest.fixture
def _no_backoff(monkeypatch):
    """Keep the retry path (three attempts per query) without its 2 s + 8 s waits."""
    monkeypatch.setattr(cc, "RETRY_BACKOFF_SECONDS", (0.0, 0.0))


def test_main_exits_zero_when_one_dataset_of_several_fails(tmp_path, monkeypatch, _no_backoff,
                                                           capsys):
    """A dropped connection on one dataset out of several must not discard the whole run — the
    dataset is named in notes.failed_datasets and the rest of the manifest still stands."""
    soma = _FakeSomaObs(_two_dataset_frame(), fail_on="'ds2'")
    monkeypatch.setattr(cs, "open_census", lambda **kwargs: _FakeCensus(soma))

    assert cc.main(["--output-dir", str(tmp_path)]) == 0
    assert soma.failures == cc.PASS2_ATTEMPTS, "the failing dataset was retried, then recorded"

    manifest = json.loads((tmp_path / f"{cc.MANIFEST_STEM}.json").read_text(encoding="utf-8"))
    failed = manifest["header"]["notes"]["failed_datasets"]
    assert set(failed) == {"ds2"}
    assert {row["dataset_id"] for row in manifest["rows"]} == {"ds1"}
    assert "FAILED" in capsys.readouterr().err


def test_main_exits_one_when_every_selected_dataset_fails(tmp_path, monkeypatch, _no_backoff,
                                                          capsys):
    soma = _FakeSomaObs(_two_dataset_frame(), fail_on="dataset_id ==")
    monkeypatch.setattr(cs, "open_census", lambda **kwargs: _FakeCensus(soma))

    assert cc.main(["--output-dir", str(tmp_path)]) == 1
    assert "every one of the 2 selected dataset(s) failed" in capsys.readouterr().err


def test_main_exits_one_when_the_coarse_filter_admits_no_dataset_at_all(tmp_path, monkeypatch,
                                                                        capsys):
    """The regression that would otherwise be a green CI run and a well-formed empty manifest.

    Zero survivors across a Census that pass 1 read tens of millions of cells from is not an empty
    result: it is a moved threshold, a renamed stratum column or a REFERENCE_LEVEL that no longer
    matches the Census's ``normal``. Nothing in the artifact distinguishes it from "no candidates
    exist", which is why the exit code has to.
    """
    thin = pd.concat([
        obs_frame(healthy_donors(n_cells=gc.MIN_CELLS - 1), dataset_id="ds1"),
        obs_frame(healthy_donors(n_cells=gc.MIN_CELLS - 1), dataset_id="ds2"),
    ], ignore_index=True)
    monkeypatch.setattr(cs, "open_census", lambda **kwargs: _FakeCensus(_FakeSomaObs(thin)))

    assert cc.main(["--output-dir", str(tmp_path)]) == 1

    err = capsys.readouterr().err
    assert "0 datasets survived the coarse filter" in err
    assert cs.REFERENCE_LEVEL in err, "the message names where to look"

    manifest = json.loads((tmp_path / f"{cc.MANIFEST_STEM}.json").read_text(encoding="utf-8"))
    assert manifest["rows"] == [], "well-formed and empty — the artifact cannot raise the alarm"
    assert manifest["header"]["notes"]["pass1"]["n_cells_seen"] == thin.shape[0]
    # ...and a run that legitimately selected nothing because it read nothing does not reach here:
    # run() refuses to write a manifest at all when pass 1 saw zero cells.
    assert cc.exit_code({"header": {"notes": {"pass1": {}, "pass2": {}}}}) == 0


def test_a_dry_run_samples_across_the_size_range_not_off_one_end():
    ordered = [f"ds{i}" for i in range(20)]  # largest first, as run() orders them
    picked = cc.spread_sample(ordered, 5)
    assert len(picked) == 5
    assert picked[0] == ordered[0] and picked[-1] == ordered[-1]
    assert cc.spread_sample(ordered, 0) == ordered
    assert cc.spread_sample(ordered[:3], 5) == ordered[:3]
    assert cc.spread_sample(ordered, 5) == picked, "deterministic: the Census is pinned"


# ---------------------------------------------------------------------------
# The manifest header, and the import boundary that keeps this suite runnable.
# ---------------------------------------------------------------------------

def test_the_header_states_the_spec_filter_not_the_per_dataset_scoping():
    """``emit_manifest`` takes its header's ``value_filter`` from ``obs.attrs``, and every pass-2
    frame carries the scoped one — true of that query, false of the run. A header that advertised
    ``... and dataset_id == '<one id>'`` would describe the candidate set incorrectly."""
    soma = _FakeSomaObs(_two_dataset_frame())
    rows, attrs, _n = cc.requery_dataset(_FakeCensus(soma), "ds1", expected_cells=240)

    manifest = cs.emit_manifest(rows, obs=cc._header_carrier(attrs))
    assert manifest["header"]["value_filter"] == cs.VALUE_FILTER
    assert "dataset_id ==" not in manifest["header"]["value_filter"]
    assert manifest["header"]["organism"] == cs.CENSUS_ORGANISM
    assert manifest["header"]["depth_column"] == "raw_sum", "pass-2 provenance is kept"


def test_the_driver_never_imports_cellxgene_census_at_module_level():
    """The extra is optional and absent on the development machine; the driver must import (and
    be testable) without it. Only ``census_select.open_census`` may reach for it, lazily."""
    sys.modules.pop("cellxgene_census", None)
    sys.modules.pop("_census_candidates", None)
    _load_driver()
    assert "cellxgene_census" not in sys.modules

    tree = ast.parse(DRIVER.read_text(encoding="utf-8"))
    top_level = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
    names = [a.name for n in top_level if isinstance(n, ast.Import) for a in n.names]
    names += [n.module for n in top_level if isinstance(n, ast.ImportFrom)]
    assert not any(str(name).startswith(("cellxgene", "tiledb")) for name in names)


def test_the_manifest_stem_says_candidates():
    """The file name is the first thing anyone reads, and this artifact is not the §1
    pre-registration of the stratum list."""
    assert cc.MANIFEST_STEM == "census_candidates"
    assert cc.DEFAULT_OUT_DIR.parent.name == "results", "gitignored: CI uploads it, git ignores it"
