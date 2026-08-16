"""The real-data stratum selector reads obs only, pins the Census, and admits nothing.

No network and no ``cellxgene-census`` anywhere in this file: every test builds the obs frame it
needs by hand, which is also the only way to exercise the pathological designs (a donor in both
conditions, an assay that perfectly separates condition) that real data supplies rarely and never on
demand. The Census access path itself is covered against a stub with the same SOMA surface.

The two things these tests are really guarding are the two that would be silent if wrong: that the
census version is a pin and not a moving alias, and that a column this module cannot compute is
written as ``pending`` rather than as a number.
"""

from __future__ import annotations

import ast
import csv
import json
import math
import sys
import types
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from pbcheck import census_select as cs  # noqa: E402
from pbcheck import design, gate_config as gc  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
SPEC = (REPO / "docs" / "PHASE0_SPEC.md").read_text(encoding="utf-8")
AMENDMENTS = (REPO / "docs" / "AMENDMENTS.md").read_text(encoding="utf-8")
SOURCE = (REPO / "src" / "pbcheck" / "census_select.py").read_text(encoding="utf-8")

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
    """An obs frame from a per-donor spec: ``donor_id``, ``disease``, ``n_cells``, overrides."""
    base = {**DEFAULT_CELL, **constants}
    rows = []
    for spec in donors:
        spec = dict(spec)
        n_cells = int(spec.pop("n_cells", 20))
        row = {**base, "dataset_id": dataset_id, "cell_type": cell_type, **spec}
        rows.extend(dict(row) for _ in range(n_cells))
    return pd.DataFrame(rows)


def healthy_donors(n_per_group=3, n_cells=20, **extra):
    donors = [{"donor_id": f"A{i}", "disease": "COVID-19", "n_cells": n_cells, **extra}
              for i in range(n_per_group)]
    donors += [{"donor_id": f"B{i}", "disease": "normal", "n_cells": n_cells, **extra}
               for i in range(n_per_group)]
    return donors


# ---------------------------------------------------------------------------
# The pin (spec §1) — the failure mode that would poison every downstream number silently.
# ---------------------------------------------------------------------------

def test_census_version_is_the_spec_pin():
    assert cs.CENSUS_VERSION == "2025-01-30"
    assert 'open_soma(census_version="2025-01-30")' in SPEC
    assert cs.SPEC_VALUE_FILTER in SPEC


def test_the_executable_filter_is_the_spec_filter_minus_the_unfilterable_organism_clause():
    """Regression from the first live run, which died on the reader constructor with
    ``SOMAError: 'Column organism does not exist in schema'`` before reading a cell.

    §1's filter text names three conditions; only two of them are obs columns. In the pinned
    Census's schema the organism is the **experiment** — ``census["census_data"]["homo_sapiens"]``,
    which ``query_obs`` has always opened — so the third clause is not a filter that is merely
    redundant, it is unparseable against that schema. The pre-registered text is kept verbatim for
    the manifest; what executes is the text minus that one conjunct, and the population selected is
    the same one either way.
    """
    assert "organism" not in cs.VALUE_FILTER, "the clause the live parser rejects"
    assert "Homo sapiens" not in cs.VALUE_FILTER
    assert cs.VALUE_FILTER == "is_primary_data == True and disease != 'na'"

    # The executable filter is the spec text with exactly one conjunct removed — nothing else was
    # relaxed, tightened or reworded while the file was open.
    dropped = "organism == 'Homo sapiens'"
    assert dropped in cs.SPEC_VALUE_FILTER
    assert " and ".join(
        part for part in cs.SPEC_VALUE_FILTER.split(" and ") if part != dropped
    ) == cs.VALUE_FILTER

    # ...and the dropped clause is realised, not abandoned.
    assert cs.CENSUS_ORGANISM == "homo_sapiens"
    assert "census['census_data']['homo_sapiens']" in cs.ORGANISM_REALIZED_BY


def test_the_mutable_aliases_are_rejected_not_used():
    """§1: "never ``latest``". The alias must be refused, and must appear nowhere as an argument."""
    for alias in cs.MUTABLE_CENSUS_ALIASES:
        with pytest.raises(ValueError, match="MUTABLE"):
            cs.open_census(census_version=alias)
    with pytest.raises(ValueError, match="AMENDMENTS"):
        cs.open_census(census_version="2024-07-01")

    tree = ast.parse(SOURCE)
    skip = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "MUTABLE_CENSUS_ALIASES" for t in node.targets
        ):
            skip.update(id(c) for c in ast.walk(node.value))
    literals = {
        node.value for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in skip
    }
    assert not literals & set(cs.MUTABLE_CENSUS_ALIASES), (
        "a mutable census alias appears as a string literal outside MUTABLE_CENSUS_ALIASES"
    )


def test_the_module_imports_without_the_census_extra():
    """``cellxgene-census`` is an optional extra; importing this module must not need it."""
    tree = ast.parse(SOURCE)
    top_level = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
    names = [a.name for n in top_level if isinstance(n, ast.Import) for a in n.names]
    names += [n.module for n in top_level if isinstance(n, ast.ImportFrom)]
    assert not any(str(n).startswith("cellxgene") for n in names)


def test_open_census_without_the_extra_says_how_to_install_it(monkeypatch):
    monkeypatch.setitem(sys.modules, "cellxgene_census", None)
    with pytest.raises(cs.CensusNotInstalled, match=r"\[census\]"):
        cs.open_census()


# ---------------------------------------------------------------------------
# query_obs against a stub with the SOMA surface (no network).
# ---------------------------------------------------------------------------

class _FakeRead:
    def __init__(self, frame):
        self._frame = frame

    def concat(self):
        return self

    def to_pandas(self):
        return self._frame


#: Columns the real Census obs schema carries that these synthetic frames do not materialise.
#: ``is_primary_data`` is one — the §1 filter tests it, and no test here needs its values.
#: ``organism`` is deliberately **not** one: in the pinned schema the organism is the experiment
#: (``census["census_data"]["homo_sapiens"]``), not a column, and a stand-in that tolerated an
#: ``organism`` clause is exactly what let a filter the live parser rejects reach production.
CENSUS_ONLY_FILTER_COLUMNS = ("is_primary_data",)


class _FakeSomaError(RuntimeError):
    """Stands in for ``tiledbsoma.SOMAError``, which is what the live reader raises."""


def _check_filter_against_schema(value_filter, schema_names) -> None:
    """Resolve every name in a value filter against the schema, as the live parser does.

    TileDB-SOMA parses the filter with ``ast`` and then resolves each name against the array's
    schema, raising ``'Column X does not exist in schema'`` when one is not there — at reader
    construction, before a single cell is read. A stand-in that skipped this step accepted filters
    the Census refuses, which is not a stand-in but a rubber stamp.
    """
    known = set(schema_names)
    for node in ast.walk(ast.parse(value_filter or "True", mode="eval")):
        if isinstance(node, ast.Name) and node.id not in known:
            raise _FakeSomaError(f"Column {node.id} does not exist in schema")


class _FakeSomaObs:
    def __init__(self, frame, *, schema_names=None):
        self._frame = frame
        self.schema = types.SimpleNamespace(
            names=list(schema_names or [*frame.columns, *CENSUS_ONLY_FILTER_COLUMNS])
        )
        self.calls = []

    def read(self, column_names=None, value_filter=None):
        self.calls.append({"column_names": list(column_names), "value_filter": value_filter})
        _check_filter_against_schema(value_filter, self.schema.names)
        return _FakeRead(self._frame[[c for c in column_names if c in self._frame.columns]].copy())


class _FakeCensus:
    def __init__(self, soma_obs):
        self._soma_obs = soma_obs

    def __getitem__(self, key):
        assert key == "census_data"
        return {cs.CENSUS_ORGANISM: types.SimpleNamespace(obs=self._soma_obs)}


def test_query_obs_asks_for_the_spec_columns_and_records_what_is_missing():
    frame = obs_frame(healthy_donors())
    frame = frame.drop(columns=["development_stage"])
    frame["library_id"] = "lib1"
    frame["raw_sum"] = 2000.0
    soma = _FakeSomaObs(frame)
    got = cs.query_obs(_FakeCensus(soma))

    call = soma.calls[0]
    assert call["value_filter"] == cs.VALUE_FILTER
    assert set(cs.OBS_COLUMNS) - {"development_stage"} <= set(call["column_names"])
    assert "library_id" in call["column_names"], "a pool/library column must be materialized (D3)"
    assert got.attrs["census_version"] == cs.CENSUS_VERSION
    assert got.attrs["missing_columns"] == ["development_stage"]
    assert got.attrs["pool_columns"] == ["library_id"]
    assert got.attrs["depth_column"] == "raw_sum"


def test_query_obs_refuses_a_schema_without_the_stratum_columns():
    frame = obs_frame(healthy_donors()).drop(columns=["donor_id"])
    with pytest.raises(ValueError, match="donor_id"):
        cs.query_obs(_FakeCensus(_FakeSomaObs(frame)))


def test_the_stand_in_rejects_a_filter_naming_a_column_the_schema_does_not_have():
    """The stand-in must fail where the Census fails, or it certifies filters that cannot run.

    This is the test that would have caught the first live dry run's death on the ground: the §1
    filter text was being sent to ``obs`` verbatim, and the old stand-in — which never resolved
    filter names against its schema — accepted it happily while tiledbsoma raised
    ``'Column organism does not exist in schema'`` at reader construction.
    """
    soma = _FakeSomaObs(obs_frame(healthy_donors()))
    with pytest.raises(_FakeSomaError, match="Column organism does not exist"):
        cs.query_obs(_FakeCensus(soma), value_filter=cs.SPEC_VALUE_FILTER)

    # The executable filter passes the same check, and `is_primary_data` resolves because the real
    # Census obs schema has it even though these synthetic frames carry no such column.
    cs.query_obs(_FakeCensus(soma), value_filter=cs.VALUE_FILTER)
    assert soma.calls[-1]["value_filter"] == cs.VALUE_FILTER


# ---------------------------------------------------------------------------
# Stratum and contrast formation (§1 Unit of analysis).
# ---------------------------------------------------------------------------

def test_each_disease_term_gets_its_own_contrast_never_any_disease():
    obs = obs_frame([
        {"donor_id": "d1", "disease": "COVID-19"},
        {"donor_id": "d2", "disease": "influenza"},
        {"donor_id": "d3", "disease": "normal"},
    ])
    contrasts = cs.build_contrasts(obs)
    assert [c.disease for c in contrasts] == ["COVID-19", "influenza"]
    assert {c.reference for c in contrasts} == {"normal"}
    assert not any("any" in c.disease.lower() for c in contrasts)


def test_a_stratum_without_normal_cells_yields_no_contrast():
    obs = obs_frame([
        {"donor_id": "d1", "disease": "COVID-19"},
        {"donor_id": "d2", "disease": "influenza"},
    ])
    assert cs.build_contrasts(obs) == []


def test_strata_are_dataset_by_cell_type():
    obs = pd.concat([
        obs_frame(healthy_donors(), cell_type="T cell"),
        obs_frame(healthy_donors(), cell_type="B cell"),
        obs_frame(healthy_donors(), dataset_id="ds2", cell_type="T cell"),
    ], ignore_index=True)
    keys = {(c.dataset_id, c.cell_type) for c in cs.build_contrasts(obs)}
    assert keys == {("ds1", "T cell"), ("ds1", "B cell"), ("ds2", "T cell")}


# ---------------------------------------------------------------------------
# Inclusion gate (§1 items 1-3 here; items 4-5 are pending by construction).
# ---------------------------------------------------------------------------

def test_thin_donors_are_dropped_not_merged():
    """§1 item 2: "donors below threshold are dropped, not merged"."""
    donors = healthy_donors(n_per_group=4, n_cells=20)
    donors[0]["n_cells"] = 9  # one cell below gate_config.MIN_CELLS
    obs = obs_frame(donors)
    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])

    assert gate.dropped_thin_donors == {"A0": 9}
    assert "A0" not in set(gate.obs["donor_id"])
    assert gate.n_donors == {"A": 3, "B": 4}
    # Dropped, not merged: the 9 cells are gone, not redistributed into a surviving donor.
    assert gate.n_cells["A"] == 60
    assert set(gate.obs.loc[gate.obs["condition"] == "COVID-19", "donor_id"].value_counts()) == {20}
    assert gate.passed


def test_the_thin_donor_threshold_comes_from_the_frozen_config():
    assert cs.apply_inclusion_gate.__defaults__ is None  # keyword-only, so read the signature
    import inspect
    default = inspect.signature(cs.apply_inclusion_gate).parameters["min_cells_per_donor"].default
    assert default == gc.MIN_CELLS == 10
    assert cs.MIN_DONORS_PER_GROUP == 3
    assert "≥ 3 distinct `donor_id` in EACH group" in SPEC


def test_fewer_than_three_donors_in_a_group_fails_the_gate():
    donors = [{"donor_id": f"A{i}", "disease": "COVID-19"} for i in range(2)]
    donors += [{"donor_id": f"B{i}", "disease": "normal"} for i in range(5)]
    obs = obs_frame(donors)
    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])
    assert not gate.passed
    assert any("group A has 2 donor(s)" in f for f in gate.failures)


def test_a_group_thinned_below_three_donors_fails_after_the_drop_not_before():
    donors = healthy_donors(n_per_group=3)
    donors[0]["n_cells"] = 4
    obs = obs_frame(donors)
    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])
    assert gate.dropped_thin_donors == {"A0": 4}
    assert not gate.passed and gate.n_donors["A"] == 2


def test_a_donor_spanning_both_conditions_is_rejected():
    """§1 item 3: donor_id must be nested within condition."""
    donors = healthy_donors(n_per_group=3)
    donors.append({"donor_id": "A0", "disease": "normal", "n_cells": 20})  # same donor, both arms
    obs = obs_frame(donors)
    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])
    assert not gate.passed
    assert any("not nested within condition" in f for f in gate.failures)


def test_one_donor_per_group_is_reported_as_collinear_with_condition():
    obs = obs_frame([
        {"donor_id": "a", "disease": "COVID-19"},
        {"donor_id": "b", "disease": "normal"},
    ])
    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])
    assert not gate.passed
    assert any("collinear with condition (n=1)" in f for f in gate.failures)
    # ...and the >= 3 rule is reported separately, so the manifest says which rule bit.
    assert any("inclusion gate item 1" in f for f in gate.failures)


def test_cells_with_a_null_donor_id_are_dropped_and_counted():
    """Regression: §1 item 3 requires donor_id PRESENT, and pandas 3.0 hides a null from the gate.

    Under ``future.infer_string`` a null survives ``.astype(str)`` instead of becoming "nan", and
    then ``value_counts``/``nunique`` skip it while ``isin(thin_donors)`` is False for it — so 15
    donor-less cells used to sail through a 3v3 gate, be counted in ``n_cells`` and travel on into
    the stratum. They must be removed explicitly and counted, like the thin-donor drop.
    """
    donors = healthy_donors(n_per_group=3, n_cells=20)
    donors.append({"donor_id": None, "disease": "COVID-19", "n_cells": 15})
    obs = obs_frame(donors)
    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])

    assert gate.dropped_missing_donor_cells == 15
    assert gate.obs.shape[0] == 120
    assert not gate.obs["donor_id"].isna().any()
    assert gate.n_cells == {"A": 60, "B": 60} and gate.n_donors == {"A": 3, "B": 3}
    assert gate.passed, "the surviving 3v3 design is still fine — the null cells simply left"

    rows = cs.screen_strata(obs)
    assert rows[0]["dropped_missing_donor_cells"] == 15
    assert rows[0]["n_cells"] == 120


def test_a_stratum_whose_donor_id_is_null_everywhere_fails_the_gate():
    obs = obs_frame([
        {"donor_id": None, "disease": "COVID-19", "n_cells": 40},
        {"donor_id": None, "disease": "normal", "n_cells": 40},
    ])
    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])
    assert not gate.passed
    assert any("null on all 80 cell(s)" in f for f in gate.failures)
    assert gate.obs.shape[0] == 0


def test_the_confound_screen_refuses_an_ungated_frame_with_null_donors():
    """A null donor would be dropped by every groupby and silently change the denominator."""
    obs = obs_frame(healthy_donors())
    obs.loc[0, "donor_id"] = None
    obs["condition"] = obs["disease"]
    with pytest.raises(ValueError, match="apply_inclusion_gate"):
        cs.confound_prescreen(obs)


def test_rows_with_a_null_stratum_key_form_no_contrast():
    obs = obs_frame(healthy_donors())
    obs.loc[0:9, "cell_type"] = None
    assert [c.cell_type for c in cs.build_contrasts(obs)] == ["T cell"]


def test_a_missing_donor_column_is_a_gate_failure_not_a_crash():
    obs = obs_frame(healthy_donors()).drop(columns=["donor_id"])
    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])
    assert gate.failures == ["donor_id_missing"]


def test_items_4_and_5_are_pending_never_passed():
    obs = obs_frame(healthy_donors())
    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])
    assert set(gate.pending) == {"integer_check", "frozen_universe_size"}
    assert "io_counts" in gate.pending["integer_check"]
    assert gate.passed, "passed() means 'everything obs can check', with 4 and 5 still pending"


def test_only_the_contrast_cells_count_toward_the_thin_donor_rule():
    """A donor's cells under a third disease term do not help it clear item 2."""
    donors = healthy_donors(n_per_group=3)
    donors.append({"donor_id": "A9", "disease": "COVID-19", "n_cells": 6})
    donors.append({"donor_id": "A9", "disease": "influenza", "n_cells": 30})
    obs = obs_frame(donors)
    covid = next(c for c in cs.build_contrasts(obs) if c.disease == "COVID-19")
    gate = cs.apply_inclusion_gate(obs, covid)
    assert gate.dropped_thin_donors == {"A9": 6}


# ---------------------------------------------------------------------------
# Confound pre-screen (§1, D3, D4) — donor level, never cell-weighted.
# ---------------------------------------------------------------------------

def _screen(donors, **kwargs):
    obs = obs_frame(donors)
    contrast = cs.build_contrasts(obs)[0]
    gate = cs.apply_inclusion_gate(obs, contrast)
    return gate, cs.confound_prescreen(gate.obs, **kwargs)


def test_the_screen_reuses_the_design_auditors_donor_level_helpers():
    """One implementation of Cramér's V in this package, not two that can drift apart."""
    assert cs._cramers_v is design._cramers_v
    assert cs._donor_level_table is design._donor_level_table
    assert cs._perfectly_separates is design._perfectly_separates


def test_perfect_separation_by_assay_excludes_the_stratum():
    donors = healthy_donors(n_per_group=3)
    for d in donors:
        d["assay"] = "10x 3' v3" if d["disease"] == "normal" else "10x 5' v1"
    _, report = _screen(donors)
    assert report.separates["assay"] and report.cramers_v["assay"] == 1.0
    assert report.exclude
    assert any("perfectly separated by 'assay'" in r for r in report.exclusion_reasons)


def test_partial_confounding_is_tagged_and_kept():
    donors = healthy_donors(n_per_group=4)
    for d in donors[:2] + donors[4:5]:   # two disease donors and one control on the other assay
        d["assay"] = "10x 5' v1"
    _, report = _screen(donors)
    v = report.cramers_v["assay"]
    assert 0 < v < cs.PERFECT_CONFOUND_V
    assert not report.exclude
    assert any(f.startswith("assay:partial_confound") for f in report.flags)
    assert "assay" in report.covariate_candidates


def test_a_covariate_level_held_by_one_donor_is_reported_not_scored():
    """One donor on its own assay is that donor, not an assay effect — but it is not hidden."""
    donors = healthy_donors(n_per_group=4)
    donors[0]["assay"] = "10x 5' v1"
    _, report = _screen(donors)
    assert report.cramers_v["assay"] is None and report.separates["assay"] is False
    assert not report.exclude
    assert any("assay:1 donor-unique level(s) set aside" in f for f in report.flags)
    assert "assay" not in report.covariate_candidates


def test_tissue_and_depth_separation_are_flagged_but_do_not_exclude():
    """§1's exclusion clause names assay/suspension/pool — three covariates, not five."""
    donors = healthy_donors(n_per_group=3)
    for d in donors:
        d["tissue_general"] = "blood" if d["disease"] == "normal" else "lung"
    _, report = _screen(donors)
    assert report.separates["tissue_general"]
    assert "tissue_general:perfect_separation" in report.flags
    assert not report.exclude
    assert cs.EXCLUDING_COVARIATES == ("assay", "suspension_type")


def test_the_association_is_measured_per_donor_not_per_cell():
    """A confound carried by small donors must not be diluted by cell counts."""
    donors = healthy_donors(n_per_group=4, n_cells=400)
    for d in donors[:2]:            # two disease donors, on their own assay, starved of cells
        d["assay"] = "10x 5' v1"
        d["n_cells"] = 10
    _, report = _screen(donors)
    obs = obs_frame(donors)
    cell_v = design._cramers_v(obs["assay"].astype(str), obs["disease"].astype(str))
    assert report.cramers_v["assay"] > 3 * cell_v


def test_a_donor_unique_covariate_does_not_count_as_a_confound():
    """A per-donor library id separates condition trivially — donor is nested in condition by
    design (§1 item 3). Excluding on it would exclude every valid case/control stratum."""
    donors = healthy_donors(n_per_group=3)
    for i, d in enumerate(donors):
        d["library_id"] = f"lib{i}"
    _, report = _screen(donors)
    assert not report.exclude
    assert any("library_id:donor-unique levels" in f for f in report.flags)
    assert report.pooled is False


def test_a_singleton_pool_coincidence_does_not_exclude_a_large_stratum():
    """Regression: the donor-unique rule must apply BEFORE the association is measured.

    Applied after (skip only when EVERY level is donor-unique), two donors who happen to share one
    library and one condition among 38 donors on unique libraries scored V = 1.0, perfect
    separation and EXCLUDE. Census library ids are near-donor-unique, so this dropped usable strata
    systematically. Only levels grouping >= 2 donors carry information beyond donor_id.
    """
    donors = []
    for i in range(20):
        donors.append({"donor_id": f"A{i}", "disease": "COVID-19",
                       "library_id": "shared_lane" if i < 2 else f"lib_A{i}"})
    for i in range(20):
        donors.append({"donor_id": f"B{i}", "disease": "normal", "library_id": f"lib_B{i}"})
    _, report = _screen(donors)

    assert not report.exclude, "a two-donor pool coincidence is not a structural confound"
    assert report.cramers_v["library_id"] is None      # one shared level: nothing to associate with
    assert report.separates["library_id"] is False
    assert report.pooled is True, "the pooling flag still fires — that is D3's job, not exclusion"
    assert any("38 donor-unique level(s) set aside" in f for f in report.flags)


def test_a_structural_pool_confound_still_excludes():
    """The other side of the same rule: pools that really group donors along condition EXCLUDE."""
    donors = []
    for i in range(10):
        donors.append({"donor_id": f"A{i}", "disease": "COVID-19",
                       "library_id": f"lane_disease_{i % 2}"})
    for i in range(10):
        donors.append({"donor_id": f"B{i}", "disease": "normal",
                       "library_id": f"lane_normal_{i % 2}"})
    _, report = _screen(donors)

    assert report.separates["library_id"] and report.cramers_v["library_id"] == 1.0
    assert report.exclude and any("'library_id'" in r for r in report.exclusion_reasons)


def test_assay_separation_survives_singleton_levels_elsewhere():
    """Setting singleton levels aside must not disarm a real assay confound in the same stratum."""
    donors = healthy_donors(n_per_group=4)
    for d in donors:
        d["assay"] = "10x 3' v3" if d["disease"] == "normal" else "10x 5' v1"
    donors.append({"donor_id": "A9", "disease": "COVID-19", "assay": "Smart-seq2", "n_cells": 20})
    _, report = _screen(donors)
    assert report.separates["assay"] and report.exclude
    assert any("1 donor-unique level(s) set aside" in f for f in report.flags)


def test_donors_sharing_a_pool_set_the_pooling_flag():
    donors = healthy_donors(n_per_group=4)
    for d in donors:
        d["pool_id"] = "poolX" if d["donor_id"] in {"A0", "B0"} else f"solo_{d['donor_id']}"
    _, report = _screen(donors)
    assert report.pooled is True
    assert report.pool_evidence["pool_id"] == {"poolX": 2}
    assert not report.exclude, "a pool shared ACROSS conditions confounds nothing"


def test_a_pool_that_separates_condition_excludes_the_stratum():
    donors = healthy_donors(n_per_group=4)
    for d in donors:
        d["pool_id"] = "poolDisease" if d["disease"] != "normal" else "poolNormal"
    _, report = _screen(donors)
    assert report.pooled is True
    assert report.exclude and any("'pool_id'" in r for r in report.exclusion_reasons)


def test_pooling_is_unresolved_not_false_when_obs_exposes_no_pool_column():
    """D3: where pooling cannot be resolved, donor-pseudobulk is a lower bound — not 'not pooled'."""
    _, report = _screen(healthy_donors())
    assert report.pooled == cs.POOLING_UNRESOLVED
    assert "pooled_flag" in report.pending
    assert "lower bound" in report.pending["pooled_flag"].lower()


def test_the_depth_bin_is_screened_when_obs_carries_counts_and_pending_when_it_does_not():
    donors = healthy_donors(n_per_group=4)
    for d in donors:                     # depth split exactly along condition
        d["raw_sum"] = 8000.0 if d["disease"] == "normal" else 1000.0
    _, report = _screen(donors)
    assert cs.DEPTH_BIN_COL in report.cramers_v
    assert report.separates[cs.DEPTH_BIN_COL]
    assert not report.exclude, "depth is handled by §7.2 depth-matching and the null, not exclusion"

    _, no_depth = _screen(healthy_donors(n_per_group=4))
    assert cs.DEPTH_BIN_COL not in no_depth.cramers_v
    assert cs.DEPTH_BIN_COL in no_depth.pending


def test_nnz_bins_depth_but_is_never_reported_as_counts_per_cell():
    """``nnz`` counts detected genes, not molecules; only ``raw_sum`` may fill a counts column."""
    donors = healthy_donors(n_per_group=4)
    for d in donors:
        d["nnz"] = 900.0 if d["disease"] == "normal" else 300.0
    obs = obs_frame(donors)
    rows = cs.screen_strata(obs)
    assert rows[0]["median_counts_per_cell_by_group"] == cs.PENDING

    gate = cs.apply_inclusion_gate(obs, cs.build_contrasts(obs)[0])
    report = cs.confound_prescreen(gate.obs)
    assert cs.DEPTH_BIN_COL in report.cramers_v


# ---------------------------------------------------------------------------
# Manifest (§1) — including what it must refuse to claim.
# ---------------------------------------------------------------------------

def _rows(donors=None, **kwargs):
    obs = obs_frame(donors or healthy_donors(n_per_group=4), **kwargs)
    return obs, cs.screen_strata(obs)


def test_the_manifest_carries_every_column_the_spec_lists():
    spec_line = next(ln for ln in SPEC.splitlines() if ln.startswith("**Manifest (auditable"))
    for token in ("dataset_id", "cell_type", "cell_type_ontology_depth", "n_cells",
                  "cells_per_donor_by_group", "median_counts_per_cell_by_group", "confound_flags",
                  "pooled_flag", "residual_df", "integer_check", "frozen_universe_size", "fitType"):
        assert token in spec_line
        assert token in cs.MANIFEST_FIELDS
    assert "n_donors_A/B" in spec_line
    assert {"n_donors_A", "n_donors_B"} <= set(cs.MANIFEST_FIELDS)
    assert "permutation count" in spec_line and "permutation_count" in cs.MANIFEST_FIELDS

    _, rows = _rows()
    assert set(rows[0]) == set(cs.MANIFEST_FIELDS)


def test_every_pending_column_says_pending_and_names_its_owner():
    _, rows = _rows()
    row = rows[0]
    assert row["gate_status"] == "candidate"
    for column in ("cell_type_ontology_depth", "integer_check", "frozen_universe_size", "fitType",
                   "sigma_donor_estimate", "envelope_min_donors_per_group", "envelope_membership"):
        assert row[column] == cs.PENDING, column
        assert column in cs.PENDING_FIELDS
    assert "io_counts" in cs.PENDING_FIELDS["integer_check"]
    assert "OPEN" in cs.PENDING_FIELDS["sigma_donor_estimate"]


def test_the_sigma_donor_mechanism_is_not_described_as_an_upper_bound():
    """Amendment 4 Part A, Correction 1 retracts Amendment 3's "upper bound" claim.

    ``sqrt(s0^2) * ln 2`` is not a ceiling on ``donor_sigma``: the ``+1`` of ``log2(CPM + 1)``
    attenuates the donor random effect along with everything else inside the logarithm, so the
    quantity understates ``donor_sigma`` whenever the technical variance is smaller than
    ``sigma^2 (1/a_bar^2 - 1)`` — which is the ordinary case on clean, wide-universe data. A string
    that still calls it an upper bound tells a reader the number errs in the safe direction when it
    does not, and that is the failure mode the correction exists to close. Both the module docstring
    and the manifest's pending-field text are pinned here, because both carried the claim.
    """
    note = cs.PENDING_FIELDS["sigma_donor_estimate"]
    assert "AUDIT QUANTITY OF UNKNOWN ERROR SIGN" in note
    assert "RETRACTS" in note and "Correction 1" in note
    assert "gates nothing" in note
    # The retraction must be visible wherever the claim was made, not only in the amendment log.
    doc = cs.__doc__
    assert "audit quantity of\n  unknown error sign" in doc
    assert "Correction 1 retracts that" in doc
    # No surviving sentence may assert the bound. The only permitted occurrences of the phrase are
    # the ones that name it as retracted.
    for text in (note, doc):
        for line in text.splitlines():
            low = line.lower()
            if "upper bound" in low:
                assert ("retract" in low or "called" in low or "also called" in low), line


def test_median_counts_per_cell_is_computed_from_obs_when_available_and_pending_otherwise():
    donors = healthy_donors(n_per_group=4)
    for d in donors:
        d["raw_sum"] = 5000.0 if d["disease"] == "normal" else 4000.0
    _, rows = _rows(donors)
    assert rows[0]["median_counts_per_cell_by_group"] == {"A": 4000.0, "B": 5000.0}

    _, plain = _rows()
    assert plain[0]["median_counts_per_cell_by_group"] == cs.PENDING


def test_the_permutation_count_is_c_of_d_choose_n_a():
    donors = [{"donor_id": f"A{i}", "disease": "COVID-19"} for i in range(4)]
    donors += [{"donor_id": f"B{i}", "disease": "normal"} for i in range(5)]
    _, rows = _rows(donors)
    row = rows[0]
    assert (row["n_donors_A"], row["n_donors_B"]) == (4, 5)
    assert row["permutation_count"] == math.comb(9, 4) == 126
    assert row["residual_df"] == 7           # D - 2, the two-parameter design the arm actually fits


def test_the_manifest_promises_no_covariate_mechanism_the_arm_does_not_have():
    """Amendment 2 Change 1: "C1 and C4 are DESeq2-specific and lapse with it".

    The shipped moderated arm fits ``~ 1 + x`` and has no covariate slot, so the manifest must not
    carry a covariate-affordability column. A partial confound is tagged and neutralised by the
    permutation null, and the tag says exactly that.
    """
    assert "C1 and C4 are DESeq2-specific and lapse with it" in AMENDMENTS
    assert not [f for f in cs.MANIFEST_FIELDS if "covariate" in f and f != "covariate_candidates"]
    assert "covariate_candidates" in cs.MANIFEST_FIELDS
    assert "lapse with it" in cs.confound_prescreen.__doc__
    assert "permutation null" in cs.confound_prescreen.__doc__


def test_no_row_is_ever_admitted_to_the_sweep():
    """Amendment 3: envelope membership needs a sigma_donor estimate, and the anchor is OPEN."""
    _, rows = _rows(healthy_donors(n_per_group=20))
    row = rows[0]
    assert row["admitted_to_sweep"] is False
    assert "sigma_donor_estimate" in row["admission_blockers"]
    assert row["envelope_membership"] == cs.PENDING
    # The one envelope number obs can supply, and it is a conditional: 20 donors/group meets the
    # envelope's demand at sigma_donor <= 0.5 (13 donors) but not at 0.7 (23 donors).
    assert row["envelope_max_sigma_supported"] == 0.5


def test_envelope_max_sigma_tracks_the_frozen_envelope():
    for entry in gc.OPERATING_ENVELOPE:
        n = int(entry["min_donors_per_group"])
        assert cs.envelope_max_sigma_supported(n) >= float(entry["sigma_donor"])
        assert cs.envelope_max_sigma_supported(n - 1) != float(entry["sigma_donor"])
    smallest = min(int(e["min_donors_per_group"]) for e in gc.OPERATING_ENVELOPE)
    assert cs.envelope_max_sigma_supported(smallest - 1) is None


def test_excluded_strata_stay_in_the_manifest_with_their_reasons():
    """D4 asks for the count and characteristics of what was excluded; survivors cannot say."""
    good = healthy_donors(n_per_group=4)
    thin = [{"donor_id": f"C{i}", "disease": "influenza", "n_cells": 20} for i in range(2)]
    obs = pd.concat([obs_frame(good), obs_frame(thin)], ignore_index=True)
    rows = cs.screen_strata(obs)
    statuses = {r["disease"]: r["gate_status"] for r in rows}
    assert statuses == {"COVID-19": "candidate", "influenza": "excluded_inclusion_gate"}

    manifest = cs.emit_manifest(rows, obs=obs)
    counts = manifest["header"]["counts"]
    assert counts["n_contrasts"] == 2
    assert counts["excluded_fraction"] == 0.5
    assert counts["by_gate_status"]["excluded_inclusion_gate"] == 1
    assert any("inclusion gate item 1" in k for k in counts["inclusion_gate_failure_reasons"])
    excluded = next(r for r in rows if r["gate_status"] == "excluded_inclusion_gate")
    assert excluded["confound_flags"] == ["not_screened:failed_inclusion_gate"]


def test_the_header_pins_the_census_version_and_the_package_versions():
    obs, rows = _rows()
    header = cs.emit_manifest(rows, obs=obs)["header"]
    assert header["census_version"] == "2025-01-30"
    # Both the filter that ran and the filter §1 pre-registers, since they are not the same string
    # and a header carrying only one of them would either misquote the protocol or misreport the
    # query. The organism clause is accounted for rather than silently missing.
    assert header["value_filter"] == cs.VALUE_FILTER
    assert header["value_filter_spec"] == cs.SPEC_VALUE_FILTER
    assert header["organism"] == cs.CENSUS_ORGANISM
    assert "not an obs filter clause" in header["organism_realized_by"]
    for package in ("cellxgene-census", "anndata", "scanpy", "numpy", "scipy", "pandas",
                    "statsmodels", "pydeseq2", "decoupler"):
        assert package in header["package_versions"]
    assert header["package_versions"]["pandas"] == pd.__version__
    assert header["config"]["pre_registered"]["min_cells"] == gc.MIN_CELLS
    assert header["operating_envelope"] == [dict(r) for r in gc.OPERATING_ENVELOPE]
    assert "CANDIDATES ONLY" in header["scope"]
    assert all(row["census_version"] == "2025-01-30" for row in rows)


def test_emit_manifest_writes_json_and_csv(tmp_path):
    obs, rows = _rows()
    cs.emit_manifest(rows, out_dir=tmp_path, obs=obs)

    written = json.loads((tmp_path / "census_strata.json").read_text(encoding="utf-8"))
    assert set(written) == {"header", "rows"}
    assert written["rows"][0]["dataset_id"] == "ds1"

    with open(tmp_path / "census_strata.csv", encoding="utf-8", newline="") as fh:
        table = list(csv.DictReader(fh))
    assert list(table[0]) == list(cs.MANIFEST_FIELDS)
    assert table[0]["census_version"] == "2025-01-30"
    assert table[0]["integer_check"] == cs.PENDING
    # dict/list cells survive the CSV as JSON rather than as str(dict).
    assert json.loads(table[0]["cells_per_donor_by_group"])["A"]["n_donors"] == 4


def test_emit_manifest_without_an_out_dir_writes_nothing(tmp_path):
    _, rows = _rows()
    manifest = cs.emit_manifest(rows)
    assert manifest["rows"] and list(tmp_path.iterdir()) == []
