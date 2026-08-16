"""The §1 stratum-list freeze, tested as the pre-registration it is.

``docs/PREREGISTRATION_STRATUM_LIST.md`` claims that a reader can re-derive every number in it from
one committed artifact, and that the 251 frozen strata follow mechanically from twelve dataset ids
with no stratum-level discretion anywhere. Both claims are only worth something if they are checked
by something other than the person who wrote them, so this file re-derives the load-bearing figures
from the raw manifest — deliberately **not** through ``freeze_stratum_list``'s own helpers where the
point is to catch that module drifting — and compares.

The failure modes worth catching here are all quiet ones:

* the source artifact is swapped, truncated or re-generated, and the pre-registered numbers silently
  come to describe different bytes;
* a dataset id is added, removed or edited after the freeze, changing the D2 denominator;
* the selection rule stops being "every candidate row of the twelve", so a stratum leaves the
  analysis set without a reason attached;
* a cells-per-donor bin turns out empty, which would make §1 (iii)'s "spanning the bins" false;
* the Layer B truncation stops matching the envelope arithmetic — the one thing that keeps the
  surviving subset from being chosen after the sigma_donor anchor lands;
* a Layer B tier goes undeclared, so the document publishes a row nothing ever checked. This is not
  hypothetical: the first version of the freeze declared only sigma 0.5 and 0.7, and the sigma 0.35
  row of §6's table — the tier at ``gate_config.POWER_EVAL_SIGMA``, where the instrument nominally
  operates — was published with a verdict that contradicted ``SPEC_DATASET_FLOOR`` in the same file;
* §6's table drifts from the artifact it describes, in either direction: the verdict column is
  ``below_spec_dataset_floor`` and is parsed back out of the document here rather than trusted;
* the two same-collection siblings stop being frozen, or stop being distinguishable from the 251;
* the emitted artifact stops being byte-reproducible, at which point "regenerate it and diff" is no
  longer an audit anyone can run.

No network, and no Census. Everything is read from ``pilot/preregistration/``.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter
from pathlib import Path

import pytest

from pbcheck import census_select as cs
from pbcheck import gate_config as gc

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "freeze_stratum_list.py"
PREREG = REPO / "pilot" / "preregistration"
DOCUMENT = REPO / "docs" / "PREREGISTRATION_STRATUM_LIST.md"


def _load_script():
    """Import the freeze driver as a module (``scripts/`` is not on the import path)."""
    if not SCRIPT.exists():  # pragma: no cover - the file is part of the commit under test
        pytest.skip("freeze_stratum_list.py not present")
    spec = importlib.util.spec_from_file_location("_freeze_stratum_list", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


fz = _load_script()


# ---------------------------------------------------------------------------
# Fixtures. The source manifest is 6.6 MB; parse it once per module.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def raw_source() -> dict:
    """The committed candidate manifest, read directly — not through ``load_source``.

    Reading it raw is the point: every re-derivation below has to be independent of the module it
    is checking, or a bug in that module would be reproduced by the test that is meant to catch it.
    """
    return json.loads(fz.SOURCE_JSON.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def frozen(raw_source) -> dict:
    """The freeze, rebuilt from the committed source (this one does go through the module)."""
    return fz.build(fz.load_source())


@pytest.fixture(scope="module")
def committed() -> dict:
    return json.loads((PREREG / f"{fz.OUT_STEM}.json").read_text(encoding="utf-8"))


def _rule(rows, dataset_ids) -> list[dict]:
    """The selection rule, written out again here rather than imported."""
    wanted = set(dataset_ids)
    return [r for r in rows if r["gate_status"] == "candidate" and r["dataset_id"] in wanted]


def _key(row) -> tuple[str, str, str]:
    return (row["dataset_id"], row["cell_type"], row["disease"])


def _min_donors(row) -> int:
    return min(row["n_donors_A"], row["n_donors_B"])


def _frozen_rows(raw_source) -> list[dict]:
    return _rule(raw_source["rows"], [d["dataset_id"] for d in fz.FROZEN_DATASETS])


def _candidates(raw_source) -> list[dict]:
    return [r for r in raw_source["rows"] if r["gate_status"] == "candidate"]


def _document_text() -> str:
    return DOCUMENT.read_text(encoding="utf-8")


def _table_cells(line: str) -> list[str]:
    """The cells of one markdown table row, with bold markers and padding stripped."""
    return [cell.replace("*", "").strip() for cell in line.strip().strip("|").split("|")]


def _leading_int(cell: str) -> int:
    match = re.match(r"\s*(\d+)", cell)
    assert match, f"no leading integer in table cell {cell!r}"
    return int(match.group(1))


def _parse_layer_b_table() -> list[dict]:
    """§6's Layer B table, read back out of the document.

    The document is the human half of the pre-registration and the artifact is the machine half;
    neither is authoritative alone, so the one number that decides §6's headline — the verdict
    against §1's dataset floor — is compared rather than trusted in either direction.
    """
    rows = []
    for line in _document_text().splitlines():
        cells = _table_cells(line)
        if len(cells) != 5 or not cells[0].startswith("≈"):
            continue
        threshold = re.match(r"≥\s*(\d+)\s*v\s*(\d+)", cells[1])
        assert threshold and threshold.group(1) == threshold.group(2), cells[1]
        rows.append({
            "sigma_donor": float(cells[0].lstrip("≈ ")),
            "min_donors_per_group": int(threshold.group(1)),
            "n_datasets": _leading_int(cells[2]),
            "n_strata": _leading_int(cells[3]),
            "verdict_below_floor": cells[4].upper() == "BELOW",
            "verdict_text": cells[4],
        })
    return rows


def _parse_counterfactual_table() -> list[dict]:
    """§6's manifest-wide table: what the 68 candidate-bearing datasets could have supported."""
    rows = []
    for line in _document_text().splitlines():
        cells = _table_cells(line)
        if len(cells) != 4 or not re.match(r"^≥\s*\d+\s*v\s*\d+$", cells[0]):
            continue
        rows.append({
            "min_donors_per_group": _leading_int(cells[0].lstrip("≥ ")),
            "n_manifest_datasets_at_tier": _leading_int(cells[1]),
            "n_datasets": _leading_int(cells[2]),
            "counterfactual_max_survivors_of_twelve": _leading_int(cells[3]),
        })
    return rows


# ---------------------------------------------------------------------------
# The source artifact
# ---------------------------------------------------------------------------


def test_source_artifacts_are_committed_with_the_pinned_hashes():
    """The pre-registered numbers describe these exact bytes and no others."""
    for path, digest, size in (
        (fz.SOURCE_JSON, fz.SOURCE_SHA256, fz.SOURCE_BYTES),
        (fz.SOURCE_CSV, fz.SOURCE_CSV_SHA256, fz.SOURCE_CSV_BYTES),
    ):
        assert path.exists(), f"{path} is not committed"
        assert path.stat().st_size == size, f"{path.name} is {path.stat().st_size} bytes"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, f"{path.name} sha256"


def test_source_header_carries_the_pinned_run_stamps(raw_source):
    header = raw_source["header"]
    assert header["generated_utc"] == fz.SOURCE_GENERATED_UTC
    assert header["census_version"] == cs.CENSUS_VERSION == "2025-01-30"
    assert len(raw_source["rows"]) == fz.SOURCE_N_ROWS == 2190
    candidates = [r for r in raw_source["rows"] if r["gate_status"] == "candidate"]
    assert len(candidates) == fz.SOURCE_N_CANDIDATES == 1197


def test_freeze_refuses_a_source_whose_bytes_differ(tmp_path, raw_source):
    """A same-length edit — the case a size check alone would wave through."""
    original = fz.SOURCE_JSON.read_bytes()
    tampered = tmp_path / fz.SOURCE_JSON.name
    # Flip one character inside a dataset_id, keeping the length identical.
    marker = b"6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3"
    assert original.count(marker) > 0
    tampered.write_bytes(original.replace(marker, b"6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc4", 1))
    assert tampered.stat().st_size == fz.SOURCE_BYTES

    with pytest.raises(fz.SourceArtifactMismatch, match="sha256"):
        fz.load_source(tampered, csv_path=None)


def test_freeze_refuses_a_truncated_source(tmp_path):
    truncated = tmp_path / fz.SOURCE_JSON.name
    truncated.write_bytes(fz.SOURCE_JSON.read_bytes()[: fz.SOURCE_BYTES // 2])
    with pytest.raises(fz.SourceArtifactMismatch, match="bytes"):
        fz.load_source(truncated, csv_path=None)


def test_freeze_refuses_a_different_run_even_if_the_hash_constant_was_edited(
    tmp_path, monkeypatch, raw_source
):
    """The hash guards the file; the header stamps guard the constant. Both are needed.

    Simulates the one edit that defeats a hash check on its own: point the pin at a different
    artifact. The header check must still fire.
    """
    other = dict(raw_source)
    other["header"] = {**raw_source["header"], "generated_utc": "2026-09-01T00:00:00+00:00"}
    path = tmp_path / "other_run.json"
    payload = json.dumps(other, ensure_ascii=False).encode("utf-8")
    path.write_bytes(payload)

    monkeypatch.setattr(fz, "SOURCE_BYTES", len(payload))
    monkeypatch.setattr(fz, "SOURCE_SHA256", hashlib.sha256(payload).hexdigest())
    with pytest.raises(fz.SourceArtifactMismatch, match="generated_utc"):
        fz.load_source(path, csv_path=None)


# ---------------------------------------------------------------------------
# The frozen set
# ---------------------------------------------------------------------------


def test_twelve_distinct_datasets_ranked_one_to_twelve():
    ids = [d["dataset_id"] for d in fz.FROZEN_DATASETS]
    assert len(ids) == 12
    assert len(set(ids)) == 12
    assert [d["rank"] for d in fz.FROZEN_DATASETS] == list(range(1, 13))
    uuid = re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$")
    assert all(uuid.match(i) for i in ids)


def test_every_frozen_dataset_holds_candidates_in_the_source(raw_source):
    present = {r["dataset_id"] for r in raw_source["rows"] if r["gate_status"] == "candidate"}
    missing = [d["dataset_id"] for d in fz.FROZEN_DATASETS if d["dataset_id"] not in present]
    assert not missing, f"frozen datasets with no candidate row: {missing}"


def test_the_frozen_set_is_exactly_the_rule_applied_to_the_manifest(raw_source, frozen):
    """No hand-picking at stratum level: the 251 are every candidate row of the twelve."""
    expected = _rule(raw_source["rows"], [d["dataset_id"] for d in fz.FROZEN_DATASETS])
    assert len(expected) == 251
    assert len(frozen["rows"]) == 251
    assert {_key(r) for r in frozen["rows"]} == {_key(r) for r in expected}


def test_no_candidate_of_a_frozen_dataset_was_left_out(raw_source, frozen):
    """The complement of the previous test, stated separately because it is the one that catches a
    quiet exclusion: a stratum dropped for being small, skewed or awkward leaves no trace."""
    included = {_key(r) for r in frozen["rows"]}
    frozen_ids = {d["dataset_id"] for d in fz.FROZEN_DATASETS}
    for row in raw_source["rows"]:
        if row["dataset_id"] in frozen_ids and row["gate_status"] == "candidate":
            assert _key(row) in included, f"candidate stratum omitted: {_key(row)}"


def test_per_dataset_counts_and_ceilings_match_the_declaration(raw_source):
    candidates = [r for r in raw_source["rows"] if r["gate_status"] == "candidate"]
    for dataset in fz.FROZEN_DATASETS:
        own = [r for r in candidates if r["dataset_id"] == dataset["dataset_id"]]
        assert len(own) == dataset["n_strata"], dataset["short"]
        assert max(_min_donors(r) for r in own) == dataset["ceiling_min_donors_per_group"], (
            dataset["short"]
        )


def test_rexach_ceiling_is_ten_not_eleven(raw_source):
    """The correction §3.2 records. The dataset's best design is A=11 vs B=10, whose ceiling
    contribution is 10; reading `max(n_donors_A)` instead gives 11 and is what went wrong."""
    rexach = "ac0c6561-7a48-4185-af6f-af799f699172"
    own = [
        r for r in raw_source["rows"]
        if r["dataset_id"] == rexach and r["gate_status"] == "candidate"
    ]
    assert max(_min_donors(r) for r in own) == 10
    assert max(r["n_donors_A"] for r in own) == 11
    declared = next(d for d in fz.FROZEN_DATASETS if d["dataset_id"] == rexach)
    assert declared["ceiling_min_donors_per_group"] == 10


def test_neither_sibling_is_in_the_frozen_set(frozen):
    """D2's denominator must not contain two datasets from one collection."""
    frozen_ids = {r["dataset_id"] for r in frozen["rows"]}
    for sibling in fz.SIBLING_DATASETS:
        assert sibling["dataset_id"] not in frozen_ids
        assert sibling["sibling_of"] in frozen_ids


def test_siblings_exist_in_the_source_and_are_therefore_a_real_exclusion(raw_source):
    """A sibling that held no candidate row would make its exclusion vacuous — and would mean the
    document names an exclusion it never had to make."""
    candidates = _candidates(raw_source)
    for sibling in fz.SIBLING_DATASETS:
        own = [r for r in candidates if r["dataset_id"] == sibling["dataset_id"]]
        assert own, sibling["dataset_id"]


def test_the_siblings_are_frozen_too_and_not_merely_named(raw_source, frozen):
    """27 runnable strata that were named, excluded from D2 and then left unlisted would have been
    selectable after the fact — the one freedom §3.1 exists to remove."""
    controls = frozen["within_collection_control_rows"]
    expected = _rule(raw_source["rows"], [s["dataset_id"] for s in fz.SIBLING_DATASETS])
    assert len(controls) == len(expected) == 27
    assert {_key(r) for r in controls} == {_key(r) for r in expected}
    assert frozen["header"]["n_within_collection_control_strata"] == 27
    per_sibling = {s["dataset_id"]: s["n_strata"] for s in fz.SIBLING_DATASETS}
    assert per_sibling == {
        "c2876b1b-06d8-4d96-a56b-5304f815b99a": 18,
        "1e5bd3b8-6a0e-4959-8d69-cafed30fe814": 9,
    }
    for sibling in fz.SIBLING_DATASETS:
        own = [r for r in expected if r["dataset_id"] == sibling["dataset_id"]]
        assert len(own) == sibling["n_strata"], sibling["short"]
        assert max(_min_donors(r) for r in own) == sibling["ceiling_min_donors_per_group"]


def test_every_row_carries_a_role_and_the_two_sets_stay_distinguishable(frozen):
    """The role travels with the row, so a row read out of the artifact still says which set it is
    in. The 251 and the 27 must never be conflated: one is the D2 denominator and one is not."""
    assert fz.ROLE_ANALYSIS_SET != fz.ROLE_WITHIN_COLLECTION_CONTROL
    assert "role" in fz.DERIVED_FIELDS
    for row in frozen["rows"]:
        assert row["role"] == fz.ROLE_ANALYSIS_SET
        assert row["dataset_rank"] in range(1, 13)
    for row in frozen["within_collection_control_rows"]:
        assert row["role"] == fz.ROLE_WITHIN_COLLECTION_CONTROL
        assert row["dataset_rank"] is None, "a control row has no place in the D2 denominator"
        assert row["expected_effect"] == fz.EXPECTED_EFFECT_NOT_APPLICABLE
        assert tuple(row) == fz.STRATUM_FIELDS
    assert not ({_key(r) for r in frozen["rows"]}
                & {_key(r) for r in frozen["within_collection_control_rows"]})


def test_control_rows_carry_no_expected_effect_label(frozen):
    """§1 (i)/(ii) is a coverage claim about the twelve. A control row labelled strong or subtle
    would be counted into a claim it is excluded from."""
    labels = {r["expected_effect"] for r in frozen["within_collection_control_rows"]}
    assert labels == {fz.EXPECTED_EFFECT_NOT_APPLICABLE}
    assert fz.EXPECTED_EFFECT_NOT_APPLICABLE not in fz.EXPECTED_EFFECT_VOCABULARY


def test_the_csv_carries_both_blocks_told_apart_by_the_role_column():
    """The CSV is the JSON's twin; one that held only the 251 would silently omit frozen content."""
    with open(PREREG / f"{fz.OUT_STEM}.csv", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 251 + 27 == 278
    counts = Counter(r["role"] for r in rows)
    assert counts == {fz.ROLE_ANALYSIS_SET: 251, fz.ROLE_WITHIN_COLLECTION_CONTROL: 27}
    assert [r["role"] for r in rows[:251]] == [fz.ROLE_ANALYSIS_SET] * 251


def test_the_document_states_the_three_rules_that_bind_the_control_set():
    text = _document_text()
    assert "within_collection_control" in text
    assert "27" in text
    for phrase in (
        "never enters the D2 denominator",
        "Promoting one to an independent dataset is an amendment",
    ):
        assert phrase in text, phrase


def test_every_frozen_row_is_a_candidate_at_the_pinned_census(frozen):
    for row in frozen["rows"]:
        assert row["gate_status"] == "candidate"
        assert row["census_version"] == cs.CENSUS_VERSION
    assert frozen["header"]["census_version"] == cs.CENSUS_VERSION


# ---------------------------------------------------------------------------
# Admission is still closed
# ---------------------------------------------------------------------------


def test_no_row_is_admitted_to_the_sweep(raw_source, frozen):
    assert not any(r["admitted_to_sweep"] for r in raw_source["rows"])
    assert not any(r["admitted_to_sweep"] for r in frozen["rows"])
    assert frozen["header"]["admission"]["admitted_to_sweep"] is False


def test_the_counts_gate_and_the_envelope_are_still_pending_on_every_row(frozen):
    """§9 items 1 and 7: the list may still shrink, and nothing here pretends otherwise."""
    for row in frozen["rows"]:
        assert row["integer_check"] == cs.PENDING
        assert row["frozen_universe_size"] == cs.PENDING
        assert row["sigma_donor_estimate"] == cs.PENDING
        assert row["envelope_membership"] == cs.PENDING


def test_pooling_is_unresolved_on_every_frozen_stratum(frozen):
    """§9 item 3: D3's 'cannot be resolved' state, so donor pseudobulk is a lower bound throughout
    and the gold-standard calibration claim is unavailable on this whole list."""
    assert {r["pooled_flag"] for r in frozen["rows"]} == {cs.POOLING_UNRESOLVED}


# ---------------------------------------------------------------------------
# D1 bins
# ---------------------------------------------------------------------------


def test_bins_are_the_declared_half_decade_sequence_and_tile_the_range():
    assert fz.CELLS_PER_DONOR_BINS == (
        (10, 30), (30, 100), (100, 300), (300, 1000), (1000, 3000), (3000, None),
    )
    lows = [low for low, _ in fz.CELLS_PER_DONOR_BINS]
    highs = [high for _, high in fz.CELLS_PER_DONOR_BINS]
    assert highs[:-1] == lows[1:], "the bins must tile with no gap and no overlap"
    assert highs[-1] is None, "the top bin must stay open"
    assert lows[0] == gc.MIN_CELLS, "the lower edge is the §1 inclusion-gate floor, not a new number"


def test_every_bin_is_occupied_by_the_frozen_set(frozen):
    """§1 (iii) requires cells-per-donor 'spanning the pre-registered bins'. An empty bin makes that
    claim false, and this is the assertion that would say so."""
    occupancy = frozen["header"]["cells_per_donor_bins"]["occupancy"]
    assert len(occupancy) == 6
    empty = [label for label, count in occupancy.items() if count == 0]
    assert not empty, f"unoccupied cells-per-donor bins: {empty}"


def test_bin_occupancy_is_recomputable_from_the_manifest(raw_source, frozen):
    medians = [
        float(r["cells_per_donor_by_group"][g]["median"])
        for r in _rule(raw_source["rows"], [d["dataset_id"] for d in fz.FROZEN_DATASETS])
        for g in ("A", "B")
    ]
    assert len(medians) == 2 * 251 == 502
    block = frozen["header"]["cells_per_donor_bins"]
    assert block["n_group_medians"] == len(medians)
    assert block["min_group_median"] == min(medians)
    assert block["max_group_median"] == max(medians)

    counts = {label: 0 for label in block["bins"]}
    for value in medians:
        for low, high in fz.CELLS_PER_DONOR_BINS:
            if value >= low and (high is None or value < high):
                counts[fz.bin_label(low, high)] += 1
                break
        else:  # pragma: no cover - guarded by the gate's >= 10 cells-per-donor rule
            pytest.fail(f"group median {value} falls in no bin")
    assert counts == block["occupancy"]


def test_a_median_below_the_gate_floor_is_an_error_not_a_bucket():
    """If the gate and the bins ever disagree about the same number, the freeze must stop."""
    with pytest.raises(fz.FrozenDeclarationMismatch):
        fz.bin_of(9.9)


def test_every_row_carries_the_bin_its_own_median_falls_in(frozen):
    for row in frozen["rows"]:
        for group in ("A", "B"):
            expected = fz.bin_of(row["cells_per_donor_by_group"][group]["median"])
            assert row[f"cells_per_donor_bin_{group}"] == expected


# ---------------------------------------------------------------------------
# Layer B
# ---------------------------------------------------------------------------


def _envelope_threshold(sigma: float) -> int:
    return next(
        int(row["min_donors_per_group"])
        for row in gc.OPERATING_ENVELOPE
        if float(row["sigma_donor"]) == sigma
    )


ENVELOPE_SIGMAS = [float(row["sigma_donor"]) for row in gc.OPERATING_ENVELOPE]


def test_layer_b_declares_every_tier_of_the_operating_envelope():
    """A tier the freeze does not declare is a table row nothing checks.

    The first version of this freeze declared 0.5 and 0.7 only, so §6's 0.2 and 0.35 rows were never
    compared against the manifest — and the 0.35 row shipped with a verdict that contradicted
    ``SPEC_DATASET_FLOOR`` in the same file. Declaring the whole envelope is what makes that
    impossible rather than unlikely.
    """
    assert [t["sigma_donor"] for t in fz.LAYER_B] == ENVELOPE_SIGMAS == [0.2, 0.35, 0.5, 0.7]


@pytest.mark.parametrize("sigma", ENVELOPE_SIGMAS)
def test_layer_b_is_exactly_the_datasets_holding_a_stratum_at_the_envelope_threshold(
    raw_source, sigma
):
    """The truncation is pre-declared; this is what stops it from being re-chosen later."""
    threshold = _envelope_threshold(sigma)
    derived = {r["dataset_id"] for r in _frozen_rows(raw_source) if _min_donors(r) >= threshold}
    declared = next(t for t in fz.LAYER_B if t["sigma_donor"] == sigma)
    assert derived == set(declared["dataset_ids"])


@pytest.mark.parametrize("sigma", ENVELOPE_SIGMAS)
def test_the_layer_b_verdict_is_computed_from_the_manifest_not_typed(raw_source, frozen, sigma):
    """``below_spec_dataset_floor`` is the verdict column of §6's table. Re-derived here from the
    raw manifest so that neither the declaration nor the emitted header can decide it alone."""
    threshold = _envelope_threshold(sigma)
    derived = {r["dataset_id"] for r in _frozen_rows(raw_source) if _min_donors(r) >= threshold}
    tier = next(t for t in frozen["header"]["layer_b_truncation"] if t["sigma_donor"] == sigma)
    assert tier["n_datasets"] == len(derived)
    assert tier["below_spec_dataset_floor"] is (len(derived) < fz.SPEC_DATASET_FLOOR)
    assert tier["n_strata"] == sum(
        1 for r in _frozen_rows(raw_source) if _min_donors(r) >= threshold
    )


def test_layer_b_thresholds_come_from_the_operating_envelope():
    """Not restated as literals: if Amendment 3's envelope ever moves, the truncation moves with it
    instead of silently describing the old one."""
    for tier in fz.LAYER_B:
        assert fz._envelope_min_donors(tier["sigma_donor"]) == _envelope_threshold(
            tier["sigma_donor"]
        )
    assert [_envelope_threshold(s) for s in ENVELOPE_SIGMAS] == [4, 8, 13, 23]


def test_three_of_the_four_tiers_fall_below_the_spec_dataset_floor(frozen):
    """§6's headline consequence, pinned so it cannot be quietly softened.

    Seven datasets is below eight. The failure starts at sigma_donor 0.35 —
    ``gate_config.POWER_EVAL_SIGMA``, the envelope boundary Amendment 3 Change 1 makes binding and
    the point at which the instrument nominally operates — not at 0.5, and the list clears §1's own
    '8-12 datasets' floor at the most optimistic tier only.
    """
    tiers = {t["sigma_donor"]: t for t in frozen["header"]["layer_b_truncation"]}
    assert [tiers[s]["n_datasets"] for s in ENVELOPE_SIGMAS] == [11, 7, 5, 3]
    assert [tiers[s]["n_strata"] for s in ENVELOPE_SIGMAS] == [227, 150, 94, 30]
    assert [tiers[s]["below_spec_dataset_floor"] for s in ENVELOPE_SIGMAS] == [
        False, True, True, True
    ]
    assert fz.SPEC_DATASET_FLOOR == 8
    assert tiers[gc.POWER_EVAL_SIGMA]["below_spec_dataset_floor"] is True


def test_the_documents_layer_b_table_agrees_with_the_artifact(frozen):
    """§6's table is parsed back out of the document and compared cell by cell.

    A pre-registration whose prose and whose artifact disagree has pre-registered neither. The
    verdict column in particular is ``below_spec_dataset_floor``, and typing it by hand is how the
    'met, with no margin' cell for seven datasets got published against a floor of eight.
    """
    parsed = _parse_layer_b_table()
    emitted = frozen["header"]["layer_b_truncation"]
    assert len(parsed) == len(emitted) == 4, f"§6's table has {len(parsed)} tier rows"
    for row, tier in zip(parsed, emitted, strict=True):
        where = f"§6 table row sigma {row['sigma_donor']}"
        assert row["sigma_donor"] == tier["sigma_donor"], where
        assert row["min_donors_per_group"] == tier["min_donors_per_group"], where
        assert row["n_datasets"] == tier["n_datasets"], where
        assert row["n_strata"] == tier["n_strata"], where
        assert row["verdict_below_floor"] is tier["below_spec_dataset_floor"], (
            f"{where}: document says {row['verdict_text']!r}, artifact says "
            f"below_spec_dataset_floor={tier['below_spec_dataset_floor']}"
        )


# ---------------------------------------------------------------------------
# What the manifest could have supported — §6's counterfactual
# ---------------------------------------------------------------------------


def test_the_manifest_wide_tier_census_is_recomputable(raw_source):
    """The numbers that refute 'only a handful clear 13 v 13 anywhere'.

    Twenty-one is not a handful, and 21 > 12 means a twelve-dataset list surviving sigma ~ 0.5 was
    assemblable from this manifest. Re-derived here from the raw rows.
    """
    candidates = _candidates(raw_source)
    assert len({r["dataset_id"] for r in candidates}) == 68
    derived = {
        threshold: len({r["dataset_id"] for r in candidates if _min_donors(r) >= threshold})
        for threshold in (4, 8, 13, 23)
    }
    assert derived == {4: 62, 8: 33, 13: 21, 23: 12}
    assert derived == fz.MANIFEST_TIER_CENSUS


def test_the_3v3_anchor_is_scarce_and_never_donor_rich(raw_source):
    """§1 (iii)'s "some exactly 3v3" is what costs a dataset slot at the hard tiers."""
    candidates = _candidates(raw_source)
    anchors = {r["dataset_id"] for r in candidates
               if r["n_donors_A"] == 3 and r["n_donors_B"] == 3}
    assert len(anchors) == fz.MANIFEST_DATASETS_WITH_EXACT_3V3 == 3
    for threshold in (13, 23):
        tier = {r["dataset_id"] for r in candidates if _min_donors(r) >= threshold}
        assert not (anchors & tier), (
            f"a dataset holding an exactly-3v3 stratum also clears {threshold}v{threshold}; §6's "
            "argument that the anchor costs a slot no longer holds"
        )
    # Rexach does clear 8v8, which is why twelve is attainable at the two easy tiers.
    assert fz.COUNTERFACTUAL_WITNESS_ANCHOR in anchors
    assert fz.COUNTERFACTUAL_WITNESS_ANCHOR in {
        r["dataset_id"] for r in candidates if _min_donors(r) >= 8
    }


@pytest.mark.parametrize("threshold", [4, 8, 13, 23])
def test_the_counterfactual_is_attained_by_an_explicit_twelve_dataset_list(raw_source, threshold):
    """Not an upper bound nobody could reach: the witness is constructed and counted.

    Rexach as the 3v3 anchor plus eleven of the twelve datasets that clear 23v23 is a
    §1 (iii)-compliant twelve, and it keeps 12 / 12 / 11 / 11 datasets across the four tiers —
    above §1's floor of eight at every one of them, where the frozen list is above it at one.
    """
    candidates = _candidates(raw_source)
    tier = {r["dataset_id"] for r in candidates if _min_donors(r) >= threshold}
    witness = set(sorted({r["dataset_id"] for r in candidates if _min_donors(r) >= 23})[:11])
    witness.add(fz.COUNTERFACTUAL_WITNESS_ANCHOR)
    assert len(witness) == 12
    attained = len(witness & tier)
    assert attained == fz.COUNTERFACTUAL_MAX_SURVIVORS[threshold]
    assert attained == fz.counterfactual_max_survivors(candidates, threshold)
    assert attained >= fz.SPEC_DATASET_FLOOR, (
        "the counterfactual list must clear §1's dataset floor at every tier, or §6's claim that "
        "the shortfall is a consequence of the selection rather than of the data is wrong"
    )


def test_the_documents_counterfactual_table_agrees_with_the_artifact(frozen):
    parsed = _parse_counterfactual_table()
    emitted = frozen["header"]["layer_b_truncation"]
    assert len(parsed) == len(emitted) == 4
    for row, tier in zip(parsed, emitted, strict=True):
        where = f"§6 counterfactual row >= {row['min_donors_per_group']}"
        assert row["min_donors_per_group"] == tier["min_donors_per_group"], where
        assert row["n_manifest_datasets_at_tier"] == tier["n_manifest_datasets_at_tier"], where
        assert row["n_datasets"] == tier["n_datasets"], where
        assert (row["counterfactual_max_survivors_of_twelve"]
                == tier["counterfactual_max_survivors_of_twelve"]), where


RETRACTED_CLAIMS = (
    "only a handful clear",
    "a property of the public data",
    "met, with no margin",
)

#: A paragraph may still contain a retracted claim — the document quotes both of them, which is how
#: a correction differs from a deletion. What it may not do is assert one.
RETRACTION_MARKERS = ("retracted", "was false", "earlier draft", "simply wrong")


@pytest.mark.parametrize("claim", RETRACTED_CLAIMS)
def test_a_retracted_claim_may_only_appear_inside_its_own_retraction(claim):
    """The two false statements §6 corrected, pinned so neither can drift back in as an assertion.

    "only a handful clear 13 v 13" was refuted by the artifact committed beside it — 21 of the 68
    candidate-bearing datasets clear it — and "met, with no margin" was the verdict given to seven
    datasets against a floor of eight. Both are quoted in the corrected text on purpose; the test
    is that every paragraph carrying one also says it is wrong.
    """
    # Normalised, because the document is hard-wrapped and a retracted phrase can straddle a line.
    paragraphs = [" ".join(p.split()) for p in _document_text().split("\n\n")]
    paragraphs = [p for p in paragraphs if claim in p]
    assert paragraphs, f"the retraction of {claim!r} has gone missing from the document"
    for paragraph in paragraphs:
        assert any(marker in paragraph for marker in RETRACTION_MARKERS), (
            f"{claim!r} appears in a paragraph that does not mark it as retracted:\n{paragraph}"
        )


# ---------------------------------------------------------------------------
# Coverage claims that the document makes and the data must support
# ---------------------------------------------------------------------------


def test_the_donor_count_axis_carries_both_ends_required_by_spec_section_1(frozen):
    """§1 (iii): 'some exactly 3v3, some >= 8v8'."""
    rows = frozen["rows"]
    exactly_3v3 = [r for r in rows if r["n_donors_A"] == 3 and r["n_donors_B"] == 3]
    at_least_8v8 = [r for r in rows if _min_donors(r) >= 8]
    assert len(exactly_3v3) == 10
    assert len(at_least_8v8) == 150
    assert len({r["dataset_id"] for r in at_least_8v8}) == 7

    # §7's partition of the 251 by A1's kill-switch boundary, pinned so the document's arithmetic
    # cannot drift from the data it describes.
    below = [r for r in rows if _min_donors(r) < 8]
    assert len(below) == 101
    assert sum(1 for r in below if _min_donors(r) == 3) == 24
    assert sum(1 for r in below if 4 <= _min_donors(r) <= 7) == 77


def test_expected_effect_labels_cover_the_frozen_arms_exactly(frozen):
    """A literature judgement is only a pre-registration if it is complete and fixed: an unlabelled
    arm could be labelled after the result, and a label with no arm is a leftover."""
    pairs = {(r["dataset_id"], r["disease"]) for r in frozen["rows"]}
    assert pairs == set(fz.EXPECTED_EFFECT)
    assert set(fz.EXPECTED_EFFECT.values()) <= set(fz.EXPECTED_EFFECT_VOCABULARY)
    for row in frozen["rows"]:
        assert row["expected_effect"] == fz.EXPECTED_EFFECT[(row["dataset_id"], row["disease"])]

    # §1 (i) and (ii) ask for 2-3 datasets on each side; the list must clear both, not one. Counted
    # over datasets rather than arms, because that is the unit §1 states and D2 clusters by.
    by_label: dict[str, set[str]] = {label: set() for label in fz.EXPECTED_EFFECT_VOCABULARY}
    for (dataset_id, _disease), label in fz.EXPECTED_EFFECT.items():
        by_label[label].add(dataset_id)
    assert len(by_label["strong"]) >= 2, "§1 (i) wants 2-3 datasets with a strong expected effect"
    assert len(by_label["subtle"]) >= 2, "§1 (ii) wants 2-3 subtle/low-effect datasets"


def test_min_donors_per_group_column_is_the_quantity_the_report_must_print(frozen):
    """§9 item 12: permutation_count flatters skewed designs, so the honest companion travels as its
    own column rather than being recomputed by whoever reads the artifact."""
    for row in frozen["rows"]:
        assert row["min_donors_per_group"] == _min_donors(row)


# ---------------------------------------------------------------------------
# Reproducibility of the artifact
# ---------------------------------------------------------------------------


def test_row_schema_carries_the_manifest_fields_verbatim_and_appends_the_derived_ones(frozen):
    assert fz.STRATUM_FIELDS[: len(cs.MANIFEST_FIELDS)] == tuple(cs.MANIFEST_FIELDS)
    assert fz.STRATUM_FIELDS[len(cs.MANIFEST_FIELDS):] == fz.DERIVED_FIELDS
    assert not set(cs.MANIFEST_FIELDS) & set(fz.DERIVED_FIELDS)
    for row in frozen["rows"]:
        assert tuple(row) == fz.STRATUM_FIELDS


def test_regeneration_is_byte_identical(tmp_path, frozen):
    first = tmp_path / "first"
    second = tmp_path / "second"
    paths_a = fz.write(frozen, first)
    paths_b = fz.write(fz.build(fz.load_source()), second)
    for a, b in zip(paths_a, paths_b, strict=True):
        assert a.read_bytes() == b.read_bytes(), f"{a.name} is not reproducible"


def test_the_committed_artifact_is_what_the_script_produces_now(tmp_path, frozen):
    """The audit the document promises — 'regenerate it and diff' — actually run."""
    fresh_json, fresh_csv = fz.write(frozen, tmp_path)
    assert fresh_json.read_bytes() == (PREREG / f"{fz.OUT_STEM}.json").read_bytes()
    assert fresh_csv.read_bytes() == (PREREG / f"{fz.OUT_STEM}.csv").read_bytes()


def test_the_artifact_carries_no_generation_timestamp(committed):
    """A timestamp of its own generation would make the artifact non-reproducible for no gain: the
    run it derives from is already stamped and the freeze date belongs to the commit."""
    header = committed["header"]
    assert "generated_utc" not in header, "the freeze must not stamp its own run time"
    assert header["source"]["generated_utc"] == fz.SOURCE_GENERATED_UTC
    assert header["frozen_date"] == fz.FREEZE_DATE == "2026-08-16"
    # The source run's stamp is the only timestamp anywhere in the artifact.
    blob = json.dumps(committed, ensure_ascii=False)
    assert blob.count(fz.SOURCE_GENERATED_UTC) == 1
    assert not re.search(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?!\+00:00)", blob)


def test_check_mode_passes_against_the_committed_tree():
    assert fz.main(["--check"]) == 0


def test_the_driver_hardcodes_no_absolute_path():
    """Everything resolves from the file's own location; a machine-specific path in a
    pre-registration driver makes it unrunnable by the reader it exists for."""
    text = SCRIPT.read_text(encoding="utf-8")
    assert not re.search(r"[A-Za-z]:[\\/]", text), "a Windows drive-letter path is hardcoded"
    assert not re.search(r"(?m)^\s*[\"'](/home/|/Users/|/mnt/)", text)


# ---------------------------------------------------------------------------
# The guards themselves
# ---------------------------------------------------------------------------


def test_a_wrong_declared_stratum_count_aborts_the_freeze(monkeypatch, raw_source):
    """The declarations are verified, not trusted: this proves the verification is live rather than
    decorative, which is the failure mode `gate_config`'s own history warns about."""
    wrong = tuple(
        {**d, "n_strata": d["n_strata"] + 1} if d["rank"] == 1 else dict(d)
        for d in fz.FROZEN_DATASETS
    )
    monkeypatch.setattr(fz, "FROZEN_DATASETS", wrong)
    with pytest.raises(fz.FrozenDeclarationMismatch, match="candidate strata"):
        fz.build(raw_source)


@pytest.mark.parametrize("index", range(4))
def test_a_wrong_layer_b_declaration_aborts_the_freeze(monkeypatch, raw_source, index):
    wrong = tuple(
        {**tier, "dataset_ids": tier["dataset_ids"][:-1]} if i == index else tier
        for i, tier in enumerate(fz.LAYER_B)
    )
    monkeypatch.setattr(fz, "LAYER_B", wrong)
    with pytest.raises(fz.FrozenDeclarationMismatch, match="Layer B"):
        fz.build(raw_source)


def test_dropping_a_layer_b_tier_aborts_the_freeze(monkeypatch, raw_source):
    """The defect this guard exists for: a tier that is not declared is a tier nobody checks, and
    §6's table would go on publishing a row for it."""
    monkeypatch.setattr(fz, "LAYER_B", fz.LAYER_B[1:])
    with pytest.raises(fz.FrozenDeclarationMismatch, match="Every tier of the envelope"):
        fz.build(raw_source)


def test_a_wrong_manifest_tier_census_aborts_the_freeze(monkeypatch, raw_source):
    """§6's counterfactual is a load-bearing retraction, so its inputs are verified like any other
    declaration rather than being prose someone typed once."""
    monkeypatch.setattr(fz, "MANIFEST_TIER_CENSUS", {**fz.MANIFEST_TIER_CENSUS, 13: 4})
    with pytest.raises(fz.FrozenDeclarationMismatch, match="candidate-bearing datasets"):
        fz.build(raw_source)


def test_a_wrong_counterfactual_maximum_aborts_the_freeze(monkeypatch, raw_source):
    monkeypatch.setattr(
        fz, "COUNTERFACTUAL_MAX_SURVIVORS", {**fz.COUNTERFACTUAL_MAX_SURVIVORS, 23: 3}
    )
    with pytest.raises(fz.FrozenDeclarationMismatch, match="could have kept"):
        fz.build(raw_source)


def test_a_wrong_sibling_stratum_count_aborts_the_freeze(monkeypatch, raw_source):
    """The controls are frozen, so their declarations are verified exactly like the twelve's."""
    wrong = tuple(
        {**s, "n_strata": s["n_strata"] + 1} if i == 0 else dict(s)
        for i, s in enumerate(fz.SIBLING_DATASETS)
    )
    monkeypatch.setattr(fz, "SIBLING_DATASETS", wrong)
    with pytest.raises(fz.FrozenDeclarationMismatch, match="candidate strata"):
        fz.build(raw_source)


def test_a_missing_source_csv_aborts_the_freeze(tmp_path):
    """An unverifiable evidence file is not evidence, so a named-and-absent CSV twin must stop the
    freeze rather than be skipped. Skipping it would turn deleting the file into a way of passing
    the hash check."""
    with pytest.raises(fz.SourceArtifactMismatch, match="CSV twin not found"):
        fz.load_source(fz.SOURCE_JSON, csv_path=tmp_path / "absent.csv")


def test_a_tampered_source_csv_aborts_the_freeze(tmp_path):
    same_length = fz.SOURCE_CSV.read_bytes().replace(b"candidate", b"candidatf", 1)
    path = tmp_path / fz.SOURCE_CSV.name
    path.write_bytes(same_length)
    assert path.stat().st_size == fz.SOURCE_CSV_BYTES
    with pytest.raises(fz.SourceArtifactMismatch, match="sha256"):
        fz.load_source(fz.SOURCE_JSON, csv_path=path)


def test_a_sibling_promoted_into_the_frozen_set_aborts_the_freeze(monkeypatch, raw_source):
    """The D2 denominator must not be able to acquire a second dataset from one collection by an
    edit that looks like a routine substitution. The list stays twelve long and correctly ranked
    here, so the only thing that can catch it is the sibling check itself."""
    promoted = tuple(
        {**d, "dataset_id": fz.SIBLING_DATASETS[0]["dataset_id"]} if d["rank"] == 2 else dict(d)
        for d in fz.FROZEN_DATASETS
    )
    monkeypatch.setattr(fz, "FROZEN_DATASETS", promoted)
    with pytest.raises(fz.FrozenDeclarationMismatch, match="independent dataset"):
        fz.build(raw_source)


def test_an_unlabelled_disease_arm_aborts_the_freeze(monkeypatch, raw_source):
    incomplete = {k: v for k, v in fz.EXPECTED_EFFECT.items()
                  if k != ("2a498ace-872a-4935-984b-1afa70fd9886", "post-COVID-19 disorder")}
    monkeypatch.setattr(fz, "EXPECTED_EFFECT", incomplete)
    with pytest.raises(fz.FrozenDeclarationMismatch, match="does not label"):
        fz.build(raw_source)


# ---------------------------------------------------------------------------
# The document and the artifact must agree
# ---------------------------------------------------------------------------


def test_the_document_names_every_frozen_dataset_and_both_siblings():
    text = DOCUMENT.read_text(encoding="utf-8")
    for dataset in fz.FROZEN_DATASETS:
        assert dataset["dataset_id"] in text, dataset["short"]
    for sibling in fz.SIBLING_DATASETS:
        assert sibling["dataset_id"] in text


def test_the_document_records_the_source_hashes_and_the_frozen_count(frozen):
    text = _document_text()
    assert fz.SOURCE_SHA256 in text
    assert fz.SOURCE_CSV_SHA256 in text
    assert fz.PROPOSAL_SHA256 in text
    assert fz.CI_RUN_ID in text
    assert str(frozen["header"]["n_frozen_strata"]) in text
    assert "docs/AMENDMENTS.md" in text or "AMENDMENTS.md" in text


def test_the_proposal_document_is_committed_with_the_pinned_hash(frozen):
    """The one act of discretion in this freeze — the choice of twelve — is justified in that file.

    It also holds the five reserves §9 item 7 forbids substituting from, and a list of candidate
    replacements that lives outside the record cannot be checked against a substitution that
    happened. Pinning it is what makes such a substitution detectable.
    """
    assert fz.PROPOSAL_MD.exists(), f"{fz.PROPOSAL_MD} is not committed"
    assert fz.PROPOSAL_MD.stat().st_size == fz.PROPOSAL_BYTES
    assert hashlib.sha256(fz.PROPOSAL_MD.read_bytes()).hexdigest() == fz.PROPOSAL_SHA256
    block = frozen["header"]["reasoning_document"]
    assert block["file"] == fz.PROPOSAL_MD.name
    assert block["sha256"] == fz.PROPOSAL_SHA256
    assert "NOT part of the binding act" in block["status"]


def test_a_missing_proposal_document_aborts_the_freeze(monkeypatch, tmp_path):
    monkeypatch.setattr(fz, "PROPOSAL_MD", tmp_path / "absent.md")
    with pytest.raises(fz.SourceArtifactMismatch, match="proposal document not found"):
        fz.load_source()


def test_the_two_open_decisions_that_could_move_the_frozen_set_say_so():
    """§9.4 and §9.5 would each change what the 251 are, and §9.5 would move Layer B. An open
    question that can silently redefine the frozen set is not an open question."""
    text = " ".join(_document_text().split())
    assert "Choosing that harmonised level is an amendment" in text
    assert "Deciding it is an amendment, because it changes what the 251 are" in text
    # §9.8 records that the shipped moderated arm fits ~ 1 + x with no covariate slot, so §9.5 may
    # not offer region-as-covariate as a live option in the same section.
    assert "no covariate slot" in text
    assert '"Or the covariate set" is struck from this item' in text


def test_the_document_does_not_claim_admission():
    """A pre-registration that reads as an admission is the one mistake that cannot be undone by a
    later correction, because the sweep would already have been run on it."""
    text = DOCUMENT.read_text(encoding="utf-8")
    assert "`admitted_to_sweep = False`" in text
    assert "admitted_to_sweep = True" not in text
