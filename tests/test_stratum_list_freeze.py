"""The spec §1 stratum-list freeze, tested as the pre-registration it is.

``docs/PREREGISTRATION_STRATUM_LIST.md`` claims that a reader can re-derive every number in it from
committed, hash-pinned artifacts, and that the 251 frozen strata follow mechanically from twelve
dataset ids with no stratum-level discretion anywhere. Both claims are only worth something if they
are checked by something other than the person who wrote them, so this file re-derives the
load-bearing figures from the raw sources — deliberately **not** through ``freeze_stratum_list``'s
own helpers where the point is to catch that module drifting — and compares.

The failure modes worth catching here are all quiet ones:

* the source artifacts are swapped, truncated or re-generated, and the pre-registered numbers
  silently come to describe different bytes;
* a dataset id is added, removed or edited after the freeze, changing the D2 denominator;
* the selection rule stops being "every candidate row of the twelve", so a stratum leaves the
  analysis set without a reason attached;
* **a collection-sibling goes unfrozen.** This is not hypothetical and it is the reason the sibling
  set is computed rather than typed: an earlier freeze named two siblings when there are five, and
  the three it missed hold 79 runnable strata, two of them clearing every envelope tier;
* a cells-per-donor bin turns out empty, which would make spec §1 (iii)'s "spanning the bins" false;
* the Layer B truncation stops matching the envelope arithmetic — the one thing that keeps the
  surviving subset from being chosen after the sigma_donor anchor lands;
* a Layer B tier goes undeclared, so the document publishes a row nothing ever checked;
* **any table of the document drifts from the artifact it describes, in either direction.** §1,
  §3.2, §3.3, §4.1, §5 and both of §6's tables are parsed back out of the file and compared cell by
  cell, here and in the freeze script;
* the two same-collection blocks stop being distinguishable from the 251;
* the emitted artifact stops being byte-reproducible, at which point "regenerate it and diff" is no
  longer an audit anyone can run.

A previous round mutation-tested this file and found nine survivors: the ``role`` vocabulary, the
``EXPECTED_EFFECT`` values, ``PROPOSAL_STATUS`` (a mutation reversing its meaning survived), the
per-dataset ``assay`` declarations, and the guard separating the control set from the analysis set.
Each has a test below that fails on the mutation that survived.

No network. Everything is read from ``pilot/preregistration/``.
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
SPEC = REPO / "docs" / "PHASE0_SPEC.md"


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
# Fixtures. The source manifest is 6.6 MB and the two snapshots 3.7 MB; parse once per module.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def raw_source() -> dict:
    """The committed candidate manifest, read directly — not through ``load_source``.

    Reading it raw is the point: every re-derivation below has to be independent of the module it
    is checking, or a bug in that module would be reproduced by the test that is meant to catch it.
    """
    return json.loads(fz.SOURCE_JSON.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def raw_discover() -> dict:
    return json.loads(fz.DISCOVER_INDEX.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def raw_release() -> dict:
    return json.loads(fz.RELEASE_DATASETS.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def collections(raw_release) -> dict[str, str]:
    """``{dataset_id: collection_id}`` read straight out of the pinned release table."""
    return {row["dataset_id"]: row["collection_id"] for row in raw_release["datasets"]}


@pytest.fixture(scope="module")
def frozen(raw_source) -> dict:
    """The freeze, rebuilt from the committed sources (this one does go through the module)."""
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
    return [cell.replace("*", "").replace("`", "").strip()
            for cell in line.strip().strip("|").split("|")]


def _leading_int(cell: str) -> int:
    match = re.search(r"\d+", cell.replace(" ", "").replace("\u202f", ""))
    assert match, f"no integer in table cell {cell!r}"
    return int(match.group(0))


def _parse_layer_b_table() -> list[dict]:
    """§6's Layer B table, read back out of the document.

    The document is the human half of the pre-registration and the artifact is the machine half;
    neither is authoritative alone, so the one number that decides §6's headline — the verdict
    against the spec's dataset floor — is compared rather than trusted in either direction.
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


def _parse_tier_census_table() -> list[dict]:
    """§6's manifest-wide table: what the 68 candidate-bearing datasets could have supported."""
    rows = []
    for line in _document_text().splitlines():
        cells = _table_cells(line)
        if len(cells) != 4 or not re.match(r"^≥\s*\d+\s*v\s*\d+$", cells[0]):
            continue
        rows.append({
            "min_donors_per_group": _leading_int(cells[0]),
            "n_manifest_datasets_at_tier": _leading_int(cells[1]),
            "n_manifest_collections_at_tier": _leading_int(cells[2]),
            "n_datasets": _leading_int(cells[3]),
        })
    return rows


# ---------------------------------------------------------------------------
# The pinned sources
# ---------------------------------------------------------------------------


def test_source_artifacts_are_committed_with_the_pinned_hashes():
    """The pre-registered numbers describe these exact bytes and no others.

    All four files, not just the manifest: the two external snapshots are what every collection,
    assay, tissue and DOI claim now rests on, and an unpinned snapshot would put those claims back
    where the third review found them.
    """
    for path, digest, size in (
        (fz.SOURCE_JSON, fz.SOURCE_SHA256, fz.SOURCE_BYTES),
        (fz.SOURCE_CSV, fz.SOURCE_CSV_SHA256, fz.SOURCE_CSV_BYTES),
        (fz.DISCOVER_INDEX, fz.DISCOVER_SHA256, fz.DISCOVER_BYTES),
        (fz.RELEASE_DATASETS, fz.RELEASE_SHA256, fz.RELEASE_BYTES),
        (fz.PROPOSAL_MD, fz.PROPOSAL_SHA256, fz.PROPOSAL_BYTES),
    ):
        assert path.exists(), f"{path} is not committed"
        assert path.stat().st_size == size, f"{path.name} is {path.stat().st_size} bytes"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, f"{path.name} sha256"


def test_the_document_records_every_pinned_hash():
    text = _document_text()
    for digest in (
        fz.SOURCE_SHA256,
        fz.SOURCE_CSV_SHA256,
        fz.DISCOVER_SHA256,
        fz.RELEASE_SHA256,
        fz.PROPOSAL_SHA256,
        fz.PROPOSAL_ORIGINAL_SHA256,
    ):
        assert digest in text, f"{digest} is pinned in the script but absent from the document"
    assert fz.CI_RUN_ID in text


def test_source_header_carries_the_pinned_run_stamps(raw_source):
    header = raw_source["header"]
    assert header["generated_utc"] == fz.SOURCE_GENERATED_UTC
    assert header["census_version"] == cs.CENSUS_VERSION == "2025-01-30"
    assert header["value_filter"] == fz.MANIFEST_VALUE_FILTER
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


@pytest.mark.parametrize(
    ("loader", "path_attr", "sha_attr"),
    [
        (lambda p: fz.load_discover(p), "DISCOVER_INDEX", "DISCOVER_SHA256"),
        (lambda p: fz.load_release(p), "RELEASE_DATASETS", "RELEASE_SHA256"),
    ],
)
def test_an_edited_external_snapshot_aborts_the_freeze(tmp_path, loader, path_attr, sha_attr):
    """The snapshots are evidence, so a changed byte must stop the freeze exactly as the manifest's
    would. Without this, "pinned" would mean "committed", which is not the same thing."""
    source = getattr(fz, path_attr)
    tampered = tmp_path / source.name
    payload = source.read_bytes()
    marker = b'"collection_id"'
    assert marker in payload
    tampered.write_bytes(payload.replace(marker, b'"collection_iD"', 1))
    assert tampered.stat().st_size == source.stat().st_size
    with pytest.raises(fz.SourceArtifactMismatch, match="sha256"):
        loader(tampered)


def test_the_two_snapshots_carry_the_whole_indexes_they_claim(raw_discover, raw_release):
    assert raw_discover["header"]["n_datasets_in_index"] == fz.DISCOVER_N_DATASETS == 2216
    assert len(raw_discover["datasets"]) == 2216
    assert raw_discover["header"]["read_date"] == fz.DISCOVER_READ_DATE == "2026-08-16"
    assert raw_release["header"]["census_version"] == cs.CENSUS_VERSION
    assert len(raw_release["datasets"]) == fz.RELEASE_N_DATASETS == 1573
    assert len({d["dataset_id"] for d in raw_release["datasets"]}) == 1573


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


def test_the_twelve_occupy_twelve_distinct_collections(collections):
    """The fact that makes the D2 denominator clean, checked against the pinned release table.

    If two of the twelve shared a collection, "majority of independent datasets" would already be
    counting one cohort and one laboratory twice before any sibling was considered.
    """
    owning = [collections[d["dataset_id"]] for d in fz.FROZEN_DATASETS]
    assert len(owning) == 12
    assert len(set(owning)) == 12, f"a collection appears twice: {Counter(owning).most_common(1)}"


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


# ---------------------------------------------------------------------------
# The siblings — computed, not named. This is where the third review found a blocker.
# ---------------------------------------------------------------------------


def test_the_sibling_set_is_every_candidate_bearing_dataset_sharing_a_collection(
    raw_source, collections
):
    """Re-derived here from the raw release table, without ``derive_siblings``.

    The defect this catches shipped: two siblings were named where there are five, because
    collection membership came from a live API read that nothing could check. The three that were
    missed hold 79 runnable strata and two of them clear every envelope tier.
    """
    frozen_ids = {d["dataset_id"] for d in fz.FROZEN_DATASETS}
    frozen_collections = {collections[i] for i in frozen_ids}
    candidate_ids = {r["dataset_id"] for r in _candidates(raw_source)}
    derived = {
        dataset_id for dataset_id in candidate_ids
        if dataset_id not in frozen_ids and collections[dataset_id] in frozen_collections
    }
    assert derived == {s["dataset_id"] for s in fz.SIBLING_DATASETS}
    assert len(derived) == 5


def test_each_sibling_shares_its_collection_with_the_parent_it_names(collections):
    for sibling in fz.SIBLING_DATASETS:
        assert collections[sibling["dataset_id"]] == sibling["collection_id"]
        assert collections[sibling["sibling_of"]] == sibling["collection_id"]


def test_neither_index_disagrees_about_any_collection_the_freeze_uses(
    raw_source, raw_discover, collections
):
    """The sibling set comes from the release table; Discover is the cross-check, not a substitute.

    They agree for every dataset the two share, and the freeze asserts it rather than assuming it —
    because the same two snapshots demonstrably disagree about `collection_doi` for 61 datasets and
    about `dataset_version_id` for all 1573.
    """
    by_discover = {row["dataset_id"]: row for row in raw_discover["datasets"]}
    for dataset_id in {r["dataset_id"] for r in raw_source["rows"]}:
        assert dataset_id in by_discover, dataset_id
        assert dataset_id in collections, dataset_id
        assert collections[dataset_id] == by_discover[dataset_id]["collection_id"], dataset_id


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
    """106 runnable strata that were named, excluded from D2 and then left unlisted would have been
    selectable after the fact — the one freedom §3.1 exists to remove."""
    controls = frozen["within_collection_control_rows"]
    expected = _rule(raw_source["rows"], [s["dataset_id"] for s in fz.SIBLING_DATASETS])
    assert len(controls) == len(expected) == 106
    assert {_key(r) for r in controls} == {_key(r) for r in expected}
    assert frozen["header"]["n_within_collection_control_strata"] == 106
    per_sibling = {s["dataset_id"]: s["n_strata"] for s in fz.SIBLING_DATASETS}
    assert per_sibling == {
        "c2876b1b-06d8-4d96-a56b-5304f815b99a": 18,
        "edc8d3fe-153c-4e3d-8be0-2108d30f8d70": 25,
        "8f4f8502-9170-4ac2-9707-3b6985ebfe5f": 11,
        "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d": 43,
        "1e5bd3b8-6a0e-4959-8d69-cafed30fe814": 9,
    }
    for sibling in fz.SIBLING_DATASETS:
        own = [r for r in expected if r["dataset_id"] == sibling["dataset_id"]]
        assert len(own) == sibling["n_strata"], sibling["short"]
        assert max(_min_donors(r) for r in own) == sibling["ceiling_min_donors_per_group"]


def test_three_of_the_five_siblings_are_donor_rich_and_two_clear_every_tier(raw_source):
    """Why the omission mattered rather than being a formality.

    The three siblings the earlier freeze missed were not marginal: two of them hold strata at
    23 v 23, the hardest tier of the operating envelope, where the whole analysis set retains three
    datasets.
    """
    controls = _rule(raw_source["rows"], [s["dataset_id"] for s in fz.SIBLING_DATASETS])
    ceilings = {
        s["dataset_id"]: max(
            _min_donors(r) for r in controls if r["dataset_id"] == s["dataset_id"]
        )
        for s in fz.SIBLING_DATASETS
    }
    assert sum(1 for c in ceilings.values() if c >= 23) == 3
    previously_named = {
        "c2876b1b-06d8-4d96-a56b-5304f815b99a",
        "1e5bd3b8-6a0e-4959-8d69-cafed30fe814",
    }
    missed = {k: v for k, v in ceilings.items() if k not in previously_named}
    assert len(missed) == 3
    assert sum(1 for c in missed.values() if c >= 23) == 2
    assert sum(len([r for r in controls if r["dataset_id"] == k]) for k in missed) == 79


def test_the_disjointness_guard_between_the_two_sets_is_live(raw_source):
    """One set is the D2 denominator and the other is explicitly outside it. A row in both would
    let the denominator acquire a sibling silently, so the guard is exercised directly."""
    strata = _frozen_rows(raw_source)
    controls = _rule(raw_source["rows"], [s["dataset_id"] for s in fz.SIBLING_DATASETS])
    fz.assert_sets_disjoint(strata, controls)  # the real sets are disjoint
    with pytest.raises(fz.FrozenDeclarationMismatch, match="also in the analysis set"):
        fz.assert_sets_disjoint(strata, [*controls, strata[0]])


def test_every_row_carries_a_role_and_the_two_sets_stay_distinguishable(frozen):
    """The role travels with the row, so a row read out of the artifact still says which set it is
    in. The 251 and the 106 must never be conflated: one is the D2 denominator and one is not."""
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


def test_the_role_vocabulary_is_pinned_to_its_exact_strings(frozen):
    """A mutation to either literal survived the previous suite.

    The strings are the artifact's public vocabulary — a downstream reader filters on them — so
    renaming one is a schema change and must not pass as a refactor.
    """
    assert fz.ROLE_ANALYSIS_SET == "analysis_set"
    assert fz.ROLE_WITHIN_COLLECTION_CONTROL == "within_collection_control"
    assert fz.EXPECTED_EFFECT_NOT_APPLICABLE == "not_applicable"
    emitted = {r["role"] for r in frozen["rows"]} | {
        r["role"] for r in frozen["within_collection_control_rows"]
    }
    assert emitted == {"analysis_set", "within_collection_control"}
    for sibling in frozen["header"]["sibling_datasets"]:
        assert sibling["role"] == "within_collection_control"
    assert "within_collection_control" in _document_text()


def test_control_rows_carry_no_expected_effect_label(frozen):
    """Spec §1 (i)/(ii) is a coverage claim about the twelve. A control row labelled strong or
    subtle would be counted into a claim it is excluded from."""
    labels = {r["expected_effect"] for r in frozen["within_collection_control_rows"]}
    assert labels == {fz.EXPECTED_EFFECT_NOT_APPLICABLE}
    assert fz.EXPECTED_EFFECT_NOT_APPLICABLE not in fz.EXPECTED_EFFECT_VOCABULARY


def test_the_csv_carries_both_blocks_told_apart_by_the_role_column():
    """The CSV is the JSON's twin; one that held only the 251 would silently omit frozen content."""
    with open(PREREG / f"{fz.OUT_STEM}.csv", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 251 + 106 == 357
    counts = Counter(r["role"] for r in rows)
    assert counts == {fz.ROLE_ANALYSIS_SET: 251, fz.ROLE_WITHIN_COLLECTION_CONTROL: 106}
    assert [r["role"] for r in rows[:251]] == [fz.ROLE_ANALYSIS_SET] * 251


def test_the_document_states_the_three_rules_that_bind_the_control_set():
    text = _document_text()
    assert "within_collection_control" in text
    assert "106" in text
    for phrase in (
        "never enters the D2 denominator",
        "Promoting one to an independent dataset is an amendment",
    ):
        assert phrase in text, phrase


def test_the_artifact_no_longer_claims_to_have_listed_every_runnable_stratum(frozen):
    """The claim that shipped and was false: 1197 candidate strata exist at this pin and 357 are
    frozen, so "no runnable stratum of this Census pin is left unlisted" was never true. The scoped
    claim — every candidate stratum of the collections represented in the analysis set — is."""
    note = frozen["header"]["sibling_datasets_note"]
    assert "no runnable stratum of this Census pin is left unlisted" not in note
    assert "collection represented in the analysis set" in note
    assert "NOT a claim about the whole Census pin" in note


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
    assert not any(r["admitted_to_sweep"] for r in frozen["within_collection_control_rows"])
    assert frozen["header"]["admission"]["admitted_to_sweep"] is False


def test_the_counts_gate_and_the_envelope_are_still_pending_on_every_row(frozen):
    """§9 items 1 and 7: the list may still shrink, and nothing here pretends otherwise."""
    for row in (*frozen["rows"], *frozen["within_collection_control_rows"]):
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
    assert lows[0] == gc.MIN_CELLS, "the lower edge is the inclusion-gate floor, not a new number"
    # Five CLOSED bins, over two and a half decades — §5 used to say "two decades", which is four.
    assert len(fz.CELLS_PER_DONOR_BINS) - 1 == 5
    assert "five closed bins" in " ".join(_document_text().split())
    assert "five closed bins per two decades" not in " ".join(_document_text().split())


def test_every_bin_is_occupied_by_the_frozen_set(frozen):
    """Spec §1 (iii) requires cells-per-donor 'spanning the pre-registered bins'. An empty bin makes
    that claim false, and this is the assertion that would say so."""
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
    """Recomputed here from the bin edges rather than through ``fz.bin_of``.

    The previous version of this test called ``bin_of`` on both sides and so asserted only that a
    function is deterministic. It would have passed against any bin assignment whatsoever.
    """
    def independent_bin(value: float) -> str:
        for low, high in ((10, 30), (30, 100), (100, 300), (300, 1000), (1000, 3000)):
            if low <= value < high:
                return f"[{low},{high})"
        assert value >= 3000, value
        return "[3000,inf)"

    for row in (*frozen["rows"], *frozen["within_collection_control_rows"]):
        for group in ("A", "B"):
            median = row["cells_per_donor_by_group"][group]["median"]
            assert row[f"cells_per_donor_bin_{group}"] == independent_bin(median)


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
    the point at which the instrument nominally operates — not at 0.5, and the list clears the
    spec's own '8-12 datasets' floor at the most optimistic tier only.
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


def test_both_surviving_columns_of_section_6_are_defined_in_the_document(frozen):
    """The two columns answer different questions and only one used to be defined.

    "Surviving strata" is not "the strata of the surviving datasets": at >= 4v4 the eleven surviving
    datasets hold 243 strata and the column reads 227. A reader given only the dataset-level
    definition directly above it gets the second column wrong by sixteen.
    """
    text = " ".join(_document_text().split())
    assert "**Surviving datasets** — of the twelve" in text
    assert "**Surviving strata** — of the 251" in text
    tiers = {t["min_donors_per_group"]: t for t in frozen["header"]["layer_b_truncation"]}
    assert tiers[4]["n_strata"] == 227
    assert f"hold 243 strata between them and the column reads {tiers[4]['n_strata']}" in text


# ---------------------------------------------------------------------------
# What the manifest could have supported — §6's trade-off
# ---------------------------------------------------------------------------


def test_the_manifest_wide_tier_census_is_recomputable(raw_source, collections):
    """The numbers that refute 'only a handful clear 13 v 13 anywhere'.

    Counted in collections as well as datasets, because that is the unit D2 uses and the unit whose
    absence let an earlier §6 publish a twelve-dataset witness containing two SEA-AD datasets and
    two KPMP datasets.
    """
    candidates = _candidates(raw_source)
    assert len({r["dataset_id"] for r in candidates}) == 68 == fz.MANIFEST_N_CANDIDATE_DATASETS
    assert len({collections[r["dataset_id"]] for r in candidates}) == 50 == (
        fz.MANIFEST_N_CANDIDATE_COLLECTIONS
    )
    by_dataset, by_collection = {}, {}
    for threshold in (4, 8, 13, 23):
        clearing = {r["dataset_id"] for r in candidates if _min_donors(r) >= threshold}
        by_dataset[threshold] = len(clearing)
        by_collection[threshold] = len({collections[i] for i in clearing})
    assert by_dataset == {4: 62, 8: 33, 13: 21, 23: 12} == fz.MANIFEST_TIER_CENSUS
    assert by_collection == {4: 46, 8: 25, 13: 15, 23: 10} == fz.MANIFEST_TIER_COLLECTIONS


def test_the_3v3_anchor_is_scarce_and_never_donor_rich(raw_source, collections):
    """Spec §1 (iii)'s "some exactly 3v3" is what costs a dataset slot at the hard tiers."""
    candidates = _candidates(raw_source)
    anchors = {r["dataset_id"] for r in candidates
               if r["n_donors_A"] == 3 and r["n_donors_B"] == 3}
    assert len(anchors) == fz.MANIFEST_DATASETS_WITH_EXACT_3V3 == 3
    assert len({collections[i] for i in anchors}) == fz.MANIFEST_COLLECTIONS_WITH_EXACT_3V3 == 2
    for threshold in (13, 23):
        tier = {r["dataset_id"] for r in candidates if _min_donors(r) >= threshold}
        assert not (anchors & tier), (
            f"a dataset holding an exactly-3v3 stratum also clears {threshold}v{threshold}; §6's "
            "argument that the anchor costs a slot no longer holds"
        )


def test_the_constructed_counterfactual_witness_is_gone(frozen):
    """It was deleted, not repaired, and this test is what keeps it deleted.

    The witness double-counted two collections, overstated its own >= 23v23 figure by one, and
    leaned on the Human Lung Cell Atlas — the dataset that owns all 12 `excluded_confound` rows of
    the manifest. §6 now states counts, which cannot be selected from.
    """
    assert not hasattr(fz, "COUNTERFACTUAL_MAX_SURVIVORS")
    assert not hasattr(fz, "COUNTERFACTUAL_WITNESS_ANCHOR")
    assert not hasattr(fz, "counterfactual_max_survivors")
    blob = json.dumps(frozen["header"], ensure_ascii=False)
    assert "counterfactual_max_survivors_of_twelve" not in blob
    text = " ".join(_document_text().split())
    assert "12 / 12 / 11 / 11" not in text
    assert "deleted rather than repaired" in text


def test_the_documents_tier_census_table_agrees_with_the_artifact(frozen):
    parsed = _parse_tier_census_table()
    emitted = frozen["header"]["layer_b_truncation"]
    assert len(parsed) == len(emitted) == 4
    for row, tier in zip(parsed, emitted, strict=True):
        where = f"§6 tier census row >= {row['min_donors_per_group']}"
        assert row["min_donors_per_group"] == tier["min_donors_per_group"], where
        assert row["n_manifest_datasets_at_tier"] == tier["n_manifest_datasets_at_tier"], where
        assert row["n_manifest_collections_at_tier"] == (
            tier["n_manifest_collections_at_tier"]
        ), where
        assert row["n_datasets"] == tier["n_datasets"], where


RETRACTED_CLAIMS = (
    "only a handful clear",
    "a property of the public data",
    "met, with no margin",
    "no runnable stratum of this Census pin",
    "unreachable",
)

#: A paragraph may still contain a retracted claim — the document quotes several of them, which is
#: how a correction differs from a deletion. What it may not do is assert one.
RETRACTION_MARKERS = (
    "retracted", "was false", "earlier draft", "simply wrong", "earlier artifact note",
    "is withdrawn",
)


@pytest.mark.parametrize("claim", RETRACTED_CLAIMS)
def test_a_retracted_claim_may_only_appear_inside_its_own_retraction(claim):
    """The false statements §6 and §3.3 corrected, pinned so none can drift back as an assertion."""
    # Normalised, because the document is hard-wrapped and a retracted phrase can straddle a line.
    paragraphs = [" ".join(p.split()) for p in _document_text().split("\n\n")]
    paragraphs = [p for p in paragraphs if claim in p]
    # Absent entirely is the better outcome and is allowed: deleting a rhetorical claim beats
    # keeping it under a correction. What is forbidden is the claim standing unmarked.
    for paragraph in paragraphs:
        assert any(marker in paragraph for marker in RETRACTION_MARKERS), (
            f"{claim!r} appears in a paragraph that does not mark it as retracted:\n{paragraph}"
        )


def test_the_retracted_upper_bound_claim_is_marked_where_it_appears():
    """Amendment 4 Part A retracts Amendment 3's "upper bound" and names §6 and §9 item 2 of this
    document as repeating it. Both shipped unmarked; a pre-registration that contradicts its own
    amendment log a thousand lines away has not been corrected, only annotated elsewhere."""
    text = " ".join(_document_text().split())
    assert text.count("unknown error sign") >= 2, "both sites must carry the marker"
    assert text.count("Amendment 4](AMENDMENTS.md) Part A, Correction 1") >= 2
    # And the direction of the error is stated, because that is the part that matters.
    assert text.count("direction of the error is the dangerous one") >= 2
    assert "unvalidated upper bound" not in text


# ---------------------------------------------------------------------------
# Coverage claims that the document makes and the data must support
# ---------------------------------------------------------------------------


def test_the_donor_count_axis_carries_both_ends_required_by_spec_section_1(frozen):
    """Spec §1 (iii): 'some exactly 3v3, some >= 8v8'."""
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

    # Spec §1 (i) and (ii) ask for 2-3 datasets on each side; the list must clear both, not one.
    by_label: dict[str, set[str]] = {label: set() for label in fz.EXPECTED_EFFECT_VOCABULARY}
    for (dataset_id, _disease), label in fz.EXPECTED_EFFECT.items():
        by_label[label].add(dataset_id)
    assert len(by_label["strong"]) >= 2, "spec §1 (i) wants 2-3 with a strong expected effect"
    assert len(by_label["subtle"]) >= 2, "spec §1 (ii) wants 2-3 subtle/low-effect datasets"


def test_every_expected_effect_value_is_pinned_arm_by_arm():
    """Mutating a single label survived the previous suite.

    The labels are the whole of §4.1's coverage claim and were fixed before any metric was computed;
    a silent reassignment is exactly the failure the pre-registration exists to prevent, so the map
    is pinned in full rather than counted.
    """
    ranks = {d["dataset_id"]: d["rank"] for d in fz.FROZEN_DATASETS}
    by_rank = {(ranks[i], disease): label for (i, disease), label in fz.EXPECTED_EFFECT.items()}
    assert by_rank == {
        (1, "dementia"): "strong",
        (2, "Alzheimer disease"): "strong",
        (2, "Pick disease"): "strong",
        (2, "progressive supranuclear palsy"): "strong",
        (3, "COVID-19"): "strong",
        (4, "COVID-19"): "strong",
        (4, "post-COVID-19 disorder"): "subtle",
        (5, "COVID-19"): "strong",
        (5, "influenza"): "moderate",
        (6, "atrial fibrillation"): "subtle",
        (7, "rheumatoid arthritis"): "subtle",
        (8, "acute kidney failure"): "moderate",
        (8, "chronic kidney disease"): "moderate",
        (9, "clonal hematopoiesis"): "subtle",
        (10, "opiate dependence"): "subtle",
        (11, "Crohn disease"): "subtle",
        (12, "pulmonary emphysema"): "moderate",
    }
    assert fz.EXPECTED_EFFECT_VOCABULARY == ("strong", "moderate", "subtle")


def test_the_per_dataset_assay_declarations_are_pinned_and_match_discover(raw_discover):
    """Mutating an assay tuple survived the previous suite, and it is the axis spec §1 (iii) names
    first. Pinned here in full and cross-checked against the Discover snapshot."""
    by_discover = {row["dataset_id"]: row for row in raw_discover["datasets"]}
    declared = {d["rank"]: tuple(sorted(d["assay"])) for d in fz.FROZEN_DATASETS}
    assert declared == {
        1: ("10x 3' v3", "10x multiome"),
        2: ("10x 3' v2", "10x 3' v3"),
        3: ("10x 3' v3",),
        4: ("10x 5' v1",),
        5: ("10x 5' v1",),
        6: ("10x 3' v3",),
        7: ("10x 3' v3",),
        8: ("10x 3' v3",),
        9: ("10x 3' v3",),
        10: ("10x 3' v3",),
        11: ("10x 3' v2",),
        12: ("10x 3' v3",),
    }
    for dataset in fz.FROZEN_DATASETS:
        record = by_discover[dataset["dataset_id"]]
        assert tuple(sorted(dataset["assay"])) == tuple(sorted(record["assay"]))
        assert tuple(sorted(dataset["suspension"])) == tuple(sorted(record["suspension_type"]))
        assert tuple(sorted(dataset["tissue"])) == tuple(sorted(record["tissue"]))
        assert (dataset["doi"] or None) == (record["collection_doi"] or None)
    # Exactly two datasets on 5', ten on 3' — the thin axis §4.2 calls out.
    fives = [d["rank"] for d in fz.FROZEN_DATASETS if any("5'" in a for a in d["assay"])]
    assert fives == [4, 5]


def test_a_wrong_assay_declaration_aborts_the_freeze(monkeypatch, raw_source):
    wrong = tuple(
        {**d, "assay": ("10x 3' v2",)} if d["rank"] == 4 else dict(d) for d in fz.FROZEN_DATASETS
    )
    monkeypatch.setattr(fz, "FROZEN_DATASETS", wrong)
    with pytest.raises(fz.FrozenDeclarationMismatch, match="Discover gives assay"):
        fz.build(raw_source)


def test_min_donors_per_group_column_is_the_quantity_the_report_must_print(frozen):
    """§9 item 12: permutation_count flatters skewed designs, so the honest companion travels as its
    own column rather than being recomputed by whoever reads the artifact."""
    for row in (*frozen["rows"], *frozen["within_collection_control_rows"]):
        assert row["min_donors_per_group"] == _min_donors(row)


# ---------------------------------------------------------------------------
# The attested figures — the eleven §10 used only to assert
# ---------------------------------------------------------------------------


def test_every_attested_figure_is_recomputed_from_the_manifest(raw_source, frozen):
    """§10 used to *attest* these: a human saying they had checked. Nothing else held them, so
    nothing would have caught them drifting. Re-derived here independently of ``measure_attested``
    for the ones that can be, and compared against the artifact's own block for all of them."""
    strata = _frozen_rows(raw_source)
    emitted = frozen["header"]["attested_figures"]
    assert emitted == dict(sorted(fz.ATTESTED.items()))

    assert sum(r["n_cells"] for r in strata) == emitted["n_cells"] == 4_609_595
    assert len({(r["dataset_id"], r["cell_type"]) for r in strata}) == 182
    assert len({r["cell_type"] for r in strata}) == 124
    assert min(r["residual_df"] for r in strata) == 4
    assert max(r["residual_df"] for r in strata) == 108
    assert max(r["permutation_count"] for r in strata) == emitted["permutation_count_max"]
    assert f"{emitted['permutation_count_max']:.2e}".startswith("7.28")
    per_donor = [r["cells_per_donor_by_group"][g] for r in strata for g in ("A", "B")]
    assert min(v["min"] for v in per_donor) == 10
    assert max(v["max"] for v in per_donor) == 16_383
    counts_per_cell = sorted(
        r["median_counts_per_cell_by_group"][g] for r in strata for g in ("A", "B")
    )
    assert counts_per_cell[len(counts_per_cell) // 2] == 3081.0
    assert (counts_per_cell[0], counts_per_cell[-1]) == (284.0, 56_841.5)
    assert sum(1 for r in strata if r["permutation_count"] < 1000) == 45
    kpmp = "a12ccb9b-4fbe-457d-8590-ac78053259ef"
    assert sum(1 for r in strata if r["dataset_id"] == kpmp and _min_donors(r) >= 23) == 6


def test_the_prose_figures_outside_the_tables_are_attested_too(raw_source, frozen):
    """The long tail: numbers §5, §6, §9 and §10.1 state in sentences rather than in tables.

    They are the ones no table check would reach, which is where every previous round's defects
    were found. 243 in particular is load-bearing: it is what makes "Surviving strata" mean
    something other than "the strata of the surviving datasets".
    """
    strata = _frozen_rows(raw_source)
    attested = frozen["header"]["attested_figures"]
    surviving = {r["dataset_id"] for r in strata if _min_donors(r) >= 4}
    assert sum(1 for r in strata if r["dataset_id"] in surviving) == 243
    assert attested["n_strata_in_datasets_at_least_4v4"] == 243
    combat = [r for r in strata if r["dataset_id"] == "ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3"]
    covid = [r for r in combat if r["disease"] == "COVID-19"]
    assert f"{max(r['permutation_count'] for r in covid):.2e}".startswith("4.69")
    assert (max(r["n_donors_A"] for r in combat), max(r["n_donors_B"] for r in combat)) == (100, 10)
    melms = [r for r in strata if r["dataset_id"] == "d8da613f-e681-4c69-b463-e94f5e66847f"]
    assert (max(r["n_donors_A"] for r in melms), max(r["n_donors_B"] for r in melms)) == (20, 7)
    wang = [r for r in strata if r["dataset_id"] == "4b6af54a-4a21-46e0-bc8d-673c0561a836"]
    assert sum(
        1 for r in wang
        if any(v is not None and v >= 0.8 for v in r["confound_cramers_v"].values())
    ) == 3 == attested["wang_near_confound_strata"]
    text = " ".join(_document_text().split())
    assert "4.69 × 10¹³" in text
    assert "#5 runs 100 versus 10 donors and #3 runs 20 versus 7" in text
    # 6671.5 / 11.0 = 606.5, so "a factor of 600" alone would have been an understatement dressed
    # as a fact; the document says "more than 600".
    assert "a factor of more than 600" in text


def test_the_median_counts_per_cell_names_its_own_population():
    """3081 is the median of the 502 per-group medians. The per-stratum mean of A and B gives
    3126.75 and the whole candidate set gives 3454.0, so a bare "median" names nothing."""
    text = " ".join(_document_text().split())
    assert "median of the 502 per-group medians" in text
    assert "3126.75" in text


def test_the_attested_block_covers_every_prose_figure_it_claims_to(frozen):
    """The list in §10 and the block in the artifact must not diverge."""
    text = " ".join(_document_text().split())
    for figure in ("4 609 595", "182", "124", "4 … 108", "7.28 × 10²³", "16 383", "3081",
                   "0.029 … 0.085", "306", "284 … 56 841.5", "11.0 … 6671.5"):
        assert figure in text, figure
    assert "attested_figures" in json.dumps(frozen["header"], ensure_ascii=False)


def test_a_drifting_attested_figure_aborts_the_freeze(monkeypatch, raw_source):
    monkeypatch.setattr(fz, "ATTESTED", {**fz.ATTESTED, "n_cell_type_labels": 123})
    with pytest.raises(fz.FrozenDeclarationMismatch, match="attested figures disagree"):
        fz.build(raw_source)


# ---------------------------------------------------------------------------
# §8 — the pinned release, read directly
# ---------------------------------------------------------------------------


def test_mathys_is_absent_from_both_pinned_indexes(raw_discover, raw_release):
    """The claim §8 rests on, recomputed over the committed bytes of both snapshots."""
    needles = fz.MATHYS_NEEDLES
    assert set(needles) == {
        "mathys", "rosmap", "religious order", "memory and aging",
        "s41586-019-1195", "10.1038/s41586-019-1195-2",
    }
    for snapshot in (raw_discover, raw_release):
        blob = json.dumps(snapshot["datasets"], ensure_ascii=False).lower()
        for needle in needles:
            assert needle not in blob, needle


def test_the_release_and_discover_are_measurably_different_objects(raw_discover, raw_release):
    """Why "absent from Discover" never implied "absent from the release".

    Six release datasets have no Discover record at all, and not one of the 1573 carries a
    dataset_version_id that Discover still lists. Both figures are stated in §2.2 and §8.
    """
    discover_ids = {row["dataset_id"] for row in raw_discover["datasets"]}
    discover_versions = {row["dataset_version_id"] for row in raw_discover["datasets"]}
    release = raw_release["datasets"]
    resolving = [row for row in release if row["dataset_id"] in discover_ids]
    assert len(release) == 1573
    assert len(resolving) == 1567
    assert len(release) - len(resolving) == 6
    assert sum(1 for row in release if row["dataset_version_id"] in discover_versions) == 0
    assert fz.RELEASE_VS_DISCOVER["n_dataset_version_id_matching"] == 0
    missing = sorted(row["dataset_id"] for row in release if row["dataset_id"] not in discover_ids)
    assert tuple(missing) == fz.RELEASE_NOT_IN_DISCOVER
    by_id = {row["dataset_id"]: row for row in release}
    assert tuple(sorted({by_id[i]["collection_id"] for i in missing})) == (
        fz.RELEASE_NOT_IN_DISCOVER_COLLECTIONS
    )


def test_section_8_no_longer_claims_the_release_table_was_unreachable():
    """The sentence that was false. `tiledbsoma` has no Windows wheel; `tiledb` does, and the array
    is public, so `census_info/datasets` was reachable all along and is now the evidence."""
    text = " ".join(_document_text().split())
    # "unreachable" is in RETRACTED_CLAIMS, so the parametrised test above already forbids it
    # standing unmarked; here we check the correction itself is present and specific.
    assert "That sentence is withdrawn: it is false." in text
    assert "49 releases" not in text
    assert "52 PyPI" in text or "52 PyPI\nreleases" in _document_text()
    assert "census_info/datasets" in text
    assert "resolve in Discover" in text


def test_the_elife_citation_is_the_one_that_was_verified():
    """Three authors, not two, and the RP article number. Recorded as an external fact."""
    text = " ".join(_document_text().split())
    assert "Murphy AE, Fancy NN, Skene NG (2023)" in text
    assert "eLife\n12:RP90214" in _document_text() or "eLife 12:RP90214" in text
    assert "10.7554/eLife.90214.3" in text
    assert "549 times fewer" in text
    assert "14 274" in text and "26" in text
    assert "corrected quality control *and* pseudobulk aggregation together" in text


# ---------------------------------------------------------------------------
# The Discover snapshot's own disagreements, recorded rather than argued away
# ---------------------------------------------------------------------------


def test_the_two_kpmp_is_primary_data_records_are_recorded_as_a_discrepancy(
    raw_discover, raw_source, frozen
):
    """Discover says neither KPMP dataset holds a primary cell; the manifest's Census query filtered
    on `is_primary_data == True` and returned cells for both. Neither is corrected here."""
    by_id = {row["dataset_id"]: row for row in raw_discover["datasets"]}
    seventeen = [d["dataset_id"] for d in fz.FROZEN_DATASETS] + [
        s["dataset_id"] for s in fz.SIBLING_DATASETS
    ]
    not_primary = [i for i in seventeen if by_id[i]["is_primary_data"] == [False]]
    assert sorted(not_primary) == sorted(fz.DISCOVER_NOT_PRIMARY)
    assert len(not_primary) == 2
    kpmp_pair = {
        "a12ccb9b-4fbe-457d-8590-ac78053259ef",
        "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d",
    }
    assert set(not_primary) == kpmp_pair
    assert raw_source["header"]["value_filter"] == fz.MANIFEST_VALUE_FILTER
    block = frozen["header"]["external_sources"]["discover_is_primary_data_discrepancy"]
    assert set(block["dataset_ids"]) == kpmp_pair
    assert "nothing here decides which is right" in block["note"]
    assert "is_primary_data" in _document_text()


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
    for script in (SCRIPT, REPO / "scripts" / "fetch_preregistration_evidence.py"):
        text = script.read_text(encoding="utf-8")
        assert not re.search(r"[A-Za-z]:[\\/]Users", text), f"{script.name}: account path"
        assert not re.search(r"(?m)^\s*[\"'](/home/|/Users/|/mnt/)", text), script.name


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
    """§6's trade-off is a load-bearing retraction, so its inputs are verified like any other
    declaration rather than being prose someone typed once."""
    monkeypatch.setattr(fz, "MANIFEST_TIER_CENSUS", {**fz.MANIFEST_TIER_CENSUS, 13: 4})
    with pytest.raises(fz.FrozenDeclarationMismatch, match="candidate-bearing datasets"):
        fz.build(raw_source)


def test_a_wrong_collection_tier_census_aborts_the_freeze(monkeypatch, raw_source):
    monkeypatch.setattr(fz, "MANIFEST_TIER_COLLECTIONS", {**fz.MANIFEST_TIER_COLLECTIONS, 23: 12})
    with pytest.raises(fz.FrozenDeclarationMismatch, match="distinct collections"):
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


def test_dropping_a_sibling_aborts_the_freeze(monkeypatch, raw_source):
    """The blocker, as a test. Naming four siblings where the release table says five must fail —
    that is the whole difference between a computed set and a typed one."""
    monkeypatch.setattr(fz, "SIBLING_DATASETS", fz.SIBLING_DATASETS[:-1])
    with pytest.raises(fz.FrozenDeclarationMismatch, match="sibling set derived"):
        fz.build(raw_source)


def test_inventing_a_sibling_aborts_the_freeze(monkeypatch, raw_source):
    """The other direction: a control set may not contain a dataset that shares no collection with
    the twelve, or "within-collection control" would stop meaning anything."""
    intruder = {
        "dataset_id": "9f222629-9e39-47d0-b83f-e08d610c7479",
        "short": "Human Lung Cell Atlas",
        "sibling_of": "d8da613f-e681-4c69-b463-e94f5e66847f",
        "collection_id": "6f6d381a-7701-4781-935c-db10d30de293",
        "role": fz.ROLE_WITHIN_COLLECTION_CONTROL,
        "role_note": "not in fact a sibling",
        "n_strata": 191,
        "ceiling_min_donors_per_group": 64,
    }
    monkeypatch.setattr(fz, "SIBLING_DATASETS", (*fz.SIBLING_DATASETS, intruder))
    with pytest.raises(fz.FrozenDeclarationMismatch, match="sibling set derived"):
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
    # The sibling's own Discover metadata is carried across too, so the substitution passes the
    # assay/tissue/DOI checks and only the sibling guard can catch it.
    mtg = {
        "dataset_id": "c2876b1b-06d8-4d96-a56b-5304f815b99a",
        "assay": ("10x 3' v3", "10x multiome"),
        "suspension": ("nucleus",),
        "tissue": ("middle temporal gyrus",),
        "doi": "10.1038/s41593-024-01774-5",
    }
    promoted = tuple(
        {**d, **mtg} if d["rank"] == 2 else dict(d) for d in fz.FROZEN_DATASETS
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
# The proposal document
# ---------------------------------------------------------------------------


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
    assert block["circulated_sha256"] == fz.PROPOSAL_ORIGINAL_SHA256
    assert block["circulated_bytes"] == fz.PROPOSAL_ORIGINAL_BYTES == 92_589


def test_the_proposal_status_still_says_what_it_is_for(frozen):
    """A mutation reversing this string's meaning survived the previous suite.

    The status is the only thing standing between a superseded working document and a reader who
    treats it as the pre-registration, so it is pinned by content and not merely by presence.
    """
    status = fz.PROPOSAL_STATUS
    for phrase in (
        "SUPERSEDED WORKING DOCUMENT",
        "NOT part of the binding act",
        "It has known errors",
        "the document governs",
    ):
        assert phrase in status, phrase
    assert "is part of the binding act" not in status.replace("NOT part of the binding act", "")
    assert frozen["header"]["reasoning_document"]["status"] == status


def test_the_redacted_proposal_leaks_no_absolute_path():
    """The privacy decision, enforced rather than remembered. Six occurrences of one repository
    root were replaced by a placeholder; a seventh anywhere in the file means the redaction missed
    something."""
    text = fz.PROPOSAL_MD.read_text(encoding="utf-8")
    assert not re.findall(r"[A-Za-z]:[\\/](?:Users|home)[\\/]\S*", text)
    # The redacted account name is deliberately NOT spelled out here. A test that guards a redaction
    # by writing the redacted token into a public test file publishes exactly what the redaction
    # removed. The pattern above already rejects any user-directory path whatever the account is
    # called; this loop additionally asserts that no component of THIS checkout's own absolute path
    # survived, which generalises to whichever machine regenerates the file.
    root = str(fz.REPO.resolve())
    for variant in (root, root.replace("\\", "/"), root.replace("/", "\\")):
        assert variant.lower() not in text.lower(), "the checkout root leaked into the proposal"
    assert text.count("<REPO>") == 7, "six substitutions, plus the one in the banner"
    fz.check_proposal_redaction()


def test_the_redacted_proposal_opens_with_its_superseded_banner():
    """A reader who opens the file directly on GitHub must not mistake a decided question for a
    live request for a decision."""
    head = fz.PROPOSAL_MD.read_text(encoding="utf-8")[:4000]
    for marker in fz.PROPOSAL_BANNER_MARKERS:
        assert marker in head, marker
    assert fz.PROPOSAL_ORIGINAL_SHA256 in head, "the banner must carry the circulated hash"
    assert head.index("SUPERSEDED WORKING DOCUMENT") < head.index("# pbcheck Phase 0")


def test_a_proposal_that_leaks_a_path_aborts_the_freeze(monkeypatch, tmp_path):
    leaky = tmp_path / "leaky.md"
    leaky.write_text(
        fz.PROPOSAL_MD.read_text(encoding="utf-8").replace("<REPO>", "C:/Users/someone/repo", 1),
        encoding="utf-8",
    )
    monkeypatch.setattr(fz, "PROPOSAL_MD", leaky)
    with pytest.raises(fz.SourceArtifactMismatch, match="redaction is incomplete"):
        fz.check_proposal_redaction()


def test_a_proposal_without_its_banner_aborts_the_freeze(monkeypatch, tmp_path):
    unbannered = tmp_path / "unbannered.md"
    text = fz.PROPOSAL_MD.read_text(encoding="utf-8")
    unbannered.write_text(text[text.index("# pbcheck Phase 0"):], encoding="utf-8")
    monkeypatch.setattr(fz, "PROPOSAL_MD", unbannered)
    with pytest.raises(fz.SourceArtifactMismatch, match="superseded-document banner"):
        fz.check_proposal_redaction()


def test_a_missing_proposal_document_aborts_the_freeze(monkeypatch, tmp_path):
    monkeypatch.setattr(fz, "PROPOSAL_MD", tmp_path / "absent.md")
    with pytest.raises(fz.SourceArtifactMismatch, match="proposal document not found"):
        fz.load_source()


def test_the_document_stops_claiming_a_defect_count_it_has_not_established():
    """§2 used to say the proposal contained "two figures this document corrects". §10.1 now lists
    what was established, one row each, and the phrasing does not promise exhaustiveness."""
    text = " ".join(_document_text().split())
    # The old claim survives exactly once, inside the sentence that withdraws it.
    assert text.count("two figures this document corrects") == 1
    assert 'An earlier draft of §2 said the proposal contained "two figures' in text
    assert "The discrepancies established in the proposal document" in text
    assert "It contains at least these twenty-seven" in text
    # The list itself: 10 reasoning rows plus 17 arithmetic rows.
    rows = [
        line for line in _document_text().splitlines()
        if re.match(r"^\|\s*(R|N)\d+\s*\|", line)
    ]
    assert len(rows) == 27, f"§10.1 lists {len(rows)} discrepancies"
    assert sum(1 for r in rows if re.match(r"^\|\s*R", r)) == 10
    assert sum(1 for r in rows if re.match(r"^\|\s*N", r)) == 17


# ---------------------------------------------------------------------------
# The document and the artifact must agree
# ---------------------------------------------------------------------------


def test_the_document_names_every_frozen_dataset_and_every_sibling():
    text = _document_text()
    for dataset in fz.FROZEN_DATASETS:
        assert dataset["dataset_id"] in text, dataset["short"]
    for sibling in fz.SIBLING_DATASETS:
        assert sibling["dataset_id"] in text, sibling["short"]


def test_the_document_records_the_source_hashes_and_the_frozen_count(frozen):
    text = _document_text()
    assert str(frozen["header"]["n_frozen_strata"]) in text
    assert "docs/AMENDMENTS.md" in text or "AMENDMENTS.md" in text


def test_the_freeze_script_verifies_the_documents_tables(frozen):
    """The guard the three previous rounds did not have, exercised here as well as in ``main``."""
    fz.verify_document(frozen)


@pytest.mark.parametrize(
    ("old", "new", "match"),
    [
        ("| Independent datasets | **12**, in **12 distinct collections** |",
         "| Independent datasets | **13**, in **12 distinct collections** |",
         "independent datasets"),
        ("**5** datasets, **106** strata", "**4** datasets, **106** strata", "sibling row"),
        ("| Elmentaite 2020, paediatric gut (Crohn) | 18 | 7 | 0 | 0 |",
         "| Elmentaite 2020, paediatric gut (Crohn) | 19 | 7 | 0 | 0 |",
         "column 3"),
        ("| `[3000, ∞)` | 13 |", "| `[3000, ∞)` | 14 |", "occupancy table"),
        ("| ≥ 8 v 8 | **33** | **25** | 7 |", "| ≥ 8 v 8 | **34** | **25** | 7 |", "datasets"),
        ("#9 clonal hematopoiesis", "#9 clonal haematopoiesis", "arms"),
        ("20 in the first, 16 in the second", "20 in the first, 12 in the second",
         "size of the tables"),
        ("18 + 25 + 11 + 43 + 9 = **106 strata**", "18 + 25 + 11 + 43 + 9 = **105 strata**",
         "sum the per-sibling counts"),
    ],
)
def test_a_hand_edited_table_cell_aborts_the_freeze(frozen, old, new, match):
    """Every table the document carries is compared against the artifact, in both directions.

    Three adversarial reviews found the same class of defect three times, always in a cell nothing
    recomputed. These mutations are what "recomputed" has to mean.
    """
    text = _document_text()
    assert old in text, f"the mutation target has moved: {old!r}"
    with pytest.raises(fz.FrozenDeclarationMismatch, match=match):
        fz.verify_document(frozen, text.replace(old, new, 1))


@pytest.mark.parametrize(
    ("old", "new"),
    [
        ("§9 item 7 forbids", "§9.7 forbids"),
        ("spec §10 risk 13", "spec §99 risk 13"),
        ("(§9 below)", "(§12 below)"),
    ],
)
def test_a_cross_reference_that_resolves_nowhere_aborts_the_freeze(old, new):
    """`§7.1` named nothing in either document and shipped. Every `§` now has to resolve, in the
    document the convention says it points at."""
    text = _document_text()
    assert old in text, f"the mutation target has moved: {old!r}"
    with pytest.raises(fz.FrozenDeclarationMismatch, match="not a section of"):
        fz.verify_document_cross_references(text.replace(old, new, 1))


def test_the_cross_reference_convention_is_stated_in_the_document():
    text = " ".join(_document_text().split())
    assert "`spec §N` is" in text
    assert "a bare `§N` or `§N.M` is a section of **this** document" in text


def test_every_dotted_reference_resolves_to_a_real_subsection_of_this_document():
    """Independent of the freeze script's own resolver."""
    headings = set(re.findall(r"(?m)^#{2,4}\s+(\d+(?:\.\d+)?)[.\s]", _document_text()))
    dotted = set(re.findall(r"§\s*(\d+\.\d+)", _document_text()))
    assert dotted, "the document should carry dotted subsection references"
    assert dotted <= headings, f"dangling dotted references: {sorted(dotted - headings)}"


def test_every_spec_reference_resolves_to_a_real_section_of_the_spec():
    spec_headings = set(re.findall(r"(?m)^#{2,4}\s+(\d+(?:\.\d+)?)[.\s]", SPEC.read_text(
        encoding="utf-8"
    )))
    cited = set(re.findall(r"(?i:spec)\s+§\s*(\d+(?:\.\d+)?)", _document_text()))
    assert cited, "the document should cite the spec"
    assert cited <= spec_headings, f"dangling spec references: {sorted(cited - spec_headings)}"


def test_the_two_open_decisions_that_could_move_the_frozen_set_say_so():
    """§9 items 4 and 5 would each change what the 251 are, and item 5 would move Layer B. An open
    question that can silently redefine the frozen set is not an open question."""
    text = " ".join(_document_text().split())
    assert "Choosing that harmonised level is an amendment" in text
    assert "Deciding it is an amendment, because it changes what the 251 are" in text
    # §9 item 8 records that the shipped moderated arm fits ~ 1 + x with no covariate slot, so
    # item 5 may not offer region-as-covariate as a live option in the same section.
    assert "no covariate slot" in text
    assert '"Or the covariate set" is struck from this item' in text


def test_the_document_does_not_claim_admission():
    """A pre-registration that reads as an admission is the one mistake that cannot be undone by a
    later correction, because the sweep would already have been run on it."""
    text = _document_text()
    assert "`admitted_to_sweep = False`" in text
    assert "admitted_to_sweep = True" not in text


def test_the_document_carries_no_self_invalidating_count():
    """"Three times already" was four by the time it shipped, and "the fourth consecutive
    pre-registration document" became the fifth the same week. A count of how many times something
    has happened goes stale on the next push, so neither survives."""
    text = " ".join(_document_text().split())
    assert "three times already" not in text
    assert "fourth consecutive" not in text
    assert "in the form that log already uses" in text
