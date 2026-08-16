"""Snapshot the two external indexes the §1 pre-registration reasons about, so they can be pinned.

``docs/PREREGISTRATION_STRATUM_LIST.md`` makes claims that the candidate manifest cannot settle —
which collection a dataset belongs to, what assay it was run on, whether Mathys 2019 is in the pinned
Census release. Until this script existed those claims rested on *live* reads of the CELLxGENE
Discover curation API, which is a moving target: a reader could not check them, and the reviews that
found defects in this document found them exactly there. Both indexes are therefore read once,
reduced to the fields the document actually uses, and committed:

* ``pilot/preregistration/discover_index_2026-08-16.json`` — the Discover curation index
  (``GET /curation/v1/datasets``), which is the only public source for a dataset's collection, assay,
  suspension type, tissue and publication DOI.
* ``pilot/preregistration/census_release_datasets_2025-01-30.json`` — the pinned Census release's own
  ``census_info/datasets`` dataframe, read straight out of the public TileDB array. This is the
  release enumerating *itself*, not an index of it.

**These snapshots are evidence, not a cache.** They are read once and committed; re-running this
script against today's upstream will produce different bytes, and that is not a regression but the
reason the bytes are pinned. ``scripts/freeze_stratum_list.py`` refuses to run unless both files hash
to the values recorded in the document, and every figure the document draws from them is recomputed
from the committed bytes on every freeze.

**What the Discover snapshot is and is not evidence of.** It was read on 2026-08-16 and describes
Discover *then*; the Census is pinned to ``2025-01-30``. The two are measurably not the same object —
see ``release_vs_discover`` in the release snapshot's header — so the Discover snapshot is evidence
about Discover, and the release snapshot is evidence about the release. The document says which of
the two each of its claims rests on.

**Why TileDB-Py and not cellxgene-census.** ``cellxgene-census`` requires ``tiledbsoma``, which
publishes no Windows wheel. ``tiledb`` (TileDB-Py) does, and ``census_info/datasets`` is an ordinary
public TileDB array under the release prefix, so the release's own dataset table is reachable here
without SOMA. TileDB-Py is a read-time dependency of this script alone: nothing in ``src/`` imports
it, no test needs it, and the freeze reads the committed JSON.

Usage (network, ~1 minute)::

    python scripts/fetch_preregistration_evidence.py --all

The reduced files are deterministic given the upstream payload — sorted by ``dataset_id``, LF line
endings, no generation timestamp of their own — so two people reading the same upstream state get the
same bytes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PREREG_DIR = REPO / "pilot" / "preregistration"

# ---------------------------------------------------------------------------
# Discover
# ---------------------------------------------------------------------------

DISCOVER_ENDPOINT = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
DISCOVER_READ_DATE = "2026-08-16"
DISCOVER_OUT = PREREG_DIR / f"discover_index_{DISCOVER_READ_DATE}.json"

#: The Discover fields the document uses, and no others. ``citation`` is retained even though the
#: document never prints it, because §8's Mathys needle search runs over it: a field a claim is
#: checked against has to be in the file the claim is checked from.
DISCOVER_FIELDS = (
    "dataset_id",
    # Retained although the document never prints it: it is the field that shows the Discover
    # index and the pinned release are not the same object (0 of 1573 release version ids appear
    # here), and a claim of that kind has to be recomputable from the committed bytes.
    "dataset_version_id",
    "collection_id",
    "collection_name",
    "collection_doi",
    "collection_doi_label",
    "citation",
    "title",
    "assay",
    "suspension_type",
    "tissue",
    "cell_count",
    "is_primary_data",
)

#: Fields Discover returns as ``[{"label": ..., "ontology_term_id": ...}, ...]``. Only the label is
#: kept: it is what the document prints and what a reader compares against the public record.
DISCOVER_ONTOLOGY_FIELDS = ("assay", "tissue")

#: §8(d)'s search terms for Mathys 2019 / ROSMAP, lower-cased. Recorded in the snapshot header with
#: the result of running them over the *whole* upstream payload, which the reduced file no longer
#: contains — see ``full_record_scan`` below.
MATHYS_NEEDLES = (
    "mathys",
    "rosmap",
    "religious order",
    "memory and aging",
    "s41586-019-1195",
    "10.1038/s41586-019-1195-2",
)

# ---------------------------------------------------------------------------
# The pinned Census release
# ---------------------------------------------------------------------------

CENSUS_VERSION = "2025-01-30"
CENSUS_S3_BUCKET = "cellxgene-census-public-us-west-2"
CENSUS_S3_REGION = "us-west-2"
CENSUS_DATASETS_URI = (
    f"s3://{CENSUS_S3_BUCKET}/cell-census/{CENSUS_VERSION}/soma/census_info/datasets"
)
RELEASE_READ_DATE = "2026-08-16"
RELEASE_OUT = PREREG_DIR / f"census_release_datasets_{CENSUS_VERSION}.json"

#: The release table's own columns, minus the two the document has no use for
#: (``soma_joinid`` is an internal row id and ``dataset_total_cell_count`` is not cited anywhere).
RELEASE_FIELDS = (
    "dataset_id",
    "dataset_version_id",
    "collection_id",
    "collection_name",
    "collection_doi",
    "collection_doi_label",
    "dataset_title",
    "dataset_h5ad_path",
    "citation",
)


def sha256_of_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _needle_hits(text: str) -> list[str]:
    lowered = text.lower()
    return [needle for needle in MATHYS_NEEDLES if needle in lowered]


def reduce_discover(records: list[dict]) -> list[dict]:
    reduced = []
    for record in records:
        row = {}
        for field in DISCOVER_FIELDS:
            value = record.get(field)
            if field in DISCOVER_ONTOLOGY_FIELDS and isinstance(value, list):
                value = sorted(
                    item["label"] if isinstance(item, dict) else item for item in value
                )
            elif isinstance(value, list) and all(isinstance(v, str) for v in value):
                value = sorted(value)
            row[field] = value
        reduced.append(row)
    reduced.sort(key=lambda r: r["dataset_id"])
    return reduced


def fetch_discover() -> dict:
    request = urllib.request.Request(DISCOVER_ENDPOINT, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=600) as response:  # noqa: S310 - pinned https URL
        raw = response.read()
    records = json.loads(raw)
    if not isinstance(records, list) or not records:
        raise RuntimeError(f"{DISCOVER_ENDPOINT} returned no records")

    # The whole-record scan, run here because the reduced file cannot support it afterwards: the
    # document must be able to say whether the needles appear in a field this snapshot dropped.
    whole_record_hits = {
        record["dataset_id"]: hits
        for record in records
        if (hits := _needle_hits(json.dumps(record, ensure_ascii=False)))
    }
    reduced = reduce_discover(records)
    retained_hits = {
        row["dataset_id"]: hits
        for row in reduced
        if (hits := _needle_hits(json.dumps(row, ensure_ascii=False)))
    }

    return {
        "header": {
            "what": "CELLxGENE Discover curation index, reduced to the fields "
                    "docs/PREREGISTRATION_STRATUM_LIST.md uses. Committed as evidence and pinned by "
                    "sha256 in §2 of that document.",
            "endpoint": DISCOVER_ENDPOINT,
            "read_date": DISCOVER_READ_DATE,
            "n_datasets_in_index": len(records),
            "n_datasets_in_snapshot": len(reduced),
            "fields_retained": list(DISCOVER_FIELDS),
            "fields_dropped_note":
                "Every other field of the upstream record is dropped. assay and tissue are reduced "
                "to their labels; the ontology term ids are not retained.",
            "upstream_payload_sha256": sha256_of_bytes(raw),
            "upstream_payload_bytes": len(raw),
            "self_sha256_note":
                "A file cannot carry its own sha256. The sha256 of THIS reduced file is recorded in "
                "§2 of docs/PREREGISTRATION_STRATUM_LIST.md and pinned in "
                "scripts/freeze_stratum_list.py, which refuses to run against different bytes. The "
                "hash above is of the raw upstream response this file was reduced from.",
            "describes":
                "Discover as of the read date. NOT the Census release: the release is pinned to "
                f"{CENSUS_VERSION} and the two indexes are measurably different objects — see "
                f"release_vs_discover in {RELEASE_OUT.name}.",
            "mathys_needles": list(MATHYS_NEEDLES),
            "full_record_scan": {
                "scope": "every field of every upstream record, before reduction",
                "n_records_scanned": len(records),
                "hits": whole_record_hits,
                "note":
                    "Run at read time against the payload whose sha256 is above. It is NOT "
                    "recomputable from this reduced file, and the document labels it as a "
                    "read-time observation rather than a machine-checked claim.",
            },
            "retained_field_scan": {
                "scope": "the retained fields of this file",
                "hits": retained_hits,
                "note": "Recomputable from these bytes; the freeze script recomputes it.",
            },
        },
        "datasets": reduced,
    }


def fetch_release(discover: dict | None = None) -> dict:
    try:
        import tiledb
    except ModuleNotFoundError as exc:  # pragma: no cover - environment-dependent
        raise RuntimeError(
            "TileDB-Py is needed to read the release's own dataset table "
            "(`pip install tiledb`). It is a read-time dependency of this script only: nothing in "
            "src/ imports it and the freeze reads the committed JSON."
        ) from exc

    context = tiledb.Ctx(tiledb.Config({
        "vfs.s3.region": CENSUS_S3_REGION,
        "vfs.s3.no_sign_request": "true",
    }))
    with tiledb.open(CENSUS_DATASETS_URI, ctx=context) as array:
        # ``.df[:]`` needs pandas < 3; the raw query path does not, and returns plain numpy arrays.
        result = array.query(use_arrow=False).multi_index[:]

    n_rows = len(result["dataset_id"])
    rows = []
    for index in range(n_rows):
        row = {}
        for field in RELEASE_FIELDS:
            value = result[field][index]
            row[field] = value.item() if hasattr(value, "item") else value
        rows.append(row)
    rows.sort(key=lambda r: r["dataset_id"])

    header = {
        "what": f"The census_info/datasets dataframe of Census release {CENSUS_VERSION}, read from "
                "the public TileDB array and reduced. This is the pinned release enumerating "
                "itself, not an index of it.",
        "uri": CENSUS_DATASETS_URI,
        "read_date": RELEASE_READ_DATE,
        "reader": "tiledb (TileDB-Py), anonymous S3 read; cellxgene-census/tiledbsoma publishes no "
                  "Windows wheel and was not used",
        "census_version": CENSUS_VERSION,
        "n_datasets": n_rows,
        "fields_retained": list(RELEASE_FIELDS),
        "fields_dropped_note":
            "soma_joinid (an internal row id) and dataset_total_cell_count (not cited) are dropped.",
        "mathys_needles": list(MATHYS_NEEDLES),
        "mathys_hits": {
            row["dataset_id"]: hits
            for row in rows
            if (hits := _needle_hits(json.dumps(row, ensure_ascii=False)))
        },
    }

    if discover is not None:
        by_discover_id = {d["dataset_id"]: d for d in discover["datasets"]}
        version_ids = {d.get("dataset_version_id") for d in discover["datasets"]}
        resolving = [r for r in rows if r["dataset_id"] in by_discover_id]
        header["release_vs_discover"] = {
            "comparand": DISCOVER_OUT.name,
            "discover_read_date": discover["header"]["read_date"],
            "n_release_datasets": n_rows,
            "n_resolving_in_discover_by_dataset_id": len(resolving),
            "n_not_resolving_in_discover": n_rows - len(resolving),
            "not_resolving_dataset_ids": sorted(
                r["dataset_id"] for r in rows if r["dataset_id"] not in by_discover_id
            ),
            "n_dataset_version_id_present_in_discover": sum(
                1 for r in rows if r["dataset_version_id"] in version_ids
            ),
            "n_collection_doi_differing": sum(
                1 for r in resolving
                if (r["collection_doi"] or None)
                != (by_discover_id[r["dataset_id"]].get("collection_doi") or None)
            ),
            "n_collection_id_differing": sum(
                1 for r in resolving
                if r["collection_id"] != by_discover_id[r["dataset_id"]]["collection_id"]
            ),
            "note":
                "'Resolves in Discover' means: the release row's dataset_id is a dataset_id in the "
                "Discover snapshot. dataset_version_id is compared too, because it is the field "
                "that shows the two indexes are not the same object.",
        }
    return {"header": header, "datasets": rows}


def render(payload: dict) -> str:
    return json.dumps(payload, indent=1, ensure_ascii=False, sort_keys=False) + "\n"


def write(payload: dict, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(render(payload))
    return path


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--discover", action="store_true", help="snapshot the Discover index")
    parser.add_argument("--release", action="store_true", help="snapshot the pinned release's table")
    parser.add_argument("--all", action="store_true", help="both")
    return parser


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)
    want_discover = args.discover or args.all
    want_release = args.release or args.all
    if not (want_discover or want_release):
        build_arg_parser().print_help()
        return 2

    discover = None
    if want_discover:
        discover = fetch_discover()
        path = write(discover, DISCOVER_OUT)
        print(f"wrote {path} ({path.stat().st_size} bytes, "
              f"{discover['header']['n_datasets_in_snapshot']} datasets)")
        print(f"  sha256 {sha256_of_bytes(path.read_bytes())}")
    if want_release:
        if discover is None and DISCOVER_OUT.exists():
            discover = json.loads(DISCOVER_OUT.read_text(encoding="utf-8"))
        release = fetch_release(discover)
        path = write(release, RELEASE_OUT)
        print(f"wrote {path} ({path.stat().st_size} bytes, "
              f"{release['header']['n_datasets']} datasets)")
        print(f"  sha256 {sha256_of_bytes(path.read_bytes())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
