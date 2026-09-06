"""Tests for the WP6 demonstration scripts and their shared helpers.

All offline: no network access, no read of a real demo dataset. ``load_and_check`` and
``pick_celltype`` are exercised on small synthetic ``AnnData`` frames; ``download`` is exercised
against local ``file://`` URLs and pre-existing files only. The frozen-list membership test reads
the committed pre-registration CSV and scans ``docs/`` and ``pilot/`` text files, both read-only.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from demo_common import DEMO_DATASETS, download, load_and_check, pick_celltype  # noqa: E402

from pbcheck import audit_schema  # noqa: E402

STRATUM_LIST_CSV = REPO_ROOT / "pilot" / "preregistration" / "stratum_list_2026-08-16.csv"

#: Text file extensions the frozen-list membership scan reads. Every file under docs/ and pilot/
#: is one of these (verified once when this test was written: csv, json, md, gitkeep).
_SCANNED_EXTENSIONS = {".md", ".csv", ".json", ".txt"}


def _synthetic_adata(n_donors: int = 4, n_cells_per_donor: int = 5) -> ad.AnnData:
    donors = []
    conditions = []
    celltypes = []
    for i in range(n_donors):
        donor = f"donor_{i}"
        condition = "ctrl" if i % 2 == 0 else "stim"
        for j in range(n_cells_per_donor):
            donors.append(donor)
            conditions.append(condition)
            # An imbalanced cell-type split so "largest level" has one clear answer.
            celltypes.append("big_type" if j < n_cells_per_donor - 1 else "small_type")
    n_obs = len(donors)
    n_var = 6
    rng = np.random.default_rng(0)
    X = rng.integers(0, 5, size=(n_obs, n_var)).astype(np.float32)
    obs = pd.DataFrame(
        {"donor": donors, "condition": conditions, "cell_type": celltypes},
        index=[f"cell_{k}" for k in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{k}" for k in range(n_var)])
    return ad.AnnData(X=X, obs=obs, var=var)


# ---------------------------------------------------------------------------
# load_and_check
# ---------------------------------------------------------------------------


def test_load_and_check_returns_adata_when_columns_present(tmp_path):
    adata = _synthetic_adata()
    path = tmp_path / "synthetic.h5ad"
    adata.write_h5ad(path)

    loaded = load_and_check(path, ("donor", "condition", "cell_type"))
    assert list(loaded.obs.columns) == list(adata.obs.columns)
    assert loaded.n_obs == adata.n_obs


def test_load_and_check_raises_with_available_columns_and_levels(tmp_path):
    adata = _synthetic_adata()
    path = tmp_path / "synthetic.h5ad"
    adata.write_h5ad(path)

    with pytest.raises(AssertionError) as excinfo:
        load_and_check(path, ("donor", "condition", "replicate"))

    message = str(excinfo.value)
    assert "replicate" in message
    assert "donor" in message and "condition" in message and "cell_type" in message
    assert "ctrl" in message and "stim" in message


# ---------------------------------------------------------------------------
# pick_celltype
# ---------------------------------------------------------------------------


def test_pick_celltype_returns_largest_level():
    adata = _synthetic_adata(n_donors=4, n_cells_per_donor=5)
    assert pick_celltype(adata, "cell_type") == "big_type"


def test_pick_celltype_raises_on_missing_column():
    adata = _synthetic_adata()
    with pytest.raises(AssertionError) as excinfo:
        pick_celltype(adata, "no_such_column")
    assert "no_such_column" in str(excinfo.value)


# ---------------------------------------------------------------------------
# download
# ---------------------------------------------------------------------------


def test_download_skips_when_dest_already_matches(tmp_path):
    dest = tmp_path / "already_here.bin"
    content = b"pbcheck demo fixture bytes"
    dest.write_bytes(content)
    sha = hashlib.sha256(content).hexdigest()

    # A URL that would raise if it were ever opened: the skip-when-matches path must never touch it.
    result = download("http://demo-common-test.invalid/should-not-be-fetched", dest, sha)

    assert result == dest
    assert dest.read_bytes() == content


def test_download_fetches_and_verifies_via_file_url(tmp_path):
    source = tmp_path / "source.bin"
    content = b"kang 2018 demo fixture, not the real file"
    source.write_bytes(content)
    sha = hashlib.sha256(content).hexdigest()

    dest = tmp_path / "downloaded.bin"
    url = source.resolve().as_uri()

    result = download(url, dest, sha)

    assert result == dest
    assert dest.read_bytes() == content
    assert not dest.with_name(dest.name + ".part").exists()


def test_download_raises_and_cleans_up_on_sha_mismatch(tmp_path):
    source = tmp_path / "source.bin"
    source.write_bytes(b"actual content")
    wrong_sha = hashlib.sha256(b"different content").hexdigest()

    dest = tmp_path / "downloaded.bin"
    url = source.resolve().as_uri()

    with pytest.raises(ValueError, match="sha256"):
        download(url, dest, wrong_sha)

    assert not dest.exists()
    assert not dest.with_name(dest.name + ".part").exists()


# ---------------------------------------------------------------------------
# scripts/demo_kang2018.py's pinned constant
# ---------------------------------------------------------------------------


def test_kang_2018_sha256_is_pinned_hex_digest():
    import demo_kang2018

    sha = demo_kang2018.KANG_2018_SHA256
    assert isinstance(sha, str)
    assert re.fullmatch(r"[0-9a-f]{64}", sha), f"KANG_2018_SHA256 is not a 64-hex sha256: {sha!r}"


# ---------------------------------------------------------------------------
# Frozen-list membership: neither demo dataset is among the 17 frozen datasets / 357 strata.
# ---------------------------------------------------------------------------


def _all_demo_tokens() -> list[str]:
    tokens: list[str] = []
    for dataset in DEMO_DATASETS.values():
        tokens.extend(dataset.search_tokens)
        tokens.extend(dataset.known_identifiers)
    return tokens


def test_membership_no_stratum_list_row_matches_a_demo_dataset():
    tokens = [t.lower() for t in _all_demo_tokens()]
    assert tokens, "DEMO_DATASETS carries no search tokens to check"

    with open(STRATUM_LIST_CSV, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    assert len(rows) == 357, f"expected 357 rows in the frozen stratum list, found {len(rows)}"

    offending = []
    for row in rows:
        for column in ("dataset_short", "dataset_id"):
            value = (row.get(column) or "").lower()
            for token in tokens:
                if token in value:
                    offending.append((column, row.get(column), token))
    assert not offending, f"frozen stratum list row(s) match a demo dataset token: {offending}"


def test_membership_docs_and_pilot_text_contain_no_demo_tokens():
    tokens = [t.lower() for t in _all_demo_tokens()]
    assert tokens, "DEMO_DATASETS carries no search tokens to check"

    offending = []
    for base in (REPO_ROOT / "docs", REPO_ROOT / "pilot"):
        for path in base.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in _SCANNED_EXTENSIONS:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore").lower()
            for token in tokens:
                if token in text:
                    offending.append((str(path.relative_to(REPO_ROOT)), token))
    assert not offending, f"docs/ or pilot/ text contains a demo dataset token: {offending}"


# ---------------------------------------------------------------------------
# Committed demo artifacts: schema-valid and under the 500 KB pre-commit cap.
# ---------------------------------------------------------------------------

_MAX_COMMITTED_FILE_BYTES = 500 * 1024


def _demo_dirs() -> list[Path]:
    demo_root = REPO_ROOT / "demo"
    if not demo_root.is_dir():
        return []
    return [p for p in demo_root.iterdir() if p.is_dir()]


@pytest.mark.parametrize("demo_dir", _demo_dirs(), ids=lambda p: p.name)
def test_committed_demo_audit_json_validates_and_is_small(demo_dir):
    audit_path = demo_dir / "pbcheck_audit.json"
    assert audit_path.exists(), f"{demo_dir} has no pbcheck_audit.json"
    payload = json.loads(audit_path.read_text(encoding="utf-8"))
    audit_schema.validate(payload)

    for path in demo_dir.iterdir():
        if path.is_file():
            size = path.stat().st_size
            assert size < _MAX_COMMITTED_FILE_BYTES, f"{path} is {size} bytes, over the 500 KB cap"


def test_demo_directory_exists_and_is_not_empty():
    demo_dirs = _demo_dirs()
    if not demo_dirs:
        pytest.skip("no demo/ directory yet")
    assert demo_dirs


# ---------------------------------------------------------------------------
# demo/README.md: carries the disclaimer required by the release plan (section 6).
# ---------------------------------------------------------------------------

_DISCLAIMER_FIRST_SENTENCE = (
    "These are demonstrations outside the pre-registered Phase 0 protocol."
)


def test_demo_readme_contains_disclaimer_first_sentence():
    demo_readme = REPO_ROOT / "demo" / "README.md"
    if not demo_readme.exists():
        pytest.skip("no demo/README.md yet")
    text = demo_readme.read_text(encoding="utf-8")
    assert _DISCLAIMER_FIRST_SENTENCE in text
