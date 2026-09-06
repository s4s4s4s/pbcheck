"""Tests for the two output documents and the writer: Markdown, HTML and ``write_outputs``.

Golden tests: for each of the four fixture payloads, the rendered Markdown and HTML are compared
byte-for-byte against a committed file under ``tests/fixtures``. Regenerate the committed files by
running this module with the environment variable ``PBCHECK_UPDATE_GOLDEN=1`` set; a pytest CLI
option cannot be registered for this instead, because ``tests/conftest.py`` belongs to another work
package during this release and this module must not add one.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

import pytest

from pbcheck.render import render_html, render_markdown, write_outputs
from pbcheck.render.sections import KeyValues, Table, build_sections

FIXTURES_DIR = Path(__file__).parent / "fixtures"
FIXTURE_NAMES = ("complete", "naive_only", "design_only", "coarse_null")

_UPDATE_GOLDEN = os.environ.get("PBCHECK_UPDATE_GOLDEN") == "1"


def _load_fixture(name: str) -> dict:
    path = FIXTURES_DIR / f"audit_payload_{name}.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _check_golden(golden_path: Path, rendered: str) -> None:
    if _UPDATE_GOLDEN:
        with open(golden_path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(rendered)
        return
    expected = golden_path.read_text(encoding="utf-8")
    assert rendered == expected, f"rendered output no longer matches {golden_path}"


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_golden_markdown(name: str) -> None:
    payload = _load_fixture(name)
    _check_golden(FIXTURES_DIR / f"audit_payload_{name}.md", render_markdown(payload))


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_golden_html(name: str) -> None:
    payload = _load_fixture(name)
    _check_golden(FIXTURES_DIR / f"audit_payload_{name}.html", render_html(payload))


def test_markdown_tables_are_pipe_tables() -> None:
    payload = _load_fixture("complete")
    rendered = render_markdown(payload)
    header_lines = [
        line for line in rendered.splitlines() if line.startswith("| ") and line.endswith(" |")
    ]
    assert header_lines, "no pipe-table rows found"

    separator_pattern = re.compile(r"^\|( ?--- ?\|)+$")
    assert [line for line in rendered.splitlines() if separator_pattern.match(line)]

    for section in build_sections(payload):
        for block in section.blocks:
            if isinstance(block, (Table, KeyValues)):
                headers = block.headers if isinstance(block, Table) else ("Field", "Value")
                assert "| " + " | ".join(headers) + " |" in rendered


def test_markdown_prints_every_read_out_sentence_in_order() -> None:
    payload = _load_fixture("complete")
    rendered = render_markdown(payload)
    positions = [rendered.index(sentence) for sentence in payload["readout"]["sentences"]]
    assert positions == sorted(positions)


def test_html_escapes_values() -> None:
    rendered = render_html(_load_fixture("complete"))
    assert "<b>" not in rendered
    assert "&lt;b&gt;" in rendered


def test_html_has_no_script_or_external_asset() -> None:
    rendered = render_html(_load_fixture("complete"))
    assert "<script" not in rendered.lower()
    assert "http://" not in rendered
    assert "https://" not in rendered
    assert "<link" not in rendered.lower()


def test_write_outputs_refuses_overwrite(tmp_path: Path) -> None:
    payload = _load_fixture("naive_only")
    write_outputs(payload, tmp_path)
    with pytest.raises(FileExistsError, match="pbcheck_audit.json"):
        write_outputs(payload, tmp_path)


def test_a_refused_write_creates_neither_file_nor_directory(tmp_path: Path) -> None:
    """Every target is checked before anything is created: a refused call leaves the filesystem
    as it found it, rather than an empty output directory the user then has to delete."""
    payload = _load_fixture("naive_only")

    out_dir = tmp_path / "out"
    write_outputs(payload, out_dir, formats=("json",))
    with pytest.raises(FileExistsError, match="pbcheck_audit.json"):
        write_outputs(payload, out_dir)
    assert not (out_dir / "pbcheck_report.md").exists()
    assert not (out_dir / "pbcheck_report.html").exists()

    never_made = tmp_path / "never" / "made"
    with pytest.raises(ValueError, match="unknown output format"):
        write_outputs(payload, never_made, formats=("json", "pdf"))
    assert not never_made.exists()
    assert not never_made.parent.exists()


def test_write_outputs_utf8(tmp_path: Path) -> None:
    payload = _load_fixture("design_only")
    out_dir = tmp_path / "etude-cafe-éà"
    written = write_outputs(payload, out_dir)

    assert set(written) == {"json", "md", "html"}
    for fmt, path in written.items():
        assert path.parent == out_dir
        content = path.read_bytes()
        assert b"\r\n" not in content
        text = content.decode("utf-8")
        if fmt == "json":
            assert json.loads(text) == payload
        elif fmt == "md":
            assert text == render_markdown(payload)
        else:
            assert text == render_html(payload)
