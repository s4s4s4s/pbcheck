"""Tests for the report section model, the Markdown/HTML renderers and ``write_outputs`` (WP2).

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

from pbcheck import audit_schema
from pbcheck.render import render_html, render_markdown, write_outputs
from pbcheck.render.sections import (
    Callout,
    KeyValues,
    Paragraph,
    Table,
    build_sections,
    prose_blocks,
    summary_lines,
)
from pbcheck.render.text import FORBIDDEN_PATTERNS

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


@pytest.fixture(params=FIXTURE_NAMES)
def fixture_payload(request: pytest.FixtureRequest) -> dict:
    return _load_fixture(request.param)


def test_fixtures_validate(fixture_payload: dict) -> None:
    audit_schema.validate(fixture_payload)


def test_render_prose_has_no_forbidden_patterns(fixture_payload: dict) -> None:
    sections = build_sections(fixture_payload)
    for text in prose_blocks(sections):
        for pattern in FORBIDDEN_PATTERNS:
            assert re.search(pattern, text) is None, (
                f"prose block matches forbidden pattern {pattern!r}: {text!r}"
            )


def test_gene_symbol_in_table_not_flagged() -> None:
    payload = _load_fixture("complete")
    sections = build_sections(payload)

    prose_texts = prose_blocks(sections)
    assert not any("GOLGA8A" in text for text in prose_texts)

    table_texts: list[str] = []
    for section in sections:
        for block in section.blocks:
            if isinstance(block, Table):
                table_texts.extend(str(cell) for row in block.rows for cell in row)
    assert any("GOLGA8A" in text for text in table_texts)

    rendered = render_markdown(payload)
    assert "GOLGA8A" in rendered


def test_all_caveats_present_by_status(fixture_payload: dict) -> None:
    rendered = render_markdown(fixture_payload)
    for caveat in fixture_payload["caveats"]:
        assert caveat["text"] in rendered, f"caveat {caveat['id']} missing from rendered report"


def test_markdown_tables_are_pipe_tables() -> None:
    payload = _load_fixture("complete")
    rendered = render_markdown(payload)
    header_lines = [line for line in rendered.splitlines() if line.startswith("| ") and line.endswith(" |")]
    assert header_lines, "no pipe-table rows found"

    separator_pattern = re.compile(r"^\|( ?--- ?\|)+$")
    separator_lines = [line for line in rendered.splitlines() if separator_pattern.match(line)]
    assert separator_lines, "no pipe-table separator rows found"

    for section in build_sections(payload):
        for block in section.blocks:
            if isinstance(block, (Table, KeyValues)):
                headers = block.headers if isinstance(block, Table) else ("Field", "Value")
                header_line = "| " + " | ".join(headers) + " |"
                assert header_line in rendered


def test_paired_floor_omitted_when_not_shown() -> None:
    payload = _load_fixture("naive_only")
    assert payload["readout"]["paired_floor_shown"] is False
    sections = build_sections(payload)

    naive_section = next(section for section in sections if section.id == "6")
    captions = [
        block.caption
        for block in naive_section.blocks
        if isinstance(block, (Table, KeyValues)) and block.caption
    ]
    assert "Paired permutation floor" not in captions
    assert "Solo permutation floor" in captions

    rendered_text = " ".join(
        block.text for block in naive_section.blocks if isinstance(block, (Paragraph, Callout))
    )
    assert "paired floor is not shown" in rendered_text


def test_summary_lines_reuse_sections(fixture_payload: dict) -> None:
    sections = build_sections(fixture_payload)
    expected: list[str] = []
    for section in sections:
        if section.id in ("1", "2"):
            expected.extend(block.text for block in section.blocks if isinstance(block, Paragraph))
    assert summary_lines(fixture_payload) == expected
    assert len(expected) > 0


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_golden_markdown(name: str) -> None:
    payload = _load_fixture(name)
    rendered = render_markdown(payload)
    _check_golden(FIXTURES_DIR / f"audit_payload_{name}.md", rendered)


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_golden_html(name: str) -> None:
    payload = _load_fixture(name)
    rendered = render_html(payload)
    _check_golden(FIXTURES_DIR / f"audit_payload_{name}.html", rendered)


def test_html_escapes_values() -> None:
    payload = _load_fixture("complete")
    rendered = render_html(payload)
    assert "<b>" not in rendered
    assert "&lt;b&gt;" in rendered


def test_html_has_no_script_or_external_asset() -> None:
    payload = _load_fixture("complete")
    rendered = render_html(payload)
    assert "<script" not in rendered.lower()
    assert "http://" not in rendered
    assert "https://" not in rendered
    assert "<link" not in rendered.lower()


def test_write_outputs_refuses_overwrite(tmp_path: Path) -> None:
    payload = _load_fixture("naive_only")
    write_outputs(payload, tmp_path)
    with pytest.raises(FileExistsError, match="pbcheck_audit.json"):
        write_outputs(payload, tmp_path)


def test_write_outputs_utf8(tmp_path: Path) -> None:
    payload = _load_fixture("design_only")
    out_dir = tmp_path / "étude-café"
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
