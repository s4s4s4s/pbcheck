"""Tests for the report section model and the Markdown renderer (WP2, part 1)."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from pbcheck import audit_schema
from pbcheck.render.markdown import render_markdown
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


def _load_fixture(name: str) -> dict:
    path = FIXTURES_DIR / f"audit_payload_{name}.json"
    return json.loads(path.read_text(encoding="utf-8"))


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
