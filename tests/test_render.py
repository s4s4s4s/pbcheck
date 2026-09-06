"""Tests for the report section model (:mod:`pbcheck.render.sections`).

What this module covers: which blocks a section carries at each status, that the read-out lines
are the payload's own sentences rather than a second derivation, the conditions the paired floor
and the ratio line hang on, and the wording gate over everything the renderer calls prose. The
Markdown and HTML documents built from these sections, and ``write_outputs``, are covered by
``tests/test_render_output.py``.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from pbcheck import audit_schema
from pbcheck.product_constants import FEW_DONORS_THRESHOLD
from pbcheck.render import render_markdown
from pbcheck.render.sections import (
    Callout,
    KeyValues,
    Paragraph,
    Table,
    build_sections,
    prose_blocks,
    summary_lines,
)
from pbcheck.render.text import STATUS_REASON_WORDS, forbidden_pattern_hits

FIXTURES_DIR = Path(__file__).parent / "fixtures"
FIXTURE_NAMES = ("complete", "naive_only", "design_only", "coarse_null")


def _load_fixture(name: str) -> dict:
    """One committed golden payload. Loaded here and in ``test_render_output.py`` alike; a shared
    helper would have to live in ``tests/conftest.py``, which belongs to another work package."""
    path = FIXTURES_DIR / f"audit_payload_{name}.json"
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(params=FIXTURE_NAMES)
def fixture_payload(request: pytest.FixtureRequest) -> dict:
    return _load_fixture(request.param)


def test_fixtures_validate(fixture_payload: dict) -> None:
    audit_schema.validate(fixture_payload)


def test_render_prose_has_no_forbidden_patterns(fixture_payload: dict) -> None:
    for prose in prose_blocks(build_sections(fixture_payload)):
        assert forbidden_pattern_hits(prose) == (), (
            f"prose block matches {forbidden_pattern_hits(prose)}: {prose!r}"
        )


def test_prose_blocks_cover_table_captions() -> None:
    """A caption is written by the renderer and reaches the reader, so the gate must scan it."""
    sections = build_sections(_load_fixture("complete"))
    captions = [
        block.caption
        for section in sections
        for block in section.blocks
        if isinstance(block, (Table, KeyValues)) and block.caption
    ]
    assert captions
    assert set(captions) <= set(prose_blocks(sections))


def test_gene_symbol_in_table_not_flagged() -> None:
    payload = _load_fixture("complete")
    sections = build_sections(payload)

    assert not any("GOLGA8A" in prose for prose in prose_blocks(sections))

    table_texts = [
        str(cell)
        for section in sections
        for block in section.blocks
        if isinstance(block, Table)
        for row in block.rows
        for cell in row
    ]
    assert any("GOLGA8A" in cell for cell in table_texts)
    assert "GOLGA8A" in render_markdown(payload)


def _markdown_prose(rendered: str) -> str:
    """The rendered Markdown with its block decoration removed.

    The renderer sets a caveat off as a quote block (a ``> `` prefix on every line) and keeps the
    line breaks its author put in (a two-space hard break at the end of each line, so the
    operating-envelope rows of note N2 stay one row per line). Both are presentation: the words
    are not touched, so a payload sentence must still appear verbatim once they are stripped.
    """
    return "\n".join(line.removeprefix("> ").removesuffix("  ") for line in rendered.split("\n"))


def test_all_caveats_present_by_status(fixture_payload: dict) -> None:
    prose = _markdown_prose(render_markdown(fixture_payload))
    for caveat in fixture_payload["caveats"]:
        assert caveat["text"] in prose, f"caveat {caveat['id']} missing from rendered report"


def test_read_out_lines_are_the_payloads_own_sentences(fixture_payload: dict) -> None:
    """The read-out lines are ``readout.sentences``, printed in order and unchanged.

    The renderer must not re-derive them: the JSON a user keeps and the report he reads would then
    be two independent renderings of the same numbers.
    """
    section = next(s for s in build_sections(fixture_payload) if s.id == "2")
    paragraphs = [block.text for block in section.blocks if isinstance(block, Paragraph)]
    sentences = fixture_payload["readout"]["sentences"]
    assert paragraphs[: len(sentences)] == sentences
    if fixture_payload["permutation_null"] is not None:
        assert sentences, "a run with a permutation null carries its read-out sentences"


def test_read_out_falls_back_to_one_note_when_no_arm_ran() -> None:
    payload = _load_fixture("design_only")
    assert payload["readout"]["sentences"] == []
    section = next(s for s in build_sections(payload) if s.id == "2")
    paragraphs = [block.text for block in section.blocks if isinstance(block, Paragraph)]
    assert len(paragraphs) == 1
    assert STATUS_REASON_WORDS[payload["status_reason"]] in paragraphs[0]


def test_ratio_line_names_the_bh_convention_of_each_arm() -> None:
    payload = _load_fixture("complete")
    assert payload["readout"]["paired_floor_shown"] is True
    assert payload["readout"]["few_donors"] is False
    ratio = _ratio_paragraph(payload)
    assert "solo BH" in ratio and "paired BH" in ratio


def test_ratio_line_omits_the_pseudobulk_arm_when_the_paired_floor_is_not_shown() -> None:
    payload = _load_fixture("coarse_null")
    assert payload["readout"]["paired_floor_shown"] is False
    ratio = _ratio_paragraph(payload)
    assert "solo BH" in ratio
    assert "paired BH" not in ratio


def test_ratio_line_says_the_ratios_are_leak_contaminated_below_the_donor_threshold() -> None:
    payload = _load_fixture("coarse_null")
    assert payload["readout"]["min_donors_per_group"] < FEW_DONORS_THRESHOLD
    assert payload["readout"]["few_donors"] is True
    assert "leak" in _ratio_paragraph(payload)


def _ratio_paragraph(payload: dict) -> str:
    """The read-out section's last paragraph, which is the ratio line."""
    section = next(s for s in build_sections(payload) if s.id == "2")
    return [block.text for block in section.blocks if isinstance(block, Paragraph)][-1]


def _naive_section_blocks(payload: dict) -> tuple[list[str], list[str]]:
    section = next(s for s in build_sections(payload) if s.id == "6")
    captions = [
        block.caption
        for block in section.blocks
        if isinstance(block, (Table, KeyValues)) and block.caption
    ]
    prose = [block.text for block in section.blocks if isinstance(block, (Paragraph, Callout))]
    return captions, prose


def test_paired_floor_omitted_when_the_arms_do_not_cover_the_same_genes() -> None:
    payload = _load_fixture("naive_only")
    assert payload["readout"]["paired_floor_shown"] is False
    captions, prose = _naive_section_blocks(payload)
    assert "Paired permutation floor" not in captions
    assert "Solo permutation floor" in captions
    n_na = payload["real_label"]["paired_bh"]["n_na_pseudobulk"]
    assert any(f"left {n_na} genes without a value" in text for text in prose)


def test_paired_floor_shown_but_never_measured_says_so_without_contradicting_itself() -> None:
    """``paired_floor_shown`` and a null ``floor_paired`` are separate conditions: the note that
    covers the first must not claim the pseudobulk arm left genes without a value."""
    payload = copy.deepcopy(_load_fixture("complete"))
    assert payload["readout"]["paired_floor_shown"] is True
    payload["permutation_null"]["naive"]["floor_paired"] = None

    captions, prose = _naive_section_blocks(payload)
    assert "Paired permutation floor" not in captions
    note = next(text for text in prose if "paired floor is not shown" in text)
    assert "recorded no paired floor" in note
    assert "without a value" not in note


def test_an_unreached_permutation_count_is_said_in_words() -> None:
    payload = copy.deepcopy(_load_fixture("complete"))
    payload["readout"]["n_perm_naive_achieved"] = None
    captions_to_items = {
        block.caption: dict(block.items)
        for section in build_sections(payload)
        for block in section.blocks
        if isinstance(block, KeyValues) and block.caption
    }
    assert captions_to_items["Permutations, naive arm"]["achieved"] == "not reached"


def test_a_paired_design_is_explained_in_full() -> None:
    """``donor_spans_conditions`` names all three parts of the reason, not only the refusal."""
    words = STATUS_REASON_WORDS["donor_spans_conditions"]
    assert "paired design" in words
    assert "paired or mixed model" in words
    assert "independent" in words

    payload = copy.deepcopy(_load_fixture("design_only"))
    payload["status_reason"] = "donor_spans_conditions"
    payload["caveats"] = [c for c in payload["caveats"] if c["id"] != "N4"]
    payload["readout"]["few_donors"] = False
    payload["readout"]["min_donors_per_group"] = payload["readout"]["few_donors_threshold"]
    audit_schema.validate(payload)
    assert words in render_markdown(payload)


def test_summary_lines_carry_the_scope_note_and_the_read_out(fixture_payload: dict) -> None:
    """The CLI summary is sections 1 and 2 verbatim, callouts included: a terminal reader gets the
    always-on scope note (N1), not the numbers alone."""
    expected = [
        block.text
        for section in build_sections(fixture_payload)
        if section.id in ("1", "2")
        for block in section.blocks
        if isinstance(block, (Paragraph, Callout))
    ]
    lines = summary_lines(fixture_payload)
    assert lines == expected
    n1 = next(c["text"] for c in fixture_payload["caveats"] if c["id"] == "N1")
    assert n1 in lines


def test_every_section_is_present_at_every_status(fixture_payload: dict) -> None:
    sections = build_sections(fixture_payload)
    assert [section.id for section in sections] == [str(i) for i in range(1, 11)]
    assert all(section.blocks for section in sections)
