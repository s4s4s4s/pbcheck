"""Tests for the report's fixed prose: the notes N1..N9, the read-out paragraph R1, the read-out
sentence templates, the glossary, and the forbidden-pattern gate they are all checked against."""

from __future__ import annotations

import collections
import inspect
import itertools
import re

import pytest

from pbcheck import gate_config, product_constants
from pbcheck.audit_schema import CAVEAT_IDS, LAMBDA_CLASS_VALUES
from pbcheck.render import text

#: The forbidden-pattern list as the wording rules state it. Pinned literally here so deleting a
#: pattern from the module breaks a test instead of silently widening what the report may say.
EXPECTED_FORBIDDEN_PATTERNS = (
    r"validit|validat",
    r"\bvalid\b",
    r"\bfindings?\b",
    r"false discover",
    r"\bpublished\b",
    r"\bproves?\b",
    r"\bGO\b",
    r"no[- ]?go",
    r"risk[ _-]score",
    r"INSTRUMENT VALID",
    r"\bTier\b",
    r"\breportab",
    r"\bconcordance\b",
    r"\bJaccard\b",
    r"\boverlap\b",
)


def _expected_envelope_rows() -> str:
    rows = []
    for row in gate_config.OPERATING_ENVELOPE:
        rows.append(
            f"sigma_donor {row['sigma_donor']}: at least {row['min_donors_per_group']} donors "
            f"per group (power at least {gate_config.POWER_TARGET} at log2FC "
            f"{gate_config.ORACLE_LOG2FC} in {gate_config.ORACLE_K} genes; grid support "
            f"'{row['grid_support']}')"
        )
    return "\n".join(rows)


def _readout_kwargs(**overrides):
    kwargs = {
        "floor_solo": "1162",
        "universe_size": 1500,
        "alpha": "0.05",
        "floor_pct": "77.5",
        "real_solo": 1301,
        "few_donors": False,
        "paired_floor_shown": True,
        "pb_real": 12,
        "pb_floor": "0",
    }
    kwargs.update(overrides)
    return kwargs


def _all_fixed_prose() -> dict[str, str]:
    """Every fixed string this module can put in a report, with dummy values interpolated."""
    dummy = collections.defaultdict(lambda: "X")
    prose = {
        f"CAVEATS[{key}]": value.format_map(dummy)
        for key, value in text.CAVEATS.items()
        if key != "R1"
    }
    prose["CAVEATS[N2] rendered"] = text.envelope_sentence()
    for few_donors, paired in itertools.product((False, True), repeat=2):
        rendered = text.readout_caveat_text(
            **_readout_kwargs(few_donors=few_donors, paired_floor_shown=paired)
        )
        prose[f"R1(few_donors={few_donors}, paired={paired})"] = rendered
    prose.update({f"SENTENCES[{k}]": v.format_map(dummy) for k, v in text.SENTENCES.items()})
    prose.update({f"GLOSSARY[{term}]": sentence for term, sentence in text.GLOSSARY})
    prose.update({f"REPORT_LINES[{k}]": v.format_map(dummy) for k, v in text.REPORT_LINES.items()})
    prose.update({f"STATUS_REASON_WORDS[{k}]": v for k, v in text.STATUS_REASON_WORDS.items()})
    return prose


# --------------------------------------------------------------------------- the envelope rows


def test_envelope_rows_render_one_row_per_envelope_point_with_its_meaning():
    rows = text.envelope_rows().split("\n")
    assert len(rows) == len(gate_config.OPERATING_ENVELOPE)
    for row, declared in zip(rows, gate_config.OPERATING_ENVELOPE, strict=True):
        assert row.startswith(f"sigma_donor {declared['sigma_donor']}: at least ")
        assert f"at least {declared['min_donors_per_group']} donors per group" in row
        assert f"power at least {gate_config.POWER_TARGET}" in row
        assert f"log2FC {gate_config.ORACLE_LOG2FC} in {gate_config.ORACLE_K} genes" in row
        assert f"grid support '{declared['grid_support']}'" in row
    assert text.envelope_rows() == _expected_envelope_rows()


def test_no_envelope_prose_calls_a_derived_or_extrapolated_row_established():
    # The word has to be absent from the sentence that encloses the rows and from the glossary
    # entry that repeats it, not only from the rows: docs/AMENDMENTS.md records two of these
    # points as derived or extrapolated rather than measured, so any sentence covering all rows
    # would assert more than the grid supports.
    assert "established" not in text.envelope_rows()
    assert "established" not in text.envelope_sentence()
    for term, sentence in text.GLOSSARY:
        assert "established" not in sentence, f"GLOSSARY[{term}] says power was established"


def test_envelope_sentence_matches_the_n2_template_filled_from_gate_config():
    expected = text.CAVEATS["N2"].format_map({
        "envelope_rows": _expected_envelope_rows(),
        "calibration_sigma": gate_config.CALIBRATION_EVAL_SIGMA,
        "oracle_donors": gate_config.ORACLE_SIM["n_donors_per_group"],
    })
    assert text.envelope_sentence() == expected


def test_module_types_no_protocol_number_by_hand():
    source = inspect.getsource(text)
    for row in gate_config.OPERATING_ENVELOPE:
        for value in (row["sigma_donor"], row["min_donors_per_group"]):
            token = re.escape(str(value))
            assert re.search(rf"\b{token}\b", source) is None, (
                f"render/text.py contains the literal envelope value {value!r}"
            )
    for value in (
        gate_config.POWER_TARGET,
        gate_config.ORACLE_K,
        gate_config.CALIBRATION_EVAL_SIGMA,
        gate_config.ORACLE_SIM["n_genes"],
        product_constants.FEW_DONORS_THRESHOLD,
    ):
        token = re.escape(str(value))
        assert re.search(rf"\b{token}\b", source) is None, (
            f"render/text.py contains the literal protocol or product value {value!r}"
        )


# --------------------------------------------------------------------------- the wording gate


def test_forbidden_patterns_are_the_pinned_literal_list():
    assert text.FORBIDDEN_PATTERNS == EXPECTED_FORBIDDEN_PATTERNS


def test_forbidden_patterns_are_compiled_once_case_insensitively():
    compiled = text.compiled_forbidden_patterns()
    assert compiled is text.compiled_forbidden_patterns()
    assert [pattern.pattern for pattern in compiled] == list(text.FORBIDDEN_PATTERNS)
    assert all(pattern.flags & re.IGNORECASE for pattern in compiled)


def test_the_valid_within_whitelist_is_gone():
    assert text.forbidden_pattern_hits("the instrument is valid within a stated envelope")
    assert text.forbidden_pattern_hits("VALIDATION of the instrument")
    assert text.forbidden_pattern_hits("no go")
    assert text.forbidden_pattern_hits("risk-score")
    assert text.forbidden_pattern_hits("the published cell")


def test_every_fixed_text_passes_the_forbidden_pattern_gate():
    for name, prose in _all_fixed_prose().items():
        assert text.forbidden_pattern_hits(prose) == (), (
            f"{name} matches a forbidden pattern: {text.forbidden_pattern_hits(prose)}"
        )


def test_status_reason_words_state_the_whole_reason_for_a_paired_design():
    """A donor measured under both conditions: the report says why that stops the run, not only
    that it does (the design is paired, no paired or mixed model exists here, and the donor
    permutation null would treat one donor's two halves as independent)."""
    words = text.STATUS_REASON_WORDS["donor_spans_conditions"]
    assert "paired design" in words
    assert "paired or mixed model" in words
    assert "does not implement" in words
    assert "independent" in words


def test_report_line_quotes_the_audited_file_name():
    """The file name is user input and reaches the header, so it is quoted and masked like a
    column name: a file called after a forbidden word cannot forge a sentence."""
    line = text.report_line("header_title", file="GO.h5ad")
    assert "'GO.h5ad'" in line
    assert text.forbidden_pattern_hits(line) == ()


def test_no_template_opens_a_quoted_span_of_its_own():
    templates = {
        **dict(text.CAVEATS),
        **dict(text.R1_CLAUSES),
        **text.SENTENCES,
        "envelope_row": text.ENVELOPE_ROW_TEMPLATE,
    }
    for name, template in templates.items():
        assert text.prose_for_pattern_check(template) == template, (
            f"template {name} carries a quote the mask would read as a delimiter, so a quoted "
            "span could hide part of pbcheck's own sentence"
        )


def test_user_values_are_quoted_and_masked_so_a_hostile_column_name_cannot_forge_prose():
    rendered = text.caveat_text("N8", col="GO")
    assert "'GO'" in rendered
    assert re.search(r"\bGO\b", rendered) is not None
    assert text.forbidden_pattern_hits(rendered) == ()

    reason = text.caveat_text("N6", reason="dtype float64 has findings")
    assert text.forbidden_pattern_hits(reason) == ()


def test_quoted_cannot_be_escaped_by_the_value_itself():
    hostile = text.quoted("GO'\nand these are findings")
    assert "\n" not in hostile
    assert hostile.count("'") == 2
    assert text.prose_for_pattern_check(f"column {hostile}") == "column ''"

    long_value = text.quoted("A" * (text.QUOTED_SPAN_LIMIT * 2))
    assert len(long_value) <= text.MAX_QUOTED_VALUE_CHARS + 2
    assert text.prose_for_pattern_check(f"column {long_value}") == "column ''"


def test_prose_for_pattern_check_masks_only_quoted_spans():
    assert text.prose_for_pattern_check("a 'b' c") == "a '' c"
    assert text.prose_for_pattern_check("no quotes here") == "no quotes here"


def test_a_possessive_apostrophe_does_not_open_a_masked_span():
    prose = "the arm's GO band and pbcheck's own rule"
    assert text.prose_for_pattern_check(prose) == prose
    assert r"\bGO\b" in text.forbidden_pattern_hits(prose)


# --------------------------------------------------------------------------- the caveat texts


def test_caveat_ids_present_are_exactly_the_schema_enum():
    assert set(text.CAVEATS.keys()) == set(CAVEAT_IDS)


def test_caveats_is_read_only():
    with pytest.raises(TypeError):
        text.CAVEATS["N1"] = "something else"


def test_caveat_text_raises_keyerror_naming_the_missing_placeholder():
    with pytest.raises(KeyError) as excinfo:
        text.caveat_text("N5")
    assert "cols" in str(excinfo.value)


def test_caveat_text_formats_with_all_placeholders_supplied():
    rendered = text.caveat_text("N9", n=3, requested=200, se=0.021)
    assert "3" in rendered
    assert "200" in rendered
    assert "0.021" in rendered


def test_n4_interpolates_the_product_threshold_and_names_its_source():
    rendered = text.caveat_text("N4")
    assert f"fewer than {product_constants.FEW_DONORS_THRESHOLD} donors" in rendered
    assert (
        "Amendment 5 Change 2 of pbcheck's protocol admits floor-based quantities outside the "
        "envelope only when every group has at least "
        f"{product_constants.FEW_DONORS_THRESHOLD} donors" in rendered
    )
    assert "min(" not in rendered


def test_n7_states_provenance_without_an_instrument_verdict():
    rendered = text.caveat_text("N7", protocol_constant_names="alpha, lambda_band")
    assert text.GATE_ARTIFACT in rendered
    assert f"sigma_donor {gate_config.CALIBRATION_EVAL_SIGMA}" in rendered
    donors = gate_config.ORACLE_SIM["n_donors_per_group"]
    assert f"{donors} against {donors} donors" in rendered
    assert f"{gate_config.ORACLE_SIM['n_genes']} genes" in rendered
    assert "that measurement is not repeated on this file" in rendered
    for claim in ("calibrat", "ready", "instrument is"):
        assert claim not in rendered.lower()


def test_n2_scopes_power_to_the_envelope_and_calibration_to_one_regime():
    rendered = text.envelope_sentence()
    assert rendered.startswith("Amendment 3 declares an operating envelope for the pseudobulk arm")
    assert "the donor count per group at which the power target is reached" in rendered
    assert "calibration was evaluated at one hard regime" in rendered
    assert "pbcheck does not estimate sigma_donor for real data" in rendered


# --------------------------------------------------------------------------- the R1 paragraph


def test_r1_names_the_bh_convention_of_every_count():
    rendered = text.readout_caveat_text(**_readout_kwargs())
    assert "corrected over the whole gene universe on its own (solo BH)" in rendered
    assert "corrected the same way, it calls 1301" in rendered
    assert "genes they have in common (paired BH)" in rendered


def test_r1_renders_the_pseudobulk_clause_only_when_the_paired_floor_is_shown():
    shown = text.readout_caveat_text(**_readout_kwargs(paired_floor_shown=True))
    hidden = text.readout_caveat_text(**_readout_kwargs(paired_floor_shown=False))
    assert "donor-pseudobulk test called 12 genes" in shown
    assert "donor-pseudobulk" not in hidden
    assert hidden.endswith("The donor is the replication unit this design supports.")


def test_r1_replaces_the_categorical_clause_below_the_donor_threshold():
    enough = text.readout_caveat_text(**_readout_kwargs(few_donors=False))
    few = text.readout_caveat_text(**_readout_kwargs(few_donors=True))
    assert "cannot be separated from that floor" in enough
    assert "cannot be separated from that floor" not in few
    assert (
        f"With fewer than {product_constants.FEW_DONORS_THRESHOLD} donors in a group" in few
    )
    assert "is not interpretable on this run" in few


def test_r1_requires_the_pseudobulk_numbers_when_the_paired_floor_is_shown():
    with pytest.raises(ValueError):
        text.readout_caveat_text(**_readout_kwargs(pb_real=None))


def test_caveat_text_refuses_to_assemble_r1():
    with pytest.raises(ValueError):
        text.caveat_text("R1")


# --------------------------------------------------------------------------- the rest


def test_lambda_class_words_are_the_schema_classes_with_neutral_phrases():
    assert tuple(text.LAMBDA_CLASS_WORDS) == LAMBDA_CLASS_VALUES
    assert dict(text.LAMBDA_CLASS_WORDS) == {
        "in_band": "inside the band",
        "above_band": "above the band",
        "below_band": "below the band",
    }
    with pytest.raises(TypeError):
        text.LAMBDA_CLASS_WORDS["in_band"] = "some other phrase"


def test_the_naive_lambda_sentence_says_the_band_is_the_pseudobulk_arms():
    sentence = text.SENTENCES["lambda_naive"]
    assert "the donor-pseudobulk arm's band" in sentence
    assert "not as the naive arm's own criterion" in sentence


def test_sentence_text_quotes_the_user_supplied_condition_levels():
    rendered = text.sentence_text(
        "donors", n_test=4, test_level="GO", n_ref=4, ref_level="ctrl", n_distinct_splits=70
    )
    assert "'GO'" in rendered
    assert text.forbidden_pattern_hits(rendered) == ()


def test_the_module_source_is_ascii_and_has_no_em_dash():
    source = inspect.getsource(text)
    assert source.isascii()
    assert "\u2014" not in source
