"""Tests for the report caveat texts (plan section 1.4, WP0)."""

from __future__ import annotations

import collections
import inspect
import re

import pytest

from pbcheck import gate_config
from pbcheck.render import text


def _expected_envelope_sentence() -> str:
    donors = " / ".join(str(row["min_donors_per_group"]) for row in gate_config.OPERATING_ENVELOPE)
    sigmas = " / ".join(str(row["sigma_donor"]) for row in gate_config.OPERATING_ENVELOPE)
    envelope_rows = f"minimum donors per group {donors} at sigma_donor {sigmas}"
    return text.CAVEATS["C2"].format_map({"envelope_rows": envelope_rows})


def test_envelope_sentence_matches_the_c2_template_filled_from_gate_config():
    assert text.envelope_sentence() == _expected_envelope_sentence()


def test_envelope_rows_matches_the_fragment_used_by_envelope_sentence():
    assert text.envelope_rows() in text.envelope_sentence()


def test_envelope_builders_type_no_hardcoded_envelope_numbers():
    # C4's "8 donors" is an unrelated fact (the product's FEW_DONORS_THRESHOLD) that happens to
    # share a digit with one envelope row, so the "no envelope number typed by hand" rule is
    # checked against the functions that build the envelope text, not the whole module: those are
    # exactly the ones the WP0 contract requires to derive every number from gate_config.
    source = inspect.getsource(text.envelope_rows) + inspect.getsource(text.envelope_sentence)
    for row in gate_config.OPERATING_ENVELOPE:
        for value in (row["sigma_donor"], row["min_donors_per_group"]):
            token = re.escape(str(value))
            assert re.search(rf"\b{token}\b", source) is None, (
                f"envelope_rows/envelope_sentence contains the literal envelope value {value!r}"
            )


def test_caveats_formatted_with_dummy_values_match_no_forbidden_pattern():
    dummy = collections.defaultdict(lambda: "X")
    for caveat_id, template in text.CAVEATS.items():
        rendered = template.format_map(dummy)
        for pattern in text.FORBIDDEN_PATTERNS:
            assert re.search(pattern, rendered) is None, (
                f"caveat {caveat_id} matches forbidden pattern {pattern!r}: {rendered!r}"
            )


def test_caveat_text_raises_keyerror_naming_the_missing_placeholder():
    with pytest.raises(KeyError) as excinfo:
        text.caveat_text("C5")
    assert "cols" in str(excinfo.value)


def test_caveat_text_formats_with_all_placeholders_supplied():
    rendered = text.caveat_text("C9", n=3, requested=200, se=0.021)
    assert "3" in rendered
    assert "200" in rendered
    assert "0.021" in rendered


def test_lambda_class_words_is_the_fixed_triple():
    assert text.LAMBDA_CLASS_WORDS == ("calibrated", "inflated", "under")


def test_caveat_ids_present_are_exactly_the_schema_enum():
    from pbcheck.audit_schema import CAVEAT_IDS

    assert set(text.CAVEATS.keys()) == set(CAVEAT_IDS)
