"""Tests for the pbcheck-audit/1 schema validator: the payload's shape, and the two rules that
make the JSON artifact as trustworthy as the rendered report (the caveat block and the readout
flags that switch conditional prose on)."""

from __future__ import annotations

import copy
from types import MappingProxyType

import pytest

from pbcheck import gate_config, io_counts, product_constants
from pbcheck.audit_schema import (
    ALWAYS_ON_CAVEAT_IDS,
    AUDIT_SCHEMA_VERSION,
    AuditSchemaError,
    empty_payload,
    validate,
)
from pbcheck.render import text


def _caveat(caveat_id: str, **values) -> dict:
    return {"id": caveat_id, "text": text.caveat_text(caveat_id, **values)}


def _floor(bh_mode: str) -> dict:
    return {
        "median_count": 1.0,
        "iqr_count": 0.0,
        "median_frac": 0.001,
        "max_count": 2.0,
        "n_perm": 200,
        "bh_mode": bh_mode,
        "mc_se": 0.01,
    }


def _permutation_null() -> dict:
    return {
        "naive": {
            "lambda": 1.0,
            "lambda_iqr": 0.1,
            "n_perm": 200,
            "b5_lambda_empirical": 1.0,
            "floor_solo": _floor("solo"),
            "floor_paired": None,
        },
        "pseudobulk": None,
        "monte_carlo": None,
        "real_split_inside_perm_range": None,
        "real_split_percentile_in_perms": None,
        "n_donors": 16,
        "n_distinct_splits": 70,
        "n_perm_naive_achieved": 200,
        "n_perm_pb_achieved": None,
        "n_perm_paired": None,
        "engine_path": "run_null",
    }


def _naive_only_real_label() -> dict:
    return {
        "naive": {
            "n_significant_paired": None,
            "n_significant_solo": 0,
            "n_tested": 0,
            "top": [],
        },
        "pseudobulk": None,
        "paired_bh": None,
        "consistent_with_permutation_path": None,
    }


def _paired_bh(n_na_pseudobulk: int = 0) -> dict:
    return {
        "n_universe": 1500,
        "n_tested_common": 1500,
        "n_na_naive": 0,
        "n_na_pseudobulk": n_na_pseudobulk,
        "n_dropped_for_fairness": 0,
        "pseudobulk_na_free": n_na_pseudobulk == 0,
    }


# --------------------------------------------------------------------------- shape


def test_empty_payload_validates():
    validate(empty_payload())


def test_empty_payload_has_the_declared_schema_version():
    assert empty_payload()["schema_version"] == AUDIT_SCHEMA_VERSION


@pytest.mark.parametrize("key", sorted(empty_payload().keys()))
def test_missing_top_level_key_is_rejected_and_named(key):
    payload = empty_payload()
    del payload[key]
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert key in str(excinfo.value)


def test_wrong_status_literal_is_rejected():
    payload = empty_payload()
    payload["status"] = "definitely_not_a_status"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "status" in str(excinfo.value)


def test_wrong_status_reason_literal_is_rejected():
    payload = empty_payload()
    payload["status_reason"] = "not_a_gate_name"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "status_reason" in str(excinfo.value)


def test_status_reason_null_is_accepted():
    payload = empty_payload()
    payload["status_reason"] = None
    validate(payload)


def test_nested_type_error_names_its_dotted_path():
    payload = empty_payload()
    payload["design"]["n_cells"] = "not an int"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("design.n_cells:")


def test_deeply_nested_type_error_names_its_dotted_path():
    payload = empty_payload()
    payload["settings"]["protocol_constants"]["alpha"] = "not a float"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("settings.protocol_constants.alpha:")


def test_nullable_block_may_be_none():
    payload = empty_payload()
    assert payload["real_label"] is None
    assert payload["permutation_null"] is None
    assert payload["counts_check"] is None
    validate(payload)


def test_wrong_type_in_a_nullable_block_is_rejected():
    payload = empty_payload()
    payload["real_label"] = 5
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("real_label:")


def test_unknown_key_is_rejected_at_every_level():
    payload = empty_payload()
    payload["extra_unknown_key"] = 1
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "extra_unknown_key" in str(excinfo.value)

    payload = empty_payload()
    payload["readout"]["bogus"] = 1
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("readout.bogus:")


def test_top_level_payload_must_be_a_dict():
    with pytest.raises(AuditSchemaError):
        validate(["not", "a", "dict"])


# --------------------------------------------------------------------------- the blocks


def test_provenance_accepts_the_frozen_constants_as_they_are_declared():
    payload = empty_payload()
    payload["provenance"]["pre_registered"] = MappingProxyType(dict(gate_config.PRE_REGISTERED))
    payload["provenance"]["operating_envelope"] = tuple(
        MappingProxyType(dict(row)) for row in gate_config.OPERATING_ENVELOPE
    )
    validate(payload)


def test_lambda_band_is_exactly_two_floats():
    payload = empty_payload()
    payload["settings"]["protocol_constants"]["lambda_band"] = list(gate_config.LAMBDA_BAND)
    validate(payload)

    for wrong in ([0.9, 1.1, 1.3], []):
        payload["settings"]["protocol_constants"]["lambda_band"] = wrong
        with pytest.raises(AuditSchemaError) as excinfo:
            validate(payload)
        assert str(excinfo.value).startswith("settings.protocol_constants.lambda_band:")


def test_counts_source_accepts_a_named_layer_and_rejects_the_bare_prefix():
    payload = empty_payload()
    payload["input"]["counts_source"] = "layers:counts"
    validate(payload)

    for wrong in ("not_a_valid_source", "layers:", "layers:   "):
        payload["input"]["counts_source"] = wrong
        with pytest.raises(AuditSchemaError) as excinfo:
            validate(payload)
        assert str(excinfo.value).startswith("input.counts_source:")


def test_counts_check_takes_an_integer_check_field_for_field():
    check = io_counts.IntegerCheck(
        passed=False,
        reason="non-integer values",
        dtype="float32",
        layer="X",
        n_cells=10,
        n_genes=20,
        n_values_checked=200,
        n_noninteger=3,
        n_negative=0,
        n_nonfinite=0,
        max_abs_value=12.5,
        examples=(1.5, 2.5),
        notes=("a note",),
    )
    payload = empty_payload()
    payload["counts_check"] = check.as_dict()
    validate(payload)


def test_universe_builder_enum_rejects_unknown_value():
    payload = empty_payload()
    payload["universe"]["builder"] = "not_a_builder"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("universe.builder:")


def test_lambda_class_enum_holds_the_neutral_band_names():
    payload = empty_payload()
    for value in ("in_band", "above_band", "below_band"):
        payload["readout"]["lambda_naive_class"] = value
        validate(payload)

    payload["readout"]["lambda_naive_class"] = "calibrated"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("readout.lambda_naive_class:")


def test_engine_path_and_floor_bh_mode_enums_reject_unknown_values():
    payload = empty_payload()
    payload["permutation_null"] = _permutation_null()
    validate(payload)

    broken = copy.deepcopy(payload)
    broken["permutation_null"]["engine_path"] = "some_other_engine"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(broken)
    assert str(excinfo.value).startswith("permutation_null.engine_path:")

    broken = copy.deepcopy(payload)
    broken["permutation_null"]["naive"]["floor_solo"]["bh_mode"] = "paired"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(broken)
    assert str(excinfo.value).startswith("permutation_null.naive.floor_solo.bh_mode:")


def test_real_label_paired_bh_is_nullable_for_the_naive_only_shape():
    payload = empty_payload()
    payload["status"] = "naive_only"
    payload["status_reason"] = "non_integer_counts"
    payload["real_label"] = _naive_only_real_label()
    payload["caveats"].append(_caveat("N6", reason="dtype float32 is not integral"))
    validate(payload)


def test_non_null_real_label_is_validated_field_by_field():
    payload = empty_payload()
    payload["real_label"] = _naive_only_real_label()
    payload["real_label"]["paired_bh"] = _paired_bh()
    validate(payload)

    broken = copy.deepcopy(payload)
    del broken["real_label"]["paired_bh"]["pseudobulk_na_free"]
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(broken)
    assert str(excinfo.value).startswith("real_label.paired_bh.pseudobulk_na_free:")


def test_settings_tool_carries_the_product_constants_the_report_documents():
    payload = empty_payload()
    for key in ("min_donors_per_group", "universe_min_total_count", "universe_min_prop"):
        assert key in payload["settings"]["tool"]


# --------------------------------------------------------------------------- the caveat block


def test_caveat_ids_are_restricted_to_the_enum():
    payload = empty_payload()
    payload["caveats"].append({"id": "C1", "text": "an id from the retired scheme"})
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "caveats" in str(excinfo.value)


def test_every_payload_carries_the_always_on_notes():
    assert ALWAYS_ON_CAVEAT_IDS == ("N1", "N2", "N3", "N7")
    for missing in ALWAYS_ON_CAVEAT_IDS:
        payload = empty_payload()
        payload["caveats"] = [
            entry for entry in payload["caveats"] if entry["id"] != missing
        ]
        with pytest.raises(AuditSchemaError) as excinfo:
            validate(payload)
        assert f"missing required caveat {missing}" in str(excinfo.value)


def test_a_caveat_may_not_appear_twice():
    payload = empty_payload()
    payload["caveats"].append(_caveat("N3"))
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "N3 appears more than once" in str(excinfo.value)


def test_n6_is_present_exactly_when_the_counts_gate_failed():
    payload = empty_payload()
    payload["status_reason"] = "non_integer_counts"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "N6" in str(excinfo.value)

    payload["caveats"].append(_caveat("N6", reason="dtype float32 is not integral"))
    validate(payload)

    payload["status_reason"] = "too_few_donors"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "N6" in str(excinfo.value)


def test_n4_is_present_exactly_when_the_readout_says_few_donors():
    payload = empty_payload()
    payload["readout"]["few_donors_threshold"] = product_constants.FEW_DONORS_THRESHOLD
    payload["readout"]["min_donors_per_group"] = 3
    payload["readout"]["few_donors"] = True
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "N4" in str(excinfo.value)

    payload["caveats"].append(_caveat("N4"))
    validate(payload)

    payload["readout"]["min_donors_per_group"] = product_constants.FEW_DONORS_THRESHOLD
    payload["readout"]["few_donors"] = False
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "N4" in str(excinfo.value)


def test_a_caveat_text_matching_a_forbidden_pattern_is_rejected():
    payload = empty_payload()
    payload["caveats"][0]["text"] = "This instrument is valid and its findings are published."
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("caveats[0].text:")


def test_a_quoted_identifier_in_a_caveat_text_is_not_a_forbidden_pattern():
    payload = empty_payload()
    payload["caveats"].append(_caveat("N8", col="GO"))
    validate(payload)


# --------------------------------------------------------------------------- the readout flags


def test_few_donors_must_agree_with_the_donor_count_and_the_threshold():
    payload = empty_payload()
    payload["readout"]["few_donors_threshold"] = product_constants.FEW_DONORS_THRESHOLD
    payload["readout"]["min_donors_per_group"] = 3
    payload["readout"]["few_donors"] = False
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("readout.few_donors:")


def test_paired_floor_shown_requires_a_paired_correction_free_of_missing_values():
    payload = empty_payload()
    payload["readout"]["paired_floor_shown"] = True
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("readout.paired_floor_shown:")

    payload["real_label"] = _naive_only_real_label()
    payload["real_label"]["paired_bh"] = _paired_bh(n_na_pseudobulk=7)
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("readout.paired_floor_shown:")

    payload["real_label"]["paired_bh"] = _paired_bh(n_na_pseudobulk=0)
    validate(payload)
