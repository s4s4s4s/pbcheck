"""Tests for the pbcheck-audit/1 schema validator (plan section WP0)."""

from __future__ import annotations

import copy

import pytest

from pbcheck.audit_schema import (
    AUDIT_SCHEMA_VERSION,
    AuditSchemaError,
    empty_payload,
    validate,
)


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


def test_caveat_ids_are_restricted_to_the_enum():
    payload = empty_payload()
    payload["caveats"] = [{"id": "C1", "text": "fine"}]
    validate(payload)

    payload["caveats"] = [{"id": "C99", "text": "not a real caveat id"}]
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert "caveats" in str(excinfo.value)


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


def test_non_null_real_label_is_validated_field_by_field():
    payload = empty_payload()
    payload["real_label"] = {
        "naive": {
            "n_significant_paired": None,
            "n_significant_solo": 0,
            "n_tested": 0,
            "top": [],
        },
        "pseudobulk": None,
        "paired_bh": {
            "n_universe": 0,
            "n_tested_common": 0,
            "n_na_naive": 0,
            "n_na_pseudobulk": 0,
            "n_dropped_for_fairness": 0,
            "pseudobulk_na_free": True,
        },
        "consistent_with_permutation_path": None,
    }
    validate(payload)

    broken = copy.deepcopy(payload)
    del broken["real_label"]["paired_bh"]["pseudobulk_na_free"]
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(broken)
    assert str(excinfo.value).startswith("real_label.paired_bh.pseudobulk_na_free:")


def test_universe_builder_enum_rejects_unknown_value():
    payload = empty_payload()
    payload["universe"]["builder"] = "not_a_builder"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("universe.builder:")


def test_counts_source_accepts_layers_prefix_and_rejects_other_strings():
    payload = empty_payload()
    payload["input"]["counts_source"] = "layers:counts"
    validate(payload)

    payload["input"]["counts_source"] = "not_a_valid_source"
    with pytest.raises(AuditSchemaError) as excinfo:
        validate(payload)
    assert str(excinfo.value).startswith("input.counts_source:")


def test_top_level_payload_must_be_a_dict():
    with pytest.raises(AuditSchemaError):
        validate(["not", "a", "dict"])
