"""The pbcheck-audit/1 payload schema and its validator.

This module defines the shape of ``pbcheck_audit.json`` and checks a payload against it: it does
not compute an audit, it tells the caller whether a payload has every key the schema requires and
no key it does not, of the right type, with the right literal values where the schema fixes them.

Two rules go past the shape, because the JSON artifact is read by machines that will never see the
rendered report and must be no easier to mislead than a reader: the caveat block has to carry the
notes the report is defined to always carry and exactly the conditional notes the payload's own
fields call for, with no text matching the report's forbidden patterns; and the readout flags that
switch conditional prose on have to agree with the fields they are derived from.

``pbcheck.audit`` fills a payload matching this schema and ``pbcheck.render`` reads it back;
neither is imported at module level. The caveat rules need the report's wording gate
(``pbcheck.render.text``), which is imported inside the two functions that use it, because
``pbcheck.render`` imports this module.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from pbcheck.product_constants import FEW_DONORS_THRESHOLD

AUDIT_SCHEMA_VERSION = "pbcheck-audit/1"


class AuditSchemaError(ValueError):
    """A payload does not match the pbcheck-audit/1 schema.

    The message names the dotted key path of the first violation found, e.g.
    ``"design.n_cells: expected int, got str"`` or ``"status_reason: missing required key"``.
    """


def _fail(path: str, message: str) -> None:
    raise AuditSchemaError(f"{path}: {message}" if path else message)


def _join(path: str, key: str) -> str:
    return f"{path}.{key}" if path else key


class _Spec:
    """Base class of the small validator DSL used to describe the schema below."""

    def check(self, value: Any, path: str) -> None:
        raise NotImplementedError


class _Str(_Spec):
    def check(self, value: Any, path: str) -> None:
        if not isinstance(value, str):
            _fail(path, f"expected str, got {type(value).__name__}")


class _Int(_Spec):
    def check(self, value: Any, path: str) -> None:
        if isinstance(value, bool) or not isinstance(value, int):
            _fail(path, f"expected int, got {type(value).__name__}")


class _Float(_Spec):
    def check(self, value: Any, path: str) -> None:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            _fail(path, f"expected float, got {type(value).__name__}")


class _Bool(_Spec):
    def check(self, value: Any, path: str) -> None:
        if not isinstance(value, bool):
            _fail(path, f"expected bool, got {type(value).__name__}")


def _is_mapping(value: Any) -> bool:
    """Any mapping, not only ``dict``: the engine hands over ``MappingProxyType`` constants."""
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    """Any sequence but a string or bytes: the engine hands over tuples as well as lists."""
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


class _AnyDict(_Spec):
    """A mapping whose shape this schema does not pin (recorded verbatim from the engine)."""

    def check(self, value: Any, path: str) -> None:
        if not _is_mapping(value):
            _fail(path, f"expected dict, got {type(value).__name__}")


class _AnyList(_Spec):
    def check(self, value: Any, path: str) -> None:
        if not _is_sequence(value):
            _fail(path, f"expected list, got {type(value).__name__}")


class _Enum(_Spec):
    def __init__(self, *values: Any) -> None:
        self.values = values

    def check(self, value: Any, path: str) -> None:
        if value not in self.values:
            _fail(path, f"expected one of {self.values!r}, got {value!r}")


class _CountsSource(_Spec):
    """``"X"`` | ``"layers:<name>"`` | ``"raw.X"``; ``null`` is handled by the caller's Nullable."""

    def check(self, value: Any, path: str) -> None:
        if not isinstance(value, str):
            _fail(path, f"expected str, got {type(value).__name__}")
            return
        if value == "X" or value == "raw.X":
            return
        if value.startswith("layers:") and value[len("layers:"):].strip():
            return
        _fail(path, f'expected "X", "layers:<name>" or "raw.X", got {value!r}')


class _Nullable(_Spec):
    def __init__(self, inner: _Spec) -> None:
        self.inner = inner

    def check(self, value: Any, path: str) -> None:
        if value is None:
            return
        self.inner.check(value, path)


class _ListOf(_Spec):
    def __init__(self, item: _Spec) -> None:
        self.item = item

    def check(self, value: Any, path: str) -> None:
        if not _is_sequence(value):
            _fail(path, f"expected list, got {type(value).__name__}")
            return
        for index, element in enumerate(value):
            self.item.check(element, f"{path}[{index}]")


class _SizedListOf(_ListOf):
    """A sequence of a fixed length, each element validated by the same spec."""

    def __init__(self, item: _Spec, length: int) -> None:
        super().__init__(item)
        self.length = length

    def check(self, value: Any, path: str) -> None:
        super().check(value, path)
        if len(value) != self.length:
            _fail(path, f"expected {self.length} items, got {len(value)}")


class _DictOf(_Spec):
    """A mapping with str keys and uniformly-typed values, with no fixed key set."""

    def __init__(self, value: _Spec) -> None:
        self.value = value

    def check(self, value: Any, path: str) -> None:
        if not _is_mapping(value):
            _fail(path, f"expected dict, got {type(value).__name__}")
            return
        for key, element in value.items():
            if not isinstance(key, str):
                _fail(path, f"expected str keys, got {type(key).__name__}")
            self.value.check(element, f"{path}.{key}")


class _DictSchema(_Spec):
    """A mapping with a fixed, required key set, each key validated by its own spec."""

    def __init__(self, fields: dict[str, _Spec]) -> None:
        self.fields = fields

    def check(self, value: Any, path: str) -> None:
        if not _is_mapping(value):
            _fail(path, f"expected dict, got {type(value).__name__}")
            return
        unknown = sorted(str(key) for key in value if key not in self.fields)
        if unknown:
            _fail(_join(path, unknown[0]), "unknown key, not in the pbcheck-audit/1 schema")
        for key, spec in self.fields.items():
            sub_path = _join(path, key)
            if key not in value:
                _fail(sub_path, "missing required key")
            spec.check(value[key], sub_path)


_STR = _Str()
_INT = _Int()
_FLOAT = _Float()
_BOOL = _Bool()
_ANY_DICT = _AnyDict()
_ANY_LIST = _AnyList()
_COUNTS_SOURCE = _CountsSource()

STATUS_VALUES = ("complete", "naive_only", "design_only")
STATUS_REASON_VALUES = (
    "design_only_requested",
    "donor_spans_conditions",
    "too_few_donors",
    "non_integer_counts",
    "universe_too_small",
    "too_few_profiles_after_thin_filter",
)
UNIVERSE_BUILDER_VALUES = ("pseudobulk_frozen", "naive_detection_fallback")
LAMBDA_CLASS_VALUES = ("in_band", "above_band", "below_band")
ENGINE_PATH_VALUES = ("run_null", "naive_null")
CAVEAT_IDS = ("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "R1")

#: Notes every payload carries, whatever its status: the scope banner, the envelope note, the
#: "none of this is a Phase 0 result" note and the provenance note. :func:`validate` requires them.
ALWAYS_ON_CAVEAT_IDS = ("N1", "N2", "N3", "N7")
RUNTIME_STAGE_KEYS = (
    "load",
    "prepare",
    "design",
    "counts",
    "pseudobulk_build",
    "universe",
    "real_label",
    "permutation_null",
    "render",
)
PACKAGE_KEYS = (
    "numpy",
    "scipy",
    "pandas",
    "anndata",
    "scanpy",
    "statsmodels",
    "decoupler",
    "pydeseq2",
)


def _floor_schema(bh_mode: str) -> _DictSchema:
    return _DictSchema({
        "median_count": _FLOAT,
        "iqr_count": _FLOAT,
        "median_frac": _FLOAT,
        "max_count": _FLOAT,
        "n_perm": _INT,
        "bh_mode": _Enum(bh_mode),
        "mc_se": _FLOAT,
    })


_INPUT_SCHEMA = _DictSchema({
    "path": _Nullable(_STR),
    "n_cells_loaded": _INT,
    "n_genes_loaded": _INT,
    "n_cells_audited": _INT,
    "n_cells_dropped_other_condition": _INT,
    "n_cells_dropped_other_celltype": _INT,
    "n_cells_dropped_missing_donor": _INT,
    "n_cells_dropped_missing_condition": _INT,
    "n_cells_dropped_missing_celltype": _INT,
    "donor_col": _STR,
    "condition_col": _STR,
    "test_level": _STR,
    "ref_level": _STR,
    "celltype_col": _Nullable(_STR),
    "celltype_value": _Nullable(_STR),
    "celltype_like_columns_found": _ListOf(_STR),
    "batch_cols": _ListOf(_STR),
    "counts_source": _Nullable(_COUNTS_SOURCE),
})

_COUNTS_CHECK_SCHEMA = _DictSchema({
    "passed": _BOOL,
    "reason": _Nullable(_STR),
    "dtype": _STR,
    "layer": _STR,
    "n_cells": _INT,
    "n_genes": _INT,
    "n_values_checked": _INT,
    "n_noninteger": _INT,
    "n_negative": _INT,
    "n_nonfinite": _INT,
    "max_abs_value": _Nullable(_FLOAT),
    "examples": _ANY_LIST,
    "notes": _ListOf(_STR),
    "rule": _STR,
})

_DESIGN_SCHEMA = _DictSchema({
    "condition_col": _STR,
    "donor_col": _STR,
    "n_cells": _INT,
    "groups": _DictOf(_INT),
    "donors_per_group": _DictOf(_INT),
    "cells_per_donor": _DictOf(_INT),
    "donor_nests_in_condition": _BOOL,
    "min_donors_per_group": _INT,
    "imbalance_ratio": _FLOAT,
    "batch_confounded": _DictOf(_FLOAT),
    "batch_separates_condition": _DictOf(_BOOL),
    "flags": _ListOf(_STR),
    "min_donors": _INT,
    "usable_for_pseudobulk": _BOOL,
})

_SETTINGS_TOOL_SCHEMA = _DictSchema({
    "n_perm_requested": _INT,
    "n_perm_pb_requested": _INT,
    "seed": _INT,
    "top_n": _INT,
    "design_only": _BOOL,
    "naive_method": _Enum("wilcoxon"),
    "naive_engine": _Enum("fast"),
    "pseudobulk_method": _Enum("moderated_ebayes"),
    "trend": _BOOL,
    "min_donors_per_group": _INT,
    "universe_min_total_count": _INT,
    "universe_min_prop": _FLOAT,
    "fallback_universe_min_prop": _FLOAT,
    "fallback_universe_min_size": _INT,
    "min_profiles_per_group_after_thin_filter": _INT,
})

_SETTINGS_PROTOCOL_CONSTANTS_SCHEMA = _DictSchema({
    "alpha": _FLOAT,
    "lambda_band": _SizedListOf(_FLOAT, 2),
    "min_universe_size": _INT,
    "min_cells": _INT,
    "min_counts": _INT,
})

_SETTINGS_SCHEMA = _DictSchema({
    "tool": _SETTINGS_TOOL_SCHEMA,
    "protocol_constants": _SETTINGS_PROTOCOL_CONSTANTS_SCHEMA,
})

_UNIVERSE_SCHEMA = _DictSchema({
    "size": _INT,
    "min_size": _INT,
    "builder": _Nullable(_Enum(*UNIVERSE_BUILDER_VALUES)),
    "builder_rule": _STR,
    "thin_donor_filter": _Nullable(_ANY_DICT),
    "profiles_per_group_after_thin_filter": _Nullable(_DictOf(_INT)),
})

_NAIVE_TOP_ITEM_SCHEMA = _DictSchema({
    "gene": _STR,
    "pval": _FLOAT,
    "padj": _FLOAT,
    "log2fc": _FLOAT,
    "pct_group": _FLOAT,
    "pct_reference": _FLOAT,
})

_PSEUDOBULK_TOP_ITEM_SCHEMA = _DictSchema({
    "gene": _STR,
    "pval": _FLOAT,
    "padj": _FLOAT,
    "log2fc": _FLOAT,
})

_REAL_LABEL_NAIVE_SCHEMA = _DictSchema({
    "n_significant_paired": _Nullable(_INT),
    "n_significant_solo": _INT,
    "n_tested": _INT,
    "top": _ListOf(_NAIVE_TOP_ITEM_SCHEMA),
})

_REAL_LABEL_PSEUDOBULK_SCHEMA = _DictSchema({
    "n_significant_paired": _INT,
    "n_tested": _INT,
    "top": _ListOf(_PSEUDOBULK_TOP_ITEM_SCHEMA),
    "moderation": _ANY_DICT,
})

_REAL_LABEL_PAIRED_BH_SCHEMA = _DictSchema({
    "n_universe": _INT,
    "n_tested_common": _INT,
    "n_na_naive": _INT,
    "n_na_pseudobulk": _INT,
    "n_dropped_for_fairness": _INT,
    "pseudobulk_na_free": _BOOL,
})

_REAL_LABEL_SCHEMA = _DictSchema({
    "naive": _REAL_LABEL_NAIVE_SCHEMA,
    "pseudobulk": _Nullable(_REAL_LABEL_PSEUDOBULK_SCHEMA),
    "paired_bh": _Nullable(_REAL_LABEL_PAIRED_BH_SCHEMA),
    "consistent_with_permutation_path": _Nullable(_BOOL),
})

_PERM_NULL_NAIVE_SCHEMA = _DictSchema({
    "lambda": _FLOAT,
    "lambda_iqr": _FLOAT,
    "n_perm": _INT,
    "b5_lambda_empirical": _FLOAT,
    "floor_solo": _floor_schema("solo"),
    "floor_paired": _Nullable(_floor_schema("paired")),
})

_PERM_NULL_PSEUDOBULK_SCHEMA = _DictSchema({
    "lambda": _FLOAT,
    "lambda_iqr": _FLOAT,
    "n_perm": _INT,
    "b5_lambda_empirical": _FLOAT,
    "floor": _floor_schema("paired"),
    "fp_rate": _FLOAT,
    "fp_rate_mc_se": _FLOAT,
})

_PERMUTATION_NULL_SCHEMA = _DictSchema({
    "naive": _PERM_NULL_NAIVE_SCHEMA,
    "pseudobulk": _Nullable(_PERM_NULL_PSEUDOBULK_SCHEMA),
    "monte_carlo": _Nullable(_ANY_DICT),
    "real_split_inside_perm_range": _Nullable(_BOOL),
    "real_split_percentile_in_perms": _Nullable(_FLOAT),
    "n_donors": _INT,
    "n_distinct_splits": _INT,
    "n_perm_naive_achieved": _INT,
    "n_perm_pb_achieved": _Nullable(_INT),
    "n_perm_paired": _Nullable(_INT),
    "engine_path": _Enum(*ENGINE_PATH_VALUES),
})

_READOUT_SCHEMA = _DictSchema({
    "lambda_naive_class": _Nullable(_Enum(*LAMBDA_CLASS_VALUES)),
    "lambda_pseudobulk_class": _Nullable(_Enum(*LAMBDA_CLASS_VALUES)),
    "naive_floor_solo": _Nullable(_DictSchema({
        "median_count": _FLOAT,
        "median_frac": _FLOAT,
        "mc_se": _FLOAT,
        "n_perm": _INT,
    })),
    "naive_real_solo": _Nullable(_INT),
    "naive_real_over_floor_solo": _Nullable(_FLOAT),
    "pseudobulk_real_over_floor": _Nullable(_FLOAT),
    "paired_floor_shown": _BOOL,
    "few_donors": _BOOL,
    "few_donors_threshold": _INT,
    "min_donors_per_group": _INT,
    "min_profiles_per_group_after_thin_filter": _Nullable(_INT),
    "n_perm_naive_requested": _INT,
    "n_perm_naive_achieved": _Nullable(_INT),
    "n_perm_pb_requested": _INT,
    "n_perm_pb_achieved": _Nullable(_INT),
    "n_distinct_splits": _Nullable(_INT),
    "sentences": _ListOf(_STR),
})

_CAVEATS_SCHEMA = _ListOf(_DictSchema({
    "id": _Enum(*CAVEAT_IDS),
    "text": _STR,
}))

_ENVELOPE_ROW_SCHEMA = _DictSchema({
    "sigma_donor": _FLOAT,
    "min_donors_per_group": _INT,
    "grid_support": _STR,
})

_PACKAGES_SCHEMA = _DictSchema({key: _STR for key in PACKAGE_KEYS})

_PROVENANCE_SCHEMA = _DictSchema({
    "pre_registered": _ANY_DICT,
    "pre_registered_source": _STR,
    "operating_envelope": _ListOf(_ENVELOPE_ROW_SCHEMA),
    "platform": _STR,
    "python": _STR,
    "packages": _PACKAGES_SCHEMA,
})

_RUNTIME_BY_STAGE_SCHEMA = _DictSchema({key: _Nullable(_FLOAT) for key in RUNTIME_STAGE_KEYS})

_TOP_LEVEL_SCHEMA = _DictSchema({
    "schema_version": _Enum(AUDIT_SCHEMA_VERSION),
    "pbcheck_version": _STR,
    "generated_utc": _STR,
    "runtime_seconds": _FLOAT,
    "runtime_by_stage_seconds": _RUNTIME_BY_STAGE_SCHEMA,
    "status": _Enum(*STATUS_VALUES),
    "status_reason": _Nullable(_Enum(*STATUS_REASON_VALUES)),
    "input": _INPUT_SCHEMA,
    "counts_check": _Nullable(_COUNTS_CHECK_SCHEMA),
    "design": _DESIGN_SCHEMA,
    "settings": _SETTINGS_SCHEMA,
    "universe": _UNIVERSE_SCHEMA,
    "real_label": _Nullable(_REAL_LABEL_SCHEMA),
    "permutation_null": _Nullable(_PERMUTATION_NULL_SCHEMA),
    "readout": _READOUT_SCHEMA,
    "caveats": _CAVEATS_SCHEMA,
    "provenance": _PROVENANCE_SCHEMA,
})


def _caveat_ids(payload: Mapping) -> list[str]:
    return [entry["id"] for entry in payload["caveats"]]


def _check_caveats(payload: Mapping) -> None:
    """The caveat block's own rules, beyond the per-entry types checked by the schema.

    Three things a payload cannot do and stay a pbcheck-audit/1 payload: carry a note twice, drop
    a note the report is defined to always carry, or carry a conditional note whose condition this
    same payload says is not met (and the reverse). The fourth rule is the wording gate: every
    caveat text is put through the report's own forbidden-pattern check, with quoted identifiers
    masked, so a text assembled outside :mod:`pbcheck.render.text` cannot reach a reader.
    """
    from pbcheck.render import text  # local: pbcheck.render imports this module's validator

    ids = _caveat_ids(payload)
    duplicates = sorted({entry for entry in ids if ids.count(entry) > 1})
    if duplicates:
        _fail("caveats", f"caveat {duplicates[0]} appears more than once")
    for required in ALWAYS_ON_CAVEAT_IDS:
        if required not in ids:
            _fail("caveats", f"missing required caveat {required}")

    expected_n6 = payload["status_reason"] == "non_integer_counts"
    if expected_n6 != ("N6" in ids):
        _fail("caveats", "N6 is present exactly when status_reason is non_integer_counts")
    expected_n4 = bool(payload["readout"]["few_donors"])
    if expected_n4 != ("N4" in ids):
        _fail("caveats", "N4 is present exactly when readout.few_donors is true")

    for index, entry in enumerate(payload["caveats"]):
        hits = text.forbidden_pattern_hits(entry["text"])
        if hits:
            _fail(
                f"caveats[{index}].text",
                f"caveat {entry['id']} matches the forbidden pattern {hits[0]!r}",
            )


def _check_readout_consistency(payload: Mapping) -> None:
    """The readout flags the report's conditional prose is switched on, against their own inputs.

    ``few_donors`` is the flag the donor-threshold note and the read-out paragraph's categorical
    clause hang on, and ``paired_floor_shown`` the flag its pseudobulk clause hangs on. Both are
    derivable from other fields of the same payload, so a payload that disagrees with itself here
    would render prose the numbers do not support.

    ``few_donors_threshold`` is checked against
    :data:`~pbcheck.product_constants.FEW_DONORS_THRESHOLD` first. It is a display copy of the
    tool's own constant, not a knob: without this check a payload could declare a threshold of its
    own choosing, keep ``few_donors`` false at any donor count and so switch the donor-threshold
    note and the read-out paragraph's categorical clause on and off from inside the payload.
    """
    readout = payload["readout"]
    if readout["few_donors_threshold"] != FEW_DONORS_THRESHOLD:
        _fail(
            "readout.few_donors_threshold",
            f"must be the tool's donor threshold {FEW_DONORS_THRESHOLD}, "
            f"got {readout['few_donors_threshold']}",
        )

    if readout["few_donors"] != (
        readout["min_donors_per_group"] < readout["few_donors_threshold"]
    ):
        _fail(
            "readout.few_donors",
            "must be true exactly when min_donors_per_group is below few_donors_threshold",
        )

    real_label = payload["real_label"]
    paired_bh = real_label["paired_bh"] if real_label is not None else None
    shown = bool(paired_bh is not None and paired_bh["n_na_pseudobulk"] == 0)
    if readout["paired_floor_shown"] and not shown:
        _fail(
            "readout.paired_floor_shown",
            "true requires real_label.paired_bh with n_na_pseudobulk == 0",
        )


def validate(payload: dict) -> None:
    """Raise :class:`AuditSchemaError` naming the dotted path of the first schema violation.

    Returns ``None`` when ``payload`` matches the pbcheck-audit/1 schema: every required key is
    present at every level and no key that the schema does not name (``null`` is accepted only
    where the schema marks a block nullable), every value has the right type, and every literal
    enum (``status``, ``status_reason``, caveat ids, ``universe.builder``,
    ``permutation_null.engine_path``, floor ``bh_mode``, ``readout`` lambda classes,
    ``schema_version`` itself) holds one of its allowed values.

    Beyond the shape, it enforces the two rules that make the JSON artifact as trustworthy as the
    rendered report: the caveat block carries the always-on notes and exactly the conditional
    notes its own fields call for, with no text matching the report's forbidden patterns
    (:func:`_check_caveats`); and the readout flags that switch conditional prose on agree with the
    fields they are derived from (:func:`_check_readout_consistency`).
    """
    if not _is_mapping(payload):
        raise AuditSchemaError(f"<root>: expected dict, got {type(payload).__name__}")
    _TOP_LEVEL_SCHEMA.check(payload, "")
    _check_caveats(payload)
    _check_readout_consistency(payload)


def _always_on_caveats(payload: Mapping) -> list[dict]:
    """The four notes every payload carries, with the placeholder payload's own values filled in.

    The texts come from :mod:`pbcheck.render.text`; nothing is written here. The import is local
    because ``pbcheck.render`` imports this module, and only this function needs the renderer.
    """
    from pbcheck.render import text  # local: pbcheck.render imports this module's validator

    names = ", ".join(payload["settings"]["protocol_constants"])
    values = {
        "N1": {"version": payload["pbcheck_version"]},
        "N2": {},
        "N3": {},
        "N7": {"protocol_constant_names": names},
    }
    return [
        {"id": caveat_id, "text": text.caveat_text(caveat_id, **values[caveat_id])}
        for caveat_id in ALWAYS_ON_CAVEAT_IDS
    ]


def empty_payload() -> dict:
    """A payload with every schema key present, every arm ``null``, ``status`` ``design_only``.

    Every key that the schema marks nullable is set to ``None``; every other key holds a
    type-correct placeholder (``""``, ``0``, ``0.0``, ``False``, ``{}`` or ``[]``), except two
    blocks that :func:`validate` constrains against each other. The caveat block carries the four
    always-on notes (:data:`ALWAYS_ON_CAVEAT_IDS`) and the donor-threshold note N4, because the
    skeleton counts no donors: ``readout.min_donors_per_group`` is ``0``, the declared threshold is
    :data:`~pbcheck.product_constants.FEW_DONORS_THRESHOLD` (the only value the validator accepts)
    and ``readout.few_donors`` is therefore true. It is meant as a starting point for building a
    real payload and as a fixture that :func:`validate` accepts.
    """
    from pbcheck.render import text  # local: pbcheck.render imports this module's validator

    input_block = {
        "path": None,
        "n_cells_loaded": 0,
        "n_genes_loaded": 0,
        "n_cells_audited": 0,
        "n_cells_dropped_other_condition": 0,
        "n_cells_dropped_other_celltype": 0,
        "n_cells_dropped_missing_donor": 0,
        "n_cells_dropped_missing_condition": 0,
        "n_cells_dropped_missing_celltype": 0,
        "donor_col": "",
        "condition_col": "",
        "test_level": "",
        "ref_level": "",
        "celltype_col": None,
        "celltype_value": None,
        "celltype_like_columns_found": [],
        "batch_cols": [],
        "counts_source": None,
    }
    design_block = {
        "condition_col": "",
        "donor_col": "",
        "n_cells": 0,
        "groups": {},
        "donors_per_group": {},
        "cells_per_donor": {},
        "donor_nests_in_condition": True,
        "min_donors_per_group": 0,
        "imbalance_ratio": 0.0,
        "batch_confounded": {},
        "batch_separates_condition": {},
        "flags": [],
        "min_donors": 3,
        "usable_for_pseudobulk": False,
    }
    settings_block = {
        "tool": {
            "n_perm_requested": 0,
            "n_perm_pb_requested": 0,
            "seed": 0,
            "top_n": 0,
            "design_only": True,
            "naive_method": "wilcoxon",
            "naive_engine": "fast",
            "pseudobulk_method": "moderated_ebayes",
            "trend": False,
            "min_donors_per_group": 0,
            "universe_min_total_count": 0,
            "universe_min_prop": 0.0,
            "fallback_universe_min_prop": 0.0,
            "fallback_universe_min_size": 0,
            "min_profiles_per_group_after_thin_filter": 0,
        },
        "protocol_constants": {
            "alpha": 0.0,
            "lambda_band": [0.0, 0.0],
            "min_universe_size": 0,
            "min_cells": 0,
            "min_counts": 0,
        },
    }
    universe_block = {
        "size": 0,
        "min_size": 0,
        "builder": None,
        "builder_rule": "",
        "thin_donor_filter": None,
        "profiles_per_group_after_thin_filter": None,
    }
    readout_block = {
        "lambda_naive_class": None,
        "lambda_pseudobulk_class": None,
        "naive_floor_solo": None,
        "naive_real_solo": None,
        "naive_real_over_floor_solo": None,
        "pseudobulk_real_over_floor": None,
        "paired_floor_shown": False,
        "few_donors": True,
        "few_donors_threshold": FEW_DONORS_THRESHOLD,
        "min_donors_per_group": 0,
        "min_profiles_per_group_after_thin_filter": None,
        "n_perm_naive_requested": 0,
        "n_perm_naive_achieved": None,
        "n_perm_pb_requested": 0,
        "n_perm_pb_achieved": None,
        "n_distinct_splits": None,
        "sentences": [],
    }
    provenance_block = {
        "pre_registered": {},
        "pre_registered_source": "",
        "operating_envelope": [],
        "platform": "",
        "python": "",
        "packages": {key: "" for key in PACKAGE_KEYS},
    }
    payload = {
        "schema_version": AUDIT_SCHEMA_VERSION,
        "pbcheck_version": "",
        "generated_utc": "",
        "runtime_seconds": 0.0,
        "runtime_by_stage_seconds": {key: None for key in RUNTIME_STAGE_KEYS},
        "status": "design_only",
        "status_reason": "design_only_requested",
        "input": input_block,
        "counts_check": None,
        "design": design_block,
        "settings": settings_block,
        "universe": universe_block,
        "real_label": None,
        "permutation_null": None,
        "readout": readout_block,
        "caveats": [],
        "provenance": provenance_block,
    }
    payload["caveats"] = [
        *_always_on_caveats(payload),
        {"id": "N4", "text": text.caveat_text("N4")},
    ]
    return payload
