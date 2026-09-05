"""Tests for :mod:`pbcheck.audit`, the audit orchestration module.

The audit is the only product module that drives the measurement engine, so these tests pin the
wiring rather than the statistics: which status a given input lands in, which arm ran, which BH
convention the headline floor carries, that the caller's object is untouched, and that the naive
half of the permutation null reproduces the engine's own numbers bit for bit rather than
approximating them. The statistics themselves are pinned by the engine's own tests
(``tests/test_permutation_mtc.py``, ``tests/test_naive_engine.py``, ``tests/test_metrics.py``).

Permutation counts are deliberately small (20 naive / 10 pseudobulk) everywhere: every assertion
here is about wiring, and the product defaults would only buy resolution no assertion reads.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from conftest import AUDIT_ORACLE_SHAPE

from pbcheck import audit, audit_schema, gate_config, gene_universe
from pbcheck.audit import AuditSettings, run_audit
from pbcheck.methods.pseudobulk import build_pseudobulk
from pbcheck.render.text import LAMBDA_CLASS_WORDS

# The oracle's own column names and levels (synthetic/oracles.py).
DONOR_COL = "donor"
CONDITION_COL = "condition"
CELLTYPE_COL = "cell_type"
CELLTYPE_VALUE = "T_cell"
TEST_LEVEL = "disease"
REF_LEVEL = "ctrl"

N_PERM = 20
N_PERM_PB = 10


def _settings(**overrides) -> AuditSettings:
    """The settings every test starts from: the oracle's columns, a small permutation budget."""
    base = dict(
        donor_col=DONOR_COL,
        condition_col=CONDITION_COL,
        test_level=TEST_LEVEL,
        ref_level=REF_LEVEL,
        celltype_col=CELLTYPE_COL,
        celltype_value=CELLTYPE_VALUE,
        n_perm=N_PERM,
        n_perm_pb=N_PERM_PB,
    )
    base.update(overrides)
    return AuditSettings(**base)


def _caveat_ids(payload: dict) -> list[str]:
    return [c["id"] for c in payload["caveats"]]


@pytest.fixture(scope="module")
def complete_run(audit_shape_oracle):
    """One ``complete`` audit of the shape oracle, plus the input as it was before the call.

    Module-scoped: the run is the expensive object in this file, and three tests read it (the
    status assertions, the read-out classes and the no-mutation check), which would otherwise pay
    for the same permutation null three times. Nothing here mutates the payload.
    """
    adata = audit_shape_oracle.adata.copy()
    before = {
        "X": np.asarray(adata.X).copy(),
        "obs": adata.obs.copy(),
        "var_names": list(adata.var_names),
        "obs_names": list(adata.obs_names),
        "layers": sorted(adata.layers.keys()),
    }
    payload = run_audit(adata, _settings())
    return {"payload": payload, "adata": adata, "before": before}


def test_audit_shape_oracle_clears_all_floors(audit_shape_oracle):
    """The fixture shape clears the three floors the audit gates on, with room to spare.

    Re-measured here rather than trusted: the comment on ``AUDIT_ORACLE_SHAPE`` states these
    numbers, and a shape that quietly stopped clearing one of them would turn every ``complete``
    test in this file into a test of a gate instead.
    """
    adata = audit_shape_oracle.adata
    donors_per_group = (
        adata.obs[[DONOR_COL, CONDITION_COL]].astype(str).drop_duplicates()[CONDITION_COL]
        .value_counts()
    )
    assert donors_per_group.min() >= audit.PRODUCT_MIN_DONORS_PER_GROUP

    pdata = build_pseudobulk(adata, donor_col=DONOR_COL, celltype_col=CELLTYPE_COL,
                             condition_col=CONDITION_COL)
    profiles_per_group = pdata.obs[CONDITION_COL].astype(str).value_counts()
    assert profiles_per_group.min() >= audit.PRODUCT_MIN_PROFILES_PER_GROUP
    assert pdata.uns["thin_donor_filter"]["n_dropped"] == 0

    universe = gene_universe.frozen_universe(
        pdata,
        min_total_count=audit.PRODUCT_UNIVERSE_MIN_TOTAL_COUNT,
        min_prop=audit.PRODUCT_UNIVERSE_MIN_PROP,
        min_size=gate_config.MIN_UNIVERSE_SIZE,
    )
    assert len(universe) == AUDIT_ORACLE_SHAPE["n_genes"]
    assert len(universe) >= gate_config.MIN_UNIVERSE_SIZE


def test_run_audit_complete_on_null_oracle(complete_run):
    """The whole path on raw counts: both arms, both nulls, a validated payload.

    ``lambda_naive_class == "inflated"`` is the instrument's own sanity condition: the oracle has a
    donor random effect and no true DE, so a per-cell test that treats cells as replicates must
    show inflation. If it did not, the audit would be reporting a number the engine's own gate
    (``scripts/synthetic_gate.py``) would refuse.
    """
    payload = complete_run["payload"]
    audit_schema.validate(payload)

    assert payload["status"] == "complete"
    assert payload["status_reason"] is None
    assert payload["input"]["counts_source"] == "X"
    assert payload["counts_check"]["passed"] is True
    assert payload["universe"]["builder"] == "pseudobulk_frozen"
    assert payload["universe"]["size"] == AUDIT_ORACLE_SHAPE["n_genes"]

    null = payload["permutation_null"]
    assert null["engine_path"] == "run_null"
    assert np.isfinite(null["naive"]["lambda"])
    assert np.isfinite(null["pseudobulk"]["lambda"])
    assert null["naive"]["floor_solo"]["bh_mode"] == "solo"
    assert null["naive"]["floor_paired"]["bh_mode"] == "paired"
    assert null["pseudobulk"]["floor"]["bh_mode"] == "paired"
    assert null["n_perm_naive_achieved"] == N_PERM
    assert null["n_perm_pb_achieved"] == N_PERM_PB
    # 4v4: C(8, 4) = 70 balanced splits, minus the true one and its complement.
    assert null["n_distinct_splits"] == 68

    readout = payload["readout"]
    assert readout["lambda_naive_class"] == "inflated"
    assert readout["lambda_pseudobulk_class"] in LAMBDA_CLASS_WORDS
    assert readout["naive_floor_solo"]["n_perm"] == N_PERM
    assert readout["naive_real_solo"] == payload["real_label"]["naive"]["n_significant_solo"]
    assert "A1" in _caveat_ids(payload)


def test_run_audit_on_no_donor_effect_oracle(audit_shape_oracle):
    """Falsification control: with ``donor_sigma = 0`` the naive arm must come back calibrated.

    Same shape, same settings, one parameter of the generator changed. A pipeline that called
    every dataset inflated would pass the test above and fail this one.
    """
    from oracles import no_donor_effect_oracle

    oracle = no_donor_effect_oracle(seed=7, **AUDIT_ORACLE_SHAPE)
    payload = run_audit(oracle.adata, _settings())

    assert payload["status"] == "complete"
    assert payload["readout"]["lambda_naive_class"] == "calibrated"


def test_design_only(small_null_adata):
    """``--design-only``: the design audit runs, nothing else does, the report still stands."""
    payload = run_audit(small_null_adata,
                        _settings(celltype_col=None, celltype_value=None, design_only=True))
    audit_schema.validate(payload)

    assert payload["status"] == "design_only"
    assert payload["status_reason"] == "design_only_requested"
    assert payload["counts_check"] is None
    assert payload["real_label"] is None
    assert payload["permutation_null"] is None
    assert payload["universe"]["builder"] is None
    assert payload["design"]["donor_nests_in_condition"] is True

    ids = _caveat_ids(payload)
    for always in ("C1", "C2", "C3", "C7"):
        assert always in ids
    assert "A1" not in ids
    # The oracle's own ``cell_type`` column matches the cell-type-like pattern, and no cell type
    # was selected, so the pooling caveat is due.
    assert "C8" in ids
    assert payload["input"]["celltype_like_columns_found"] == [CELLTYPE_COL]


def test_paired_design_is_gated(small_null_adata):
    """A donor measured under both conditions stops the DE arms before any matrix is touched."""
    obs = small_null_adata.obs
    spanning = str(obs[DONOR_COL].astype(str).iloc[0])
    cells = np.flatnonzero(obs[DONOR_COL].astype(str).to_numpy() == spanning)
    condition = obs[CONDITION_COL].astype(str).to_numpy()
    condition[cells[: len(cells) // 2]] = TEST_LEVEL
    condition[cells[len(cells) // 2:]] = REF_LEVEL
    small_null_adata.obs[CONDITION_COL] = condition

    payload = run_audit(small_null_adata, _settings(celltype_col=None, celltype_value=None))
    audit_schema.validate(payload)

    assert payload["status"] == "design_only"
    assert payload["status_reason"] == "donor_spans_conditions"
    assert payload["design"]["donor_nests_in_condition"] is False
    assert payload["counts_check"] is None
    assert payload["real_label"] is None
    assert payload["permutation_null"] is None


def test_too_few_donors_gated(small_null_adata):
    """2v2 is below the donor gate, whatever the counts look like."""
    donor_of = small_null_adata.obs[DONOR_COL].astype(str)
    condition_of = small_null_adata.obs[CONDITION_COL].astype(str)
    keep_donors = []
    for level in (TEST_LEVEL, REF_LEVEL):
        keep_donors += sorted(set(donor_of[condition_of == level]))[:2]
    subset = small_null_adata[donor_of.isin(keep_donors).to_numpy()].copy()

    payload = run_audit(subset, _settings(celltype_col=None, celltype_value=None))
    audit_schema.validate(payload)

    assert payload["status"] == "design_only"
    assert payload["status_reason"] == "too_few_donors"
    assert payload["readout"]["min_donors_per_group"] == 2
    assert payload["design"]["min_donors"] == audit.PRODUCT_MIN_DONORS_PER_GROUP
    assert payload["permutation_null"] is None


def test_non_integer_x_drops_pseudobulk_and_runs_naive(audit_shape_oracle):
    """A matrix that is not raw counts: the pseudobulk arm is dropped, never rounded.

    The naive arm still runs, on the matrix as found, and the payload says so in three places at
    once: the status reason, the universe builder and the C6 caveat. The action sentence A1, which
    compares the two arms, is not rendered because one of them does not exist here.
    """
    adata = audit_shape_oracle.adata.copy()
    adata.X = np.asarray(adata.X, dtype=float) / 2.0

    payload = run_audit(adata, _settings())
    audit_schema.validate(payload)

    assert payload["status"] == "naive_only"
    assert payload["status_reason"] == "non_integer_counts"
    assert payload["input"]["counts_source"] is None
    assert payload["counts_check"]["passed"] is False
    assert payload["universe"]["builder"] == "naive_detection_fallback"
    assert payload["universe"]["size"] >= audit.FALLBACK_UNIVERSE_MIN_SIZE
    assert payload["universe"]["builder_rule"] == audit.FALLBACK_UNIVERSE_RULE

    null = payload["permutation_null"]
    assert null["engine_path"] == "naive_null"
    assert null["pseudobulk"] is None
    assert null["n_perm_pb_achieved"] is None
    assert null["naive"]["floor_paired"] is None
    assert payload["real_label"]["pseudobulk"] is None
    assert payload["readout"]["lambda_pseudobulk_class"] is None
    assert payload["readout"]["paired_floor_shown"] is False

    ids = _caveat_ids(payload)
    assert "C6" in ids
    assert "A1" not in ids
    c6 = next(c["text"] for c in payload["caveats"] if c["id"] == "C6")
    assert payload["counts_check"]["reason"] in c6


def test_naive_null_matches_run_null_bitwise(audit_shape_oracle):
    """``audit.naive_null`` is the engine's naive half, not a second implementation of it.

    Same data, same seed, same requested counts: the per-permutation p-value matrix and the solo
    #DEG series must be identical value for value, as ``tests/test_naive_engine.py`` demands of
    the two naive engines. Anything less would mean the two statuses that use ``naive_null``
    report a different statistic from the one ``complete`` reports.
    """
    from pbcheck.permutation import run_null

    adata = audit_shape_oracle.adata.copy()
    settings = _settings()
    pdata = build_pseudobulk(adata, donor_col=DONOR_COL, celltype_col=CELLTYPE_COL,
                             condition_col=CONDITION_COL)
    universe = gene_universe.frozen_universe(
        pdata,
        min_total_count=audit.PRODUCT_UNIVERSE_MIN_TOTAL_COUNT,
        min_prop=audit.PRODUCT_UNIVERSE_MIN_PROP,
        min_size=gate_config.MIN_UNIVERSE_SIZE,
    )

    ours = audit.naive_null(adata, universe, settings)
    theirs = run_null(
        adata, universe,
        donor_col=DONOR_COL, condition_col=CONDITION_COL,
        test_level=TEST_LEVEL, ref_level=REF_LEVEL, celltype_col=CELLTYPE_COL,
        n_perm=settings.n_perm, n_perm_pb=settings.n_perm_pb, fdr=settings.alpha,
        seed=settings.seed, naive_engine="fast",
    )

    assert ours["n_perm_naive"] == theirs["n_perm_naive"]
    assert np.array_equal(ours["naive_pvals"], theirs["naive_pvals"], equal_nan=True)
    assert np.array_equal(ours["naive_ndeg_solo"].counts, theirs["naive_ndeg_solo"].counts)
    assert ours["naive_ndeg_solo"].bh_mode == theirs["naive_ndeg_solo"].bh_mode
    assert ours["n_donors"] == theirs["n_donors"]
    assert ours["real_split_inside_perm_range"] == theirs["real_split_inside_perm_range"]
    assert ours["real_split_percentile_in_perms"] == theirs["real_split_percentile_in_perms"]


def test_counts_from_layers_bitwise_equal_to_counts_in_x(audit_shape_oracle):
    """Where the counts live must not change a single number the audit reports.

    Left: a file whose ``X`` is normalised and whose raw counts sit in ``layers["counts"]``.
    Right: the same counts in ``X``. The audit must find the layer, say so in ``counts_source``
    and in the check's own ``layer`` label, and produce the same read-out and the same real-label
    p-values as the run that read ``X``.
    """
    counts = np.asarray(audit_shape_oracle.adata.X, dtype=float)

    from_x = audit_shape_oracle.adata.copy()
    from_x.X = counts.copy()

    from_layer = audit_shape_oracle.adata.copy()
    from_layer.layers["counts"] = counts.copy()
    depth = counts.sum(axis=1, keepdims=True)
    from_layer.X = 1e4 * counts / np.where(depth == 0, 1.0, depth)

    payload_x = run_audit(from_x, _settings())
    payload_layer = run_audit(from_layer, _settings())

    assert payload_x["input"]["counts_source"] == "X"
    assert payload_layer["input"]["counts_source"] == "layers:counts"
    assert payload_layer["counts_check"]["layer"] == "layers:counts"
    assert payload_layer["status"] == payload_x["status"] == "complete"

    assert payload_layer["readout"] == payload_x["readout"]
    for arm, columns in (("naive", ("pval", "padj")), ("pseudobulk", ("pval", "padj"))):
        left = payload_x["real_label"][arm]["top"]
        right = payload_layer["real_label"][arm]["top"]
        assert [g["gene"] for g in left] == [g["gene"] for g in right]
        for column in columns:
            assert np.array_equal(
                np.array([g[column] for g in left]),
                np.array([g[column] for g in right]),
                equal_nan=True,
            )
    assert (payload_layer["permutation_null"]["naive"]["lambda"]
            == payload_x["permutation_null"]["naive"]["lambda"])
    assert (payload_layer["permutation_null"]["pseudobulk"]["lambda"]
            == payload_x["permutation_null"]["pseudobulk"]["lambda"])


def test_run_audit_does_not_mutate_input(complete_run):
    """The caller's object comes back exactly as it went in, down to the ``.obs`` columns.

    The audit adds a constant stratum column, drops unused categories and subsets cells; every one
    of those happens on the copy ``prepare_stratum`` makes, and a regression here would silently
    edit a user's in-memory data.
    """
    adata = complete_run["adata"]
    before = complete_run["before"]

    assert np.array_equal(np.asarray(adata.X), before["X"])
    assert adata.obs.equals(before["obs"])
    assert list(adata.obs.columns) == list(before["obs"].columns)
    assert audit.STRATUM_COL not in adata.obs.columns
    assert list(adata.var_names) == before["var_names"]
    assert list(adata.obs_names) == before["obs_names"]
    assert sorted(adata.layers.keys()) == before["layers"]


def test_no_protocol_constant_from_gate_config_for_perms():
    """The permutation budget is the product's own decision, not a pre-registered constant.

    ``audit.PRODUCT_N_PERM is gate_config.N_PERM`` says nothing about small ints, so the check is
    on the source: the module must not name ``gate_config.N_PERM`` at all. The protocol constants
    it is allowed to read are the ones the payload lists under ``protocol_constants``.
    """
    source = Path(audit.__file__).read_text(encoding="utf-8")
    assert "gate_config.N_PERM" not in source
    assert isinstance(audit.PRODUCT_N_PERM, int)
    assert isinstance(audit.PRODUCT_N_PERM_PB, int)
    assert AuditSettings(donor_col=DONOR_COL, condition_col=CONDITION_COL,
                         test_level=TEST_LEVEL, ref_level=REF_LEVEL).n_perm == audit.PRODUCT_N_PERM
