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

import inspect
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import AUDIT_ORACLE_SHAPE

from pbcheck import audit, audit_schema, gate_config, gene_universe, io_counts
from pbcheck.audit import AuditSettings, run_audit
from pbcheck.methods.pseudobulk import build_pseudobulk
from pbcheck.render import sections, text
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


#: What may legitimately sit outside every stage: the payload's own bookkeeping (the provenance
#: block reads installed package versions, the caveat texts are rendered, the schema is walked).
#: It is a fixed cost that does not grow with the data, hence an absolute allowance; the fraction
#: on top only absorbs scheduling noise. The defect these bounds exist to catch charged about a
#: quarter of a product-sized run to no stage at all.
STAGE_TABLE_OVERHEAD_SECONDS = 0.2
STAGE_TABLE_OVERHEAD_FRACTION = 0.05


def _assert_stage_table_accounts_for_the_run(payload: dict) -> None:
    """The per-stage table must account for the run's wall clock bar the payload's bookkeeping.

    The table is what the report publishes as where the time went, so a piece of work performed
    outside every ``audit._stage`` block would understate it without changing any shape.
    """
    runtime = payload["runtime_seconds"]
    accounted = sum(v for v in payload["runtime_by_stage_seconds"].values() if v is not None)
    assert accounted <= runtime
    assert runtime - accounted <= (STAGE_TABLE_OVERHEAD_SECONDS
                                   + STAGE_TABLE_OVERHEAD_FRACTION * runtime)


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
        "var": adata.var.copy(),
        "layers": {name: np.asarray(matrix).copy()
                   for name, matrix in adata.layers.items()},
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

    ``lambda_naive_class == "above_band"`` is the instrument's own sanity condition: the oracle has
    a donor random effect and no true DE, so a per-cell test that treats cells as replicates must
    land above the donor-pseudobulk arm's band. If it did not, the audit would be reporting a
    number the engine's own gate (``scripts/synthetic_gate.py``) would refuse.
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
    assert readout["lambda_naive_class"] == "above_band"
    assert readout["lambda_pseudobulk_class"] in LAMBDA_CLASS_WORDS
    assert readout["naive_floor_solo"]["n_perm"] == N_PERM
    assert readout["naive_real_solo"] == payload["real_label"]["naive"]["n_significant_solo"]
    assert "R1" in _caveat_ids(payload)


def test_run_audit_on_no_donor_effect_oracle(audit_shape_oracle):
    """Falsification control: with ``donor_sigma = 0`` the naive arm must land inside the band.

    Same shape, same settings, one parameter of the generator changed. A pipeline that put every
    dataset above the band would pass the test above and fail this one.
    """
    from oracles import no_donor_effect_oracle

    oracle = no_donor_effect_oracle(seed=7, **AUDIT_ORACLE_SHAPE)
    payload = run_audit(oracle.adata, _settings())

    assert payload["status"] == "complete"
    assert payload["readout"]["lambda_naive_class"] == "in_band"


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
    assert payload["readout"]["sentences"] == []
    for always in ("N1", "N2", "N3", "N7"):
        assert always in ids
    assert "R1" not in ids
    # The oracle's own ``cell_type`` column matches the cell-type-like pattern, and no cell type
    # was selected, so the pooling caveat is due.
    assert "N8" in ids
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
    once: the status reason, the universe builder and the N6 caveat. The action sentence R1, which
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
    assert "N6" in ids
    assert "R1" not in ids
    n6 = next(c["text"] for c in payload["caveats"] if c["id"] == "N6")
    # The reason is the user's own file talking, so the note carries it as a quoted identifier,
    # shortened by ``text.quoted`` when it is longer than one masked span may be.
    assert text.quoted(payload["counts_check"]["reason"]) in n6
    assert payload["counts_check"]["reason"].startswith("non-integer values in X")

    # The naive-only path summarises its own null after the engine call; that work is charged to
    # the permutation stage rather than falling outside the table.
    _assert_stage_table_accounts_for_the_run(payload)


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
    assert adata.var.equals(before["var"])
    assert list(adata.var.columns) == list(before["var"].columns)
    assert sorted(adata.layers.keys()) == sorted(before["layers"])
    for name, matrix in before["layers"].items():
        assert np.array_equal(np.asarray(adata.layers[name]), matrix, equal_nan=True)


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


# ---------------------------------------------------------------------------
# Part 2: the rest of the work package's test list.
# ---------------------------------------------------------------------------

#: Permutation counts of the ``slow`` end-to-end run below. Written out here rather than read from
#: ``gate_config``: the audit's permutation budget is a product decision, and a test that borrowed
#: the gate script's number would tie the two together where the module deliberately does not.
SLOW_N_PERM = 200
SLOW_N_PERM_PB = 200


def _pseudobulk_universe(adata, celltype_col: str = CELLTYPE_COL) -> list[str]:
    """The frozen universe the audit would build for ``adata``, from the audit's own parameters."""
    pdata = build_pseudobulk(adata, donor_col=DONOR_COL, celltype_col=celltype_col,
                             condition_col=CONDITION_COL)
    return gene_universe.frozen_universe(
        pdata,
        min_total_count=audit.PRODUCT_UNIVERSE_MIN_TOTAL_COUNT,
        min_prop=audit.PRODUCT_UNIVERSE_MIN_PROP,
        min_size=gate_config.MIN_UNIVERSE_SIZE,
    )


def _top_column(payload: dict, arm: str, column: str) -> np.ndarray:
    """One column of an arm's top-gene table, as an array, for bitwise comparisons."""
    return np.array([gene[column] for gene in payload["real_label"][arm]["top"]])


@pytest.fixture(scope="module")
def engine_null(audit_shape_oracle):
    """One ``permutation.run_null`` on the shape oracle at this file's permutation budget.

    Module-scoped for the same reason as ``complete_run``: it is the engine-side counterpart the
    Monte-Carlo SE test compares against, and nothing mutates it.
    """
    from pbcheck.permutation import run_null

    adata = audit_shape_oracle.adata.copy()
    settings = _settings()
    universe = _pseudobulk_universe(adata)
    result = run_null(
        adata, universe,
        donor_col=DONOR_COL, condition_col=CONDITION_COL,
        test_level=TEST_LEVEL, ref_level=REF_LEVEL, celltype_col=CELLTYPE_COL,
        n_perm=settings.n_perm, n_perm_pb=settings.n_perm_pb, fdr=settings.alpha,
        seed=settings.seed, naive_engine="fast",
    )
    return {"result": result, "universe": universe}


def test_positive_oracle_real_label_counts():
    """With true DE in the data the pseudobulk arm must find some of it.

    The counterpart of the null-oracle tests above: those pin that the audit reports inflation
    where the truth is "nothing", this one pins that a run whose truth is "something" does not come
    back empty. Without it every assertion in this file would still pass on a pipeline whose
    pseudobulk arm returned zero significant genes for any input.
    """
    from oracles import positive_oracle

    oracle = positive_oracle(n_de=100, seed=3, log2fc=2.0, **AUDIT_ORACLE_SHAPE)
    payload = run_audit(oracle.adata, _settings())
    audit_schema.validate(payload)

    assert payload["status"] == "complete"
    assert payload["real_label"]["pseudobulk"]["n_significant_paired"] > 0
    assert payload["readout"]["pseudobulk_real_over_floor"] > 1.0


def test_counts_layer_flag_and_raw_fallback(audit_shape_oracle):
    """The three ways the caller can point the audit at a matrix that is not ``X``.

    ``--counts-layer NAME`` names a layer outright; ``.raw`` is the last fallback and may carry a
    different ``var`` than the file itself, which the working object must follow rather than the
    file's; a layer name the file does not have is a caller error, not a silent fallback to ``X``.
    """
    counts = np.asarray(audit_shape_oracle.adata.X, dtype=float)
    depth = counts.sum(axis=1, keepdims=True)
    normalised = 1e4 * counts / np.where(depth == 0, 1.0, depth)

    named_layer = audit_shape_oracle.adata.copy()
    named_layer.layers["my_counts"] = counts.copy()
    named_layer.X = normalised.copy()
    payload = run_audit(named_layer, _settings(counts_layer="my_counts"))
    audit_schema.validate(payload)
    assert payload["status"] == "complete"
    assert payload["input"]["counts_source"] == "layers:my_counts"
    assert payload["counts_check"]["layer"] == "layers:my_counts"

    # ``.raw`` with a wider ``var`` than the file: two all-zero genes that exist only in ``.raw``.
    raw_counts = np.hstack([counts, np.zeros((counts.shape[0], 2))])
    raw_var = pd.DataFrame(
        index=list(audit_shape_oracle.adata.var_names) + ["raw_only_0", "raw_only_1"]
    )
    raw_var["only_in_raw"] = True
    from_raw = audit_shape_oracle.adata.copy()
    from_raw.X = normalised.copy()
    from_raw.raw = ad.AnnData(X=raw_counts, obs=from_raw.obs.copy(), var=raw_var)
    payload_raw = run_audit(from_raw, _settings())
    audit_schema.validate(payload_raw)
    assert payload_raw["status"] == "complete"
    assert payload_raw["input"]["counts_source"] == "raw.X"
    assert payload_raw["counts_check"]["layer"] == "raw.X"
    # The two genes that exist only in ``.raw`` are all-zero, so they cannot enter the universe;
    # the audit still read the wider matrix, not the file's ``X``.
    assert payload_raw["universe"]["size"] == AUDIT_ORACLE_SHAPE["n_genes"]

    with pytest.raises(audit.AuditInputError, match="no layer 'not_here'"):
        run_audit(audit_shape_oracle.adata.copy(), _settings(counts_layer="not_here"))


def test_too_few_profiles_after_thin_filter(audit_shape_oracle):
    """Two ctrl donors thinned below ``gate_config.MIN_CELLS`` leave that group 2 profiles.

    The gate is on the profiles that survive aggregation, not on the donors the design audit
    counted: three donors that contribute a profile each are a different thing from three donors
    two of which are dropped for having 5 cells. The naive arm still runs, on every cell.
    """
    adata = audit_shape_oracle.adata
    donor = adata.obs[DONOR_COL].astype(str).to_numpy()
    condition = adata.obs[CONDITION_COL].astype(str).to_numpy()
    thinned = sorted(set(donor[condition == REF_LEVEL]))[:2]

    keep = np.ones(adata.n_obs, dtype=bool)
    for name in thinned:
        cells = np.flatnonzero(donor == name)
        keep[cells[5:]] = False
    work = adata[keep].copy()

    payload = run_audit(work, _settings())
    audit_schema.validate(payload)

    assert payload["status"] == "naive_only"
    assert payload["status_reason"] == "too_few_profiles_after_thin_filter"
    assert payload["universe"]["builder"] == "pseudobulk_frozen"
    profiles = payload["universe"]["profiles_per_group_after_thin_filter"]
    assert profiles[REF_LEVEL] == 2
    assert profiles[TEST_LEVEL] == AUDIT_ORACLE_SHAPE["n_donors_per_group"]
    assert payload["readout"]["min_profiles_per_group_after_thin_filter"] == 2
    assert payload["universe"]["thin_donor_filter"]["n_dropped"] == 2

    assert payload["permutation_null"]["engine_path"] == "naive_null"
    assert payload["permutation_null"]["pseudobulk"] is None
    assert payload["real_label"]["pseudobulk"] is None
    assert payload["readout"]["paired_floor_shown"] is False
    # The naive arm saw the thinned donors too: the thin-donor filter is an aggregation rule.
    assert payload["design"]["donors_per_group"][REF_LEVEL] == \
        AUDIT_ORACLE_SHAPE["n_donors_per_group"]
    assert "R1" not in _caveat_ids(payload)


def test_unused_condition_level_does_not_distort_design(audit_shape_oracle):
    """A third condition level is dropped from the stratum and counted, and changes nothing else.

    One donor per group is relabelled to a level the audit was not asked about, so the audited
    stratum is a balanced 3v3 and the design numbers must read as such: an implementation that
    audited the design before subsetting would report four donors per group and a cell-count
    imbalance that the audited data does not have.

    ``design_only`` because every assertion here is on the ``input`` and ``design`` blocks, both
    built before any matrix is touched; the DE arms are covered by the tests above.
    """
    adata = audit_shape_oracle.adata.copy()
    donor = adata.obs[DONOR_COL].astype(str).to_numpy()
    condition = adata.obs[CONDITION_COL].astype(str).to_numpy()
    relabelled = [sorted(set(donor[condition == level]))[0] for level in (TEST_LEVEL, REF_LEVEL)]
    condition[np.isin(donor, relabelled)] = "other_condition"
    adata.obs[CONDITION_COL] = pd.Categorical(condition)

    payload = run_audit(adata, _settings(design_only=True))
    audit_schema.validate(payload)

    n_dropped = 2 * AUDIT_ORACLE_SHAPE["n_cells_per_donor"]
    assert payload["input"]["n_cells_dropped_other_condition"] == n_dropped
    assert payload["input"]["n_cells_audited"] == adata.n_obs - n_dropped

    design = payload["design"]
    assert set(design["groups"]) == {TEST_LEVEL, REF_LEVEL}
    assert design["donors_per_group"] == {
        TEST_LEVEL: AUDIT_ORACLE_SHAPE["n_donors_per_group"] - 1,
        REF_LEVEL: AUDIT_ORACLE_SHAPE["n_donors_per_group"] - 1,
    }
    assert design["imbalance_ratio"] == 1.0
    assert not [flag for flag in design["flags"] if "imbalance" in flag]


def test_missing_obs_values_dropped_and_counted(small_null_adata):
    """Cells with no donor or no condition are dropped, each for its own recorded reason."""
    obs = small_null_adata.obs
    donor = np.array(obs[DONOR_COL].astype(object).to_numpy(), dtype=object)
    condition = np.array(obs[CONDITION_COL].astype(object).to_numpy(), dtype=object)
    donor[:10] = None
    condition[10:15] = None
    small_null_adata.obs[DONOR_COL] = donor
    small_null_adata.obs[CONDITION_COL] = condition

    payload = run_audit(small_null_adata,
                        _settings(celltype_col=None, celltype_value=None, design_only=True))
    audit_schema.validate(payload)

    block = payload["input"]
    assert block["n_cells_dropped_missing_donor"] == 10
    assert block["n_cells_dropped_missing_condition"] == 5
    assert block["n_cells_dropped_missing_celltype"] == 0
    assert block["n_cells_audited"] == block["n_cells_loaded"] - 15
    assert payload["design"]["n_cells"] == block["n_cells_audited"]


def test_too_many_missing_raises(small_null_adata):
    """Past :data:`audit.MAX_MISSING_OBS_FRACTION` the file is refused instead of half-audited."""
    condition = np.array(small_null_adata.obs[CONDITION_COL].astype(object).to_numpy(),
                         dtype=object)
    n_missing = int(np.ceil(0.6 * condition.size))
    condition[:n_missing] = None
    small_null_adata.obs[CONDITION_COL] = condition

    with pytest.raises(audit.AuditInputError, match="missing a value"):
        run_audit(small_null_adata, _settings(celltype_col=None, celltype_value=None))


def test_duplicate_var_names_raise_audit_input_error(small_null_adata):
    """Duplicated gene names are refused: every downstream index is by gene name."""
    names = list(small_null_adata.var_names)
    names[1] = names[0]
    small_null_adata.var_names = names

    with pytest.raises(audit.AuditInputError, match="var_names_make_unique"):
        run_audit(small_null_adata, _settings(celltype_col=None, celltype_value=None,
                                              design_only=True))


def test_celltype_subset_and_constant_column(audit_shape_oracle):
    """Pooling and an explicit single-cell-type selection are the same run, and say so differently.

    The oracle carries one cell type, so selecting it explicitly and pooling everything must give
    the same numbers gene for gene; what differs is the bookkeeping. Pooling adds the constant
    stratum column so the engine's ``celltype_col`` is never ``None``, records the cell-type-like
    column it found and raises the pooling caveat; the explicit run does neither.
    """
    settings_pooled = _settings(celltype_col=None, celltype_value=None)
    work, _ = audit.prepare_stratum(audit_shape_oracle.adata.copy(), settings_pooled)
    assert list(work.obs[audit.STRATUM_COL].astype(str).unique()) == [audit.STRATUM_VALUE]

    pooled = run_audit(audit_shape_oracle.adata.copy(), settings_pooled)
    explicit = run_audit(audit_shape_oracle.adata.copy(), _settings())
    audit_schema.validate(pooled)

    assert pooled["status"] == explicit["status"] == "complete"
    assert pooled["readout"] == explicit["readout"]
    for arm in ("naive", "pseudobulk"):
        assert ([gene["gene"] for gene in pooled["real_label"][arm]["top"]]
                == [gene["gene"] for gene in explicit["real_label"][arm]["top"]])
        for column in ("pval", "padj"):
            assert np.array_equal(_top_column(pooled, arm, column),
                                  _top_column(explicit, arm, column), equal_nan=True)

    assert pooled["input"]["celltype_col"] is None
    assert pooled["input"]["celltype_like_columns_found"] == [CELLTYPE_COL]
    assert "N8" in _caveat_ids(pooled)
    assert explicit["input"]["celltype_value"] == CELLTYPE_VALUE
    assert "N8" not in _caveat_ids(explicit)

    # The mirror of the "column without a value" error: a value with no column to read it from
    # would otherwise be copied into the payload as a subset that was never taken.
    with pytest.raises(audit.AuditInputError, match="without a cell-type column"):
        run_audit(audit_shape_oracle.adata.copy(),
                  _settings(celltype_col=None, celltype_value=CELLTYPE_VALUE))
    with pytest.raises(audit.AuditInputError, match="without a cell-type value"):
        run_audit(audit_shape_oracle.adata.copy(),
                  _settings(celltype_col=CELLTYPE_COL, celltype_value=None))


def test_universe_too_small_becomes_design_only(small_null_adata):
    """80 genes cannot make a universe of 200, and no arm is run on the ones that are left."""
    payload = run_audit(small_null_adata, _settings(celltype_col=None, celltype_value=None))
    audit_schema.validate(payload)

    assert payload["status"] == "design_only"
    assert payload["status_reason"] == "universe_too_small"
    assert payload["universe"]["builder"] == "pseudobulk_frozen"
    assert 0 < payload["universe"]["size"] < gate_config.MIN_UNIVERSE_SIZE
    assert payload["universe"]["min_size"] == gate_config.MIN_UNIVERSE_SIZE
    assert payload["real_label"] is None
    assert payload["permutation_null"] is None
    assert "R1" not in _caveat_ids(payload)


def test_paired_bh_consistent_with_run_null(complete_run):
    """The real-label paired correction and the permutation path corrected the same gene set.

    ``run_null`` runs its own real-label pass; if its common tested set differed from the one the
    payload reports, the floors would be a null for a different set of genes than the real-label
    counts they are compared against, and nothing in the payload's shape would show it.
    """
    payload = complete_run["payload"]
    assert payload["real_label"]["consistent_with_permutation_path"] is True
    assert payload["real_label"]["paired_bh"]["n_universe"] == payload["universe"]["size"]
    assert (payload["real_label"]["paired_bh"]["n_tested_common"]
            <= payload["real_label"]["paired_bh"]["n_universe"])


def test_paired_floor_hidden_when_pseudobulk_nans(audit_shape_oracle, complete_run):
    """The paired floor is shown only while the pseudobulk arm tested every gene of the universe.

    The construction the plan names for producing a NaN pseudobulk p-value (a gene that is zero in
    every pseudobulk profile of one group) does not produce one, and this test measures that rather
    than assuming it: ``moderated.log_cpm`` is log2(CPM + 1), so a profile matrix of finite counts
    gives finite ``Y``, a finite residual variance, and ``moderated_t`` returns NaN only for a
    non-finite ``s2``. The read-out rule is therefore pinned directly, on a real payload whose
    paired bookkeeping is set to the state the rule is about.
    """
    zeroed = audit_shape_oracle.adata.copy()
    matrix = np.asarray(zeroed.X, dtype=float)
    is_ref = zeroed.obs[CONDITION_COL].astype(str).to_numpy() == REF_LEVEL
    universe = _pseudobulk_universe(zeroed)
    silenced = universe[0]
    matrix[is_ref, list(zeroed.var_names).index(silenced)] = 0.0
    zeroed.X = matrix

    payload = run_audit(zeroed, _settings())
    assert silenced in _pseudobulk_universe(zeroed)
    assert payload["real_label"]["paired_bh"]["n_na_pseudobulk"] == 0
    assert payload["readout"]["paired_floor_shown"] is True

    assert payload["permutation_null"]["naive"]["floor_solo"]["bh_mode"] == "solo"
    assert payload["permutation_null"]["naive"]["floor_paired"]["bh_mode"] == "paired"
    assert payload["permutation_null"]["pseudobulk"]["floor"]["bh_mode"] == "paired"
    assert complete_run["payload"]["readout"]["paired_floor_shown"] is True
    # The other side of the rule (a real pseudobulk NaN switching the paired floor off) is
    # exercised end to end in
    # ``test_paired_floor_hidden_when_the_pseudobulk_arm_returns_a_nan``.


def test_solo_floor_mc_se_matches_engine_formula(engine_null, complete_run):
    """The audit's Monte-Carlo SE is the engine's estimator, applied to one more series.

    ``run_null`` computes the SE for the paired series only, and the headline floor of this report
    is the solo one, so the formula lives in ``audit._floor_mc_se``. Pinned by running it on the
    engine's own paired series and demanding the engine's own number, bit for bit.
    """
    result = engine_null["result"]
    n_genes = len(engine_null["universe"])
    monte_carlo = result["monte_carlo"]

    assert (audit._floor_mc_se(result["naive_ndeg_paired"].counts)
            == monte_carlo["naive_floor_mc_se"])
    assert audit._floor_mc_se(result["pb_ndeg"].counts) == monte_carlo["pb_floor_mc_se"]

    solo = audit._floor_block(result["naive_ndeg_solo"], n_genes)
    assert solo["bh_mode"] == "solo"
    assert solo["mc_se"] == audit._floor_mc_se(result["naive_ndeg_solo"].counts)
    assert np.isfinite(solo["mc_se"])

    payload = complete_run["payload"]
    assert (payload["permutation_null"]["naive"]["floor_paired"]["mc_se"]
            == payload["permutation_null"]["monte_carlo"]["naive_floor_mc_se"])
    assert (payload["readout"]["naive_floor_solo"]["mc_se"]
            == payload["permutation_null"]["naive"]["floor_solo"]["mc_se"])


def test_readout_classes_match_gate_bands():
    """The three class names are the band's own three cases, and they name a position only.

    The names describe where a lambda falls against ``gate_config.LAMBDA_BAND``, the
    donor-pseudobulk arm's band, and award no property to the arm they describe; the renderer's
    phrase table is keyed by exactly these three names.
    """
    low, high = gate_config.LAMBDA_BAND
    assert (low, high) == (0.9, 1.1)

    assert audit._lambda_class(0.85) == "below_band"
    assert audit._lambda_class(1.0) == "in_band"
    assert audit._lambda_class(1.15) == "above_band"
    assert audit._lambda_class(low) == "in_band"
    assert audit._lambda_class(high) == "in_band"
    assert audit._lambda_class(None) is None
    assert audit._lambda_class(float("nan")) is None
    assert set(LAMBDA_CLASS_WORDS) == {"in_band", "above_band", "below_band"}


def test_achieved_perm_counts_and_coarse_caveat(audit_shape_oracle):
    """3v3 has 18 distinct donor splits, whatever the caller requested, and the report says so.

    The payload's achieved counts come from the engine's return, never from the request: at 3v3
    ``build_perms`` enumerates the 18 balanced splits that are neither the true one nor its
    complement, and a report that printed the 1000 that were asked for would overstate the
    resolution of every floor on the page by a factor of 55.
    """
    adata = audit_shape_oracle.adata
    donor = adata.obs[DONOR_COL].astype(str)
    condition = adata.obs[CONDITION_COL].astype(str)
    keep_donors = []
    for level in (TEST_LEVEL, REF_LEVEL):
        keep_donors += sorted(set(donor[condition == level]))[:3]
    work = adata[donor.isin(keep_donors).to_numpy()].copy()

    payload = run_audit(work, _settings(n_perm=1000, n_perm_pb=N_PERM_PB))
    audit_schema.validate(payload)

    assert payload["status"] == "complete"
    readout = payload["readout"]
    assert readout["n_distinct_splits"] == 18
    assert readout["n_perm_naive_requested"] == 1000
    assert readout["n_perm_naive_achieved"] == 18
    assert readout["n_perm_naive_achieved"] < audit.COARSE_NULL_THRESHOLD
    assert readout["n_perm_pb_achieved"] <= N_PERM_PB
    assert readout["naive_floor_solo"]["n_perm"] == 18

    caveat = next(c["text"] for c in payload["caveats"] if c["id"] == "N9")
    assert "18" in caveat and "1000" in caveat


def test_same_test_and_ref_level_is_refused_as_input(audit_shape_oracle):
    """One level named as both sides of the contrast is refused before anything is computed.

    Every per-level presence check passes for such settings (the level does exist), and the donor
    gate reads the same group twice, so it sees a healthy donor count for a stratum that has one
    group. Without this check the settings reached the DE arms and failed inside scanpy after the
    whole run had been paid for.
    """
    for extra in ({}, {"design_only": True}):
        with pytest.raises(audit.AuditInputError, match="both the test and the reference level"):
            run_audit(audit_shape_oracle.adata.copy(),
                      _settings(test_level=REF_LEVEL, ref_level=REF_LEVEL, **extra))

    settings = _settings(test_level=TEST_LEVEL, ref_level=TEST_LEVEL)
    with pytest.raises(audit.AuditInputError) as caught:
        audit.prepare_stratum(audit_shape_oracle.adata.copy(), settings)
    message = str(caught.value)
    assert CONDITION_COL in message
    assert TEST_LEVEL in message and REF_LEVEL in message  # the available values are listed


def test_universe_filter_parameters_are_the_engine_defaults():
    """The product's universe filter is the engine's own rule, bound to it rather than re-typed.

    ``audit`` passes both parameters explicitly so the payload cannot report a filter the call did
    not use, which makes a silent drift from the engine's default possible; this pins the two
    equal at both ends, against the engine constant and against the signature default the audit
    would inherit if it stopped passing them.
    """
    assert audit.PRODUCT_UNIVERSE_MIN_TOTAL_COUNT == io_counts.UNIVERSE_MIN_TOTAL_COUNT
    assert audit.PRODUCT_UNIVERSE_MIN_PROP == io_counts.UNIVERSE_MIN_PROP

    defaults = inspect.signature(gene_universe.frozen_universe).parameters
    assert defaults["min_total_count"].default == audit.PRODUCT_UNIVERSE_MIN_TOTAL_COUNT
    assert defaults["min_prop"].default == audit.PRODUCT_UNIVERSE_MIN_PROP


def test_runtime_by_stage_accounts_for_the_whole_run(complete_run):
    """The published per-stage table is not allowed to lose a piece of the run.

    The summarisation of the permutation null (the empirical permutation p-values behind every
    lambda, computed over the whole permutation matrix) runs after the engine call returns, and
    charging it to no stage understated the table by about a quarter of the run at product
    permutation counts.
    """
    payload = complete_run["payload"]
    stages = payload["runtime_by_stage_seconds"]
    assert stages["permutation_null"] is not None
    assert stages["real_label"] is not None
    _assert_stage_table_accounts_for_the_run(payload)


@pytest.mark.slow
def test_run_audit_complete_at_product_counts_on_gate_shape():
    """The whole path at the gate's own oracle shape and 200/200 permutations.

    The fast tests above run 20 permutations on 600 genes, which is enough to pin wiring and
    nothing else. This one runs the shape the instrument was measured on (8v8 donors, 1500
    genes) at a permutation budget whose floor is readable, and demands the instrument's own
    sanity condition on that floor: on a synthetic null with donor structure the per-cell arm's
    permutation floor must be at least ``gate_config.INSTRUMENT_NAIVE_FLOOR_FRAC_MIN`` of the
    universe. A wiring that reported the paired floor, the wrong universe or a floor from a
    different series would pass every fast test in this file and fail here.
    """
    from oracles import null_oracle

    oracle = null_oracle(seed=1, **dict(gate_config.ORACLE_SIM))
    payload = run_audit(oracle.adata, _settings(n_perm=SLOW_N_PERM, n_perm_pb=SLOW_N_PERM_PB))
    audit_schema.validate(payload)

    assert payload["status"] == "complete"
    assert payload["universe"]["size"] == gate_config.ORACLE_SIM["n_genes"]
    assert payload["permutation_null"]["n_perm_naive_achieved"] == SLOW_N_PERM
    assert payload["permutation_null"]["n_perm_pb_achieved"] == SLOW_N_PERM_PB

    readout = payload["readout"]
    assert readout["lambda_naive_class"] == "above_band"
    assert (payload["permutation_null"]["naive"]["lambda"]
            >= gate_config.INSTRUMENT_LAMBDA_NAIVE_MIN)
    assert readout["naive_floor_solo"]["median_frac"] >= \
        gate_config.INSTRUMENT_NAIVE_FLOOR_FRAC_MIN
    assert readout["n_perm_naive_achieved"] == SLOW_N_PERM
    assert "N4" not in _caveat_ids(payload)
    assert "N9" not in _caveat_ids(payload)


def test_readout_sentences_are_written_by_the_audit(complete_run):
    """The read-out lines are payload fields: filled here, printed verbatim by the renderer.

    Each line is checked against the payload value it interpolates, so a renderer that prints
    ``readout.sentences`` cannot show a number the payload does not carry. The naive lambda line
    is also where the band is attributed: the band is the donor-pseudobulk arm's, shown to
    describe the naive number and not as the naive arm's own criterion.
    """
    payload = complete_run["payload"]
    readout = payload["readout"]
    naive_null = payload["permutation_null"]["naive"]
    lines = readout["sentences"]

    assert len(lines) == 5
    assert all("{" not in line and "}" not in line for line in lines)

    floor_line, lambda_line, pseudobulk_line, donor_line, profile_line = lines
    assert audit.format_scalar(naive_null["floor_solo"]["median_count"]) in floor_line
    assert audit.format_scalar(readout["naive_real_solo"]) in floor_line
    assert "solo BH" in floor_line
    # The coarse clause and note N9 answer the same condition: the achieved permutation count.
    coarse = readout["n_perm_naive_achieved"] < audit.COARSE_NULL_THRESHOLD
    assert (audit.COARSE_NULL_NOTE in floor_line) is coarse
    assert ("N9" in _caveat_ids(payload)) is coarse

    assert audit.format_scalar(naive_null["lambda"]) in lambda_line
    assert LAMBDA_CLASS_WORDS[readout["lambda_naive_class"]] in lambda_line
    assert str(gate_config.LAMBDA_BAND[0]) in lambda_line
    assert "the donor-pseudobulk arm's band" in lambda_line

    pb_null = payload["permutation_null"]["pseudobulk"]
    assert audit.format_scalar(pb_null["lambda"]) in pseudobulk_line
    assert audit.format_scalar(pb_null["fp_rate"]) in pseudobulk_line
    assert "paired BH" in pseudobulk_line

    assert text.quoted(TEST_LEVEL) in donor_line
    assert text.quoted(REF_LEVEL) in donor_line
    assert str(payload["permutation_null"]["n_distinct_splits"]) in donor_line

    profiles = payload["universe"]["profiles_per_group_after_thin_filter"]
    assert str(profiles[TEST_LEVEL]) in profile_line
    assert str(gate_config.MIN_CELLS) in profile_line

    for line in lines:
        assert text.forbidden_pattern_hits(line) == ()


def test_readout_sentence_says_the_pseudobulk_arm_did_not_run(audit_shape_oracle):
    """With no pseudobulk arm the read-out says so in words, and says why."""
    normalised = audit_shape_oracle.adata.copy()
    normalised.X = np.asarray(normalised.X, dtype=float) / 2.0

    payload = run_audit(normalised, _settings())
    lines = payload["readout"]["sentences"]

    assert payload["status"] == "naive_only"
    assert payload["permutation_null"]["pseudobulk"] is None
    assert lines[2] == text.sentence_text(
        "pseudobulk_not_run",
        reason_text=sections.STATUS_REASON_WORDS[payload["status_reason"]],
    )


def test_format_scalar_is_the_reports_own_number_format():
    """The sentences and the report's tables format a number the same way; two formatters would
    show one payload value differently in two places of the same report."""
    for value in (None, True, False, 0.5, 1.0 / 3.0, 1234567.0, 12, "0.9"):
        assert audit.format_scalar(value) == sections._fmt(value)
    assert audit.format_scalar(1.0 / 3.0) == "0.3333"


def test_few_donors_flag_note_and_readout_clause(complete_run):
    """Below the donor threshold: the flag, the note and the read-out's replaced clause.

    The oracle shape is under ``product_constants.FEW_DONORS_THRESHOLD`` donors per group, so the
    read-out paragraph must not carry the categorical clause: change 2 of the fifth amendment
    (docs/AMENDMENTS.md) admits the floor-based ratio outside the envelope only when every group
    has at least that many donors.
    """
    payload = complete_run["payload"]
    readout = payload["readout"]

    assert readout["min_donors_per_group"] == AUDIT_ORACLE_SHAPE["n_donors_per_group"]
    assert readout["few_donors_threshold"] == audit.FEW_DONORS_THRESHOLD
    assert readout["few_donors"] is True
    assert "N4" in _caveat_ids(payload)

    r1 = next(c["text"] for c in payload["caveats"] if c["id"] == "R1")
    assert text.R1_CLAUSES["separable"] not in r1
    assert text.R1_CLAUSES["leak_contaminated"].format(
        threshold=audit.FEW_DONORS_THRESHOLD) in r1

    # The pseudobulk clause is rendered exactly when the paired floor is shown, with the numbers
    # of the arm it names and the BH convention they were produced under.
    pb_clause = text.R1_CLAUSES["pseudobulk"].format(
        pb_real=payload["real_label"]["pseudobulk"]["n_significant_paired"],
        pb_floor=f"{payload['permutation_null']['pseudobulk']['floor']['median_count']:.0f}",
    )
    assert (pb_clause in r1) is readout["paired_floor_shown"]


def test_paired_floor_hidden_when_the_pseudobulk_arm_returns_a_nan(audit_shape_oracle,
                                                                   monkeypatch):
    """A real NaN in the pseudobulk arm switches the paired floor and its clause off.

    The construction the plan named does not produce a NaN on this engine (see
    ``test_paired_floor_hidden_when_pseudobulk_nans``); the only documented source is a non-finite
    residual variance in ``moderated_t``. That state is injected at the engine boundary the audit
    calls, so the paired bookkeeping, the read-out flag and the read-out paragraph are reached the
    way a real NaN would reach them, rather than by editing the finished payload.
    """
    real_pseudobulk_de = audit.pseudobulk_de

    def one_nan_gene(*args, **kwargs):
        result = real_pseudobulk_de(*args, **kwargs)
        result.table.loc[result.table.index[0], ["pval", "padj"]] = np.nan
        return result

    monkeypatch.setattr(audit, "pseudobulk_de", one_nan_gene)
    payload = run_audit(audit_shape_oracle.adata.copy(), _settings())
    audit_schema.validate(payload)

    assert payload["status"] == "complete"
    assert payload["real_label"]["paired_bh"]["n_na_pseudobulk"] == 1
    assert payload["real_label"]["paired_bh"]["pseudobulk_na_free"] is False
    assert payload["readout"]["paired_floor_shown"] is False

    r1 = next(c["text"] for c in payload["caveats"] if c["id"] == "R1")
    assert "paired BH)." not in r1
    assert text.R1_CLAUSES["pseudobulk"].split("{")[0] not in r1


def test_audit_h5ad_reads_the_file_and_records_it(small_null_adata, tmp_path):
    """``audit_h5ad`` audits what it read, and the payload says which file and how long the read
    took: the path, the load stage and a total runtime that includes it."""
    path = tmp_path / "stratum.h5ad"
    small_null_adata.write_h5ad(path)

    settings = _settings(celltype_col=None, celltype_value=None, design_only=True)
    payload = audit.audit_h5ad(path, settings)
    audit_schema.validate(payload)

    in_memory = run_audit(small_null_adata, settings)
    assert payload["input"]["path"] == str(path)
    assert payload["input"]["n_cells_loaded"] == in_memory["input"]["n_cells_loaded"]
    assert payload["design"]["donors_per_group"] == in_memory["design"]["donors_per_group"]
    assert payload["status"] == in_memory["status"] == "design_only"
    assert in_memory["input"]["path"] is None

    load_seconds = payload["runtime_by_stage_seconds"]["load"]
    assert load_seconds is not None and load_seconds > 0
    assert payload["runtime_seconds"] >= load_seconds
    assert audit.audit_h5ad(str(path), settings)["input"]["path"] == str(path)
