"""The counts gate reads X, drops what is not raw counts, and corrects nothing silently.

No network and no ``cellxgene-census`` anywhere in this file. Every test builds its own AnnData,
which is also the only way to exercise the cases that matter and that real data supplies rarely on
demand: a matrix that has been normalised, one carrying a NaN, a stratum whose universe collapses
below C5, a load that disagrees with the obs snapshot it was planned from. The Census fetch path
itself is covered against a stub module with the same surface.

The three things these tests really guard are the three that would be silent if wrong: that
**integral values in a float dtype pass** (Census raw is float32 and dropping it would throw away
most of the Census), that a **fractional** matrix is dropped *with its reason* rather than rounded,
and that a manifest row is never mutated in place — a filled row must be a new object, or the
committed artifact and the in-memory one drift apart with nothing to show for it.
"""

from __future__ import annotations

import ast
import copy
import csv
import json
import sys
import types
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from scipy import sparse as sp  # noqa: E402

from pbcheck import census_select as cs  # noqa: E402
from pbcheck import gate_config as gc  # noqa: E402
from pbcheck import io_counts as io  # noqa: E402
from pbcheck.methods import pseudobulk as pb  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
SPEC = (REPO / "docs" / "PHASE0_SPEC.md").read_text(encoding="utf-8")
SOURCE = (REPO / "src" / "pbcheck" / "io_counts.py").read_text(encoding="utf-8")

DISEASE = "COVID-19"
REFERENCE = "normal"


def stratum_adata(
    *,
    n_donors_per_group=3,
    n_cells=12,
    n_genes=300,
    seed=0,
    sparse_x=True,
    dataset_id="ds1",
    cell_type="T cell",
    cells_per_donor=None,
    depth_by_group=None,
    dtype=np.float32,
):
    """A small stratum with donor structure and raw integer counts in ``.X``.

    ``cells_per_donor`` overrides the per-donor cell count by donor id (to build a donor that the
    thin filter must drop); ``depth_by_group`` multiplies one group's counts (to build a
    sequencing-depth confound) and makes the counts **exact** rather than Poisson, so that the donor
    depth medians take one value per group — census_select bins depth into tertiles, and six noisy
    medians would put one donor of each group in a shared middle bin for reasons that have nothing
    to do with the confound being tested.
    """
    rng = np.random.default_rng(seed)
    donors = (
        [(f"A{i}", DISEASE) for i in range(n_donors_per_group)]
        + [(f"B{i}", REFERENCE) for i in range(n_donors_per_group)]
    )
    blocks, rows = [], []
    for donor, disease in donors:
        n = (cells_per_donor or {}).get(donor, n_cells)
        scale = (depth_by_group or {}).get(disease, 1.0)
        blocks.append(
            np.full((n, n_genes), 6.0 * scale) if depth_by_group
            else rng.poisson(6.0 * scale, size=(n, n_genes))
        )
        rows.extend(
            {
                "dataset_id": dataset_id,
                "cell_type": cell_type,
                "donor_id": donor,
                "disease": disease,
                "condition": disease,
                "assay": "10x 3' v3",
                "suspension_type": "cell",
                "tissue_general": "blood",
            }
            for _ in range(n)
        )
    X = np.vstack(blocks).astype(dtype)
    obs = pd.DataFrame(rows)
    obs.index = [f"c{i:05d}" for i in range(obs.shape[0])]
    adata = ad.AnnData(X=sp.csr_matrix(X) if sparse_x else X, obs=obs)
    adata.var_names = [f"g{i:04d}" for i in range(n_genes)]
    return adata


def manifest_row(adata, disease=DISEASE):
    """The ``census_select`` row this stratum would have produced from obs alone."""
    rows = cs.screen_strata(adata.obs.drop(columns=["condition"]))
    return next(r for r in rows if r["disease"] == disease)


class NoDensify(sp.csr_matrix):
    """A sparse matrix that refuses to become dense — the point of the sparse path."""

    def toarray(self, *a, **kw):  # pragma: no cover - the test asserts it is never reached
        raise AssertionError("X was densified")

    def todense(self, *a, **kw):  # pragma: no cover - same
        raise AssertionError("X was densified")


# ---------------------------------------------------------------------------
# Inclusion-gate item 4 (§1; §10 risk 6) — integrality, on the values and not on the dtype.
# ---------------------------------------------------------------------------

def test_integral_values_in_a_float_dtype_pass():
    """Census raw is float32 carrying 0.0/1.0/7.0. Those ARE integer counts."""
    X = sp.csr_matrix(np.array([[0, 3, 0], [5, 0, 11]], dtype=np.float32))
    check = io.check_integer_counts(X)
    assert check.passed and check.reason is None
    assert check.dtype.startswith("float32")
    assert check.manifest_value == io.INTEGER_CHECK_PASSED == "passed"
    assert check.n_values_checked == 3, "only the stored entries are scanned"
    assert "Raw integer counts confirmed at load" in SPEC


def test_an_integer_dtype_passes_without_a_floor_test():
    check = io.check_integer_counts(np.array([[0, 2], [4, 7]], dtype=np.int64))
    assert check.passed and check.n_noninteger == 0
    assert check.n_values_checked == 4


def test_fractional_values_fail_with_a_reason_and_examples():
    """§1 item 4 / §10 risk 6: a normalised or log1p matrix is DROPPED, never rounded."""
    X = sp.csr_matrix(np.array([[0.0, 1.7], [2.5, 0.0]], dtype=np.float32))
    check = io.check_integer_counts(X)
    assert not check.passed
    assert check.n_noninteger == 2
    assert "non-integer values in X" in check.reason
    assert "normalised and/or log-transformed" in check.reason
    assert check.manifest_value.startswith(io.INTEGER_CHECK_FAILED_PREFIX)
    assert len(check.examples) == 2
    with pytest.raises(io.NonIntegerCounts, match="DROPPED"):
        io.assert_integer_counts(X)


def test_the_failure_says_drop_and_never_round():
    X = sp.csr_matrix(np.array([[1.5]], dtype=np.float32))
    with pytest.raises(io.NonIntegerCounts) as exc:
        io.assert_integer_counts(X)
    assert "Do not round" in str(exc.value)
    assert "round" not in io.check_integer_counts(X).manifest_value


def test_non_finite_and_negative_values_are_separate_failures():
    nan = io.check_integer_counts(np.array([[1.0, np.nan]], dtype=np.float32))
    assert not nan.passed and nan.n_nonfinite == 1 and "non-finite" in nan.reason

    neg = io.check_integer_counts(np.array([[1.0, -2.0]], dtype=np.float32))
    assert not neg.passed and neg.n_negative == 1 and "negative" in neg.reason
    assert neg.n_noninteger == 0, "-2.0 is integral; it fails for being negative, not fractional"


def test_a_sparse_matrix_is_never_densified():
    dense = np.zeros((400, 500), dtype=np.float32)
    dense[3, 4] = 7.0
    X = NoDensify(sp.csr_matrix(dense))
    assert io.check_integer_counts(X).passed
    assert io.counts_per_cell(X)[3] == 7.0


def test_implicit_zeros_are_not_checked_and_an_all_zero_x_says_so():
    X = sp.csr_matrix((5, 9), dtype=np.float32)
    check = io.check_integer_counts(X)
    assert check.passed and check.n_values_checked == 0
    assert any("structural zero" in n for n in check.notes)


def test_an_empty_x_is_a_failure_not_a_pass():
    check = io.check_integer_counts(sp.csr_matrix((0, 0), dtype=np.float32))
    assert not check.passed and "empty" in check.reason


def test_a_dense_matrix_is_scanned_in_chunks_with_the_same_verdict():
    dense = np.arange(60, dtype=np.float64).reshape(12, 5)
    assert io.check_integer_counts(dense, chunk_rows=5).passed
    dense[7, 2] = 0.25
    bad = io.check_integer_counts(dense, chunk_rows=5)
    assert not bad.passed and bad.n_noninteger == 1


def test_a_lazy_or_backed_x_raises_rather_than_being_densified():
    class Lazy:
        shape = (10, 10)

    with pytest.raises(TypeError, match="will not densify"):
        io.check_integer_counts(Lazy())


def test_an_anndata_may_be_passed_directly():
    adata = stratum_adata(n_genes=20, n_cells=10)
    assert io.check_integer_counts(adata).passed
    assert io.check_integer_counts(adata).n_cells == adata.n_obs


# ---------------------------------------------------------------------------
# Depth from X: the per-cell counts census_select could not read out of obs.
# ---------------------------------------------------------------------------

def test_counts_per_cell_matches_the_row_sums_on_both_layouts():
    adata = stratum_adata(n_genes=30, n_cells=10, seed=3)
    sparse_sums = io.counts_per_cell(adata)
    dense = adata.copy()
    dense.X = np.asarray(adata.X.todense())
    assert np.allclose(sparse_sums, io.counts_per_cell(dense, chunk_rows=7))
    assert np.allclose(sparse_sums, np.asarray(adata.X.sum(axis=1)).ravel())


def test_median_counts_per_cell_is_keyed_by_the_census_select_group_tags():
    adata = stratum_adata(n_genes=40, n_cells=10, depth_by_group={DISEASE: 3.0}, seed=5)
    medians = io.median_counts_per_cell_by_group(
        adata, test_level=DISEASE, ref_level=REFERENCE
    )
    assert set(medians) == {"A", "B"}
    assert medians["A"] > 2 * medians["B"]


# ---------------------------------------------------------------------------
# attach_obs (§9 item 2).
# ---------------------------------------------------------------------------

def test_attach_obs_derives_condition_and_the_depth_column():
    adata = stratum_adata(n_genes=20, n_cells=10)
    del adata.obs["condition"]
    out = io.attach_obs(adata, inplace=False)
    assert cs.CONDITION_COL in out.obs.columns
    assert set(out.obs[cs.CONDITION_COL]) == {DISEASE, REFERENCE}
    assert cs.COUNTS_COLUMN in out.obs.columns
    assert np.allclose(out.obs[cs.COUNTS_COLUMN].to_numpy(), io.counts_per_cell(out))
    assert cs.CONDITION_COL not in adata.obs.columns, "inplace=False must not touch the input"


def test_attach_obs_drops_donorless_cells_and_records_how_many():
    """§1 item 3 opens with "donor_id present"; a cell without one has no replication unit."""
    adata = stratum_adata(n_genes=20, n_cells=10)
    adata.obs["donor_id"] = adata.obs["donor_id"].astype(object)
    adata.obs.iloc[:4, adata.obs.columns.get_loc("donor_id")] = None
    out = io.attach_obs(adata)
    assert out.n_obs == adata.n_obs - 4
    assert out.uns["pbcheck_obs"]["n_missing_donor_cells"] == 4
    assert out.uns["pbcheck_obs"]["missing_donor_cells_dropped"] is True


def test_attach_obs_refuses_a_third_condition_level_instead_of_subsetting():
    adata = stratum_adata(n_genes=20, n_cells=10)
    adata.obs.iloc[:3, adata.obs.columns.get_loc("condition")] = "influenza"
    contrast = cs.Contrast("ds1", "T cell", DISEASE, REFERENCE)
    with pytest.raises(ValueError, match="will not subset silently"):
        io.attach_obs(adata, contrast, inplace=False)


# ---------------------------------------------------------------------------
# Inclusion-gate item 5 / C5, and §3's post-aggregation donor count.
# ---------------------------------------------------------------------------

def test_the_frozen_universe_is_sized_from_real_counts_and_clears_c5():
    adata = stratum_adata(n_genes=300, n_cells=12, seed=7)
    check = io.frozen_universe_check(adata, test_level=DISEASE, ref_level=REFERENCE)
    assert check.passed and check.size == 300
    assert check.min_size == gc.MIN_UNIVERSE_SIZE == 200
    assert check.profiles_by_group == {"A": 3, "B": 3}
    assert check.thin_donor_filter["n_dropped"] == 0


def test_a_small_universe_is_a_c5_skip_and_still_reports_its_size():
    """§1 item 5: below the minimum gene count the stratum is SKIPPED — with the number kept."""
    adata = stratum_adata(n_genes=50, n_cells=12, seed=8)
    check = io.frozen_universe_check(adata, test_level=DISEASE, ref_level=REFERENCE)
    assert not check.passed
    assert check.skip_statuses == (io.STATUS_SKIPPED_UNIVERSE_TOO_SMALL,)
    assert check.size == 50, "the measured size is kept: 'too small' without a number is not audit"
    assert "C5" in check.reasons[0] and "SKIP" in check.reasons[0]


def test_the_thin_donor_filter_is_the_arms_own_at_the_pre_registered_thresholds():
    """Amendment 2 Change 7: dropped, not merged — and at one set of thresholds, not two."""
    assert io.gc.MIN_CELLS == pb.MIN_CELLS == 10
    assert io.gc.MIN_COUNTS == pb.MIN_COUNTS == 1000
    assert "MIN_CELLS" in SOURCE and "min_cells: int = gc.MIN_CELLS" in SOURCE

    adata = stratum_adata(n_genes=300, n_cells=12, cells_per_donor={"A0": 4}, seed=9)
    check = io.frozen_universe_check(adata, test_level=DISEASE, ref_level=REFERENCE)
    assert check.thin_donor_filter["n_dropped"] == 1
    assert "dropped, not merged" in check.thin_donor_filter["semantics"]
    assert check.profiles_by_group["A"] == 2
    assert io.STATUS_SKIPPED_THIN_PSEUDOBULK in check.skip_statuses
    assert any("§3" in r for r in check.reasons)


# ---------------------------------------------------------------------------
# The §1 pre-screen's depth bin, which needs counts.
# ---------------------------------------------------------------------------

def test_the_depth_bin_is_screened_from_x_and_never_excludes():
    adata = stratum_adata(n_genes=40, n_cells=12, depth_by_group={DISEASE: 5.0}, seed=11)
    screen = io.depth_screen(adata)
    assert screen.bin_column == cs.DEPTH_BIN_COL
    assert screen.separates, "depth split exactly along condition is a perfect separation"
    assert screen.source == "row sums of the loaded X"
    assert cs.EXCLUDING_COVARIATES == ("assay", "suspension_type"), "depth is not an excluder"


def test_the_depth_screen_reuses_census_selects_own_prescreen():
    """One implementation of the donor-level V and of the singleton-level rule, not two."""
    assert "cs.confound_prescreen(" in SOURCE
    for reimplementation in ("from pbcheck.design import", "_perfectly_separates", "qcut"):
        assert reimplementation not in SOURCE, reimplementation


# ---------------------------------------------------------------------------
# Filling the manifest row — the pending columns census_select owed to this module.
# ---------------------------------------------------------------------------

def test_every_pending_column_this_module_owes_is_filled():
    adata = stratum_adata(n_genes=300, n_cells=12, seed=13)
    row = manifest_row(adata)
    assert row["integer_check"] == cs.PENDING
    assert row["frozen_universe_size"] == cs.PENDING
    assert row["median_counts_per_cell_by_group"] == cs.PENDING
    assert cs.DEPTH_BIN_COL not in row["confound_cramers_v"]

    out = io.update_manifest_row(row, adata)
    assert out["integer_check"] == "passed"
    assert out["frozen_universe_size"] == 300
    assert set(out["median_counts_per_cell_by_group"]) == {"A", "B"}
    assert cs.DEPTH_BIN_COL in out["confound_cramers_v"]
    assert out["gate_status"] == "candidate"
    assert out["counts_provenance"]["layer"] == io.RAW_LAYER == "raw"
    assert "X is loaded from the Census **raw** layer" in SPEC


def test_the_columns_this_module_does_not_own_stay_pending():
    """fitType belongs to the arm (C5 / Amendment 2 Change 4); ontology depth to controls (D5)."""
    adata = stratum_adata(n_genes=300, n_cells=12, seed=14)
    out = io.update_manifest_row(manifest_row(adata), adata)
    assert out["fitType"] == cs.PENDING
    assert out["cell_type_ontology_depth"] == cs.PENDING
    assert out["sigma_donor_estimate"] == cs.PENDING
    assert out["envelope_membership"] == cs.PENDING


def test_no_row_is_admitted_by_filling_items_4_and_5():
    """Amendment 3 Change 1: sigma_donor and envelope membership are still OPEN."""
    adata = stratum_adata(n_donors_per_group=20, n_genes=300, n_cells=12, seed=15)
    out = io.update_manifest_row(manifest_row(adata), adata)
    assert out["admitted_to_sweep"] is False
    assert "integer_check" not in out["admission_blockers"]
    assert "frozen_universe_size" not in out["admission_blockers"]
    assert out["admission_blockers"] == ["sigma_donor_estimate", "envelope_membership"]


def test_the_row_is_not_mutated():
    adata = stratum_adata(n_genes=300, n_cells=12, seed=16)
    row = manifest_row(adata)
    before = copy.deepcopy(row)
    out = io.update_manifest_row(row, adata)

    assert row == before, "the input row must survive byte for byte"
    assert out is not row
    assert out["confound_flags"] is not row["confound_flags"]
    assert out["admission_blockers"] is not row["admission_blockers"]
    assert out["gate_failures"] is not row["gate_failures"]
    assert out["frozen_universe_size"] != row["frozen_universe_size"]


def test_a_non_integer_stratum_is_dropped_with_its_reason_and_nothing_is_derived_from_it():
    adata = stratum_adata(n_genes=300, n_cells=12, seed=17)
    row = manifest_row(adata)
    adata.X = adata.X.multiply(0.5).tocsr()  # a CPM-like rescale: no longer counts

    out = io.update_manifest_row(row, adata)
    assert out["gate_status"] == io.STATUS_DROPPED_NON_INTEGER
    assert out["integer_check"].startswith("failed: ")
    assert any("inclusion gate item 4" in f for f in out["gate_failures"])
    # Nothing computed off the magnitudes of a matrix we have just refused.
    assert out["frozen_universe_size"] == io.NOT_COMPUTED
    assert out["median_counts_per_cell_by_group"] == io.NOT_COMPUTED
    assert cs.DEPTH_BIN_COL not in out["confound_cramers_v"]
    assert io.NOT_COMPUTED != cs.PENDING, "'not measured on purpose' is not 'owed by someone'"
    assert io.STATUS_DROPPED_NON_INTEGER in out["admission_blockers"]


def test_a_small_universe_marks_the_row_skipped_per_c5():
    adata = stratum_adata(n_genes=50, n_cells=12, seed=18)
    out = io.update_manifest_row(manifest_row(adata), adata)
    assert out["gate_status"] == io.STATUS_SKIPPED_UNIVERSE_TOO_SMALL
    assert out["frozen_universe_size"] == 50
    assert out["frozen_universe_detail"]["min_size"] == gc.MIN_UNIVERSE_SIZE
    assert any("C5" in f for f in out["gate_failures"])


def test_an_exclusion_already_on_the_row_is_not_overwritten():
    adata = stratum_adata(n_genes=50, n_cells=12, seed=19)
    row = manifest_row(adata)
    row["gate_status"] = "excluded_confound"
    out = io.update_manifest_row(row, adata)
    assert out["gate_status"] == "excluded_confound", "the obs-side exclusion happened first"
    assert io.STATUS_SKIPPED_UNIVERSE_TOO_SMALL in out["admission_blockers"]


def test_a_snapshot_median_from_obs_is_kept_and_the_x_value_recorded_beside_it():
    adata = stratum_adata(n_genes=300, n_cells=12, seed=20)
    obs = adata.obs.drop(columns=["condition"]).copy()
    obs[cs.COUNTS_COLUMN] = io.counts_per_cell(adata)
    row = next(r for r in cs.screen_strata(obs) if r["disease"] == DISEASE)
    assert isinstance(row["median_counts_per_cell_by_group"], dict)

    out = io.update_manifest_row(row, adata)
    assert out["median_counts_per_cell_by_group"] == row["median_counts_per_cell_by_group"]
    assert "median_counts_per_cell_by_group_from_X" in out["counts_provenance"]
    assert out["load_vs_obs_snapshot"]["agrees"] is True


def test_update_manifest_row_requires_attach_obs_rather_than_deriving_columns_itself():
    adata = stratum_adata(n_genes=40, n_cells=12, seed=21)
    row = manifest_row(adata)
    del adata.obs["condition"]
    with pytest.raises(ValueError, match="attach_obs"):
        io.update_manifest_row(row, adata)


# ---------------------------------------------------------------------------
# The obs snapshot vs the load — flagged, never quietly adopted.
# ---------------------------------------------------------------------------

def test_a_donor_missing_from_the_load_is_flagged_and_the_row_keeps_the_snapshots_numbers():
    adata = stratum_adata(n_genes=300, n_cells=12, seed=22)
    row = manifest_row(adata)
    short = adata[adata.obs["donor_id"] != "A0"].copy()

    out = io.update_manifest_row(row, short)
    report = out["load_vs_obs_snapshot"]
    assert report["agrees"] is False
    assert any(f.startswith("n_donors_A:") for f in report["flags"])
    assert any(f.startswith("n_cells:") for f in report["flags"])
    assert out["n_donors_A"] == row["n_donors_A"] == 3, "no silent correction"
    assert out["n_cells"] == row["n_cells"]


def test_the_gate_result_form_compares_donor_ids_not_just_counts():
    adata = stratum_adata(n_genes=40, n_cells=12, seed=23)
    obs = adata.obs.drop(columns=["condition"])
    contrast = next(c for c in cs.build_contrasts(obs) if c.disease == DISEASE)
    gate = cs.apply_inclusion_gate(obs, contrast)

    swapped = adata.copy()
    swapped.obs["donor_id"] = swapped.obs["donor_id"].replace({"A0": "A9"})
    report = io.reconcile_with_obs_snapshot(gate, swapped)
    assert report.agrees is False
    assert report.donors_missing["A"] == ["A0"]
    assert report.donors_unexpected["A"] == ["A9"]


def test_a_donor_the_obs_gate_dropped_as_thin_is_not_a_surprise_in_the_load():
    """A load filters on dataset/cell type/disease, not on donor: thin donors come back."""
    adata = stratum_adata(n_genes=300, n_cells=12, cells_per_donor={"A0": 4}, seed=24)
    obs = adata.obs.drop(columns=["condition"])
    contrast = next(c for c in cs.build_contrasts(obs) if c.disease == DISEASE)
    gate = cs.apply_inclusion_gate(obs, contrast)
    assert gate.dropped_thin_donors == {"A0": 4}

    report = io.reconcile_with_obs_snapshot(gate, adata)
    assert report.donors_dropped_by_gate["A"] == ["A0"]
    assert "A" not in report.donors_unexpected
    assert not any("donors_A" in f for f in report.flags)


def test_no_magnitude_of_a_refused_matrix_appears_inside_the_reconciliation():
    """Regression: the "nothing derived from a failed matrix" rule binds at every depth of the row.

    The top-level column was already ``not_computed``, while
    ``load_vs_obs_snapshot['loaded']['median_counts_per_cell_by_group']`` still carried real medians
    computed off the very matrix the row had just refused — the same violation, one dict level down.
    """
    adata = stratum_adata(n_genes=300, n_cells=12, seed=40)
    row = manifest_row(adata)
    adata.X = adata.X.multiply(0.5).tocsr()

    out = io.update_manifest_row(row, adata)
    loaded = out["load_vs_obs_snapshot"]["loaded"]
    assert out["integer_check"].startswith("failed: ")
    assert loaded["median_counts_per_cell_by_group"] == io.NOT_COMPUTED
    assert not [v for v in loaded.values() if isinstance(v, float)]
    assert "integer-counts gate" in loaded["magnitudes_note"]
    assert not any(f.startswith("median_counts_per_cell_") for f in
                   out["load_vs_obs_snapshot"]["flags"])
    # ...and the donor/cell comparison still ran: those numbers do not depend on the values.
    assert loaded["n_cells"] == adata.n_obs
    assert loaded["n_donors"] == {"A": 3, "B": 3}


def test_a_non_integer_rescale_is_not_re_reported_as_a_depth_discrepancy():
    """Regression: the depth cross-check must not restate the integrality failure under a new name.

    With the snapshot's ``raw_sum`` medians present, a ×0.5 rescale used to surface as
    ``agrees=False`` with "obs snapshot 1809, X row sums 904.5" — a fabricated *depth* disagreement
    that is the non-integrality already reported by ``integer_check``, counted twice.
    """
    adata = stratum_adata(n_genes=300, n_cells=12, seed=41)
    obs = adata.obs.drop(columns=["condition"]).copy()
    obs[cs.COUNTS_COLUMN] = io.counts_per_cell(adata)
    row = next(r for r in cs.screen_strata(obs) if r["disease"] == DISEASE)
    assert isinstance(row["median_counts_per_cell_by_group"], dict)

    adata.X = adata.X.multiply(0.5).tocsr()
    out = io.update_manifest_row(row, adata)
    report = out["load_vs_obs_snapshot"]
    assert report["flags"] == []
    assert report["agrees"] is True, "the only thing wrong with this stratum is its integrality"
    assert out["integer_check"].startswith("failed: ")
    assert out["median_counts_per_cell_by_group"] == row["median_counts_per_cell_by_group"]
    assert "median_counts_per_cell_by_group_from_X" not in out["counts_provenance"]


def test_a_donor_no_longer_thin_in_the_load_is_flagged_not_netted_out():
    """Regression: carrying the snapshot's thin-donor drop across is only legal while it still holds.

    A donor dropped at 4 cells that arrives with 400 used to be netted out unconditionally, so the
    largest possible snapshot-vs-load divergence produced ``agrees=True, flags=()`` and left a trace
    only in a count nothing flagged on.
    """
    snapshot_adata = stratum_adata(n_genes=40, n_cells=12, cells_per_donor={"A0": 4}, seed=42)
    row = manifest_row(snapshot_adata)
    assert row["dropped_thin_donors"] == {"A0": 4}

    load = stratum_adata(n_genes=40, n_cells=12, cells_per_donor={"A0": 400}, seed=42)
    report = io.reconcile_with_obs_snapshot(row, load)

    assert report.agrees is False
    assert report.thin_donors_no_longer_thin == {"A": {"A0": 400}}
    flag = next(f for f in report.flags if "no longer" in f or "dropped as thin" in f)
    assert "'A0' was dropped as thin by the obs snapshot (4 cell(s))" in flag
    assert f"400 cell(s) in the load (>= {gc.MIN_CELLS})" in flag
    assert report.as_dict()["donors_dropped_as_thin_but_not_thin_in_the_load"] == {"A": {"A0": 400}}
    # The net comparison still uses the snapshot's own decision, so it is not double-reported.
    assert report.loaded["n_cells"] == report.expected["n_cells"] == 60
    assert report.loaded["n_cells_including_gate_dropped_donors"] == 460


def test_a_donor_still_thin_in_the_load_stays_unflagged():
    """The other side of the same rule: 4 cells then, 5 now — still thin, still not a surprise."""
    snapshot_adata = stratum_adata(n_genes=40, n_cells=12, cells_per_donor={"A0": 4}, seed=43)
    row = manifest_row(snapshot_adata)
    load = stratum_adata(n_genes=40, n_cells=12, cells_per_donor={"A0": 5}, seed=43)

    report = io.reconcile_with_obs_snapshot(row, load)
    assert report.agrees is True
    assert report.thin_donors_no_longer_thin == {}
    assert report.donors_dropped_by_gate == {"A": ["A0"]}


def test_a_depth_that_contradicts_the_snapshot_is_flagged():
    adata = stratum_adata(n_genes=100, n_cells=12, seed=25)
    obs = adata.obs.drop(columns=["condition"]).copy()
    obs[cs.COUNTS_COLUMN] = io.counts_per_cell(adata) * 10  # the snapshot saw a different matrix
    row = next(r for r in cs.screen_strata(obs) if r["disease"] == DISEASE)

    report = io.reconcile_with_obs_snapshot(row, adata)
    assert report.agrees is False
    assert any(f.startswith("median_counts_per_cell_") for f in report.flags)


def test_a_row_that_came_back_through_the_csv_still_reconciles():
    """Manifest cells are JSON in the CSV; a round-tripped row must not silently lose its numbers."""
    adata = stratum_adata(n_genes=100, n_cells=12, seed=26)
    row = manifest_row(adata)
    round_tripped = dict(row)
    for key in ("cells_per_donor_by_group", "dropped_thin_donors"):
        round_tripped[key] = json.dumps(row[key])

    report = io.reconcile_with_obs_snapshot(round_tripped, adata)
    assert report.expected["n_cells_by_group"] == {"A": 36, "B": 36}
    assert report.agrees is True


# ---------------------------------------------------------------------------
# The extended manifest.
# ---------------------------------------------------------------------------

def test_the_extended_manifest_is_a_superset_and_the_csv_keeps_the_added_columns(tmp_path):
    adata = stratum_adata(n_genes=300, n_cells=12, seed=27)
    obs = adata.obs.drop(columns=["condition"])
    rows = [io.update_manifest_row(r, adata) for r in cs.screen_strata(obs)]

    assert io.MANIFEST_FIELDS[:len(cs.MANIFEST_FIELDS)] == tuple(cs.MANIFEST_FIELDS)
    assert set(rows[0]) == set(io.MANIFEST_FIELDS)

    manifest = io.emit_manifest(rows, out_dir=tmp_path, obs=obs)
    assert manifest["header"]["census_version"] == cs.CENSUS_VERSION
    assert manifest["header"]["counts_gate"]["n_integer_check_failed"] == 0

    with open(tmp_path / "census_strata.csv", encoding="utf-8", newline="") as fh:
        table = list(csv.DictReader(fh))
    assert list(table[0]) == list(io.MANIFEST_FIELDS)
    assert json.loads(table[0]["integer_check_detail"])["passed"] is True
    assert table[0]["integer_check"] == "passed"


def test_the_counts_gate_accounting_says_it_is_not_the_obs_excluded_fraction():
    adata = stratum_adata(n_genes=50, n_cells=12, seed=28)
    obs = adata.obs.drop(columns=["condition"])
    rows = [io.update_manifest_row(r, adata) for r in cs.screen_strata(obs)]
    header = io.emit_manifest(rows)["header"]

    gate = header["counts_gate"]
    assert gate["by_counts_gate_status"] == {io.STATUS_SKIPPED_UNIVERSE_TOO_SMALL: 1}
    assert gate["counts_excluded_fraction"] == 1.0
    assert header["counts"]["excluded_fraction"] == 0.0, "the obs side excluded nothing"
    assert "NOT included" in gate["note"]


# ---------------------------------------------------------------------------
# The Census fetch wrapper: lazy import, pinned version, raw layer. No network.
# ---------------------------------------------------------------------------

def test_the_module_imports_without_the_census_extra():
    tree = ast.parse(SOURCE)
    top_level = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
    names = [a.name for n in top_level if isinstance(n, ast.Import) for a in n.names]
    names += [n.module for n in top_level if isinstance(n, ast.ImportFrom)]
    assert not any(str(n).startswith("cellxgene") for n in names)


def test_the_census_version_is_imported_not_restated():
    assert io.CENSUS_VERSION == cs.CENSUS_VERSION == "2025-01-30"
    assert "2025-01-30" not in SOURCE, "the pin lives in census_select; a second copy can drift"


def test_loading_without_the_extra_says_how_to_install_it(monkeypatch):
    monkeypatch.setitem(sys.modules, "cellxgene_census", None)
    with pytest.raises(cs.CensusNotInstalled, match=r"\[census\]"):
        io.load_stratum("ds1", "T cell")


def test_the_heavy_pseudobulk_dependencies_are_imported_lazily():
    """Importing the integrality gate must not pull decoupler and pydeseq2 in with it."""
    tree = ast.parse(SOURCE)
    top_level = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
    modules = [n.module for n in top_level if isinstance(n, ast.ImportFrom)]
    assert "pbcheck.methods.pseudobulk" not in modules
    assert "from pbcheck.methods.pseudobulk import build_pseudobulk" in SOURCE


class _FakeCensus:
    def __init__(self):
        self.closed = False

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.closed = True
        return False


def _fake_census_module(adata, *, modern=True):
    calls = {}
    census = _FakeCensus()

    def open_soma(census_version=None):
        calls["open_soma"] = {"census_version": census_version}
        return census

    if modern:
        def get_anndata(census, organism, measurement_name="RNA", *, X_name="raw",
                        obs_value_filter=None, obs_column_names=None):
            calls["get_anndata"] = dict(
                census=census, organism=organism, measurement_name=measurement_name,
                X_name=X_name, obs_value_filter=obs_value_filter,
                obs_column_names=obs_column_names,
            )
            return adata
    else:
        def get_anndata(census, organism, measurement_name="RNA", *, X_name="raw",
                        obs_value_filter=None, column_names=None):
            calls["get_anndata"] = dict(
                census=census, organism=organism, measurement_name=measurement_name,
                X_name=X_name, obs_value_filter=obs_value_filter, column_names=column_names,
            )
            return adata

    module = types.SimpleNamespace(open_soma=open_soma, get_anndata=get_anndata)
    return module, calls, census


def test_load_stratum_asks_for_the_raw_layer_at_the_pinned_version(monkeypatch):
    adata = stratum_adata(n_genes=40, n_cells=12, seed=29)
    del adata.obs["condition"]
    module, calls, census = _fake_census_module(adata)
    monkeypatch.setitem(sys.modules, "cellxgene_census", module)

    got = io.load_stratum("ds1", "T cell", diseases=(DISEASE, REFERENCE))

    assert calls["open_soma"]["census_version"] == cs.CENSUS_VERSION
    call = calls["get_anndata"]
    assert call["X_name"] == "raw", "§1: X is loaded from the Census raw layer"
    assert call["measurement_name"] == "RNA"
    assert call["organism"] == cs.CENSUS_ORGANISM
    assert cs.VALUE_FILTER in call["obs_value_filter"]
    assert "dataset_id == 'ds1'" in call["obs_value_filter"]
    assert "cell_type == 'T cell'" in call["obs_value_filter"]
    assert "disease in ['COVID-19', 'normal']" in call["obs_value_filter"]
    assert set(cs.OBS_COLUMNS) <= set(call["obs_column_names"])

    assert census.closed, "a census this function opened is a census it closes"
    assert got.uns["pbcheck_load"]["census_version"] == cs.CENSUS_VERSION
    assert got.uns["pbcheck_integer_check"]["passed"] is True
    assert cs.CONDITION_COL in got.obs.columns, "attach_obs ran"


def test_load_stratum_adapts_to_the_older_column_names_argument(monkeypatch):
    adata = stratum_adata(n_genes=20, n_cells=12, seed=30)
    module, calls, _ = _fake_census_module(adata, modern=False)
    monkeypatch.setitem(sys.modules, "cellxgene_census", module)

    io.load_stratum("ds1", "T cell")
    assert set(cs.OBS_COLUMNS) <= set(calls["get_anndata"]["column_names"]["obs"])


def test_a_census_handle_passed_in_is_not_closed(monkeypatch):
    adata = stratum_adata(n_genes=20, n_cells=12, seed=31)
    module, calls, _ = _fake_census_module(adata)
    monkeypatch.setitem(sys.modules, "cellxgene_census", module)

    handle = _FakeCensus()
    io.load_stratum("ds1", "T cell", census=handle)
    assert handle.closed is False
    assert "open_soma" not in calls
    assert calls["get_anndata"]["census"] is handle


def test_a_non_integer_load_is_recorded_not_raised(monkeypatch):
    """§1 item 4's answer is a DROP in the manifest; a loader that raised would abort the sweep."""
    adata = stratum_adata(n_genes=20, n_cells=12, seed=32)
    adata.X = adata.X.multiply(0.1).tocsr()
    module, _, _ = _fake_census_module(adata)
    monkeypatch.setitem(sys.modules, "cellxgene_census", module)

    got = io.load_stratum("ds1", "T cell")
    assert got.uns["pbcheck_integer_check"]["passed"] is False
    assert "non-integer" in got.uns["pbcheck_integer_check"]["reason"]


def test_a_value_that_cannot_be_quoted_is_refused_rather_than_escaped():
    with pytest.raises(ValueError, match="quote or backslash"):
        io.stratum_value_filter("ds1", "Müller's cell")


def test_load_contrast_carries_both_levels_of_the_contrast(monkeypatch):
    adata = stratum_adata(n_genes=20, n_cells=12, seed=33)
    del adata.obs["condition"]
    module, calls, _ = _fake_census_module(adata)
    monkeypatch.setitem(sys.modules, "cellxgene_census", module)

    contrast = cs.Contrast("ds1", "T cell", DISEASE, REFERENCE)
    io.load_contrast(contrast)
    assert "disease in ['COVID-19', 'normal']" in calls["get_anndata"]["obs_value_filter"]


def test_a_contrast_from_another_stratum_is_refused(monkeypatch):
    module, _, _ = _fake_census_module(stratum_adata(n_genes=10, n_cells=12, seed=34))
    monkeypatch.setitem(sys.modules, "cellxgene_census", module)
    other = cs.Contrast("ds2", "B cell", DISEASE, REFERENCE)
    with pytest.raises(ValueError, match="not the"):
        io.load_stratum("ds1", "T cell", contrast=other)
