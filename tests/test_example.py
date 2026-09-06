"""Tests for :mod:`pbcheck.example`, the offline quickstart generator.

These pin the contract the CLI's ``pbcheck example`` command and the quickstart docs
rely on: the shape produced is exactly what :func:`pbcheck.audit.run_audit` needs (raw
integer counts, a donor column, a two-level condition column, a cell-type column), the
draw is deterministic for a given seed and differs across seeds, and running the audit
on the small shape reaches ``status == "complete"`` with the naive arm flagged as
inflated, the one property this generator exists to demonstrate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pbcheck.audit import AuditSettings, run_audit
from pbcheck.example import REFERENCE_SHAPE, SMALL_SHAPE, example_adata


def test_shape_and_obs_columns():
    adata = example_adata(seed=0)
    assert adata.n_obs == 2 * SMALL_SHAPE["n_donors_per_group"] * SMALL_SHAPE["n_cells_per_donor"]
    assert adata.n_vars == SMALL_SHAPE["n_genes"]
    for col in ("donor", "condition", "cell_type"):
        assert col in adata.obs.columns
    assert set(adata.obs["condition"].unique()) == {"ctrl", "disease"}
    assert set(adata.obs["cell_type"].unique()) == {"example"}
    assert adata.obs["donor"].nunique() == 2 * SMALL_SHAPE["n_donors_per_group"]


def test_x_is_raw_integer_counts():
    adata = example_adata(seed=0)
    x = np.asarray(adata.X)
    assert np.issubdtype(x.dtype, np.integer)
    assert (x >= 0).all()
    # Not degenerate: some genes are actually expressed.
    assert x.sum() > 0


def test_deterministic_per_seed():
    a = example_adata(seed=3)
    b = example_adata(seed=3)
    np.testing.assert_array_equal(np.asarray(a.X), np.asarray(b.X))
    pd.testing.assert_series_equal(
        a.obs["donor"].astype(str).reset_index(drop=True),
        b.obs["donor"].astype(str).reset_index(drop=True),
    )


def test_different_across_seeds():
    a = example_adata(seed=0)
    b = example_adata(seed=1)
    assert not np.array_equal(np.asarray(a.X), np.asarray(b.X))


def test_reference_shape_constant_matches_plan():
    assert REFERENCE_SHAPE == {"n_genes": 8000, "n_donors_per_group": 8, "n_cells_per_donor": 625}


def test_small_shape_audit_reaches_complete_and_is_inflated():
    adata = example_adata(seed=0)
    settings = AuditSettings(
        donor_col="donor",
        condition_col="condition",
        test_level="disease",
        ref_level="ctrl",
        n_perm=20,
        n_perm_pb=10,
    )
    payload = run_audit(adata, settings)
    readout = payload["readout"]
    assert payload["status"] == "complete"
    assert readout["lambda_naive_class"] == "above_band"
    # The three floors this generator exists to clear: the naive solo floor,
    # the naive paired floor (only rendered when the pseudobulk arm NaNs no
    # gene), and the pseudobulk floor itself.
    naive_floor_solo = readout["naive_floor_solo"]
    assert naive_floor_solo is not None
    assert naive_floor_solo["median_count"] is not None
    assert readout["paired_floor_shown"] is True
    assert readout["pseudobulk_real_over_floor"] is not None
