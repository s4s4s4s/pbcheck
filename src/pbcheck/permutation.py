"""Donor-permutation null (spec §4).

The unit of permutation is the **donor**, never the cell. One permutation reassigns the
condition-label vector across donors while holding the group sizes fixed; every cell keeps its real
donor id and its real counts, so within-donor correlation and donor structure stay intact and only
the donor<->condition association is broken. Under this sharp null the correct answer is "no
association", so the naive test's #DEG is an upper bound on its pseudoreplication false-positive floor
and each arm's p-values reveal its calibration.

For speed, donor pseudobulk profiles are aggregated **once**: they are invariant under relabeling, so
each permutation only re-runs DESeq2 on a relabeled metadata vector. The naive per-cell test does
depend on labels and is re-run each permutation (it is cheap).
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd

from pbcheck.methods.naive import naive_de
from pbcheck.methods.pseudobulk import build_pseudobulk, deseq_from_pdata
from pbcheck import mtc


def _true_map(adata, donor_col, condition_col):
    dc_map = (
        adata.obs[[donor_col, condition_col]].astype(str).drop_duplicates()
        .set_index(donor_col)[condition_col]
    )
    return dc_map  # Series: donor -> condition


def build_perms(donors, true_test_set, *, n_perm=1000, seed=0, max_enumerate=200):
    """Balanced donor-label permutations, excluding the identity and its exact complement.

    ``true_test_set`` is the set of donors that are really in the test group. Each returned
    permutation is a set of donors assigned to the test group (size = |true_test_set|).

    Enumerate outright when ``total <= max_enumerate`` (small D, cheap to list exhaustively). When
    D is larger but the number of *distinct available* sets (``total`` minus the identity and its
    complement) still does not exceed ``n_perm``, also enumerate rather than rejection-sample:
    rejection sampling draws WITHOUT replacement from a pool no bigger than the request, so the
    `while len(perms) < n_perm` loop can never fill and hangs forever (e.g. 5v5, C(10,5) = 252
    available sets, default n_perm = 1000 — this used to hang indefinitely before this fix). Only
    when more distinct sets exist than requested is rejection sampling used, where it terminates.
    """
    donors = list(donors)
    n_test = len(true_test_set)
    true_test_set = frozenset(true_test_set)
    complement = frozenset(donors) - true_test_set

    from math import comb
    total = comb(len(donors), n_test)
    available = total - (2 if complement != true_test_set else 1)
    perms: list[frozenset] = []
    if total <= max_enumerate or available <= n_perm:
        for combo in combinations(donors, n_test):
            s = frozenset(combo)
            if s == true_test_set or s == complement:
                continue
            perms.append(s)
    else:
        rng = np.random.default_rng(seed)
        seen = {true_test_set, complement}
        arr = np.array(donors, dtype=object)
        while len(perms) < n_perm:
            pick = frozenset(rng.choice(arr, size=n_test, replace=False).tolist())
            if pick in seen:
                continue
            seen.add(pick)
            perms.append(pick)
    return perms


def _labels_for(donors_index, test_set, test_level, ref_level):
    return np.where(pd.Index(donors_index).isin(test_set), test_level, ref_level)


def run_null(
    adata,
    universe,
    *,
    donor_col="donor",
    condition_col="condition",
    test_level="disease",
    ref_level="ctrl",
    celltype_col="cell_type",
    n_perm=40,
    n_perm_pb=20,
    fdr=0.05,
    seed=0,
    n_cpus=4,
):
    """Run both arms under the donor-permutation null.

    Returns a dict with per-permutation p-value matrices (aligned to ``universe``), post-BH #DEG
    arrays, per-permutation per-group cell totals, and the donor/perm bookkeeping.
    """
    tmap = _true_map(adata, donor_col, condition_col)
    donors = list(tmap.index)
    true_test = set(tmap.index[tmap == test_level])

    perms = build_perms(donors, true_test, n_perm=max(n_perm, n_perm_pb), seed=seed)
    perms_naive = perms[:n_perm]
    perms_pb = perms[:n_perm_pb]

    G = len(universe)
    uni_index = pd.Index(universe, name="gene")

    # Pre-aggregate donor pseudobulk ONCE (invariant under relabeling).
    pdata = build_pseudobulk(adata, donor_col=donor_col, celltype_col=celltype_col,
                             condition_col=condition_col)
    cells_per_donor = adata.obs.groupby(donor_col, observed=True).size()

    # ---- naive arm ----
    naive_pvals = np.full((len(perms_naive), G), np.nan)
    naive_ndeg = np.zeros(len(perms_naive), dtype=int)
    cell_totals = []
    a = adata.copy()
    for i, tset in enumerate(perms_naive):
        a.obs["_perm"] = _labels_for(a.obs[donor_col].astype(str), tset, test_level, ref_level)
        res = naive_de(a, condition_col="_perm", test_level=test_level, ref_level=ref_level,
                       genes=universe)
        res = mtc.bh_over_universe(res, universe, alpha=fdr)
        naive_pvals[i] = res.table["pval"].reindex(uni_index).to_numpy()
        naive_ndeg[i] = res.n_significant(fdr=fdr)
        n_test_cells = int(cells_per_donor[list(tset)].sum())
        cell_totals.append((n_test_cells, int(cells_per_donor.sum()) - n_test_cells))

    # ---- pseudobulk arm ----
    pb_pvals = np.full((len(perms_pb), G), np.nan)
    pb_ndeg = np.zeros(len(perms_pb), dtype=int)
    pb_cell_totals = []
    donor_names = pd.Index(pdata.obs[donor_col].astype(str)) if donor_col in pdata.obs else pdata.obs_names
    for i, tset in enumerate(perms_pb):
        n_test_cells = int(cells_per_donor[list(tset)].sum())
        pb_cell_totals.append((n_test_cells, int(cells_per_donor.sum()) - n_test_cells))
        cond_vals = _labels_for(donor_names, tset, test_level, ref_level)
        res = deseq_from_pdata(pdata, condition_col=condition_col, test_level=test_level,
                               ref_level=ref_level, universe=universe,
                               condition_values=cond_vals, fdr=fdr, n_cpus=n_cpus)
        res = mtc.bh_over_universe(res, universe, alpha=fdr)
        pb_pvals[i] = res.table["pval"].reindex(uni_index).to_numpy()
        pb_ndeg[i] = res.n_significant(fdr=fdr)

    # ---- A2: is the real split comparable to the permutations it is judged against? ----
    #
    # Full cell-count stratification (restricting the permutation set to splits whose per-group
    # cell totals match the real one) is not reachable at the donor counts real strata carry: at
    # 3v3 there are only C(6,3) - 2 = 18 permutations to begin with, and filtering on cell count
    # would leave single digits. So the spec's own alternative is used instead — log the per-group
    # cell totals for BOTH arms and check explicitly that the real split sits inside the
    # permutation distribution. If it does not, the floor is being compared against permutations
    # of a systematically different size, and the comparison is size-confounded rather than a
    # clean statement about the replication unit.
    real_test_cells = int(cells_per_donor[list(true_test)].sum())
    perm_test_cells = np.array([t for t, _ in cell_totals], dtype=float)
    inside = bool(perm_test_cells.size and
                  perm_test_cells.min() <= real_test_cells <= perm_test_cells.max())
    pct = (float((perm_test_cells <= real_test_cells).mean()) if perm_test_cells.size
           else float("nan"))

    # ---- Monte-Carlo error (spec §4) ----
    #
    # The floor and the calibration claim are estimates over a finite permutation set. Reporting
    # them without their sampling error invites reading a gap that the resolution cannot support:
    # a false-positive rate from 20 permutations has SE ~0.11 and cannot separate 0.05 from 0.15.
    n_n, n_p = len(perms_naive), len(perms_pb)
    pb_fp_rate = float(np.mean(pb_ndeg >= 1)) if n_p else float("nan")
    mc = {
        "naive_floor_median": float(np.median(naive_ndeg)) if n_n else float("nan"),
        "naive_floor_mc_se": float(np.std(naive_ndeg, ddof=1) / np.sqrt(n_n)) if n_n > 1 else float("nan"),
        "pb_floor_median": float(np.median(pb_ndeg)) if n_p else float("nan"),
        "pb_floor_mc_se": float(np.std(pb_ndeg, ddof=1) / np.sqrt(n_p)) if n_p > 1 else float("nan"),
        "pb_fp_rate": pb_fp_rate,
        "pb_fp_rate_mc_se": (float(np.sqrt(pb_fp_rate * (1 - pb_fp_rate) / n_p)) if n_p
                             else float("nan")),
    }
    # The calibration claim is only as strong as the gap between the arms relative to this error.
    gap = mc["naive_floor_median"] - mc["pb_floor_median"]
    mc["floor_gap"] = float(gap)
    mc["floor_gap_over_mc_se"] = (float(gap / mc["naive_floor_mc_se"])
                                  if mc["naive_floor_mc_se"] and np.isfinite(mc["naive_floor_mc_se"])
                                  and mc["naive_floor_mc_se"] > 0 else float("inf"))

    return {
        "naive_pvals": naive_pvals,
        "pb_pvals": pb_pvals,
        "naive_ndeg": naive_ndeg,
        "pb_ndeg": pb_ndeg,
        "cell_totals": cell_totals,
        "pb_cell_totals": pb_cell_totals,
        "real_test_cells": real_test_cells,
        "real_split_inside_perm_range": inside,
        "real_split_percentile_in_perms": pct,
        "monte_carlo": mc,
        "n_perm_naive": len(perms_naive),
        "n_perm_pb": len(perms_pb),
        "n_donors": len(donors),
        "G": G,
    }
