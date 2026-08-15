"""Unit tests for the permutation constructor and the shared BH — fast, no DE runs."""

import numpy as np
import pandas as pd

from pbcheck.permutation import build_perms, labels_for, _labels_for
from pbcheck import mtc
from pbcheck.methods.de import DEResult


def test_perms_exclude_identity_and_complement_small_D():
    donors = list("abcdef")
    true_test = {"a", "b", "c"}
    perms = build_perms(donors, true_test, n_perm=100)
    # C(6,3) = 20, minus the identity and its complement = 18.
    assert len(perms) == 18
    assert true_test not in perms
    assert (set(donors) - true_test) not in perms
    # All permutations keep the group size fixed.
    assert all(len(p) == 3 for p in perms)


def test_perms_5v5_with_default_n_perm_terminates():
    """Regression: build_perms used to hang forever on a 5v5 design with default n_perm.

    C(10,5) = 252, so only 250 distinct non-trivial label-sets exist (252 minus the identity and
    its complement) -- fewer than the default n_perm=1000. The old code always rejection-sampled
    once total > max_enumerate (default 200), so `while len(perms) < n_perm` could never reach
    1000 distinct draws from a pool of 250 and looped forever. The fix falls back to full
    enumeration whenever the available distinct-set count is <= n_perm, so this must now return
    promptly with exactly the 250 available sets.
    """
    donors = [f"d{i}" for i in range(10)]
    true_test = set(donors[:5])
    perms = build_perms(donors, true_test)  # default n_perm=1000, max_enumerate=200
    assert len(perms) == 250
    assert true_test not in perms
    assert (set(donors) - true_test) not in perms
    assert all(len(p) == 5 for p in perms)
    assert len(set(perms)) == 250  # distinct


def test_perms_sampled_when_large_D():
    donors = [f"d{i}" for i in range(16)]  # C(16,8) huge -> sampling path
    true_test = set(donors[:8])
    perms = build_perms(donors, true_test, n_perm=30, seed=0)
    assert len(perms) == 30
    assert all(len(p) == 8 for p in perms)
    assert len(set(perms)) == 30  # distinct


def test_labels_for_maps_test_and_ref():
    labels = labels_for(pd.Index(["a", "b", "c"]), {"a", "c"}, "disease", "ctrl")
    assert list(labels) == ["disease", "ctrl", "disease"]


def test_private_labels_for_alias_still_resolves():
    """``scripts/pb_calibration_probe.py`` imports the private name and is left untouched."""
    assert _labels_for is labels_for


def test_bh_over_universe_length_and_missing():
    tab = pd.DataFrame(
        {"pval": [0.001, 0.5], "padj": [np.nan, np.nan], "log2fc": [3.0, 0.1]},
        index=pd.Index(["g1", "g2"], name="gene"),
    )
    res = DEResult(method="x", table=tab, contrast=("condition", "disease", "ctrl"))
    universe = ["g1", "g2", "g3"]  # g3 untested
    out = mtc.bh_over_universe(res, universe)
    assert len(out.table) == 3
    assert np.isnan(out.table.loc["g3", "padj"])  # untested -> NaN, never significant
    assert "g3" not in out.significant()
    assert "g1" in out.significant()  # strong signal survives BH over 2 tested genes


def _res(method, genes, pvals):
    return DEResult(
        method=method,
        table=pd.DataFrame(
            {"pval": np.asarray(pvals, float), "padj": np.nan, "log2fc": 0.0},
            index=pd.Index(genes, name="gene"),
        ),
        contrast=("condition", "disease", "ctrl"),
    )


def test_bh_both_arms_equalises_the_denominator():
    """One NaN in one arm must not give the arms different BH denominators.

    This is the fairness property the whole comparison rests on: if the arms are corrected over
    different m, the pseudobulk arm becomes systematically more liberal per gene than the arm it
    is the control for. The naive Wilcoxon arm never NaNs; DESeq2 does on degenerate strata.
    """
    universe = [f"g{i}" for i in range(5)]
    p = [0.001, 0.01, 0.02, 0.03, 0.04]
    naive = _res("naive", universe, p)
    pb = _res("pseudobulk", universe, [np.nan] + p[1:])  # one gene untested on one side only

    paired = mtc.bh_both_arms(naive, pb, universe, alpha=0.05)

    assert paired.n_na_pseudobulk == 1 and paired.n_na_naive == 0
    assert paired.dropped_for_fairness == ["g0"]
    assert paired.n_tested_common == 4

    # Same m on both sides => identical padj on the genes both arms tested.
    a = paired.naive.table["padj"].to_numpy()[1:]
    b = paired.pseudobulk.table["padj"].to_numpy()[1:]
    assert np.allclose(a, b)
    # The gene only one arm tested is excluded from BOTH, never counted in one arm alone.
    assert np.isnan(paired.naive.table["padj"].to_numpy()[0])
    assert np.isnan(paired.pseudobulk.table["padj"].to_numpy()[0])

    # And the old per-arm path really would have diverged, which is why this function exists.
    solo_naive = mtc.bh_over_universe(naive, universe, alpha=0.05).table["padj"].to_numpy()[1:]
    assert not np.allclose(solo_naive, b)


def test_bh_both_arms_strict_mode_raises_on_disagreement():
    universe = ["g0", "g1", "g2"]
    naive = _res("naive", universe, [0.01, 0.02, 0.03])
    pb = _res("pseudobulk", universe, [np.nan, 0.02, 0.03])
    try:
        mtc.bh_both_arms(naive, pb, universe, strict=True)
    except ValueError as e:
        assert "disagree on testable genes" in str(e)
    else:
        raise AssertionError("strict mode must raise when the arms disagree")


def test_bh_complete_null_reference_is_median_zero_not_alpha_times_G():
    """Spec B1: under the complete null BH rejects nothing in most replicates.

    Pinned by a test because the deleted `alpha * G` benchmark (~750 rejections at G = 15000) is
    exactly what a later reader is likely to reintroduce when writing the GO/NO-GO logic.
    """
    ref = mtc.bh_complete_null_reference(alpha=0.05)
    rng = np.random.default_rng(0)
    G, reps = 1500, 400
    counts = [(mtc.bh(rng.uniform(size=G), alpha=0.05) < 0.05).sum() for _ in range(reps)]

    assert np.median(counts) == ref["median_rejections"] == 0
    assert np.mean(counts) < 1.0, f"mean rejections {np.mean(counts)} should be of order alpha"
    assert np.mean(counts) < 0.05 * G / 100  # nowhere near alpha*G = 75


def test_monte_carlo_error_and_split_check_are_reported():
    """Spec §4 requires the MC error and that the real split be checked against the permutations.

    Both were absent: the 'SE << floor gap' calibration check could only ever pass by assertion,
    and the per-permutation cell totals were logged for one arm and never consumed by anything.
    """
    import warnings
    warnings.filterwarnings("ignore")
    from oracles import null_oracle
    from pbcheck.permutation import run_null
    from pbcheck.gene_universe import frozen_universe
    from pbcheck.methods.pseudobulk import build_pseudobulk

    o = null_oracle(n_donors_per_group=4, n_cells_per_donor=60, n_genes=200, seed=11)
    uni = frozen_universe(build_pseudobulk(o.adata))
    r = run_null(o.adata, uni, n_perm=10, n_perm_pb=6, fdr=0.05)

    # Both arms log their per-group cell totals now, not just the naive one.
    assert len(r["cell_totals"]) == r["n_perm_naive"]
    assert len(r["pb_cell_totals"]) == r["n_perm_pb"]

    # The real split is actually compared against the permutation distribution.
    assert r["real_split_inside_perm_range"] is True
    assert 0.0 <= r["real_split_percentile_in_perms"] <= 1.0

    mc = r["monte_carlo"]
    for k in ("naive_floor_mc_se", "pb_fp_rate", "pb_fp_rate_mc_se", "floor_gap"):
        assert k in mc
    # A binomial rate over n draws must carry sqrt(p(1-p)/n), so a handful of permutations
    # cannot resolve alpha from a few times alpha -- the point of reporting it at all.
    p, n = mc["pb_fp_rate"], r["n_perm_pb"]
    assert np.isclose(mc["pb_fp_rate_mc_se"], np.sqrt(p * (1 - p) / n))
    # And the headline gap must be large relative to its own sampling error.
    assert mc["floor_gap_over_mc_se"] > 5


def test_run_null_carries_the_paired_bh_bookkeeping():
    """Amendment 2 Change 2: cross-arm numbers come from one common tested set, and it is visible.

    The erratum being closed is that ``bh_both_arms`` existed, was tested, and no production caller
    used it. A test that only checked #DEG would not have caught that -- so this asserts the paired
    bookkeeping is actually present in the output, which is only possible if the paired path ran.
    """
    import warnings
    warnings.filterwarnings("ignore")
    from oracles import null_oracle
    from pbcheck.permutation import run_null
    from pbcheck.gene_universe import frozen_universe
    from pbcheck.methods.pseudobulk import build_pseudobulk

    o = null_oracle(n_donors_per_group=4, n_cells_per_donor=60, n_genes=200, seed=12)
    uni = frozen_universe(build_pseudobulk(o.adata), min_size=50)
    r = run_null(o.adata, uni, n_perm=8, n_perm_pb=5, fdr=0.05)

    assert r["monte_carlo"]["bh_mode"].startswith("paired")
    assert len(r["paired_bh"]) == r["n_perm_paired"] == 5
    for bk in r["paired_bh"] + [r["paired_bh_real"]]:
        assert set(bk) == {"n_tested_common", "n_na_naive", "n_na_pseudobulk",
                           "n_dropped_for_fairness"}
        # The common set is what BH saw on BOTH sides; it can never exceed the universe.
        assert 0 < bk["n_tested_common"] <= len(uni)
        assert bk["n_dropped_for_fairness"] == len(uni) - bk["n_tested_common"]

    # Real-label p-values are carried so the B5 empirical construction has a reference (spec §4).
    assert r["naive_pvals_real"].shape == (len(uni),)
    assert r["pb_pvals_real"].shape == (len(uni),)
    # The realised prior of the moderated arm rides along (the fitType analog, spec C5).
    assert "prior_df_d0" in r["pb_moderation"]

    # Naive-only permutations beyond the paired set are reported separately, not pooled into the
    # cross-arm comparison.
    assert len(r["naive_ndeg_solo"]) == r["n_perm_naive"] == 8
