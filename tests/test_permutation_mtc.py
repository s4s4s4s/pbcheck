"""Unit tests for the permutation constructor and the shared BH — fast, no DE runs."""

import numpy as np
import pandas as pd

from pbcheck.permutation import build_perms, _labels_for
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


def test_perms_sampled_when_large_D():
    donors = [f"d{i}" for i in range(16)]  # C(16,8) huge -> sampling path
    true_test = set(donors[:8])
    perms = build_perms(donors, true_test, n_perm=30, seed=0)
    assert len(perms) == 30
    assert all(len(p) == 8 for p in perms)
    assert len(set(perms)) == 30  # distinct


def test_labels_for_maps_test_and_ref():
    labels = _labels_for(pd.Index(["a", "b", "c"]), {"a", "c"}, "disease", "ctrl")
    assert list(labels) == ["disease", "ctrl", "disease"]


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
