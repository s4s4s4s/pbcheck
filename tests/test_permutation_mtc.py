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
