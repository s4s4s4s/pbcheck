"""The fast naive engine is scanpy's statistic, not an approximation of it — held to that here.

:mod:`pbcheck.methods.naive_engine` replaces one ``rank_genes_groups`` run per donor-permutation
with one ranking pass per stratum. That is only admissible if it reproduces the pre-registered
statistic (spec §2, ``method='wilcoxon'`` with ``tie_correct=True``) exactly, so these tests
compare it against :func:`pbcheck.methods.naive.naive_de` — the slow path itself, still callable —
rather than against the formula the fast path implements.

The tolerance is pinned where the arithmetic actually lands, not at a comfortable ``atol``:

* **Wilcoxon p-values: bitwise.** Ranks are exact multiples of 0.5 and the largest rank sum a
  frozen stratum can produce is ~1.5e11, so every partial sum is exactly representable and
  float64 addition is associative over them. There is no rounding to differ over, and the tests
  assert equality of the raw bits — a ``rtol`` here would hide a real divergence.
* **t-test_overestim_var and log2fc: floating-point level.** Group means and variances are
  reassociated sums (per donor, then per group, instead of one pass over the group), and
  reassociating a float sum is not bit-exact in general. ``rtol=1e-12`` on p-values is roughly
  four orders of magnitude tighter than the measured disagreement, which is zero on every case in
  this file; the looser bound is what the arithmetic guarantees, not what it delivers.

Also pinned: that the two engines agree end-to-end through :func:`pbcheck.permutation.run_null`,
that sparse and dense inputs give the same answer, and that the gene chunking — the thing that
keeps a 547665-cell stratum off the heap — cannot move a number.
"""

from __future__ import annotations

import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from hypothesis import assume, given, settings, strategies as st
from scipy import sparse

from pbcheck.methods.naive import naive_de
from pbcheck.methods.naive_engine import (
    DEFAULT_CHUNK_BYTES,
    NaiveRelabelEngine,
    _average_ranks,
)

warnings.filterwarnings("ignore")

_SETTINGS = settings(max_examples=30, deadline=None)

#: Gene shapes worth generating, because each one exercises a different branch of the ranking.
#: ``zero`` collapses the tie correction to 0 and drives scanpy's ``0/0 -> z = 0`` repair;
#: ``flat`` is constant in the raw counts (which per-cell normalisation then spreads out again, so
#: the genuinely-constant case is built explicitly in
#: :func:`test_degenerate_genes_are_reproduced_not_repaired`); ``binary`` and ``few`` make giant
#: tie blocks; ``wide`` makes almost none.
_GENE_KINDS = ("zero", "flat", "binary", "few", "wide")


def _build_counts(kinds, n_cells, rng):
    # The leading column is always strictly positive so that no cell ends up with a zero library
    # size. A zero-count cell is a scanpy warning path, not a ranking one, and real strata have
    # none -- the counts gate drops them long before this arm sees them.
    X = np.zeros((n_cells, len(kinds) + 1), dtype=np.float32)
    X[:, 0] = rng.integers(1, 60, size=n_cells)
    for j, kind in enumerate(kinds, start=1):
        if kind == "zero":
            continue
        elif kind == "flat":
            X[:, j] = float(rng.integers(1, 5))
        elif kind == "binary":
            X[:, j] = rng.integers(0, 2, size=n_cells)
        elif kind == "few":
            X[:, j] = rng.integers(0, 4, size=n_cells)
        else:
            X[:, j] = rng.integers(0, 60, size=n_cells)
    return X


def _make_adata(cells_per_donor, kinds, seed, *, sparse_x=False, X=None):
    rng = np.random.default_rng(seed)
    n_cells = int(sum(cells_per_donor))
    X = _build_counts(kinds, n_cells, rng) if X is None else X
    donors = np.repeat([f"d{i}" for i in range(len(cells_per_donor))], cells_per_donor)
    obs = pd.DataFrame({"donor": donors}, index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"g{j}" for j in range(X.shape[1])])
    a = ad.AnnData(X=sparse.csr_matrix(X) if sparse_x else X, obs=obs, var=var)
    return a


def _label(a, test_donors, *, test_level="disease", ref_level="ctrl"):
    a.obs["condition"] = np.where(a.obs["donor"].isin(list(test_donors)), test_level, ref_level)
    return a


def _reference(a, test_donors, method):
    """The slow path: scanpy, re-run from scratch, on exactly the same AnnData."""
    _label(a, test_donors)
    return naive_de(a, condition_col="condition", test_level="disease", ref_level="ctrl",
                    genes=list(a.var_names), method=method)


def _split(cells_per_donor, n_test):
    """A donor split whose two sides both clear scanpy's two-cell minimum, or ``None``."""
    donors = [f"d{i}" for i in range(len(cells_per_donor))]
    test = set(donors[:n_test])
    n_a = int(sum(c for d, c in zip(donors, cells_per_donor) if d in test))
    n_b = int(sum(cells_per_donor)) - n_a
    return test if n_a >= 2 and n_b >= 2 else None


def _assert_same_result(fast, slow, method):
    got = fast.table
    ref = slow.table.reindex(got.index)
    p_fast = got["pval"].to_numpy(dtype=float)
    p_slow = ref["pval"].to_numpy(dtype=float)

    if method == "wilcoxon":
        assert np.array_equal(p_fast, p_slow), (
            "wilcoxon p-values must be BITWISE identical -- ranks are exact half-integers and "
            f"their sums are exactly representable, so there is nothing to round: max ulp "
            f"{int(np.max(np.abs(p_fast.view(np.int64) - p_slow.view(np.int64))))}"
        )
    else:
        np.testing.assert_allclose(p_fast, p_slow, rtol=1e-12, atol=0.0)
    # Whatever the tolerance, no gene may change side of any plausible threshold.
    for alpha in (0.05, 0.01, 1e-4):
        assert np.array_equal(p_fast < alpha, p_slow < alpha)

    # Expressed fractions are integer counts over an integer denominator: exact, both ways.
    for col in ("pct_group", "pct_reference"):
        assert np.array_equal(got[col].to_numpy(), ref[col].to_numpy()), col
    # log2fc rides along in the table and is read by nothing in the decision path; it is a
    # float32 cast of reassociated float64 means, so it is checked at float32 resolution.
    np.testing.assert_allclose(got["log2fc"].to_numpy(dtype=float),
                               ref["log2fc"].to_numpy(dtype=float), rtol=1e-6, atol=1e-6)
    assert fast.method == slow.method
    assert fast.contrast == slow.contrast


# ---------------------------------------------------------------------------
# The vendored kernels against scanpy's own.
# ---------------------------------------------------------------------------

@given(
    n_cells=st.integers(min_value=2, max_value=90),
    kinds=st.lists(st.sampled_from(_GENE_KINDS), min_size=1, max_size=5),
    seed=st.integers(min_value=0, max_value=9999),
)
@_SETTINGS
def test_average_ranks_and_tie_term_match_scanpys_kernels_bitwise(n_cells, kinds, seed):
    """``_average_ranks`` is scanpy's ``rankdata`` + ``_tiecorrect``, rewritten off numba.

    Compared against the private names directly: that duplication is the one place the fast path
    restates scanpy's code instead of calling it, so it is the one place a scanpy change could move
    the statistic without any other test noticing. If the private names move this skips and the
    end-to-end comparisons below remain the binding check.
    """
    sc_kernels = pytest.importorskip("scanpy.tools._rank_genes_groups")
    if not hasattr(sc_kernels, "rankdata") or not hasattr(sc_kernels, "_tiecorrect"):
        pytest.skip("scanpy moved rankdata/_tiecorrect; the end-to-end tests still bind")

    block = _build_counts(kinds, n_cells, np.random.default_rng(seed))
    ranks, tie_sum = _average_ranks(block)
    size = np.float64(n_cells)
    tie_coef = 1.0 - tie_sum / (size ** 3 - size) if n_cells >= 2 else np.ones(len(kinds))

    ref_ranks = sc_kernels.rankdata(block)
    assert np.array_equal(ranks, ref_ranks)
    assert np.array_equal(tie_coef, sc_kernels._tiecorrect(ref_ranks))


# ---------------------------------------------------------------------------
# The engine against the slow path.
# ---------------------------------------------------------------------------

@given(
    cells_per_donor=st.lists(st.integers(min_value=1, max_value=9), min_size=2, max_size=6),
    kinds=st.lists(st.sampled_from(_GENE_KINDS), min_size=1, max_size=5),
    n_test=st.integers(min_value=1, max_value=5),
    seed=st.integers(min_value=0, max_value=9999),
)
@_SETTINGS
def test_wilcoxon_matches_the_slow_path_bitwise(cells_per_donor, kinds, n_test, seed):
    """Across donor counts, group sizes, tie structures and unequal cells per donor."""
    assume(1 <= n_test < len(cells_per_donor))
    test_donors = _split(cells_per_donor, n_test)
    assume(test_donors is not None)

    a = _make_adata(cells_per_donor, kinds, seed)
    slow = _reference(a, test_donors, "wilcoxon")
    fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                         method="wilcoxon").test(test_donors)
    _assert_same_result(fast, slow, "wilcoxon")


@given(
    cells_per_donor=st.lists(st.integers(min_value=1, max_value=9), min_size=2, max_size=6),
    kinds=st.lists(st.sampled_from(_GENE_KINDS), min_size=1, max_size=5),
    n_test=st.integers(min_value=1, max_value=5),
    seed=st.integers(min_value=0, max_value=9999),
)
@_SETTINGS
def test_ttest_overestim_var_matches_the_slow_path(cells_per_donor, kinds, n_test, seed):
    """Spec §2's robustness variant, reconstructed from n_d, sum x and sum x^2.

    Admissible only because scanpy's variance comes from ``E[x^2] - E[x]^2``
    (``fast_array_utils.stats.mean_var``, both the dense and the sparse kernel), which the per-donor
    sufficient statistics reproduce. A two-pass variance would be better conditioned and would be a
    *different* number from the pre-registered one, so it is not used.
    """
    assume(1 <= n_test < len(cells_per_donor))
    test_donors = _split(cells_per_donor, n_test)
    assume(test_donors is not None)

    a = _make_adata(cells_per_donor, kinds, seed)
    slow = _reference(a, test_donors, "t-test_overestim_var")
    fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                         method="t-test_overestim_var").test(test_donors)
    _assert_same_result(fast, slow, "t-test_overestim_var")


@pytest.mark.parametrize("method", ["wilcoxon", "t-test_overestim_var"])
def test_degenerate_genes_are_reproduced_not_repaired(method):
    """All-zero and constant genes: tie coefficient 0, std_dev 0, and scanpy's ``0/0 -> z = 0``.

    The engine must land on p = 1 there for the same reason scanpy does. Getting this wrong in the
    optimistic direction turns every silent gene in a stratum into a rejection, which under a
    permutation null is precisely the false-positive floor being measured.

    "Constant" has to be built deliberately: a constant *raw count* is not constant after per-cell
    normalisation, because each cell is divided by its own library size. Here every cell carries
    the same total (g1 + g2 + g3 = 9), so g1's normalised value is identical in all of them.
    """
    n_cells = 16
    varying = np.arange(n_cells, dtype=np.float32) % 7
    X = np.stack([
        np.zeros(n_cells, dtype=np.float32),          # g0: all zero
        np.full(n_cells, 3.0, dtype=np.float32),      # g1: constant after normalisation too
        varying,                                      # g2
        6.0 - varying,                                # g3: keeps every library size at 9
    ], axis=1)
    a = _make_adata([4, 4, 4, 4], [], seed=3, X=X)
    test_donors = {"d0", "d1"}
    slow = _reference(a, test_donors, method)
    fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                         method=method).test(test_donors)
    _assert_same_result(fast, slow, method)
    assert fast.table.loc["g0", "pval"] == 1.0      # all-zero
    assert fast.table.loc["g1", "pval"] == 1.0      # constant


@pytest.mark.parametrize("method", ["wilcoxon", "t-test_overestim_var"])
def test_single_donor_dominated_stratum(method):
    """One donor carrying almost every cell — the shape that makes the naive arm's floor explode."""
    a = _make_adata([200, 3, 2, 4], ["binary", "few", "wide", "zero", "flat"], seed=5)
    for test_donors in ({"d0"}, {"d1", "d2"}, {"d1", "d2", "d3"}):
        slow = _reference(a, test_donors, method)
        fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                             method=method).test(test_donors)
        _assert_same_result(fast, slow, method)


@pytest.mark.parametrize("method", ["wilcoxon", "t-test_overestim_var"])
def test_sparse_input_gives_the_same_answer_as_dense(method):
    """Two claims, because they are different claims.

    Against scanpy on the *same* storage format the engine must be exact — scanpy's own numbers
    differ slightly between formats (its docs say so, and its sparse and dense mean/variance
    kernels reduce in different orders), so comparing engine-on-sparse to scanpy-on-dense would
    measure scanpy against itself. Across formats the engine must agree with itself to
    floating-point.
    """
    kinds = ["zero", "binary", "few", "wide", "flat"]
    dense = _make_adata([6, 5, 4, 7], kinds, seed=11)
    spars = _make_adata([6, 5, 4, 7], kinds, seed=11, sparse_x=True)
    assert sparse.issparse(spars.X) and not sparse.issparse(dense.X)
    test_donors = {"d0", "d3"}

    for a in (dense, spars):
        slow = _reference(a, test_donors, method)
        fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                             method=method).test(test_donors)
        _assert_same_result(fast, slow, method)

    p_dense = NaiveRelabelEngine.from_adata(dense, donor_col="donor", method=method) \
        .test(test_donors).table["pval"].to_numpy()
    p_sparse = NaiveRelabelEngine.from_adata(spars, donor_col="donor", method=method) \
        .test(test_donors).table["pval"].to_numpy()
    np.testing.assert_allclose(p_sparse, p_dense, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("chunk_bytes", [1, 64, 4096, DEFAULT_CHUNK_BYTES])
def test_gene_chunking_cannot_move_a_number(chunk_bytes):
    """Chunking is over genes and a rank is a whole-column quantity, so width must be inert.

    ``chunk_bytes=1`` forces one gene per chunk, which is the configuration a 547665-cell stratum
    degrades to and the one no other test would reach.
    """
    kinds = ["zero", "binary", "few", "wide", "flat", "binary"]
    a = _make_adata([5, 6, 4, 5], kinds, seed=17, sparse_x=True)
    test_donors = {"d0", "d2"}
    slow = _reference(a, test_donors, "wilcoxon")
    fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                         method="wilcoxon",
                                         chunk_bytes=chunk_bytes).test(test_donors)
    _assert_same_result(fast, slow, "wilcoxon")


def test_gene_restriction_happens_before_normalisation():
    """The universe restriction changes every value, so it must be applied where naive_de does.

    Normalising over all genes and then subsetting gives different library sizes from subsetting
    and then normalising. ``naive_de`` does the latter; if the engine did the former this test's
    two results would disagree while every all-genes test still passed.
    """
    a = _make_adata([5, 5, 5, 5], ["wide", "few", "wide", "binary", "wide"], seed=23)
    subset = ["g0", "g1", "g3"]
    test_donors = {"d0", "d1"}
    _label(a, test_donors)
    slow = naive_de(a, condition_col="condition", test_level="disease", ref_level="ctrl",
                    genes=subset)
    fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=subset,
                                         method="wilcoxon").test(test_donors)
    _assert_same_result(fast, slow, "wilcoxon")
    # ...and the restricted answer really is different from the all-genes one, or the test above
    # would pass for the wrong reason.
    all_genes = NaiveRelabelEngine.from_adata(a, donor_col="donor",
                                              method="wilcoxon").test(test_donors)
    assert not np.array_equal(all_genes.table.loc[subset, "pval"].to_numpy(),
                              fast.table["pval"].to_numpy())


def test_engine_refuses_labelings_the_slow_path_refuses():
    """An optimisation must not answer where the statistic is undefined."""
    a = _make_adata([4, 4, 4], ["wide", "binary"], seed=29)
    eng = NaiveRelabelEngine.from_adata(a, donor_col="donor", method="wilcoxon")
    with pytest.raises(ValueError, match="must split the donors"):
        eng.test(set())
    with pytest.raises(ValueError, match="must split the donors"):
        eng.test({"d0", "d1", "d2"})

    thin = _make_adata([1, 6, 6], ["wide", "binary"], seed=31)
    thin_eng = NaiveRelabelEngine.from_adata(thin, donor_col="donor", method="wilcoxon")
    with pytest.raises(ValueError, match="fewer than 2 cells"):
        thin_eng.test({"d0"})

    with pytest.raises(ValueError, match="method must be one of"):
        NaiveRelabelEngine.from_adata(a, donor_col="donor", method="t-test")


@pytest.mark.slow
def test_run_null_agrees_between_the_two_engines():
    """The end-to-end claim: swapping the engine changes the cost and nothing else.

    This is the form of the claim that matters, because ``run_null``'s outputs are what the gate
    and the sweep read. Everything downstream of the naive p-values -- both BH conventions, the
    Monte-Carlo block, the floor -- has to come out identical, not merely close.
    """
    from oracles import null_oracle

    from pbcheck.gene_universe import frozen_universe
    from pbcheck.methods.pseudobulk import build_pseudobulk
    from pbcheck.permutation import run_null

    o = null_oracle(n_donors_per_group=4, n_cells_per_donor=50, n_genes=150, seed=41)
    uni = frozen_universe(build_pseudobulk(o.adata), min_size=50)
    kw = dict(n_perm=7, n_perm_pb=4, fdr=0.05, seed=2)
    fast = run_null(o.adata, uni, naive_engine="fast", **kw)
    slow = run_null(o.adata, uni, naive_engine="scanpy", **kw)

    assert fast["naive_engine"] == "fast" and slow["naive_engine"] == "scanpy"
    assert fast["naive_engine_setup"]["n_genes"] == len(uni)
    assert slow["naive_engine_setup"] is None

    for key in ("naive_pvals", "naive_pvals_real", "pb_pvals", "pb_pvals_real"):
        assert np.array_equal(fast[key], slow[key], equal_nan=True), key
    for key in ("naive_ndeg_paired", "naive_ndeg_solo", "pb_ndeg"):
        assert fast[key].bh_mode == slow[key].bh_mode
        assert np.array_equal(fast[key].counts, slow[key].counts), key
    assert fast["monte_carlo"] == slow["monte_carlo"]
    assert fast["paired_bh"] == slow["paired_bh"]

    with pytest.raises(ValueError, match="naive_engine must be"):
        run_null(o.adata, uni, naive_engine="turbo", **kw)
