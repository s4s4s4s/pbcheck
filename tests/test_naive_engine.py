"""The fast naive engine is scanpy's statistic, not an approximation of it — held to that here.

:mod:`pbcheck.methods.naive_engine` replaces one ``rank_genes_groups`` run per donor-permutation
with one ranking pass per stratum. That is only admissible if it reproduces the pre-registered
statistic (spec §2, ``method='wilcoxon'`` with ``tie_correct=True``) exactly, so these tests
compare it against :func:`pbcheck.methods.naive.naive_de` — the slow path itself, still callable —
rather than against the formula the fast path implements.

The tolerance is derived from the arithmetic, not chosen to make the suite pass. The two methods
sit in different numerical regimes and are asserted differently for that reason.

**Wilcoxon — bitwise, everywhere.** Ranks are exact multiples of 0.5 and the largest rank sum a
frozen stratum can produce is ~1.5e11, so every partial sum is exactly representable, float64
addition is associative over them, and there is nothing to round. These tests compare raw bits; an
``rtol`` here would hide a real divergence. This is the pre-registered primary test (spec §2).

**t-test_overestim_var — its p-value is not a machine-independent quantity, and nothing here
pretends otherwise.** scanpy takes the group variance from :func:`fast_array_utils.stats.mean_var`,
which evaluates ``E[x^2] - E[x]^2`` and squares **in the stored dtype** — float32 for a normalised
count matrix. With ``kappa = E[x^2] / var`` the condition number of that subtraction,

    relative error of the variance  ~  (eps32/2 + n * eps64) * kappa,    eps32/2 = 5.96e-08

so the float32 squaring dominates the float64 accumulation by eight orders of magnitude. On
log1p-normalised data kappa is large by construction: the values are logs of a normalised count, so
they cluster tightly about a mean of order 9 while the within-group variance can be 1e-13. Measured
on these generators kappa reaches 1e6 routinely and 1.9e14 once a group's cells collapse onto
adjacent float32 values (a two-gene cell saturates at log1p(1e4) = 9.21). At kappa = 1.9e14 the
computed variance is 4.79e-06 against a true 4.55e-13 — seven orders of magnitude wrong, **zero
significant digits left** — and two equally valid evaluations of the identical expression then
differ by percent in the p-value. That was measured with scanpy on both sides; it is a property of
the pre-registered formula, not of this optimisation.

The t-test is therefore asserted where the claim is well conditioned and the assertion can be tight:

1. :func:`test_group_sufficient_statistics_are_correctly_rounded` — the engine's per-group ``n``,
   ``sum x`` and ``sum x**2`` against :func:`math.fsum`, which is exactly rounded. Sums of
   like-signed terms, perfectly conditioned; this is where a real reduction bug would live.
2. :func:`test_ttest_formula_is_scanpys_source_bitwise` — the engine's p-values against a
   transcription of scanpy's ``_RankGenes.t_test`` fed identical statistics, at **zero** tolerance
   and with **unequal group sizes**, which is what makes the ``nobs2 = ns_group`` overestim
   substitution observable at all.
3. Only then the end-to-end comparison against scanpy, bounded by the rounding envelope the
   formula's own arithmetic admits (:func:`_ttest_rounding_envelope`). It collapses towards 1e-15
   where kappa is small and widens only where the formula has run out of digits. It is deliberately
   not the only t-test assertion: on balanced groups an envelope cannot see a wrong ``nobs2`` —
   measured, 0 of 1077 injected ``nobs2`` bugs caught by the envelope, all of them caught by (2).

Also pinned: that the two engines agree end-to-end through :func:`pbcheck.permutation.run_null`,
that sparse and dense inputs give the same answer, and that the gene chunking — the thing that
keeps a 547665-cell stratum off the heap — cannot move a number.

**Domain.** The property generators draw only strata the spec would admit: at least
``pbcheck.design.audit_design``'s ``min_donors`` = 3 donors per group (spec §1) and at least
``pbcheck.gate_config.MIN_CELLS`` = 10 cells per donor (spec §1 item 2, §3). That is the domain the
sweep runs in, and importing the two constants keeps generator and gate in step.

That restriction is **not** what makes these tests machine-independent and must not be read as the
fix: kappa does not improve with n — measured worst case 1.9e14 at 2 cells per group and still
9.9e13 at 120 — so a wider domain would only change how often a badly conditioned gene is drawn,
never whether one can be. The conditioning is handled by asserting on the right layer, above.
:func:`test_one_pass_variance_has_no_digits_on_saturated_genes` keeps the degenerate regime pinned
rather than merely excluded from the generator.
"""

from __future__ import annotations

import inspect
import math
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from hypothesis import given, settings, strategies as st
from scipy import sparse
from scipy.stats import ttest_ind_from_stats

from pbcheck.design import audit_design
from pbcheck.gate_config import MIN_CELLS
from pbcheck.methods.naive import naive_de
from pbcheck.methods.naive_engine import (
    DEFAULT_CHUNK_BYTES,
    NaiveRelabelEngine,
    _average_ranks,
)

warnings.filterwarnings("ignore")

#: ``derandomize=True`` is load-bearing, not tidiness. Without it hypothesis seeds itself from
#: entropy on every run, so the set of examples explored differs between this machine and each CI
#: runner, and a green local run is not evidence that CI will be green — which is exactly how the
#: t-test conditioning below reached ``main`` unnoticed: 7200 local configurations found a maximum
#: engine-vs-scanpy difference of 0.0, and CI's first run found 2.7e-03. Derandomised, every
#: machine explores the identical sequence, so local green and CI green mean the same thing. The
#: cost is that repeated runs no longer discover new inputs over time; the example budget is raised
#: to compensate, and the search is cheap here because the generated strata are tiny.
_SETTINGS = settings(max_examples=120, deadline=None, derandomize=True)

#: Spec §1's inclusion gate, read from the code that enforces it rather than restated as literals,
#: so that a change to either constant moves the generators with it.
_MIN_DONORS_PER_GROUP = inspect.signature(audit_design).parameters["min_donors"].default
_MIN_CELLS_PER_DONOR = MIN_CELLS

_EPS64 = float(np.finfo(np.float64).eps)
_EPS32 = float(np.finfo(np.float32).eps)

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


def _normalised(a):
    """The float32 matrix both paths actually rank / sum, and the test group mask."""
    b = a.copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    values = b.X.toarray() if sparse.issparse(b.X) else np.asarray(b.X)
    return values, np.asarray(a.obs["condition"] == "disease")


def _ttest_rounding_envelope(values, mask, n_a, n_b):
    """Per gene, the ``[p_lo, p_hi]`` scanpy's own variance formula can land in.

    The centre is the best available float64 evaluation of the formula: ``sum x`` and ``sum x**2``
    taken with :func:`math.fsum` (exactly rounded) over squares formed in **float64**. The half
    width is

        delta = (eps32/2 + n * eps64) * E[x**2] * n/(n-1)

    with two distinct jobs, and the first is the one that matters:

    * ``eps32/2`` covers **which dtype the squares are formed in**, which is not a fixed property
      of the statistic. ``fast_array_utils.power(x, 2)`` is called with no ``dtype``, so it squares
      in the stored dtype — float32 for a normalised count matrix — but whether a given platform's
      stack actually presents the matrix as float32 at that point is not something either
      implementation controls. Evidence that this is real and not hypothetical: on the machine this
      was developed on, engine and scanpy agree bitwise on 7200 configurations, while CI's scanpy
      returned the value exact rational arithmetic gives (0.338437873214868 where this machine
      gives 0.337523825174102). Both are legitimate evaluations; the envelope has to contain both.
    * ``n * eps64`` covers the float64 accumulation order — a per-donor reduction here, one pass
      over the group in scanpy. This term is eight orders of magnitude smaller.

    Both are then multiplied by the condition number of ``E[x**2] - E[x]**2`` implicitly: the
    subtraction is what turns an absolute error in ``E[x**2]`` into a relative error in the
    variance, so the interval self-widens exactly where the formula has lost digits and collapses
    towards 1e-15 where it has not.

    ``nobs2`` is ``n_a``: ``t-test_overestim_var`` substitutes the *test* group's size for the
    reference group's, so the envelope makes the same substitution or it brackets a different
    statistic.

    Assumes p is extremal at the interval ends in each variance. It is monotone in each through
    ``|t|``; the Welch–Satterthwaite ``df`` also moves, but its effect on p is orders of magnitude
    below the interval width, and the unperturbed value is evaluated too so the centre is covered.
    """
    lo = np.empty(values.shape[1])
    hi = np.empty(values.shape[1])
    for j in range(values.shape[1]):
        stats = []
        for sel, n in ((mask, n_a), (~mask, n_b)):
            x = values[sel, j].astype(np.float64)
            mean = math.fsum(x) / n
            e2 = math.fsum(x * x) / n
            var = (e2 - mean ** 2) * (n / (n - 1))
            delta = (_EPS32 / 2 + n * _EPS64) * e2 * (n / (n - 1))
            stats.append((mean, var, delta))
        (m_a, v_a, d_a), (m_b, v_b, d_b) = stats
        corners = []
        for va in (max(v_a - d_a, 0.0), v_a, v_a + d_a):
            for vb in (max(v_b - d_b, 0.0), v_b, v_b + d_b):
                with np.errstate(invalid="ignore", divide="ignore"):
                    _, p = ttest_ind_from_stats(m_a, np.sqrt(va), n_a, m_b, np.sqrt(vb), n_a,
                                                equal_var=False)
                corners.append(1.0 if np.isnan(p) else float(p))
        lo[j], hi[j] = min(corners), max(corners)
    return lo, hi


def _assert_same_result(fast, slow, method, adata=None):
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
        # Nothing may change side of any plausible threshold. Only asserted for the arm where the
        # p-value is machine-independent; see the module docstring for why the t-test's is not.
        for alpha in (0.05, 0.01, 1e-4):
            assert np.array_equal(p_fast < alpha, p_slow < alpha)
    else:
        assert adata is not None, "the t-test comparison needs the data to derive its bound"
        values, mask = _normalised(adata)
        n_a = int(mask.sum())
        lo, hi = _ttest_rounding_envelope(values, mask, n_a, int(values.shape[0] - n_a))
        width = hi - lo
        for name, p in (("engine", p_fast), ("scanpy", p_slow)):
            outside = np.where((p < lo - 1e-15) | (p > hi + 1e-15))[0]
            assert outside.size == 0, (
                f"{name} left the rounding envelope of scanpy's own variance formula at genes "
                f"{outside.tolist()}: p={p[outside]} vs [{lo[outside]}, {hi[outside]}] -- that is "
                "a difference the formula's arithmetic does not account for, i.e. a real bug"
            )
        gap = np.abs(p_fast - p_slow)
        assert np.all(gap <= width + 1e-15), (
            f"engine and scanpy differ by more than the formula's own rounding freedom: "
            f"max gap {gap.max():.3e} against envelope width {width[gap.argmax()]:.3e}"
        )

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
# The engine against the slow path, over the domain spec §1 admits.
# ---------------------------------------------------------------------------

#: Donors per group and cells per donor, drawn only where the spec would admit the stratum.
#: ``max`` values are kept small so 120 examples stay in milliseconds, not because larger is
#: untested -- the fixed-example tests below cover a 200-cell donor and a 6-donor stratum.
_ADMISSIBLE = dict(
    donors_per_group=st.integers(min_value=_MIN_DONORS_PER_GROUP,
                                 max_value=_MIN_DONORS_PER_GROUP + 2),
    cells=st.integers(min_value=_MIN_CELLS_PER_DONOR, max_value=_MIN_CELLS_PER_DONOR + 8),
    kinds=st.lists(st.sampled_from(_GENE_KINDS), min_size=1, max_size=5),
    seed=st.integers(min_value=0, max_value=9999),
)


def _admissible_stratum(donors_per_group, cells, seed):
    """Cells per donor for a stratum inside the gate, with the per-donor counts deliberately
    unequal (real strata never have equal donors, and an equal-``n`` shape hides a whole class of
    reduction bug)."""
    rng = np.random.default_rng(seed)
    n_donors = 2 * donors_per_group
    per_donor = cells + rng.integers(0, 5, size=n_donors)
    return [int(c) for c in per_donor], set(f"d{i}" for i in range(donors_per_group))


@given(**_ADMISSIBLE)
@_SETTINGS
def test_wilcoxon_matches_the_slow_path_bitwise(donors_per_group, cells, kinds, seed):
    """Across donor counts, group sizes, tie structures and unequal cells per donor."""
    cells_per_donor, test_donors = _admissible_stratum(donors_per_group, cells, seed)
    a = _make_adata(cells_per_donor, kinds, seed)
    slow = _reference(a, test_donors, "wilcoxon")
    fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                         method="wilcoxon").test(test_donors)
    _assert_same_result(fast, slow, "wilcoxon", a)


@given(**_ADMISSIBLE)
@_SETTINGS
def test_ttest_overestim_var_matches_the_slow_path(donors_per_group, cells, kinds, seed):
    """Spec §2's robustness variant, reconstructed from n_d, sum x and sum x^2.

    Admissible only because scanpy's variance comes from ``E[x^2] - E[x]^2``
    (``fast_array_utils.stats.mean_var``, both the dense and the sparse kernel), which the per-donor
    sufficient statistics reproduce. A two-pass variance would be better conditioned and would be a
    *different* number from the pre-registered one, so it is not used.

    Bounded by the formula's own rounding freedom rather than by a fixed ``rtol`` -- see the module
    docstring, and the two tests below that carry the tight half of the claim.
    """
    cells_per_donor, test_donors = _admissible_stratum(donors_per_group, cells, seed)
    a = _make_adata(cells_per_donor, kinds, seed)
    slow = _reference(a, test_donors, "t-test_overestim_var")
    fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                         method="t-test_overestim_var").test(test_donors)
    _assert_same_result(fast, slow, "t-test_overestim_var", a)


@given(**_ADMISSIBLE)
@_SETTINGS
def test_group_sufficient_statistics_are_correctly_rounded(donors_per_group, cells, kinds, seed):
    """The tight half of the t-test claim: the engine's ``n``, ``sum x``, ``sum x**2`` are right.

    Checked against :func:`math.fsum`, which is exactly rounded, so this is the correctness of the
    reduction itself and not a comparison of two roundings. These are sums of like-signed terms and
    are perfectly conditioned -- a genuine bug in the per-donor segmentation, the donor ordering or
    the float32-then-float64 squaring shows up here at once, where the p-value would smear it
    through a subtraction that can have no digits left.
    """
    cells_per_donor, test_donors = _admissible_stratum(donors_per_group, cells, seed)
    a = _make_adata(cells_per_donor, kinds, seed)
    _label(a, test_donors)
    eng = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                        method="t-test_overestim_var")
    values, mask = _normalised(a)

    for sel, want_mask in ((mask, eng.donors.isin(list(test_donors))),
                           (~mask, ~eng.donors.isin(list(test_donors)))):
        n = int(sel.sum())
        assert int(eng.cells_per_donor[want_mask].sum()) == n
        got_sum = eng.donor_sum[want_mask].sum(axis=0)
        got_sq = eng.donor_sumsq[want_mask].sum(axis=0)
        for j in range(values.shape[1]):
            col = values[sel, j]
            exact_sum = math.fsum(col.astype(np.float64))
            # scanpy squares in the STORED dtype before accumulating; the reference has to do the
            # same or it would be measuring a different quantity.
            exact_sq = math.fsum((col * col).astype(np.float64))
            for label, got, exact in (("sum x", got_sum[j], exact_sum),
                                      ("sum x**2", got_sq[j], exact_sq)):
                bound = n * _EPS64 * max(abs(exact), 1e-300)
                assert abs(got - exact) <= bound, (
                    f"{label} gene {j}: {got!r} vs exactly-rounded {exact!r}, "
                    f"error {abs(got - exact):.3e} > n*eps bound {bound:.3e}"
                )


def _scanpy_ttest_from_stats(mean_a, var_a, n_a, mean_b, var_b):
    """scanpy ``_RankGenes.t_test('t-test_overestim_var')``, transcribed line for line.

    From ``scanpy/tools/_rank_genes_groups.py``: ``ns_rest = ns_group`` for the overestim variant
    (the reference group's own size never reaches the test), ``equal_var=False`` for Welch, then
    ``scores[isnan] = 0`` and ``pvals[isnan] = 1``. Written out here so that changing any of those
    four decisions in the engine turns into a red test at **zero** tolerance -- an envelope over
    balanced groups cannot see a wrong ``nobs2``, and this can.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        scores, pvals = ttest_ind_from_stats(
            mean1=mean_a, std1=np.sqrt(var_a), nobs1=n_a,
            mean2=mean_b, std2=np.sqrt(var_b), nobs2=n_a,
            equal_var=False,
        )
    scores = np.where(np.isnan(scores), 0.0, scores)
    return scores, np.where(np.isnan(pvals), 1.0, pvals)


@pytest.mark.parametrize("cells_per_donor,n_test", [
    # Deliberately UNBALANCED: with n_a == n_b the overestim substitution nobs2 = n_a is invisible,
    # so a balanced-only test would pass with the reference group's size in its place.
    ([12, 13, 14, 40, 41, 42], 3),
    ([10, 11, 12, 13, 14, 60], 3),
    ([30, 31, 32, 10, 11, 12], 3),
])
def test_ttest_formula_is_scanpys_source_bitwise(cells_per_donor, n_test):
    """Given identical statistics, the engine's t-test is scanpy's, bit for bit.

    This is the assertion that carries spec fidelity for the robustness variant: ddof, Welch, the
    NaN repairs, and above all ``nobs2 = ns_group``. It is exact because both sides are handed the
    same variances -- the ill-conditioned subtraction has already happened and is common to both,
    so nothing here is at the mercy of a summation order.
    """
    test_donors = set(f"d{i}" for i in range(n_test))
    a = _make_adata(cells_per_donor, ["wide", "few", "flat", "zero"], seed=101)
    _label(a, test_donors)
    eng = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                        method="t-test_overestim_var")
    got = eng.test(test_donors)

    mask = eng.donors.isin(list(test_donors))
    n_a = int(eng.cells_per_donor[mask].sum())
    n_b = int(eng.cells_per_donor[~mask].sum())
    assert n_a != n_b, "this test is worthless on balanced groups"

    def stats(sel, n):
        mean = eng.donor_sum[sel].sum(axis=0) / n
        var = (eng.donor_sumsq[sel].sum(axis=0) / n - mean ** 2) * (n / (n - 1))
        return mean, var

    m_a, v_a = stats(mask, n_a)
    m_b, v_b = stats(~mask, n_b)
    _, expected = _scanpy_ttest_from_stats(m_a, v_a, n_a, m_b, v_b)
    assert np.array_equal(got.table["pval"].to_numpy(), expected)

    # ...and the substitution really is load-bearing on this data, so the test can fail.
    with np.errstate(invalid="ignore", divide="ignore"):
        _, wrong = ttest_ind_from_stats(m_a, np.sqrt(v_a), n_a, m_b, np.sqrt(v_b), n_b,
                                        equal_var=False)
    wrong = np.where(np.isnan(wrong), 1.0, wrong)
    assert not np.allclose(expected, wrong, rtol=1e-6), (
        "using the reference group's nobs would give the same answer here, so this configuration "
        "cannot detect the overestim substitution -- pick a more unbalanced one"
    )


def test_one_pass_variance_has_no_digits_on_saturated_genes():
    """Pin the ill-conditioning itself, rather than only excluding it from the generators.

    Two cells whose normalised values land one float32 ulp apart give
    ``kappa = E[x**2]/var ~ 1.9e14``. scanpy's ``E[x**2] - E[x]**2`` then returns a variance seven
    orders of magnitude too large, so the pre-registered t-test's p-value on such a gene is
    arithmetic noise -- for scanpy, for the engine, and for any other implementation of the same
    expression. Asserted here so that the next person to see a percent-level cross-machine
    disagreement on this arm finds the explanation in a test instead of rediscovering it in CI.

    The Wilcoxon arm, which is the pre-registered primary, is immune: it never subtracts two nearly
    equal large numbers.
    """
    base = np.float32(9.2103405)
    x32 = np.array([base, np.nextafter(base, np.float32(1e9))], dtype=np.float32)
    x = x32.astype(np.float64)

    var_two_pass = float(((x - x.mean()) ** 2).sum() / (len(x) - 1))
    e2 = float((x32 * x32).astype(np.float64).sum() / len(x))
    var_one_pass = (e2 - (x.sum() / len(x)) ** 2) * (len(x) / (len(x) - 1))
    kappa = e2 / var_two_pass

    assert kappa > 1e13, kappa
    assert _EPS32 / 2 * kappa > 1, "float32 squaring must be the dominant error term here"
    assert var_one_pass / var_two_pass > 1e6, (
        f"the one-pass formula should be orders of magnitude off here: "
        f"{var_one_pass:.3e} vs {var_two_pass:.3e}"
    )
    # And the float64 accumulation is NOT the culprit -- squaring in float64 recovers the answer.
    var_f64_squares = (float((x * x).sum() / len(x)) - (x.sum() / len(x)) ** 2) * 2
    assert abs(var_f64_squares - var_two_pass) / var_two_pass < 1e-3


def test_the_generators_stay_inside_the_inclusion_gate():
    """The admissible domain is read from the gate, not restated, so the two cannot drift apart."""
    assert _MIN_DONORS_PER_GROUP == 3 and _MIN_CELLS_PER_DONOR == 10
    cells_per_donor, test_donors = _admissible_stratum(_MIN_DONORS_PER_GROUP,
                                                       _MIN_CELLS_PER_DONOR, seed=0)
    assert len(test_donors) >= _MIN_DONORS_PER_GROUP
    assert len(cells_per_donor) - len(test_donors) >= _MIN_DONORS_PER_GROUP
    assert min(cells_per_donor) >= _MIN_CELLS_PER_DONOR


@pytest.mark.parametrize("method", ["wilcoxon", "t-test_overestim_var"])
def test_at_the_inclusion_gates_own_boundary(method):
    """Exactly 3 donors per group and exactly 10 cells per donor — the smallest admitted stratum.

    Fixed, not generated: the boundary of the domain is where an off-by-one in the generator would
    stop covering the sweep's own worst case, so it is pinned independently of the strategy.
    """
    cells_per_donor = [_MIN_CELLS_PER_DONOR] * (2 * _MIN_DONORS_PER_GROUP)
    test_donors = set(f"d{i}" for i in range(_MIN_DONORS_PER_GROUP))
    a = _make_adata(cells_per_donor, ["zero", "flat", "binary", "few", "wide"], seed=7)
    slow = _reference(a, test_donors, method)
    fast = NaiveRelabelEngine.from_adata(a, donor_col="donor", genes=list(a.var_names),
                                         method=method).test(test_donors)
    _assert_same_result(fast, slow, method, a)


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
    _assert_same_result(fast, slow, method, a)
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
        _assert_same_result(fast, slow, method, a)


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
        _assert_same_result(fast, slow, method, a)

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
    _assert_same_result(fast, slow, "wilcoxon", a)


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
