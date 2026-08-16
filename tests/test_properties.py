"""Property-based tests (hypothesis) — the two invariants worth generalizing beyond fixed examples.

Both are pinned by example-based regression tests elsewhere (``test_permutation_mtc.py``); these
generalize across a range of inputs rather than a handful of hand-picked ones, which is exactly
where the historical bugs in this code lived (``build_perms`` used to hang on a *specific* 5v5
input; a monotonicity break in BH would show up on *some* p-value array, not necessarily the ones
anyone thought to write down). ``max_examples`` is kept modest so this stays fast enough to run on
every commit rather than being relegated to ``slow``.
"""

from __future__ import annotations

from math import comb

import numpy as np
from hypothesis import given, settings, strategies as st

from pbcheck.permutation import build_perms
from pbcheck import mtc

#: ``derandomize=True`` for the same reason ``tests/test_naive_engine.py`` sets it: without it
#: hypothesis reseeds from entropy every run, so this machine and each CI runner explore different
#: examples and a green local run says nothing about CI. That is not hypothetical here — it is how
#: a numerical disagreement in the naive engine passed locally across 7200 configurations and failed
#: on CI's first attempt. Derandomised, every machine walks the identical sequence, so local green
#: and CI green mean the same thing. The cost is that repeated runs no longer turn up new inputs
#: over time; the example budget is doubled to buy some of that back in a single pass, which both
#: of these cheap properties can afford.
_SETTINGS = settings(max_examples=100, deadline=None, derandomize=True)


# ---------------------------------------------------------------------------
# build_perms invariants (spec §4: the donor-permutation null's construction).
# ---------------------------------------------------------------------------

@given(
    n_per_group=st.integers(min_value=3, max_value=8),
    n_perm_request=st.integers(min_value=5, max_value=60),
    seed=st.integers(min_value=0, max_value=10_000),
)
@_SETTINGS
def test_build_perms_invariants(n_per_group, n_perm_request, seed):
    """Across varying donor/group sizes: right group sizes, no identity/complement, right count.

    ``build_perms`` documents its own count contract (enumerate outright when the total
    combination count is small, or when fewer distinct non-trivial label-sets exist than
    requested; otherwise reject-sample ``n_perm`` distinct sets) — that contract is what is
    checked here, computed independently from the donor count rather than copied from the
    function's internals, so a change to the *behaviour* (not just its internal bookkeeping)
    would be caught. Termination itself is implicit: this test hanging is exactly the historical
    bug (a 5v5 design under the old code) that the two regression tests in
    test_permutation_mtc.py pin a single instance of.
    """
    donors = [f"d{i:03d}" for i in range(2 * n_per_group)]
    true_test = set(donors[:n_per_group])
    complement = set(donors) - true_test

    perms = build_perms(donors, true_test, n_perm=n_perm_request, seed=seed)

    # Group size is fixed by construction, always.
    assert all(len(p) == n_per_group for p in perms)
    # Never the real labeling or its exact complement.
    assert frozenset(true_test) not in perms
    assert frozenset(complement) not in perms
    # Every returned labeling is distinct.
    assert len(set(perms)) == len(perms)

    total = comb(len(donors), n_per_group)
    available = total - (2 if complement != true_test else 1)
    max_enumerate = 200  # build_perms' own default
    expected = available if (total <= max_enumerate or available <= n_perm_request) else n_perm_request
    assert len(perms) == expected


# ---------------------------------------------------------------------------
# BH monotonicity (spec §5: the shared multiple-testing correction).
# ---------------------------------------------------------------------------

@given(
    pvals=st.lists(
        st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
        min_size=1, max_size=300,
    ),
    a=st.floats(min_value=1e-4, max_value=0.999, allow_nan=False),
    b=st.floats(min_value=1e-4, max_value=0.999, allow_nan=False),
)
@_SETTINGS
def test_bh_rejections_monotone_in_alpha(pvals, a, b):
    """More hypotheses are rejected at a looser (larger) alpha, never fewer, on any p-value array.

    A real invariant of the Benjamini-Hochberg step-up procedure, not a tautology about
    threshold comparison: a bug that swapped in ``multipletests(...)[0]`` (the alpha-dependent
    reject boolean) somewhere BH's *adjusted p-values* belong, or that flipped a comparison
    direction, would break this on some generated array even though it might pass every
    hand-written example.
    """
    alpha1, alpha2 = sorted((a, b))
    p = np.asarray(pvals, dtype=float)

    n_rejected_1 = int((mtc.bh(p, alpha=alpha1) < alpha1).sum())
    n_rejected_2 = int((mtc.bh(p, alpha=alpha2) < alpha2).sum())

    assert n_rejected_1 <= n_rejected_2
