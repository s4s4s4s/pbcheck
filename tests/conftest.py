import sys
from pathlib import Path

import pytest

# Make the synthetic-oracle generators importable in tests (they live outside the package on purpose:
# they are the correctness spec, not shipped runtime code).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "synthetic"))

from oracles import null_oracle  # noqa: E402

#: Shape shared by the tests refactored onto ``small_null_adata`` below (see
#: tests/test_methods.py): small enough to stay in milliseconds, with donor structure but no
#: true DE, and >= 3 donors/group so it clears the design auditor's own inclusion threshold.
_SMALL_SHAPE = dict(n_donors_per_group=3, n_cells_per_donor=40, n_genes=80)


@pytest.fixture(scope="session")
def _small_null_oracle():
    """The one expensive ``simulate()`` draw behind :func:`small_null_adata`, built once per
    session. Tests must not mutate this directly — use ``small_null_adata``, which copies it."""
    return null_oracle(seed=99, **_SMALL_SHAPE)


@pytest.fixture
def small_null_adata(_small_null_oracle):
    """A small, fixed-seed ``AnnData`` with donor structure and no true DE (``synthetic.oracles``,
    the correctness spec — see the module docstring there).

    Function-scoped so each test gets its own object and may safely mutate ``.obs`` in place
    (add a column, relabel, etc.) without leaking into other tests, while the underlying
    ``simulate()`` call — the actual cost — runs only once per test session via the
    session-scoped ``_small_null_oracle`` it copies.

    For a scenario that needs specific sizing, true DE, or ``donor_sigma`` this fixture does not
    provide, call ``oracles.null_oracle`` / ``positive_oracle`` / ``no_donor_effect_oracle``
    directly as the existing tests do; this fixture is scaffolding for tests that only need *some*
    small donor-structured dataset, not a replacement for every oracle call in the suite.
    """
    return _small_null_oracle.adata.copy()


@pytest.fixture
def donor_design_df(small_null_adata):
    """One row per donor: ``donor`` id and its ``condition``, derived from ``small_null_adata``.

    The shape several design-auditor tests build by hand when they need to reason at the donor
    level rather than the cell level (see the module docstring of ``pbcheck.design`` on why that
    distinction is the whole point of the auditor).
    """
    return (
        small_null_adata.obs[["donor", "condition"]]
        .astype(str)
        .drop_duplicates()
        .reset_index(drop=True)
    )
