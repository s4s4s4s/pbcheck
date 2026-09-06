import subprocess
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


#: Shape of :func:`audit_shape_oracle`: the smallest oracle draw that clears all three floors the
#: audit gates on at once. Measured on this shape: ``build_pseudobulk`` keeps 8 of 8 donor profiles
#: (4 per group, above ``audit.PRODUCT_MIN_PROFILES_PER_GROUP``), ``frozen_universe`` keeps all 600
#: genes (well above ``gate_config.MIN_UNIVERSE_SIZE`` = 200), and 4 donors per group clears
#: ``audit.PRODUCT_MIN_DONORS_PER_GROUP``. ``tests/test_audit.py`` re-measures all three so the
#: numbers in this comment cannot rot silently.
AUDIT_ORACLE_SHAPE = dict(n_genes=600, n_donors_per_group=4, n_cells_per_donor=60)


@pytest.fixture(scope="session")
def audit_shape_oracle():
    """The one oracle draw the audit tests share: no true DE, donor effect present, 480 cells.

    Session-scoped because ``simulate()`` at this shape is the expensive part of the audit tests
    and the draw is immutable; every test that needs to touch the data copies ``.adata`` first,
    as ``small_null_adata`` does for the smaller shape above.
    """
    return null_oracle(seed=7, **AUDIT_ORACLE_SHAPE)


#: The branch the release checks diff against. A CI checkout of a branch, a pull request or a
#: tag holds it only as the remote-tracking ref (tests.yml and release.yml fetch the full history
#: so that ref exists); a checkout with neither is too shallow to diff, which is an error in the
#: checkout, not a reason to skip a check.
_BASE_REF_CANDIDATES = ("main", "origin/main")


@pytest.fixture(scope="session")
def base_ref() -> str:
    """``main`` where the checkout has it, else ``origin/main``; fails if neither resolves."""
    root = Path(__file__).resolve().parents[1]
    for ref in _BASE_REF_CANDIDATES:
        cp = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--verify", "--quiet", ref],
            capture_output=True,
            text=True,
        )
        if cp.returncode == 0:
            return ref
    raise AssertionError(
        "neither main nor origin/main resolves in this checkout; the release checks need the base "
        "branch (tests.yml and release.yml check out with fetch-depth 0)"
    )
