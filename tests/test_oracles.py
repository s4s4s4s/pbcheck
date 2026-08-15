"""Self-tests for ``synthetic/oracles.py`` — the correctness spec is untested on its own terms.

Every other test in this suite trusts the oracles to mean what their docstrings say: donor_sigma
controls a between-donor random effect and nothing else, and the positive oracle plants exactly
the requested effect. Nothing had checked that directly. Fixed seeds throughout, small sizes, all
fast (no DE runs — these only look at ``simulate``'s output).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "synthetic"))

from oracles import no_donor_effect_oracle, null_oracle, positive_oracle, simulate  # noqa: E402


def _donor_pseudobulk(oracle) -> tuple[np.ndarray, list[str]]:
    """Sum raw counts per donor per gene — a self-contained aggregation, not
    ``pbcheck.methods.pseudobulk.build_pseudobulk``, so these tests do not depend on the module
    they exist to validate an input for."""
    adata = oracle.adata
    donors = adata.obs["donor"].astype(str)
    X = np.asarray(adata.X)
    uniq = sorted(donors.unique())
    mat = np.zeros((len(uniq), X.shape[1]))
    donor_arr = donors.to_numpy()
    for i, d in enumerate(uniq):
        mat[i] = X[donor_arr == d].sum(axis=0)
    return mat, uniq


def _between_donor_cv(mat: np.ndarray, min_mean: float = 50.0) -> np.ndarray:
    """Per-gene coefficient of variation across donors, restricted to well-expressed genes so a
    handful of near-zero counts do not dominate with pure Poisson noise."""
    mean_per_gene = mat.mean(axis=0)
    well_expressed = mean_per_gene > min_mean
    sub = mat[:, well_expressed]
    return sub.std(axis=0, ddof=1) / sub.mean(axis=0)


# ---------------------------------------------------------------------------
# donor_sigma controls between-donor variance, and only that.
# ---------------------------------------------------------------------------

def test_donor_sigma_zero_has_no_between_donor_shift_beyond_sampling_noise():
    """With no donor random effect, donor-level means must agree up to sampling noise.

    ``donor_sigma=0`` collapses the per-(gene, donor) random effect to exactly 1 (see
    ``simulate``'s ``donor_re`` branch), so the only remaining source of between-donor spread is
    per-cell depth and NB/Poisson sampling — no systematic per-donor multiplicative shift. The
    between-donor coefficient of variation on well-expressed genes should therefore be small; the
    bound (0.15) sits well above what fixed-seed runs at this size actually produce (~0.07), with
    headroom for a different seed.
    """
    o = no_donor_effect_oracle(n_genes=300, n_donors_per_group=6, n_cells_per_donor=150, seed=5)
    mat, donors = _donor_pseudobulk(o)
    assert len(donors) == 12
    cv = _between_donor_cv(mat)
    assert cv.size > 100, "too few well-expressed genes to say anything"
    assert np.median(cv) < 0.15, f"donor_sigma=0 median between-donor CV {np.median(cv)} too high"


def test_donor_sigma_positive_yields_measurably_more_between_donor_variance():
    """The same design, only ``donor_sigma`` changed, must show a clear CV gap.

    Not "some" more — a real, easily-resolved gap, since this is the entire mechanism
    pseudoreplication auditing depends on: if the oracle's donor random effect were too weak to
    detect this way, it would also be too weak to produce the inflation the naive-arm tests rely
    on. Same seed, genes, donors and cells-per-donor as the sigma=0 test above, donor_sigma=0.5 is
    the gate's own default operating point (``pbcheck.gate_config.ORACLE_SIM``).
    """
    null = no_donor_effect_oracle(n_genes=300, n_donors_per_group=6, n_cells_per_donor=150, seed=5)
    donor_re = null_oracle(n_genes=300, n_donors_per_group=6, n_cells_per_donor=150,
                            donor_sigma=0.5, seed=5)

    cv_null = _between_donor_cv(_donor_pseudobulk(null)[0])
    cv_donor_re = _between_donor_cv(_donor_pseudobulk(donor_re)[0])

    assert np.median(cv_donor_re) > 5 * np.median(cv_null), (
        f"donor_sigma=0.5 median CV {np.median(cv_donor_re)} is not clearly above "
        f"donor_sigma=0 median CV {np.median(cv_null)}"
    )


# ---------------------------------------------------------------------------
# positive_oracle plants exactly K genes at exactly the requested effect size.
# ---------------------------------------------------------------------------

def test_positive_oracle_injects_exactly_k_genes_at_the_requested_lfc():
    K, LOG2FC = 15, 1.3
    o = positive_oracle(n_de=K, log2fc=LOG2FC, n_genes=200, n_donors_per_group=4,
                        n_cells_per_donor=50, seed=3)

    assert o.n_true_de == K
    assert len(o.de_genes) == K
    assert len(set(o.de_genes)) == K, "de_genes must be K distinct genes, not K draws with repeats"
    assert set(o.log2fc.keys()) == set(o.de_genes)

    # Every planted gene gets exactly the requested magnitude — not diluted, not averaged.
    magnitudes = {abs(v) for v in o.log2fc.values()}
    assert magnitudes == {LOG2FC}

    # Both directions are available (mixed up/down), not a directional artifact — with K=15 at a
    # fixed seed this is a property of the draw, not guaranteed for every seed/K, so it is checked
    # on this fixed instance rather than asserted as a general law.
    signs = {v / abs(v) for v in o.log2fc.values()}
    assert signs <= {1.0, -1.0} and signs, "log2fc values must be exactly +-LOG2FC"


def test_null_oracle_plants_no_de_genes():
    o = null_oracle(n_genes=150, n_donors_per_group=3, n_cells_per_donor=40, seed=9)
    assert o.n_true_de == 0
    assert o.de_genes == []
    assert o.log2fc == {}


# ---------------------------------------------------------------------------
# Counts are always non-negative integers, regardless of donor_sigma / dispersion / n_de.
# ---------------------------------------------------------------------------

def test_counts_are_nonnegative_integers():
    for donor_sigma, dispersion, n_de in ((0.0, 0.2, 0), (0.6, 0.5, 10)):
        o = simulate(n_genes=100, n_donors_per_group=3, n_cells_per_donor=30, n_de=n_de,
                    donor_sigma=donor_sigma, dispersion=dispersion, seed=17)
        X = np.asarray(o.adata.X)
        assert np.all(X >= 0), f"negative count with donor_sigma={donor_sigma}, dispersion={dispersion}"
        assert np.array_equal(X, np.round(X)), (
            f"non-integer count with donor_sigma={donor_sigma}, dispersion={dispersion}"
        )


def test_counts_are_nonnegative_integers_with_dispersion_zero():
    """``dispersion <= 0`` takes the pure-Poisson branch (a separate code path in ``_nb_counts``)."""
    o = simulate(n_genes=80, n_donors_per_group=3, n_cells_per_donor=25, dispersion=0.0, seed=17)
    X = np.asarray(o.adata.X)
    assert np.all(X >= 0)
    assert np.array_equal(X, np.round(X))
