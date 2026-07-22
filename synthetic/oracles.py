"""Synthetic single-cell oracles with known ground truth.

These generators are the *correctness spec* for the pbcheck engine. A pseudoreplication
auditor is only trustworthy if it behaves correctly on data where we know the answer:

  * ``simulate(...)`` builds an ``AnnData`` of raw counts with a realistic **donor random
    effect** (the between-donor variability that pseudoreplication ignores) and an optional
    **condition effect** shared consistently across the donors of a group.

Two canonical oracles are built on top of it:

  * NULL      — donor random effect present, **no** true condition effect (n_de = 0).
                A correct method must call ~0 DE; a naive per-cell test will call many
                false positives. The gap is the inflation pbcheck exists to measure.
  * POSITIVE  — donor random effect present **and** ``n_de`` genes with a real,
                donor-consistent fold change. A correct method must recover those genes
                (sensitivity) without exploding elsewhere (specificity).

The single most important control is ``donor_sigma=0``: with no donor random effect the
naive and pseudobulk analyses must **agree** (inflation ~0). That falsifies the alternative
explanation that any naive-vs-pseudobulk gap is merely a power/sample-size artifact rather
than pseudoreplication. Every calibration run exercises it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import anndata as ad


@dataclass
class Oracle:
    """A simulated dataset plus its ground truth."""

    adata: ad.AnnData
    de_genes: list[str] = field(default_factory=list)  # genes with a true condition effect
    log2fc: dict[str, float] = field(default_factory=dict)  # true log2 fold change per DE gene
    params: dict = field(default_factory=dict)

    @property
    def n_true_de(self) -> int:
        return len(self.de_genes)


def _nb_counts(mean: np.ndarray, dispersion: float, rng: np.random.Generator) -> np.ndarray:
    """Negative-binomial counts via a Gamma-Poisson mixture.

    Parameterization matches DESeq2/edgeR: Var = mean + dispersion * mean**2.
    dispersion (a.k.a. alpha) -> Gamma shape = 1/dispersion.
    """
    mean = np.clip(mean, 1e-9, None)
    if dispersion <= 0:
        return rng.poisson(mean).astype(np.int64)
    shape = 1.0 / dispersion
    gamma = rng.gamma(shape=shape, scale=mean / shape)
    return rng.poisson(gamma).astype(np.int64)


def simulate(
    n_genes: int = 2000,
    n_donors_per_group: int = 4,
    n_cells_per_donor: int = 300,
    n_de: int = 0,
    *,
    donor_sigma: float = 0.35,
    log2fc: float = 1.0,
    dispersion: float = 0.2,
    mean_log_mu: float = 1.0,
    mean_log_sigma: float = 1.2,
    depth_log_sigma: float = 0.3,
    one_cell_type: bool = True,
    seed: int = 0,
) -> Oracle:
    """Simulate a two-group scRNA-seq dataset with a donor random effect.

    Model for the expected count of gene g in cell c of donor d (group k):

        E[count] = mu_g * donor_re[g, d] * cond_fc[g, k] * size_factor[c]

    where ``mu_g`` is a per-gene baseline (log-normal across genes), ``donor_re`` is a
    per-(gene, donor) log-normal random effect with sd ``donor_sigma`` (the source of
    pseudoreplication), ``cond_fc`` is 1 except for the ``n_de`` true-DE genes in the
    "disease" group, and ``size_factor`` is a per-cell depth multiplier. Counts are NB.

    Set ``donor_sigma=0`` for the falsification control (naive and pseudobulk must agree).
    """
    rng = np.random.default_rng(seed)

    n_donors = 2 * n_donors_per_group
    donors = [f"d{i:02d}" for i in range(n_donors)]
    groups = np.array(["ctrl"] * n_donors_per_group + ["disease"] * n_donors_per_group)

    genes = [f"g{j:04d}" for j in range(n_genes)]

    # Per-gene baseline mean (log-normal across genes).
    mu_g = np.exp(rng.normal(mean_log_mu, mean_log_sigma, size=n_genes))

    # Condition fold change: 1 everywhere except the true-DE genes in the disease group.
    cond_fc = np.ones(n_genes)
    de_idx = rng.choice(n_genes, size=n_de, replace=False) if n_de > 0 else np.array([], dtype=int)
    # Mix up- and down-regulation for the true DE set.
    signs = rng.choice([-1.0, 1.0], size=len(de_idx))
    cond_fc[de_idx] = 2.0 ** (signs * log2fc)
    de_genes = [genes[i] for i in de_idx]
    true_log2fc = {genes[i]: float(signs[k] * log2fc) for k, i in enumerate(de_idx)}

    # Per-(gene, donor) random effect: log-normal, mean-preserving (center at -sigma^2/2).
    if donor_sigma > 0:
        donor_re = np.exp(
            rng.normal(-0.5 * donor_sigma**2, donor_sigma, size=(n_genes, n_donors))
        )
    else:
        donor_re = np.ones((n_genes, n_donors))

    blocks_X = []
    obs_donor: list[str] = []
    obs_group: list[str] = []
    for d_idx, donor in enumerate(donors):
        grp = groups[d_idx]
        fc = cond_fc if grp == "disease" else np.ones(n_genes)
        # Expected per-gene mean for this donor's cells (before per-cell depth).
        donor_mean = mu_g * donor_re[:, d_idx] * fc  # (n_genes,)
        depth = np.exp(rng.normal(0.0, depth_log_sigma, size=n_cells_per_donor))  # (n_cells,)
        # (n_cells, n_genes) expected means.
        lam = depth[:, None] * donor_mean[None, :]
        counts = _nb_counts(lam, dispersion, rng)
        blocks_X.append(counts)
        obs_donor.extend([donor] * n_cells_per_donor)
        obs_group.extend([grp] * n_cells_per_donor)

    X = np.vstack(blocks_X).astype(np.float32)
    obs = pd.DataFrame(
        {
            "donor": pd.Categorical(obs_donor),
            "condition": pd.Categorical(obs_group, categories=["ctrl", "disease"]),
            "cell_type": pd.Categorical(["T_cell"] * X.shape[0]) if one_cell_type else None,
        }
    )
    if not one_cell_type:
        obs = obs.drop(columns=["cell_type"])
    adata = ad.AnnData(X=X, obs=obs)
    adata.var_names = genes
    adata.obs_names = [f"c{i:06d}" for i in range(X.shape[0])]

    return Oracle(
        adata=adata,
        de_genes=de_genes,
        log2fc=true_log2fc,
        params=dict(
            n_genes=n_genes,
            n_donors_per_group=n_donors_per_group,
            n_cells_per_donor=n_cells_per_donor,
            n_de=n_de,
            donor_sigma=donor_sigma,
            log2fc=log2fc,
            dispersion=dispersion,
            seed=seed,
        ),
    )


def null_oracle(seed: int = 0, **kw) -> Oracle:
    """No true DE, donor random effect present -> a correct method calls ~0 DE."""
    kw.setdefault("n_de", 0)
    return simulate(seed=seed, **kw)


def positive_oracle(n_de: int = 100, seed: int = 0, **kw) -> Oracle:
    """``n_de`` genes with a real, donor-consistent fold change (default |log2FC|=1)."""
    kw.setdefault("log2fc", 1.0)
    return simulate(n_de=n_de, seed=seed, **kw)


def no_donor_effect_oracle(seed: int = 0, **kw) -> Oracle:
    """Falsification control: donor_sigma=0 -> naive and pseudobulk must agree."""
    kw["donor_sigma"] = 0.0
    kw.setdefault("n_de", 0)
    return simulate(seed=seed, **kw)


if __name__ == "__main__":
    o = null_oracle(seed=1)
    print("NULL oracle:", o.adata)
    print("donors:", o.adata.obs.donor.nunique(), "| groups:",
          o.adata.obs.condition.value_counts().to_dict(), "| true DE:", o.n_true_de)
    print("count matrix dtype/min/max:", o.adata.X.dtype, o.adata.X.min(), o.adata.X.max())
