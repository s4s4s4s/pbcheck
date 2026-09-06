"""A small, deterministic example generator for the ``pbcheck`` quickstart.

``example_adata`` is product code for a five-minute demo: it builds a synthetic
``AnnData`` with the columns :func:`pbcheck.audit.run_audit` needs (a donor column, a
two-level condition column, a cell-type column) and raw integer counts, so a new user
can run ``pbcheck example`` immediately followed by ``pbcheck audit`` without hunting
for a real dataset first.

It is **not** the Phase 0 oracle in :mod:`synthetic.oracles`. That module is the
correctness spec the pbcheck engine is validated against (known ground truth, tunable
true-DE genes, used by the test suite). This module exists only to hand a first-time
user something to point the CLI at; it makes no claim about statistical realism beyond
"has a donor random effect, so the naive analysis is inflated and the pseudobulk one is
not", which is the one property a quickstart needs to demonstrate.

The counts come from a per-donor gamma-Poisson (negative-binomial) mixture: a per-gene
baseline mean is perturbed by a per-(gene, donor) log-normal offset with standard
deviation ``donor_sigma`` (the donor random effect responsible for pseudoreplication)
and a per-cell log-normal depth offset, then Poisson counts are drawn from that mean. No
condition effect is injected, so any genes the naive per-cell test calls significant are
false positives driven by donor pseudoreplication, exactly what the audit is built to
flag.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

#: The shape ``pbcheck example --shape small`` (and this module's own defaults) build.
SMALL_SHAPE = {"n_genes": 600, "n_donors_per_group": 4, "n_cells_per_donor": 60}

#: The shape ``pbcheck example --shape reference`` builds: closer to a real dataset's
#: size, at the cost of a slower draw.
REFERENCE_SHAPE = {"n_genes": 8000, "n_donors_per_group": 8, "n_cells_per_donor": 625}

_DISPERSION = 0.2
_MEAN_LOG_MU = 1.0
_MEAN_LOG_SIGMA = 1.0
_DEPTH_LOG_SIGMA = 0.25
_CELL_TYPE = "example"


def example_adata(
    seed: int = 0,
    *,
    n_genes: int = SMALL_SHAPE["n_genes"],
    n_donors_per_group: int = SMALL_SHAPE["n_donors_per_group"],
    n_cells_per_donor: int = SMALL_SHAPE["n_cells_per_donor"],
    donor_sigma: float = 0.5,
) -> AnnData:
    """Build a small synthetic ``AnnData`` for the ``pbcheck`` quickstart.

    Two donor groups (``ctrl`` and ``disease``, ``n_donors_per_group`` donors each) of
    ``n_cells_per_donor`` cells each, raw integer counts over ``n_genes`` genes, and a
    per-(gene, donor) log-normal random effect of standard deviation ``donor_sigma`` that
    makes cells from the same donor correlated, the source of the naive arm's inflation.
    No condition effect is injected: this is a demo of pseudoreplication, not of a
    real signal.

    ``obs`` columns: ``donor`` (categorical, e.g. ``d00``..), ``condition``
    (categorical, ``ctrl``/``disease``), ``cell_type`` (categorical, constant
    ``"example"``, present so callers can exercise ``--celltype``/``--celltype-value``
    even though this generator only has one cell type).

    Deterministic for a given ``seed`` (and the other arguments); different seeds draw
    different data. Not the Phase 0 oracle: see the module docstring.
    """
    rng = np.random.default_rng(seed)

    n_donors = 2 * n_donors_per_group
    donor_ids = [f"d{i:02d}" for i in range(n_donors)]
    donor_conditions = ["ctrl"] * n_donors_per_group + ["disease"] * n_donors_per_group
    gene_ids = [f"g{j:04d}" for j in range(n_genes)]

    # Per-gene baseline mean, log-normal across genes.
    mu_g = np.exp(rng.normal(_MEAN_LOG_MU, _MEAN_LOG_SIGMA, size=n_genes))

    # Per-(gene, donor) random effect: log-normal, mean-preserving, the donor
    # pseudoreplication signal. donor_sigma=0 would make donors interchangeable.
    if donor_sigma > 0:
        donor_re = np.exp(
            rng.normal(-0.5 * donor_sigma**2, donor_sigma, size=(n_genes, n_donors))
        )
    else:
        donor_re = np.ones((n_genes, n_donors))

    shape = 1.0 / _DISPERSION

    blocks_x: list[np.ndarray] = []
    obs_donor: list[str] = []
    obs_condition: list[str] = []
    obs_names: list[str] = []

    for d, (donor, condition) in enumerate(zip(donor_ids, donor_conditions)):
        # Per-cell depth offset, log-normal, mean-preserving.
        depth = np.exp(rng.normal(-0.5 * _DEPTH_LOG_SIGMA**2, _DEPTH_LOG_SIGMA, size=n_cells_per_donor))
        mean = mu_g[None, :] * donor_re[:, d][None, :] * depth[:, None]
        mean = np.clip(mean, 1e-9, None)
        gamma = rng.gamma(shape=shape, scale=mean / shape)
        counts = rng.poisson(gamma).astype(np.int32)
        blocks_x.append(counts)
        obs_donor.extend([donor] * n_cells_per_donor)
        obs_condition.extend([condition] * n_cells_per_donor)
        obs_names.extend(f"{donor}_c{c:04d}" for c in range(n_cells_per_donor))

    x = np.concatenate(blocks_x, axis=0)

    obs = pd.DataFrame(
        {
            "donor": pd.Categorical(obs_donor),
            "condition": pd.Categorical(obs_condition, categories=["ctrl", "disease"]),
            "cell_type": pd.Categorical([_CELL_TYPE] * x.shape[0]),
        },
        index=pd.Index(obs_names, name=None),
    )
    var = pd.DataFrame(index=pd.Index(gene_ids, name=None))

    return AnnData(X=x, obs=obs, var=var)
