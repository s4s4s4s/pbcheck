"""The frozen universe is the one thing both arms and every permutation must agree on."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from pbcheck.gene_universe import UniverseTooSmall, frozen_universe, require_min_size


def _pdata(counts: np.ndarray) -> AnnData:
    n_donors, n_genes = counts.shape
    return AnnData(
        X=counts.astype(float),
        obs=pd.DataFrame(index=[f"d{i}" for i in range(n_donors)]),
        var=pd.DataFrame(index=[f"g{j}" for j in range(n_genes)]),
    )


def test_universe_is_label_agnostic_and_sorted():
    """Correction A3: the filter must not depend on the condition labels.

    It never sees them — the function takes no label argument — so the check that matters is that
    the result depends only on the counts, and is deterministic and ordered.
    """
    rng = np.random.default_rng(0)
    counts = rng.poisson(20, size=(8, 300))
    a = frozen_universe(_pdata(counts))
    b = frozen_universe(_pdata(counts))
    assert a == b == sorted(a)
    # Permuting the donor rows is exactly what a relabelling does to this filter: nothing.
    assert frozen_universe(_pdata(counts[::-1])) == a


def test_low_expression_genes_are_dropped():
    counts = np.zeros((6, 4))
    counts[:, 0] = 100          # expressed everywhere, high total -> keep
    counts[:5, 1] = 10          # detected in 5/6 donors, total 50 -> keep
    counts[0, 2] = 3            # detected in 1/6 donors, total 3 -> drop
    counts[:, 3] = 0            # never detected -> drop
    keep = frozen_universe(_pdata(counts), min_total_count=15, min_prop=0.5)
    assert keep == ["g0", "g1"]


def test_min_size_gate_is_enforced_when_given_and_off_by_default():
    """Regression: min_size used to be accepted, documented, and then ignored.

    An argument that looks like a gate but silently does nothing is worse than no argument — the
    caller believes a degenerate stratum was skipped when it was not.
    """
    counts = np.zeros((6, 4))
    counts[:, 0] = 100
    tiny = _pdata(counts)

    # Default: no gate, caller decides.
    assert frozen_universe(tiny) == ["g0"]

    # Given: actually enforced.
    with pytest.raises(UniverseTooSmall, match="too sparse"):
        frozen_universe(tiny, min_size=200)


def test_require_min_size_passes_through_when_large_enough():
    genes = [f"g{i}" for i in range(250)]
    assert require_min_size(genes, 200) is genes
