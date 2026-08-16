"""One ranking pass per stratum, then every donor relabeling is an addition.

The naive arm (:mod:`pbcheck.methods.naive`) is re-run from scratch on every donor-permutation:
copy the AnnData, renormalise, re-log, re-rank all cells for all genes. Measured at ~1.5e-7 s per
cell per gene, that is ~2870 core-hours for the naive arm alone over the 251 frozen strata at spec
§4's ``n_perm`` = 1000 — which makes the pre-registered sweep infeasible, not slow.

**It is exactly recoverable in one pass**, because the donor-permutation null relabels donors and
never moves a cell:

* per-cell normalisation and ``log1p`` are functions of the cell alone, so the value matrix is
  identical under every labeling;
* the ranks of a gene over the **pooled** cells of both groups are therefore identical too — scanpy
  ranks the union of the two group masks, and under a donor relabeling that union is the same set
  of cells;
* the tie correction depends only on the pooled tie structure, so it is invariant as well;
* and the Wilcoxon rank-sum statistic of the test group is the **sum of the per-donor rank sums**
  of the donors assigned to it, because a rank sum is a sum over cells and the donors partition the
  cells.

So a stratum needs one ranking pass; after it every permutation is a ``n_donors x n_genes`` masked
row-sum. The same collapse holds for spec §2's robustness variant ``t-test_overestim_var``: scanpy
takes its group means and variances from :func:`fast_array_utils.stats.mean_var`, which computes
``E[x^2] - E[x]^2`` (checked in the installed 1.12.2 stack, dense and sparse paths alike), so the
per-donor sufficient statistics ``n_d``, ``sum x``, ``sum x^2`` reconstruct them.

WHAT "EXACT" MEANS HERE, precisely, because the whole value of this module is that it is not an
approximation of a pre-registered statistic:

* **Ranks — bitwise.** scanpy's ``rankdata`` assigns ``0.5 * (block_start + block_last + 2)``, an
  exact multiple of 0.5. :func:`_average_ranks` computes the same integers by the same rule.
* **Rank sums — bitwise, in any summation order.** Every rank is a multiple of 0.5 and the largest
  attainable sum is ``n^2/2`` — 1.5e11 for the biggest frozen stratum, far under 2**53. Every
  partial sum is therefore exactly representable, so float64 addition is exact and associative
  here, and per-donor-then-per-group gives the identical bits to scanpy's single group-wise sum.
* **Tie coefficient — same expression, same operand order**, including the strictly sequential
  reduction (scanpy runs ``_tiecorrect`` under numba, whose ``.sum()`` accumulates left to right).
* **Means and variances — the same algebra, and NOT reliably the same bits.** These feed ``log2fc``
  (both methods) and the whole ``t-test_overestim_var`` statistic, and the operation that loses the
  bits is named precisely: :func:`fast_array_utils.stats.mean_var` computes the variance as
  ``E[x**2] - E[x]**2``, forming the squares in the array's **stored dtype** — float32 for a
  normalised count matrix — before accumulating them in float64.

  Write ``kappa = E[x**2] / var``. The relative error of that variance is then about
  ``(eps32/2 + n*eps64) * kappa``, in which ``eps32/2 = 5.96e-08`` dominates by eight orders of
  magnitude. On log1p-normalised data ``kappa`` is large by construction — the values are logs of a
  normalised count, so they sit tightly around a mean of order 9 while the within-group variance
  can be 1e-13 — and it does **not** shrink as cells are added (measured: 1.9e14 at 2 cells per
  group, still 9.9e13 at 120). At ``kappa`` = 1.9e14 the formula returns a variance seven orders of
  magnitude too large, and the p-value built on it is arithmetic noise.

  That is a property of the pre-registered statistic, not of this module — it was measured with
  scanpy on both sides of the comparison. The consequence to plan around: **an
  ``t-test_overestim_var`` p-value is not reproducible across machines** on genes whose group
  values collapse onto neighbouring float32 values, whichever engine computes it. Observed as
  exactly that: engine and scanpy agree bitwise on 7200 configurations here, while CI's scanpy
  returned the value exact rational arithmetic gives. What this module guarantees, and what
  ``tests/test_naive_engine.py`` asserts tightly, is one level down — the per-donor ``n``,
  ``sum x`` and ``sum x**2`` are correctly rounded, and the formula applied to them is scanpy's,
  bit for bit. A better-conditioned two-pass variance here would be more accurate and would be a
  **different number** from the pre-registered one, so it is not used.

  The Wilcoxon arm — the pre-registered primary (spec §2) — is untouched by any of this. It never
  subtracts two nearly equal large numbers, which is why it is bitwise on every platform.

Nothing here is a new statistic. The slow path stays callable and the tests run the two against
each other on the same inputs, which is the only form of the claim a reviewer can check.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import norm, ttest_ind_from_stats

from pbcheck.methods.de import DEResult

#: Methods the engine reproduces. ``wilcoxon`` is the pre-registered naive test (spec §2, with
#: ``tie_correct=True``); ``t-test_overestim_var`` is §2's robustness variant. ``t-test`` is
#: deliberately absent: it differs from the overestim variant only in ``nobs2``, but it is not in
#: the spec and an untested third path is a liability, not a feature.
SUPPORTED_METHODS = ("wilcoxon", "t-test_overestim_var")

#: Working-memory budget for one gene chunk of the setup pass, in bytes. The ranking builds about
#: eight ``n_cells x chunk`` temporaries (sorted values, order, two boundary masks, block starts and
#: ends, the ranks, the tie accumulator), so the per-element cost is budgeted at 64 bytes rather
#: than at the 8 an ``float64`` array alone would suggest. Chunking is over GENES and never over
#: cells: a rank is defined over the whole pooled column and cannot be split.
DEFAULT_CHUNK_BYTES = 256 * 1024 * 1024
_BYTES_PER_ELEMENT = 64

#: Above this many cells the exactly-representable-integer argument for the ranks would need
#: float64 arithmetic that a float32 input would not supply. 2**22 = 4194304 cells is an order of
#: magnitude past the largest frozen stratum (547665); the guard exists so the claim above is
#: enforced rather than asserted.
_MAX_CELLS_FOR_NATIVE_DTYPE_RANKS = 2 ** 22


def _average_ranks(block: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Column-wise average ranks and the per-column tie term, as scanpy computes them.

    Returns ``(ranks, tie_sum)`` where ``ranks`` is float64 ``(n_cells, n_genes_in_block)`` and
    ``tie_sum[j]`` is ``sum_over_tie_blocks(t**3 - t)`` for column ``j``.

    This is scanpy 1.12.2's ``scanpy.tools._rank_genes_groups.rankdata`` and ``_tiecorrect``
    written against numpy instead of numba, kept here rather than imported because those are
    private names. ``tests/test_naive_engine.py`` asserts the two agree **bitwise**, so the
    duplication is pinned rather than hoped for, and a scanpy change that moves the statistic turns
    into a red test instead of a silent divergence.

    The tie term is accumulated with :func:`numpy.add.accumulate` down the column rather than with
    ``.sum()``: numba reduces sequentially and numpy's ``.sum()`` pairwise, and for the tie blocks
    of a real stratum ``t**3`` exceeds 2**53 and the two orders need not agree to the last bit.
    Non-block positions contribute an exact ``0.0``, which a sequential accumulation ignores.
    """
    n = block.shape[0]
    # Tie order is irrelevant -- every member of a tie block is given the same average rank -- so
    # the default (unstable) sort is used, exactly as scanpy's ``rankdata`` does.
    order = np.argsort(block, axis=0)
    srt = np.take_along_axis(block, order, axis=0)

    # A position starts a tie block if its value differs from the one below it in sorted order;
    # it ends one if it differs from the one above.
    differs = srt[1:] != srt[:-1]
    is_first = np.empty(srt.shape, dtype=bool)
    is_first[0] = True
    is_first[1:] = differs
    is_last = np.empty(srt.shape, dtype=bool)
    is_last[-1] = True
    is_last[:-1] = differs
    del differs, srt

    # int32 positions halve the working set; ``starts + lasts + 2`` peaks at 2n + 2, and the guard
    # in :data:`_MAX_CELLS_FOR_NATIVE_DTYPE_RANKS` keeps n three orders of magnitude below int32.
    pos = np.arange(n, dtype=np.int32)[:, None]
    starts = np.maximum.accumulate(np.where(is_first, pos, -1), axis=0)
    lasts = np.minimum.accumulate(np.where(is_last, pos, n)[::-1], axis=0)[::-1]

    # scanpy: ranked = 0.5 * (count[dense] + count[dense - 1] + 1), i.e. 0.5 * (end_exclusive +
    # start + 1) with end_exclusive = last + 1. Integer inside, so the float64 result is exact.
    ranked_sorted = 0.5 * (starts + lasts + 2).astype(np.float64)
    ranks = np.empty(ranked_sorted.shape, dtype=np.float64)
    np.put_along_axis(ranks, order, ranked_sorted, axis=0)
    del ranked_sorted, order

    cnt = (lasts - starts + 1).astype(np.float64)
    del starts, lasts
    tie = np.where(is_first, cnt ** 3 - cnt, 0.0)
    np.add.accumulate(tie, axis=0, out=tie)
    return ranks, tie[-1].copy()


def _segment_sums(values: np.ndarray, starts: np.ndarray) -> np.ndarray:
    """Sum contiguous row segments. ``starts`` must be strictly increasing and start at 0."""
    return np.add.reduceat(values, starts, axis=0)


@dataclass(frozen=True)
class NaiveRelabelEngine:
    """Everything a stratum's naive test needs, reduced to per-donor sufficient statistics.

    Build once with :meth:`from_adata`, then call :meth:`test` for each donor labeling. The result
    is a :class:`~pbcheck.methods.de.DEResult` with the same columns and the same numbers as
    :func:`pbcheck.methods.naive.naive_de` on the same AnnData with the same labels — see the
    module docstring for what "the same numbers" means at the bit level, and
    ``tests/test_naive_engine.py`` for the property tests that hold it there.

    One difference from ``naive_de`` and it is deliberate: the table comes back in **gene order**
    (the ``genes`` the engine was built over), not in scanpy's score-descending order. Every
    consumer in the tree reindexes onto the frozen universe before reading anything, so the order
    carries no information; reproducing an ordering that ``argpartition`` leaves unspecified for
    tied scores would be reproducing an implementation detail rather than the statistic.
    """

    genes: pd.Index
    donors: pd.Index
    method: str
    n_cells: int
    cells_per_donor: np.ndarray            # (D,) int64
    donor_sum: np.ndarray                  # (D, G) float64, over log1p(normalised) values
    donor_nnz: np.ndarray                  # (D, G) int64
    donor_rank_sums: np.ndarray | None = None   # (D, G) float64 — wilcoxon only
    tie_coef: np.ndarray | None = None          # (G,)  float64 — wilcoxon only
    donor_sumsq: np.ndarray | None = None       # (D, G) float64 — t-test only
    setup: dict = field(default_factory=dict)   # provenance for the run artifact

    # -- construction ------------------------------------------------------------------------

    @classmethod
    def from_adata(
        cls,
        adata,
        *,
        donor_col: str = "donor",
        genes=None,
        method: str = "wilcoxon",
        target_sum: float = 1e4,
        chunk_bytes: int = DEFAULT_CHUNK_BYTES,
    ) -> NaiveRelabelEngine:
        """Run the one-off pass: normalise, log, rank in gene chunks, reduce to per-donor rows.

        ``adata``, ``genes`` and ``target_sum`` mean exactly what they mean in
        :func:`pbcheck.methods.naive.naive_de`, and the gene restriction is applied **before**
        normalisation there and here — the library size a cell is scaled by depends on which genes
        are in the universe, so doing it in the other order would silently change every value.

        Memory: the value matrix is never densified as a whole. A sparse input is converted to CSC
        once so that a gene chunk is a cheap column slice (peak: two copies of the sparse matrix
        during the conversion), and only one ``n_cells x chunk`` block is dense at a time. The
        per-donor stores are ``n_donors x n_genes``, which is megabytes.
        """
        if method not in SUPPORTED_METHODS:
            raise ValueError(f"method must be one of {SUPPORTED_METHODS}, got {method!r}")

        a = adata if genes is None else adata[:, list(genes)]
        a = a.copy()
        sc.pp.normalize_total(a, target_sum=target_sum)
        sc.pp.log1p(a)

        gene_index = pd.Index(a.var_names, name="gene")
        donor_labels = np.asarray(a.obs[donor_col].astype(str))
        donors = pd.Index(np.unique(donor_labels), name=donor_col)
        donor_of_cell = donors.get_indexer(donor_labels)
        if (donor_of_cell < 0).any():  # pragma: no cover — get_indexer over its own uniques
            raise ValueError(f"cells with an unmappable {donor_col!r}")

        # Sort cells by donor so each donor is a contiguous row segment and the per-donor
        # reductions are one reduceat instead of D masked sums. Ranks are per-cell values, so
        # permuting rows before ranking cannot change any of them.
        cell_order = np.argsort(donor_of_cell, kind="stable")
        cells_per_donor = np.bincount(donor_of_cell, minlength=len(donors)).astype(np.int64)
        if (cells_per_donor == 0).any():  # pragma: no cover — donors come from the data
            raise ValueError("a donor with no cells cannot be a permutation unit")
        seg_starts = np.concatenate(([0], np.cumsum(cells_per_donor)[:-1]))

        X = a.X
        if issparse(X):
            X = X.tocsc()
            X.eliminate_zeros()
        X = X[cell_order]
        del a

        n_cells, n_genes = int(X.shape[0]), len(gene_index)
        n_donors = len(donors)
        want_ranks = method == "wilcoxon"
        want_sumsq = method == "t-test_overestim_var"

        donor_sum = np.zeros((n_donors, n_genes), dtype=np.float64)
        donor_nnz = np.zeros((n_donors, n_genes), dtype=np.int64)
        donor_rank_sums = np.zeros((n_donors, n_genes), dtype=np.float64) if want_ranks else None
        tie_coef = np.zeros(n_genes, dtype=np.float64) if want_ranks else None
        donor_sumsq = np.zeros((n_donors, n_genes), dtype=np.float64) if want_sumsq else None

        width = max(1, min(n_genes, chunk_bytes // max(n_cells * _BYTES_PER_ELEMENT, 1)))
        size = np.float64(n_cells)
        tie_denominator = size ** 3 - size

        for left in range(0, n_genes, width):
            right = min(left + width, n_genes)
            block = X[:, left:right]
            block = block.toarray() if issparse(block) else np.ascontiguousarray(block)
            if not np.isfinite(block).all():
                raise ValueError(
                    "non-finite value in the normalised matrix: ranks and the tie structure are "
                    "undefined for NaN and scanpy's own ranking would order it arbitrarily"
                )

            # scanpy squares in the STORED dtype and only then accumulates in float64
            # (fast_array_utils' power(x, 2) takes no dtype; mean(..., dtype=float64) does the
            # accumulating). Reproduce that order, not the more accurate one.
            block64 = block.astype(np.float64)
            donor_sum[:, left:right] = _segment_sums(block64, seg_starts)
            if want_sumsq:
                donor_sumsq[:, left:right] = _segment_sums((block * block).astype(np.float64),
                                                           seg_starts)
            donor_nnz[:, left:right] = _segment_sums((block != 0).astype(np.int64), seg_starts)
            del block64

            if want_ranks:
                if block.dtype != np.float64 and n_cells > _MAX_CELLS_FOR_NATIVE_DTYPE_RANKS:
                    block = block.astype(np.float64)
                ranks, tie_sum = _average_ranks(block)
                donor_rank_sums[:, left:right] = _segment_sums(ranks, seg_starts)
                tie_coef[left:right] = (1.0 if n_cells < 2
                                        else 1.0 - tie_sum / tie_denominator)
                del ranks
            del block

        return cls(
            genes=gene_index,
            donors=donors,
            method=method,
            n_cells=n_cells,
            cells_per_donor=cells_per_donor,
            donor_sum=donor_sum,
            donor_nnz=donor_nnz,
            donor_rank_sums=donor_rank_sums,
            tie_coef=tie_coef,
            donor_sumsq=donor_sumsq,
            setup={
                "n_cells": n_cells,
                "n_genes": n_genes,
                "n_donors": n_donors,
                "gene_chunk_width": int(width),
                "target_sum": float(target_sum),
                "method": method,
            },
        )

    # -- per-labeling evaluation -------------------------------------------------------------

    def _group_mask(self, test_donors) -> np.ndarray:
        mask = self.donors.isin(list(test_donors))
        n_test = int(mask.sum())
        if n_test == 0 or n_test == len(self.donors):
            raise ValueError(
                f"a labeling must split the donors: {n_test} of {len(self.donors)} in the test "
                "group leaves the other group empty and no two-sample test is defined"
            )
        return mask

    def test(
        self,
        test_donors,
        *,
        condition_col: str = "condition",
        test_level: str = "disease",
        ref_level: str = "ctrl",
    ) -> DEResult:
        """The naive test for one donor labeling, as a :class:`DEResult` over ``self.genes``."""
        mask = self._group_mask(test_donors)
        n_a = int(self.cells_per_donor[mask].sum())
        n_b = self.n_cells - n_a
        if n_a < 2 or n_b < 2:
            raise ValueError(
                f"group sizes {n_a} and {n_b}: scanpy refuses a group with fewer than 2 cells, "
                "so the engine must not silently answer where the slow path raises"
            )

        mean_group = self.donor_sum[mask].sum(axis=0) / n_a
        mean_ref = self.donor_sum[~mask].sum(axis=0) / n_b

        if self.method == "wilcoxon":
            scores, pvals = self._wilcoxon(mask, n_a, n_b)
        else:
            scores, pvals = self._ttest_overestim_var(mask, n_a, n_b, mean_group, mean_ref)
        del scores  # scanpy's z / t is not part of the DEResult contract

        # scanpy: foldchanges = (expm1(mean_group) + 1e-9) / (expm1(mean_rest) + 1e-9), then log2,
        # stored float32 (its recarray dtype, which naive_de's dataframe inherits).
        with np.errstate(over="ignore", invalid="ignore"):
            fold = (np.expm1(mean_group) + 1e-9) / (np.expm1(mean_ref) + 1e-9)
            log2fc = np.log2(fold).astype(np.float32)

        table = pd.DataFrame(
            {
                "pval": pvals,
                "padj": np.nan,          # raw; BH deferred to mtc so both arms match (spec §5)
                "log2fc": log2fc,
                "pct_group": self.donor_nnz[mask].sum(axis=0) / n_a,
                "pct_reference": self.donor_nnz[~mask].sum(axis=0) / n_b,
            },
            index=self.genes,
        )
        return DEResult(
            method=f"naive[{self.method}]",
            table=table,
            contrast=(condition_col, test_level, ref_level),
        )

    def _wilcoxon(self, mask, n_a, n_b):
        """scanpy's tie-corrected rank-sum z, assembled from the per-donor rank sums.

        Mirrors ``_RankGenes.wilcoxon``'s reference branch operand for operand::

            std_dev = sqrt(tc_coef * n_active * m_active * (n_active + m_active + 1) / 12.0)
            scores  = (rank_sum - n_active * ((n_active + m_active + 1) / 2.0)) / std_dev
            scores[isnan(scores)] = 0
            pvals   = 2 * norm.sf(|scores|)

        The ``0/0 -> 0`` branch is scanpy's handling of a gene where every cell ties (the tie
        coefficient is then exactly 0, so is ``std_dev``, and the rank sum is exactly its own
        expectation): the p-value becomes 1 rather than NaN. It is reproduced, not repaired --
        this module's job is to be the same statistic, and an all-zero gene must not become a
        rejection here that it is not in the slow path.
        """
        rank_sum = self.donor_rank_sums[mask].sum(axis=0)
        std_dev = np.sqrt(self.tie_coef * n_a * n_b * (n_a + n_b + 1) / 12.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            scores = (rank_sum - (n_a * ((n_a + n_b + 1) / 2.0))) / std_dev
        scores = np.where(np.isnan(scores), 0.0, scores)
        return scores, 2 * norm.sf(np.abs(scores))

    def _ttest_overestim_var(self, mask, n_a, n_b, mean_group, mean_ref):
        """Welch's t on the group means/variances, with the reference group's ``nobs`` faked.

        scanpy's ``t-test_overestim_var`` passes ``nobs2 = ns_group`` -- the *test* group's size --
        which is the "hack for overestimating the variance for small groups" its source calls it.
        That is the statistic spec §2 names as the robustness variant, so it is reproduced exactly,
        including that ``n_b`` never reaches the t-test.

        The variances come from the per-donor ``sum x`` and ``sum x^2`` through the same
        ``E[x^2] - E[x]^2`` that :func:`fast_array_utils.stats.mean_var` uses, with the same
        ``correction=1`` rescaling. **Read the module docstring's fourth bullet before relying on
        a p-value from this path**: that subtraction has condition number ``E[x^2]/var``, which on
        log1p data reaches 1e14 and does not improve with more cells, so on saturated genes the
        pre-registered formula has no significant digits left and its p-value is not reproducible
        across machines — for scanpy as much as for this engine. Inheriting scanpy's formula
        inherits its conditioning exactly, which is the point; a better-conditioned two-pass
        variance would be a *different* number from the pre-registered one.
        """
        var_group = self._group_variance(mask, mean_group, n_a)
        var_ref = self._group_variance(~mask, mean_ref, n_b)
        # sqrt of a variance that cancellation drove a few ulp below zero yields NaN, which is
        # exactly the path scanpy lands in and repairs below with scores -> 0, pvals -> 1.
        with np.errstate(invalid="ignore", divide="ignore"):
            scores, pvals = ttest_ind_from_stats(
                mean1=mean_group, std1=np.sqrt(var_group), nobs1=n_a,
                mean2=mean_ref, std2=np.sqrt(var_ref), nobs2=n_a,   # the overestim hack
                equal_var=False,
            )
        scores = np.where(np.isnan(scores), 0.0, scores)
        pvals = np.where(np.isnan(pvals), 1.0, pvals)
        return scores, pvals

    def _group_variance(self, mask, mean, n):
        var = self.donor_sumsq[mask].sum(axis=0) / n - mean ** 2
        if n != 1:
            var = var * (n / (n - 1))
        return var
