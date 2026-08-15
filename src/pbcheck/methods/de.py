"""Common representation of a differential-expression result.

Both the naive per-cell test and the correct pseudobulk test return a :class:`DEResult`, so the
inflation metric can compare them on the same footing (same gene universe, same significance rule).
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass
class DEResult:
    """A DE result for one contrast.

    ``table`` is indexed by gene and has columns ``pval``, ``padj``, ``log2fc``. Genes that could
    not be tested (filtered out, or NA padj from independent filtering / Cook's outliers) are simply
    absent or carry NaN in ``padj`` and are never counted as significant.
    """

    method: str
    table: pd.DataFrame
    contrast: tuple[str, str, str]  # (factor, test_level, ref_level)

    def __post_init__(self) -> None:
        missing = {"pval", "padj", "log2fc"} - set(self.table.columns)
        if missing:
            raise ValueError(f"DEResult.table missing columns: {sorted(missing)}")

    @property
    def genes_tested(self) -> pd.Index:
        """Genes with a non-NA adjusted p-value (i.e. actually testable)."""
        return self.table.index[self.table["padj"].notna()]

    def significant(self, fdr: float = 0.05, lfc: float | None = None) -> set[str]:
        """Genes called DE at ``padj < fdr`` and (optionally) ``|log2fc| > lfc``.

        Per the synthetic calibration (docs/PILOT_FINDINGS.md) the primary metric uses ``lfc=None``:
        a strict fold-change filter masks pseudoreplication false positives.
        """
        t = self.table
        sig = t["padj"] < fdr
        if lfc is not None:
            sig &= t["log2fc"].abs() > lfc
        return set(t.index[sig.fillna(False)])

    def n_significant(self, fdr: float = 0.05, lfc: float | None = None) -> int:
        return len(self.significant(fdr=fdr, lfc=lfc))
