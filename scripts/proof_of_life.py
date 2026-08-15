"""Proof of life: does naive per-cell DE inflate where pseudobulk does not?

Runs the pbcheck comparison on synthetic oracles with known ground truth, BEFORE the
engine is formalized. This validates three things at once: (1) the phenomenon pbcheck
audits is real and large, (2) the synthetic oracle generator behaves, (3) the decoupler
2.x / PyDESeq2 0.5.4 API path works end-to-end.

HISTORICAL -- predates the engine and the frozen protocol; not comparable to it.
This script predates ``pbcheck.mtc``, ``pbcheck.methods`` and ``docs/PHASE0_SPEC.md``,
and uses statistics the current protocol forbids:

  * it reads scanpy's own ``pvals_adj`` for the naive arm directly. ``pbcheck.mtc``
    requires one paired BH over a gene universe shared with the pseudobulk arm (spec
    §5); the two arms here are NOT corrected over the same tested set.
  * the pseudobulk arm here has no ``poscounts`` size factors and no Cook's-outlier
    handling (spec C1/C2), unlike ``pbcheck.methods.pseudobulk``.
  * it applies a strict ``|log2FC| > 1`` filter on top of FDR, which
    ``docs/PILOT_FINDINGS.md`` (2026-07-19) found MASKS the pseudoreplication
    phenomenon (drops the naive false-positive count from ~1500 to 0-34).

Its numbers are historical evidence of *why* the engine was built the way it was, not
numbers comparable to ``scripts/synthetic_gate.py`` or anything else in ``pbcheck``.
Importing this module emits a warning for exactly that reason.

Expected:
  NULL + donor effect     -> naive calls many false DE; pseudobulk ~0. Large inflation.
  NULL + NO donor effect  -> naive and pseudobulk agree (~0). Inflation collapses.  [control]
  POSITIVE + donor effect -> pseudobulk recovers the true genes; naive calls them + extras.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

warnings.warn(
    "scripts/proof_of_life.py predates pbcheck's engine and protocol: it consumes "
    "scanpy's own pvals_adj instead of the paired BH pbcheck.mtc requires, has no "
    "poscounts/Cook's-outlier handling, and applies the |log2FC| > 1 filter that "
    "docs/PILOT_FINDINGS.md found masks the phenomenon. Its numbers are historical, "
    "not comparable to scripts/synthetic_gate.py or the engine. See the module "
    "docstring.",
    FutureWarning,
    stacklevel=2,
)

warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import decoupler as dc  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.ds import DeseqStats  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "synthetic"))
from oracles import null_oracle, positive_oracle, no_donor_effect_oracle  # noqa: E402

FDR = 0.05
LFC = 1.0


def naive_de(adata, gene_universe=None) -> set[str]:
    """What a careless analyst does: per-cell Wilcoxon across cells, cells as replicates."""
    a = adata.copy()
    if gene_universe is not None:
        a = a[:, list(gene_universe)].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.tl.rank_genes_groups(a, groupby="condition", groups=["disease"], reference="ctrl",
                            method="wilcoxon")
    df = sc.get.rank_genes_groups_df(a, group="disease")
    sig = df[(df["pvals_adj"] < FDR) & (df["logfoldchanges"].abs() > LFC)]
    return set(sig["names"])


def pseudobulk_de(adata) -> tuple[set[str], list[str]]:
    """The correct analysis: aggregate to one profile per donor, DESeq2 across donors."""
    pdata = dc.pp.pseudobulk(adata, sample_col="donor", groups_col="cell_type", mode="sum")
    dc.pp.filter_by_expr(pdata, group="condition", min_count=10, min_total_count=15)
    genes_tested = list(pdata.var_names)

    counts = pd.DataFrame(np.asarray(pdata.X).round().astype(int),
                          index=pdata.obs_names, columns=pdata.var_names)
    metadata = pdata.obs[["condition"]].copy()
    metadata["condition"] = metadata["condition"].astype(str)

    dds = DeseqDataSet(counts=counts, metadata=metadata, design="~condition",
                       ref_level=["condition", "ctrl"], quiet=True, n_cpus=4)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "disease", "ctrl"], alpha=FDR, quiet=True)
    st.summary()
    res = st.results_df.dropna(subset=["padj"])
    sig = res[(res["padj"] < FDR) & (res["log2FoldChange"].abs() > LFC)]
    return set(sig.index), genes_tested


def evaluate(name: str, oracle) -> dict:
    a = oracle.adata
    pb_hits, tested = pseudobulk_de(a)
    naive_hits = naive_de(a, gene_universe=tested)  # same gene universe = fair comparison
    truth = set(oracle.de_genes)

    ratio = (len(naive_hits) / len(pb_hits)) if pb_hits else float("inf")
    row = {
        "scenario": name,
        "genes_tested": len(tested),
        "true_DE": len(truth),
        "naive_calls": len(naive_hits),
        "pseudobulk_calls": len(pb_hits),
        "inflation_ratio": ratio,
    }
    if truth:
        row["pb_sensitivity"] = round(len(pb_hits & truth) / len(truth), 3)
        row["pb_false_pos"] = len(pb_hits - truth)
        row["naive_false_pos"] = len(naive_hits - truth)
    return row


def main() -> None:
    kw = dict(n_genes=2000, n_donors_per_group=4, n_cells_per_donor=300, dispersion=0.2)
    rows = []
    print(">>> NULL + donor effect (expect naive >> pseudobulk)")
    rows.append(evaluate("null_donor_effect", null_oracle(seed=1, donor_sigma=0.35, **kw)))
    print(">>> NULL + NO donor effect  [falsification control: expect ratio ~1]")
    rows.append(evaluate("null_no_donor_effect", no_donor_effect_oracle(seed=1, **kw)))
    print(">>> POSITIVE + donor effect (expect pseudobulk recovers true genes)")
    rows.append(evaluate("positive_donor_effect", positive_oracle(n_de=100, seed=1,
                                                                  donor_sigma=0.35, **kw)))

    df = pd.DataFrame(rows)
    pd.set_option("display.width", 200, "display.max_columns", 20)
    print("\n" + "=" * 78)
    print(df.to_string(index=False))
    print("=" * 78)


if __name__ == "__main__":
    main()
