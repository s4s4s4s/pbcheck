# Environment notes (verified 2026-07-19)

Machine: Windows 11, Python 3.12.10, 20 CPUs. venv at `.venv/` (Git Bash: `.venv/Scripts/python.exe`).
`requirements.lock` holds the exact frozen versions. No conda; no system C compiler (all deps installed from wheels).

## Installed & smoke-tested (single process, pilot import order)
numpy 2.4.6 · scipy 1.18.0 · pandas 3.0.3 · anndata · scanpy 1.12.2 · statsmodels 0.14.6 ·
pydeseq2 0.5.4 · decoupler 2.2.0 · jinja2 3.1.6 · matplotlib 3.11.1 · tqdm · pytest 9.1.1 · hypothesis

`cellxgene-census` is intentionally NOT yet installed — it is only needed for the real-data audit
(`pip install -e ".[census]"`). Phase-0 engine validation runs entirely on synthetic oracles, no downloads.

## API drift to respect (these differ from older tutorials / the design docs)

**decoupler 2.x (2.2.0)** — the old `dc.get_pseudobulk(...)` is GONE. Use:
- `dc.pp.pseudobulk(adata, sample_col, groups_col, layer=None, raw=False, mode='sum', ...) -> AnnData`
  aggregates to one profile per (sample_col value × groups_col value). For our audit:
  `sample_col="donor"`, `groups_col="cell_type"`, `mode="sum"` on **raw integer counts**.
- Gene filtering (edgeR filterByExpr port): `dc.pp.filter_by_expr(pdata, group=..., min_count=..., min_total_count=...)`
  and `dc.pp.filter_by_prop(...)`. Also `dc.pp.filter_samples(...)`.

**PyDESeq2 0.5.4**:
- `DeseqDataSet(counts=<DataFrame samples×genes, raw ints>, metadata=<DataFrame>, design="~condition", ...)`
  (formula string via `design=`; `min_replicates=7` default triggers Cook's refit).
- `DeseqStats(dds, contrast=["condition", "disease", "ctrl"], alpha=0.05)`, then `.summary()` -> `.results_df`.

**scanpy 1.12.2** naive per-cell DE: `sc.tl.rank_genes_groups(adata, groupby="condition", method="wilcoxon")`,
results in `adata.uns["rank_genes_groups"]` (pull with `sc.get.rank_genes_groups_df`).

## Gotchas seen
- First import of `anndata` after a fresh install threw a transient `_proxy` DLL load error once, then
  never again — a first-load race, not a real breakage. Full stack is stable in one process.
- numpy 2.4 / pandas 3.0 are bleeding-edge; keep `requirements.lock` and pin if anything downstream breaks.
