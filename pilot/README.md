# Phase 0 pilot — real-data sweep

The measurement instrument is built and runs end to end on synthetic oracles with known ground truth
(`../scripts/synthetic_gate.py`, `../docs/PILOT_FINDINGS.md`). What remains for the actual GO/NO-GO gate is to
run the **same code** over real public datasets and see whether the naive-vs-pseudobulk inflation holds up
where the confounds are messy. Methodology: [`../docs/PHASE0_SPEC.md`](../docs/PHASE0_SPEC.md).

> [!important]
> **The instrument is not yet validated.** At the parameters the spec pre-registers, the pseudobulk arm does
> not currently meet its own binding validity gate (spec §8(a)/§8(c)): `lambda_pseudobulk` sits outside the
> pre-registered [0.9, 1.1] band, the permutation-null false-positive rate exceeds alpha, and sensitivity at
> the pre-registered effect size (log2FC = 1.0, K = 200) falls short of the required 0.60. Because the
> pseudobulk arm is the denominator of every inflation number, no result here may be read as a finding until
> that is resolved. Diagnosis is in progress; any change to the frozen protocol will be recorded, dated, in
> `../docs/AMENDMENTS.md` before it is applied.

## What is done (this repo, runs offline)

- `pbcheck.gene_universe` — label-agnostic frozen universe (spec §3, A3)
- `pbcheck.methods.naive` / `pbcheck.methods.pseudobulk` — the two DE arms (poscounts, cooks_filter=False)
- `pbcheck.mtc` — identical BH over the frozen universe for both arms (§5)
- `pbcheck.permutation` — donor-permutation null, both arms (§4)
- `pbcheck.metrics` — genomic-inflation λ, permutation floor, signal-above-floor (§6)
- `pbcheck.design` — obs-only design auditor (§7)
- `synthetic/oracles.py` + `scripts/synthetic_gate.py` — validation on known truth (§8 b/c)

### Known gaps in what is listed above

Several spec corrections are named in module docstrings but are **not** implemented; they are listed here so
that nothing reads as done when it is not:

- **A2** — permutations are not cell-count stratified; per-permutation cell totals are logged but never consumed.
- **B5** — λ is computed from analytic Wilcoxon / DESeq2 p-values, not from empirical permutation-null p-values.
- **C3** — the native independent-filtered `padj` cross-check is never computed.
- **C5** — the minimum-frozen-universe SKIP gate is never called (`min_size` is accepted and ignored).
- Thin-donor filtering (`min_cells` / `min_counts`) does not run: decoupler 2.x has no such parameters.
- `naive.py` does not pass `tie_correct=True` / `pts=True`, which spec §2 pins verbatim.
- `design.py` computes the confounding Cramér's V per cell rather than per donor, which can miss a genuine
  donor-level confound.

## What remains (needs `pip install -e ".[census]"`)

Per spec §9, in order:

1. **`census_select.py`** — open a *pinned* Census version, query obs, apply the inclusion gate
   (≥3 donors/group, ≥10 cells/donor, integer counts), confound pre-screen (Cramér's V + perfect
   separation + pooling flag), and **pre-register** the stratum list before any metric is computed.
2. **`io_counts.py`** — load a `(dataset_id × cell_type)` stratum, assert integer raw counts.
3. **`controls.py`** — cells-per-donor sweep (the primary conditioning axis, D1), donors-per-group
   sweep, depth-match downsampling.
4. **`decision.py`** — the pre-registered GO/NO-GO rule, clustered by dataset (D2), pseudobulk
   validity gate first.
5. **`report.py`** — jinja2 HTML: per-stratum/per-dataset tables, null-distribution plots,
   floor-vs-cells-per-donor curves, λ, oracle pass/fail, provenance manifest.

## First pass (spec §1)

8–12 datasets chosen to **span** the outcome space (strong-effect, subtle, varied assay/tissue/donor
count, cells-per-donor across pre-registered bins) — not cherry-picked wins. Real anchor: **Mathys
2019** AD snRNA-seq, to reproduce Murphy & Skene 2023 qualitatively (§8 d). Pre-register the stratum
list first.

## Decision rule (short form)

GO if, on a majority of **independent datasets** at matched cells-per-donor: pseudobulk is calibrated
(perm-null rejects ~0) **and** powered (synthetic positive), naive λ ≥ 2, and the naive permutation
floor is far above the BH complete-null expectation (~0). NO-GO if pseudobulk fails its validity gate,
or naive inflation vanishes once cells-per-donor is matched (i.e. it was a depth/cell-count artifact).
