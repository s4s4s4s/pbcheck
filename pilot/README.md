# Phase 0 pilot — real-data sweep

The measurement instrument is built and runs end to end on synthetic oracles with known ground truth
(`../scripts/synthetic_gate.py`, `../docs/PILOT_FINDINGS.md`). What remains for the actual GO/NO-GO gate is to
run the **same code** over real public datasets and see whether the naive-vs-pseudobulk inflation holds up
where the confounds are messy. Methodology: [`../docs/PHASE0_SPEC.md`](../docs/PHASE0_SPEC.md).

> [!important]
> **The instrument is still not validated, but the reason changed.** Amendment 1's finding — that the
> pseudobulk arm (then DESeq2-Wald) fails calibration — is **superseded** by
> [Amendment 2](../docs/AMENDMENTS.md): the arm's test was replaced with moderated eBayes, and under it
> `lambda_pseudobulk` sits inside the pre-registered [0.9, 1.1] band and the permutation-null false-positive
> rate meets its target (0.05 against alpha = 0.05, though only 2/40 permutations — Monte-Carlo SE 0.034
> means this run cannot yet distinguish 0.05 from ~0.12). Both calibration criteria of the binding validity
> gate (spec §8(a)) are therefore met, for the first time. **Power is not**: sensitivity at the
> pre-registered effect size (log2FC = 1.0, K = 200; §8(c)) is 0.35 against the required 0.60, and the
> committed test-selection grid shows no test — moderated or not — clears 0.60 at this `sigma_donor`. That
> points to a pending decision about §8(c) itself, not an implementation gap. Because the pseudobulk arm is
> the denominator of every inflation number, no result here may be read as a finding until the gate passes
> in full. Any further change to the frozen protocol will be recorded, dated, in `../docs/AMENDMENTS.md`
> before it is applied.

## What is done (this repo, runs offline)

- `pbcheck.gene_universe` — label-agnostic frozen universe (spec §3, A3), minimum size enforced (C5)
- `pbcheck.methods.naive` — the naive per-cell arm, with the settings §2 pins verbatim
- `pbcheck.methods.moderated` / `pbcheck.methods.pseudobulk` — the pseudobulk arm: **moderated eBayes**
  (Amendment 2 Change 1) with the pre-registered thin-donor filter (Change 7); `deseq_from_pdata` retained,
  superseded, to reproduce Amendment-1-era numbers
- `pbcheck.mtc` — one **paired** BH over one common tested set for both arms (§5, Amendment 2 Change 2)
- `pbcheck.permutation` — donor-permutation null, both arms (§4)
- `pbcheck.metrics` — genomic-inflation λ, permutation floor, signal-above-floor (§6)
- `pbcheck.gate_config` — the pre-registered thresholds, in one place, tagged by provenance
- `pbcheck.design` — obs-only design auditor (§7)
- `synthetic/oracles.py` + `scripts/synthetic_gate.py` — validation on known truth (§8 b/c)

### Known gaps in what is listed above

Listed so that nothing reads as done when it is not. Statuses below are as of **Amendment 2 (2026-08-15)**;
each item names the amendment entry that settles it.

- **A2 — partially addressed, full stratification deferred to Phase 1** (Amendment 2 Change 6). Permutations
  are *not* cell-count stratified. What runs is a range check: per-permutation per-group cell totals are
  logged for both arms and the real split is verified to sit inside the permutation distribution. Rationale
  for the deferral, including why it does not gate the pseudobulk validity claim, is in the amendment.
- **B5 — implemented, and its literal construction measured and rejected** (Amendment 2 Change 3). The
  binding λ is computed from each arm's own p-values under the permutation null. B5's "rank of the observed
  statistic within the donor-permutation null" is implemented as `metrics.empirical_perm_pvalues` and
  reported as a **machinery check only**: under the sharp null it is uniform by construction (measured:
  λ = 0.93 where the binding λ was 26.08, and identical on a held-out permuted labeling), so it carries no
  calibration information.
- **C3 — superseded / retired** (Amendment 2 Change 5). The native independent-filtered `padj` cross-check
  is DESeq2-specific and dies with DESeq2. Replaced by prior-strength disclosure (`d0`, shrinkage factor) and
  a null p-value uniformity check.
- **C5 — implemented** (Amendment 2 Change 4). The minimum-frozen-universe gate is enforced and the gate
  script passes `min_size`; the realised prior (`d0`, shrinkage, trend flag, residual df) is persisted as the
  moderated analog of DESeq2's `fitType`.
- **Thin-donor filtering (`min_cells` / `min_counts`) — implemented** (Amendment 2 Change 7). decoupler 2.x
  has no such parameters, so it is applied manually after aggregation on decoupler's own `psbulk_cells` /
  `psbulk_counts`, with the spec's semantics: dropped, not merged.
- `naive.py` **does** pass `tie_correct=True` / `pts=True` as of R0 (`c13a21e`), asserted by capturing the
  actual call in `tests/test_methods.py`. Remaining gap: on scanpy 1.12.2 only `pct_nz_group` is returned
  when an explicit `reference` is given, so the min-expression sensitivity check (§6/B5) has one side of the
  expressed fractions available, not both.
- **`design.py`'s Cramér's V — resolved** (`aecf8c5`, *"Measure design confounding per donor, not per cell"*).
  Batch/assay/suspension/pool confounding is now measured on the deduplicated donor-level table, not
  weighted by cell count; a design where one batch is entirely one condition now scores V = 0.577
  (was 0.128, cell-weighted, which raised no flag).
- **The `sigma_donor` of the synthetic oracle is a free knob, unanchored to any real
  mean-dispersion / donor-variance trend (§8(b)).** Both amendments close on this as the outstanding
  threat to every power number computed here. **Still open, and it is the important one.**

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
