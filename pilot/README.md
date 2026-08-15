# Phase 0 pilot — real-data sweep

The measurement instrument is built and runs end to end on synthetic oracles with known ground truth
(`../scripts/synthetic_gate.py`, `../docs/PILOT_FINDINGS.md`). What remains for the actual GO/NO-GO gate is to
run the **same code** over real public datasets and see whether the naive-vs-pseudobulk inflation holds up
where the confounds are messy. Methodology: [`../docs/PHASE0_SPEC.md`](../docs/PHASE0_SPEC.md).

> [!important]
> **The instrument now passes its validity gate — within a stated operating envelope, and not outside it.**
> The history matters and is kept: Amendment 1 found the pseudobulk arm (then DESeq2-Wald) miscalibrated;
> [Amendment 2](../docs/AMENDMENTS.md) replaced its test with moderated eBayes, which fixed calibration but
> left power at 0.35 against the required 0.60; [Amendment 3](../docs/AMENDMENTS.md) re-scoped **where**
> §8(c)'s threshold binds, from an arbitrary simulator setting to the boundary of a declared envelope, and
> raised the permutation count 40 → 200 before re-running. The 2026-08-15 run
> ([`gate/synthetic_gate_2026-08-15.json`](gate/synthetic_gate_2026-08-15.json), 511.8 s) passes all six
> criteria: `lambda_pseudobulk` 1.01 in the pre-registered [0.9, 1.1] band and permutation-null FP rate
> 0.035 (7/200, MC SE 0.013) — both at the hard regime `sigma_donor` = 0.5 — and sensitivity 0.86 against
> the required 0.60 at the **unchanged** pre-registered effect size (log2FC = 1.0, K = 200), evaluated at
> the envelope boundary `sigma_donor` = 0.35 with 8 donors/group.
>
> **Read the scope before reading the verdict.** Power 0.60 at `sigma_donor` = 0.5 with 8 donors/group is
> **still unmet (0.35) and is not claimed**; the envelope requires minimum donors/group of 4 / 8 / 13 / 23
> at `sigma_donor` 0.2 / 0.35 / 0.5 / 0.7, so a stratum near 0.5 needs ≥ ~13 donors per group or must be
> excluded. That is a *narrower* claim than before, not a relaxed one. And `sigma_donor` itself has never
> been anchored to real data, so **whether any real stratum falls inside the envelope is unknown** —
> Amendment 3 supplies the per-stratum estimation mechanism and explicitly not the anchor. Stratum
> inclusion in the real sweep therefore now requires a per-stratum `sigma_donor` estimate and envelope
> membership, on top of everything §1 already pins. Because the pseudobulk arm is the denominator of every
> inflation number, no result here may be read as a finding, and every number in this repo is still
> synthetic. Any further change to the frozen protocol will be recorded, dated, in
> `../docs/AMENDMENTS.md` before it is applied.

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
- `pbcheck.census_select` — obs-only stratum selection (§1): pinned census version, `(dataset_id × cell_type)`
  strata with one `disease vs normal` contrast per term, the inclusion-gate items obs can decide (thin donors
  dropped not merged, ≥ 3 donors/group, donor nested in condition), the donor-level confound pre-screen with
  the pooling flag, and the candidate manifest (JSON + CSV). **Candidates, not a pre-registration** — see the
  gap note below
- `pbcheck.io_counts` — the counts gate (§9 item 2): raw-integrality assertion on the values and not the
  dtype (item 4 / §10 risk 6 — a fractional matrix is **dropped with its reason, never rounded**; sparse `X`
  is never densified), the frozen universe sized through the arm's own aggregation with C5 enforced (item 5)
  and §3's ≥ 3 post-aggregation profiles reported, the remaining pending manifest columns filled, and a
  **discrepancy flag** where the load disagrees with the obs snapshot. The core takes an `AnnData` and no
  network; the Census fetch is one lazily-imported wrapper on the pinned version's raw layer
- `scripts/census_candidates.py` + `.github/workflows/census-candidates.yml` — the driver that runs §9 item 1
  over the **whole** pinned Census, and the manual CI job that dispatches it. Two passes, because one query
  does not fit a 16 GB runner: a streaming pass folds ~50-65 M cells into per-`(dataset, cell type, disease,
  donor)` counts without materialising them and decides only *which* datasets are worth reading; a second
  pass re-reads each of those and hands the per-cell frame to `census_select.screen_strata` unchanged. It
  adds no method and moves no threshold. Its tests run offline like everything else listed here — **the run
  itself has not been made**, see below
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
  actual call in `tests/test_methods.py`. The follow-on gap — on scanpy 1.12.2 `rank_genes_groups_df`
  returns only `pct_nz_group` when an explicit `reference` is given — **is now closed**: the reference
  fraction was never lost, only dropped by the dataframe helper. `uns['rank_genes_groups']['pts']` is a
  genes x groups frame carrying every level of the groupby, so `naive.py` reads it directly and reindexes
  from `var_names` order onto the ranked gene order, emitting both `pct_group` and `pct_reference`. The
  min-expression sensitivity check (§6/B5) now has both sides. This is an implementation fix, not a
  protocol change — the spec already required both fractions — so it carries no amendment. The dataframe
  columns remain as a fallback for a scanpy version that stops populating `uns['pts']`.
- **`design.py`'s Cramér's V — resolved** (`aecf8c5`, *"Measure design confounding per donor, not per cell"*).
  Batch/assay/suspension/pool confounding is now measured on the deduplicated donor-level table, not
  weighted by cell count; a design where one batch is entirely one condition now scores V = 0.577
  (was 0.128, cell-weighted, which raised no flag).
- **The `sigma_donor` of the synthetic oracle is a free knob, unanchored to any real
  mean-dispersion / donor-variance trend (§8(b)).** Both amendments close on this as the outstanding
  threat to every power number computed here. **Still open, and it is the important one.**

## What remains

Per spec §9, in order. Items 1 and 2 (`census_select`, `io_counts`) are built and listed above; the
`[census]` extra is needed only to *fetch* — every check in either module runs on a locally supplied
`AnnData`, and their tests use no network at all.

Before any of the modules below: **the candidate run itself has not been made.** The driver and its CI job
exist; nothing has yet read the real Census. The first dispatch is the `census candidates` workflow with
`dry_run: true` — pass 1 streams the whole Census either way, and the cells/second it reports is the only
honest input to the full run's budget against the runner's hard 6-hour kill. Its output is an artifact
(candidate manifest + log), uploaded and not committed. What comes back is a list of **candidates**: the
choice of 8–12 datasets from it, and the freezing of that choice, is the separate human act described at the
end of this section.

1. **`controls.py`** — cells-per-donor sweep (the primary conditioning axis, D1), donors-per-group
   sweep, depth-match downsampling, cell-type annotation ontology depth (D5 — still a `pending`
   column of the selection manifest, and the one `io_counts` explicitly does not fill: depth is a
   property of the CL graph, not of `X`).
2. **`decision.py`** — the pre-registered GO/NO-GO rule, clustered by dataset (D2), pseudobulk
   validity gate first.
3. **`report.py`** — jinja2 HTML: per-stratum/per-dataset tables, null-distribution plots,
   floor-vs-cells-per-donor curves, λ, oracle pass/fail, provenance manifest.

Not a module, and not automatable: **the §1 pre-registration of the stratum list itself** — choosing
8–12 datasets that *span* the outcome space and freezing that list by commit before any metric is
computed. `census_select` proposes candidates and stamps every row `admitted_to_sweep = False`, and
`io_counts` — which settles two of that row's blockers — leaves it `False` too; selecting from them
is a judgement, and admission additionally needs the per-stratum `sigma_donor` estimate and envelope
membership Amendment 3 leaves open.

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
