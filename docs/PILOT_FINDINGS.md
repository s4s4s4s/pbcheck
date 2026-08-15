# Pilot findings — synthetic calibration (2026-07-19)

First empirical run of the pbcheck comparison on synthetic oracles with **known zero true DE**
(`synthetic/oracles.py`), before formalizing the engine. Naive per-cell DE = scanpy
`rank_genes_groups` Wilcoxon on `condition` across all cells. Correct DE = decoupler pseudobulk
per donor → PyDESeq2 `~condition`. 4 vs 4 donors, 2000 genes, dispersion 0.2, seed 7.

## Result: on a null with donor structure, naive DE flags most of the genome; pseudobulk flags ~none

| donor_sigma | cells/donor | naive @FDR<0.05 | naive @FDR&\|LFC\|>1 | pseudobulk @FDR<0.05 |
|:-----------:|:-----------:|:---------------:|:-------------------:|:--------------------:|
| 0.2 | 300  | 1093 | 0   | 0 |
| 0.2 | 1000 | 1503 | 0   | 0 |
| 0.4 | 300  | 1496 | 34  | 1 |
| 0.4 | 1000 | 1722 | 32  | 0 |
| 0.6 | 300  | 1639 | 214 | 0 |
| 0.6 | 1000 | 1786 | 218 | 1 |
| 0.8 | 300  | 1680 | 463 | 2 |
| 0.8 | 1000 | 1820 | 452 | 7 |

Truth = 0 DE genes in every row. Naive calls **1093–1820 of 2000 genes** (55–91%); pseudobulk calls **0–7**.

## Load-bearing consequences for the engine

1. **The primary inflation metric is the FDR-only significant-gene count, not FDR+LFC.** A strict
   `|log2FC|>1` filter *masks* the phenomenon (drops naive from ~1500 to 0–34 at low donor variance).
   Report inflation primarily at FDR<0.05 alone; keep FDR+LFC as a secondary, more conservative view.
2. **More cells per donor increases naive false positives** (1093→1503 at σ=0.2). This is the
   signature of pseudoreplication — pseudo-replicates inflate apparent significance without adding
   information — and is itself a diagnostic pbcheck can surface.
3. **The effect is robust and large even at modest donor variance** (σ=0.2 already gives >1000 false
   genes). Pseudobulk holds 0–7 across the whole grid: correct FDR control.

## Falsification control passes

With `donor_sigma=0` (no donor random effect) the naive count collapses to **0** (see
`scripts/proof_of_life.py`, scenario `null_no_donor_effect`). So the inflation is attributable to
donor structure specifically, not to the naive test's larger sample size / power. This is the control
that stops a spurious GO: the gap only appears when pseudoreplication is present.

## Caveats (why this is calibration, not the gate)

- Synthetic data, one cell type, balanced design, our own generative model — the real GO/NO-GO gate
  is the same measurement over 20–50 **real** CELLxGENE datasets, where confounds are messier.
- Pseudobulk sensitivity on the POSITIVE oracle was modest (~35% at \|log2FC\|=1, 4v4 donors, only 2
  false positives) — expected for n=4/group; not a concern for the *null* inflation measurement, but
  the metric must not over-credit pseudobulk with power it lacks.
- These runs use FDR/LFC thresholds and a single dispersion; the real sweep spans a grid.

---

## Instrument validation with the verified methodology (2026-07-19, later same day)

> [!caution]
> **SUPERSEDED 2026-07-22 — this section's verdict is retracted.** The pseudobulk arm does **not** pass the
> spec's binding validity gate at its pre-registered thresholds. `λ_pseudobulk = 1.25` is outside the
> pre-registered band [0.9, 1.1]; the permutation-null false-positive rate is 0.35–0.50 against α = 0.05; and
> power at the pre-registered oracle (log2FC = 1.0, K = 200 — **not** the 1.5 / 150 used below) is 0.47–0.57
> against a required 0.60. Item 1 of "A statistics correction I made by hand" is itself **wrong** and is
> retracted. Diagnosis, evidence and the corrected statement are in
> [`AMENDMENTS.md`](AMENDMENTS.md) — Amendment 1. The numbers in the table below reproduce exactly; it is
> their *interpretation* as a passing gate that does not stand.

After the `pbcheck-phase0-design` workflow produced `PHASE0_SPEC.md`, the engine was rebuilt to the
verified methodology (label-agnostic frozen universe, poscounts size factors, `cooks_filter=False`,
identical BH over both arms, donor-permutation null, and the **genomic-inflation factor λ** as the
primary size-robust metric). `scripts/synthetic_gate.py` runs the whole instrument on known truth
(8 vs 8 donors, 1500 genes, σ_donor=0.5) and reports:

| quantity | value | reads as |
|---|---|---|
| λ_naive (perm-null p-values) | **54.8** | grossly inflated (was 50.9 before `tie_correct=True` was passed, as spec §2 pins) |
| λ_pseudobulk | 1.25 | ⚠️ STALE READING — **outside** the pre-registered band [0.9, 1.1]; arm is anti-conservative (perm-null FP 0.35–0.50), not "~calibrated" |
| naive perm-floor (median #DEG) | **1166 / 1500 (77.8%)** | pseudoreplication floor is most of the genome (1138 / 75.9% before tie correction) |
| pseudobulk perm-floor (median #DEG) | **0** (mean 0.60) | correct FDR control under the null |
| pseudobulk power (positive oracle, \|log2FC\|=1.5) | **0.83** | ⚠️ STALE — measured on an *easier* oracle than the spec pre-registers; at log2FC=1.0 / K=200 it is 0.47–0.57, **below** the binding 0.60 |

~~All six instrument checks pass → **INSTRUMENT VALID, ready for the real CELLxGENE sweep.**~~

**RETRACTED.** The gate passed only because two pre-registered criteria were never wired into it and a third
was evaluated at a substituted effect size. With the criteria restored, the instrument reports
**INSTRUMENT NEEDS ATTENTION**. See [`AMENDMENTS.md`](AMENDMENTS.md).

### A statistics correction I made by hand (the mandate: don't vibecode the stats)

1. ~~**FWER vs FDR.** The spec's calibration criterion "fraction of permutations with ≥1 rejection ≤ α"
   is a *family-wise* bound; BH controls the *false-discovery rate*, not P(≥1 rejection).~~
   **RETRACTED 2026-07-22 — this was wrong, and it was the load-bearing error.** Under the **complete
   null** every rejection is false, so V = R, FDP = 1{R ≥ 1}, and FDR = E[FDP] = P(R ≥ 1). FWER and FDR
   *coincide* there, BH bounds it, and the spec's original criterion stood. Retiring it on this argument
   removed the one check the arm would otherwise have had to meet; the measured FP rate is 0.35–0.50
   against α = 0.05. That it appeared under a heading advertising hand-verification is the lesson worth
   keeping. See [`AMENDMENTS.md`](AMENDMENTS.md).
2. **Positive-oracle "empirical FDR" is not method error.** It sits at ~0.10–0.15 across all regimes
   and rises with σ_donor — because the simulator's realized donor random effects create genuine
   between-group differences in non-injected genes that DESeq2 rightly detects, which the injected-only
   "truth" miscounts as false positives. So FDR calibration is judged by the permutation null; the
   positive oracle judges **power** only, with a loose non-catastrophic sanity bound.
3. ~~**λ is the naive-inflation flag, not a pseudobulk verdict.** λ_pb≈1.25 at n=8 with *zero* null
   rejections is expected small-sample DESeq2 behavior, not miscalibration.~~
   **RETRACTED 2026-07-22.** It is not a small-sample effect: λ_pb does not decay from n = 4 to n = 32,
   and on counts drawn from DESeq2's *own* exact-NB model the Wald FP rate is 0.425 at 8v8, reaching
   nominal only near 48 donors per group. A donor-level Welch t-test on the identical matrix is calibrated
   at every n ≥ 4. It is miscalibration.

~~Net: the primary measurement (naive inflation) is robust and validated; the pseudobulk validity gate
(calibrated + powered) is satisfied where the spec says it should be (≥8v8, detectable effect).~~

**Net, as of 2026-07-22:** the naive-arm measurement stands and reproduces exactly. The pseudobulk validity
gate is **not** satisfied — and since that arm is the denominator of every inflation number, no figure here
may be read as a finding until the gate is met. See [`AMENDMENTS.md`](AMENDMENTS.md), Amendment 1.

---

## Addendum (2026-08-15) — the pseudobulk arm's test was replaced; calibration now passes, power still does not

This is an addendum, not a rewrite: everything above stands as the record of what was found and retracted at
the time. What follows is what changed since, under [`AMENDMENTS.md`](AMENDMENTS.md) Amendment 2.

Amendment 1 (above) left the pseudobulk arm's replacement test open, pending a comparison against the
moderated methods it had not tested. That comparison is the committed 146-cell test-selection grid,
`pilot/testsel/summary.{csv,json}` (`72dec7b`), read by `scripts/analyze_test_selection.py`: six donor-level
tests (Welch t, pooled t, Wilcoxon, `ebayes`, `ebayes_trend`, `voom`) plus `wald`, swept over `sigma_het`,
`sigma_donor`/donor-count tiers, and three generative arms. All six donor-level tests are calibrated on the
grid; `wald` is not. Among the calibrated tests, the three moderated ones (`ebayes`, `ebayes_trend`, `voom`)
separate cleanly on worst-case power (minimax regret), are within noise of each other, and `ebayes` is chosen
— simplest estimator, no measurable power cost against `ebayes_trend`, ~10× cheaper than `voom` for an
equal-or-worse power. Full arithmetic is in Amendment 2, Change 1.

Re-running the gate (`scripts/synthetic_gate.py`) under the moderated arm, with the paired BH now actually
wired into `run_null` (Amendment 2, Change 2 — an erratum: it had been specified and tested since before
Amendment 1 but no production caller used it), gives:

| quantity | this run (2026-08-15) | Amendment-1-era (DESeq2-Wald) |
|---|---|---|
| λ_pseudobulk | **1.02** (in-band [0.9, 1.1]) | 1.25 (out of band) |
| pseudobulk perm-null FP rate | **0.05** at α = 0.05 — marginal: 2/40 permutations, Monte-Carlo SE 0.034, cannot yet separate 0.05 from ~0.12 | 0.35–0.50 |
| pseudobulk power (log2FC = 1.0, K = 200, the pre-registered oracle) | **0.35** | 0.47–0.57 |
| required power (§8(c)) | ≥ 0.60 | ≥ 0.60 |

Calibration moved from clearly failing to meeting its criteria (with the FP-rate reading marked marginal
rather than clean — more permutations are needed to resolve it). Power did **not** move to passing, and this
is not a surprise: it is exactly what the grid predicts. At `sigma_donor = 0.5` — the gate's operating point
— no test in the grid, moderated or not, reaches power 0.60 anywhere; the moderated tests only clear 0.60 at
`sigma_donor = 0.35` (smallest sufficient n = 8 donors/group). See Amendment 2's power-frontier table.

This sharpens, rather than resolves, the open question both amendments have carried from the start: whether
real strata's `sigma_donor` sits near 0.35 (where the pre-registered oracle is achievable) or near 0.5–0.7
(where Amendment 1 already noted that §8(c) itself, not the test, would need amending). That decision is
**pending** — `sigma_donor` remains an unanchored free knob of the simulator (§8(b)) — and is not taken by
this addendum. Until it is, the gate correctly reports `INSTRUMENT NEEDS ATTENTION`: calibrated, not yet
powered, at the pre-registered operating point.
