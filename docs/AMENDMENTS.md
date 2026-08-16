# pbcheck — Amendment Log

Pre-registration amendments to [`PHASE0_SPEC.md`](PHASE0_SPEC.md). The spec is frozen; it is changed
only through a numbered entry here, written **before** the change is applied. Each entry is dated, names
the exact spec sections it changes, states what changes and why, and — bindingly — records **what data was
already visible** when the amendment was written.

---

## Amendment 1 (2026-07-22) — Pseudobulk arm fails its binding validity gate; FP-rate criterion restored; a statistics error in PILOT_FINDINGS retracted

### Data visible at the time of this amendment (full disclosure)

At the time of writing we have seen **only synthetic-oracle results from our own generative model**
(`synthetic/oracles.py`: a log-normal donor random effect on an NB count model). We have **not run a single
real CELLxGENE Census dataset**; oracle (d), the Mathys 2019 real anchor (spec §8(d)), is untouched. **No
stratum list has been pre-registered.**

The evidence base is a diagnostic sweep (`scripts/pb_calibration_probe.py`) over
test ∈ {DESeq2-Wald, hand-rolled LRT, donor-level Welch t on log2(CPM+1), donor-level Wilcoxon},
`sigma_donor` ∈ {0, 0.2, 0.35, 0.5, 0.7}, donors/group ∈ {4, 6, 8, 12, 16, 20, 24, 32, 48}, seeds {7, 8, 9}
with independent replications {41, 53, 67, 101, 202, 303}, plus a generator arm that varies **only** the donor
random-effect distribution (log-normal / variance-matched gamma / exact-NB) at matched mean and variance.

All numbers below are computed, not asserted. Permutation-null figures for the donor-level tests are
calibrated **by construction** (donor profiles are fixed and only labels move) and are therefore **not** used
as evidence; only fresh-null cells (independently simulated data, real labels) and exact-NB generator cells
carry evidential weight.

### What we found (the trigger)

At the operating point of the repo's own instrument (8v8 donors, 1500 genes, `sigma_donor` = 0.5,
dispersion 0.2), the pseudobulk arm **fails the binding validity gate** of decision rule item 1 / §8(a):

- `lambda_pseudobulk` = 1.21–1.25, outside the pre-registered band [0.9, 1.1], reproducible across 20/20
  permutations and every seed — a systematic offset, not noise.
- Empirical permutation-null **FP rate = 0.35–0.50** at alpha = 0.05 (exact binomial P ≈ 1e-17 … 1e-23),
  roughly 7–10× nominal; the p-value tail is anti-conservative monotonically (P(p ≤ .001) ≈ 2–3× nominal).
- Power at the **pre-registered** oracle parameters (§8(c): log2FC = 1.0, K = 200) is 0.47–0.57 at 8v8,
  below the binding ≥ 0.60.

`scripts/synthetic_gate.py` previously reported PASS only because it (i) substituted an easier oracle
(log2FC = 1.5, K = 150), (ii) never wired the lambda-band or FP-rate criteria into its checks — the FP-rate
variable was computed and discarded — and (iii) relied on the retracted argument below.

### Diagnosis (why the arm fails)

**The FP failure is a property of the DESeq2 NB-GLM at small donor counts, not an artifact of our simulator.**
On counts drawn from DESeq2's *own* exact-NB model (no cells, no aggregation), the Wald FP rate is still
**0.425 at 8v8** (exact binomial P = 2.6e-23), and does not reach nominal until roughly **48 donors per group**
(FP = 0.425 / 0.25 / 0.20 / 0.15 / 0.10 at n = 8 / 12 / 20 / 32 / 48). At `sigma_donor` = 0 the FP rate is
nominal, so the failure requires donor-level overdispersion — it is not merely a small-sample effect.

On the **identical** aggregated donor matrix, a donor-level Welch t-test on log2(CPM + 1) is calibrated
(lambda ≈ 0.99, FP 0.03–0.05) at every tested n ≥ 4 and every sigma 0.2–0.7, on fresh-null data. Swapping Wald
for a likelihood-ratio test does **not** help — the LRT reproduces Wald numerically (corr(−log10 p) = 0.99993),
which localizes the fault to the NB-GLM's dispersion estimation rather than to the Wald statistic.
(Note for implementers: PyDESeq2 0.5.4 exposes no LRT; the comparison used a hand-rolled `nbinomLRT`
validated against `dds` internals to 0.0 absolute error.)

The lambda **magnitude** carries an additional component from our pre-registered log-normal donor random
effect, which DESeq2's dispersion model under-fits: under the same Wald test, lambda is in-band on
well-specified arms (gamma 1.040, exact-NB 1.035) but 1.198 on the log-normal arm. A closed form,
lambda ≈ alpha_true / alpha_MLE, reproduces the whole lambda-vs-sigma curve to < 0.03. This is **not** removed
by recalibrating the mean-dispersion trend (§8(b)), because the simulator's mean-variance relation is already
exactly the NB-quadratic form DESeq2 assumes (Var/Var_theory = 0.97–1.02). It is also cured by the test change.

### Change 1 — the FP-rate criterion is restored as binding (decision rule item 1; §8(a))

The empirical permutation-null FP rate (fraction of permutations with ≥ 1 BH rejection ≤ alpha) is
**re-instated as a binding calibration criterion**, co-equal with `lambda_pseudobulk` ∈ [0.9, 1.1]. It was
correct as originally written and must not be replaced by a median-rejection surrogate.

### Change 2 — the power oracle is held at pre-registered parameters (§8(c))

The power oracle is evaluated **only** at log2FC = 1.0, K = 200. Substituting an easier effect size in gate
code is prohibited.

### Change 3 — the gate must fail when the criteria fail

`scripts/synthetic_gate.py` now encodes both restored criteria as pass/fail at the spec's own thresholds. On
the current DESeq2-Wald arm the gate therefore reports **INSTRUMENT NEEDS ATTENTION**, which is the truthful
status of the instrument as of this date.

### Correction to `PILOT_FINDINGS.md` (retraction of a statistics error)

`PILOT_FINDINGS.md` (2026-07-19), section *"A statistics correction I made by hand"*, item 1, claimed that the
criterion *"fraction of permutations with ≥ 1 rejection ≤ alpha"* is a family-wise (FWER) bound that BH does
not control, and used this to retire the FP-rate criterion. **This is wrong and is retracted.**

Under the **complete null** every rejection is false, so V = R, FDP = 1{R ≥ 1}, and FDR = E[FDP] = P(R ≥ 1).
BH controls FDR ≤ alpha, hence P(R ≥ 1) ≤ alpha. FWER and FDR **coincide** under the complete null, and the
spec's original criterion stood. The observed 0.35–0.50 is genuine anti-conservativity of DESeq2-Wald, not an
artifact of BH dependence — the donor-level t-test, using the identical BH over the identical universe on the
identical data, is calibrated at 0.03–0.05.

### Still open — deliberately NOT decided in this amendment

The replacement test is **not** fixed here. A donor-level Welch t-test is demonstrably calibrated, but
`limma-voom` and `edgeR`-QLF — the actual field-standard pseudobulk methods — were **not tested**, and a
moderated empirical-Bayes method may be both calibrated *and* more powerful than a flat Welch t. Fixing the
test before that comparison would leave power on the table. That comparison is the next diagnostic.

### What this does NOT settle (carried to the stratum-list pre-registration)

The real per-stratum `sigma_donor` is **unknown and unanchored** — it is currently a free knob of the
simulator. Under a calibrated donor-level test the binding constraint flips from calibration to **power**, and
the minimum donors per group for power 0.60 at log2FC = 1.0 / K = 200 is a steep function of `sigma_donor`:
n* = 4 / 8 / 13 / 23 at sigma = 0.2 / 0.35 / 0.5 / 0.7 (derived, then validated numerically, |error| < 0.033).

Therefore `sigma_donor` **must** be pinned to a real empirical mean-dispersion / donor-variance trend (§8(b))
**before** any minimum-donor stratum-inclusion rule is pre-registered; otherwise that threshold can be wrong by
a factor of 4–5 and the pre-registration is worthless. A rule of the form "≥ 12 donors per group" is **not**
supported by anything measured here and would certify strata whose FP rate is 3–9× alpha.

If real strata carry `sigma_donor` ≈ 0.5–0.7, then log2FC = 1.0 / K = 200 is unachievable at realistic donor
counts by **any** test, and §8(c) itself — not only the choice of test — must be amended.

*Author attests: the synthetic evidence above is all that was available; no real data informed this amendment.*

---

## Amendment 2 (2026-08-15) — The pseudobulk arm's test becomes moderated eBayes; the paired BH is wired (erratum); B5 and C5 restored, C3 superseded, A2 deferred, the thin-donor filter implemented

Amendment 1 established that the pseudobulk arm fails its binding validity gate under DESeq2-Wald and
**deliberately left the replacement test open**, pending a comparison against the moderated methods
(limma-voom, eBayes) it had not tested. That comparison has since been run in full and committed. This
amendment closes the question, and settles five further items that were claimed in module docstrings but
not implemented.

Spec sections touched: §3 (pseudobulk DE model), §6 (λ), §9 items 5/6/7/8, decision rule item 1 / §8(a)
(the criteria are unchanged; the arm they are applied to changes), §1 inclusion-gate item 2, and
corrections A2, B5, C3, C5. Correction **C2's fairness requirement and §5 are NOT changed** — Change 2 is
an implementation erratum against an unchanged spec.

### Data visible at the time of this amendment (full disclosure)

**The selection data has been seen in full, before the selection rule was written down.** This is the
opposite of the pre-registration ideal and it is stated plainly rather than glossed: the grid was run,
committed, and read; only then was this amendment written. That is precisely why the rule below is
stated as a **deterministic criterion mechanically applied to already-committed data**, evaluated by a
script anyone can re-run, rather than as a judgement made cell-by-cell. A reader who distrusts the
choice can recompute it — the inputs are frozen in git and the criterion is arithmetic.

Visible:

1. **The full 146-cell test-selection grid**, `pilot/testsel/summary.{csv,json}`, committed at
   **`72dec7b`** (2026-07-23, *"grid summary: full run on the PC (146 cells, failed=0)"*). Its axes:
   test ∈ {ttest (Welch), pooled_t, wilcoxon, ebayes, ebayes_trend, voom, wald} × `sigma_het` ∈
   {0.0, 0.4, 0.8, 1.2} × (`sigma_donor`, `n_donors`) tiers × generative arm ∈ {directnb, gamma,
   lognormal}. 24 cells for each of the six donor-level tests, 2 for `wald`, `failed = 0`. Null
   calibration is measured on **fresh** independent simulations with truth = 0 (150–3000 replicates per
   cell), never on the donor-label permutation null — which, as Amendment 1 recorded, is calibrated by
   construction for permutation-invariant tests and therefore carries no evidential weight here.
2. **The analyzer's verdict**, `scripts/analyze_test_selection.py`, reproduced while writing this entry.
3. Everything that was visible for Amendment 1.
4. **No real data.** Oracle (d), the Mathys 2019 anchor (§8(d)), remains untouched. No CELLxGENE stratum
   has been run and no stratum list has been pre-registered. Every number below is synthetic.

### What we found

**(i) The shrinkage hypothesis is NOT supported.** The grid was built to refute a specific worry: that
variance-borrowing methods become anti-conservative as gene-level variance heterogeneity rises (as the
realised prior degrees of freedom `d0` fall). They do not. At the decisive column (`gen_arm=directnb`,
`sigma_donor=0.5`, 8v8, fresh null), FP rate against α = 0.05 across `sigma_het` = 0.0 / 0.4 / 0.8 / 1.2:

| test | | het 0.0 | het 0.4 | het 0.8 | het 1.2 | d0 at het 1.2 |
|---|---|---|---|---|---|---|
| ttest (Welch) | flat | 0.030 | 0.027 | 0.029 | 0.029 | — |
| pooled_t | flat | 0.047 | 0.044 | 0.042 | 0.045 | 0.000 (unmoderated by construction) |
| wilcoxon | flat | 0.000 | 0.000 | 0.000 | 0.000 | — |
| ebayes | shrink | 0.031 | 0.033 | 0.046 | **0.046** | **2.130** |
| ebayes_trend | shrink | 0.032 | 0.032 | 0.045 | **0.048** | **2.243** |
| voom | shrink | 0.030 | 0.032 | 0.044 | **0.041** | **2.226** |

At het = 1.2 the prior is genuinely weak (`d0` ≈ 2.1–2.2 against d = 14 residual df, i.e. shrinkage
factor `d0/(d0+d)` ≈ 0.13 — the moderation is doing real work and is not a no-op), and all three
moderated tests remain calibrated: `ebayes` FP = 0.046 (exact-binomial P(≥ observed | rate = α) = 0.853),
λ = 0.998; `ebayes_trend` FP = 0.048 (P = 0.733), λ = 0.997; `voom` FP = 0.041 (P = 0.993), λ = 0.999.
Every one is inside the pre-registered λ band [0.9, 1.1] and none is significantly anti-conservative.

**(ii) DESeq2-Wald fails on the grid too, independently of Amendment 1's sweep.** Its two grid cells
(directnb, `sigma_donor` 0.5, 8v8) give FP = 0.26 (binomial P = 1.1e-17) at het 0 and FP = 0.56
(P = 5.8e-68) at het 0.8. Amendment 1's diagnosis reproduces on a differently-seeded, independently-run
grid.

**(iii) Only the moderated tests reach the pre-registered power target anywhere in the grid.** Smallest
donors-per-group attaining sensitivity ≥ 0.60 at the pre-registered oracle (log2FC = 1.0, K = 200), over
*calibrated* cells only, at `sigma_het` = 0:

| test | σ_donor = 0.35 | σ_donor = 0.5 | σ_donor = 0.7 |
|---|---|---|---|
| ttest / pooled_t / wilcoxon | > grid | > grid | > grid |
| ebayes / ebayes_trend / voom | **8** | > grid | > grid |

The unshrunken donor-level tests — including the Welch t that Amendment 1 showed to be *calibrated* —
never reach 0.60 anywhere in the grid. Calibration alone was never the criterion; §8(c) is binding too.
Amendment 1's refusal to fix the test on calibration evidence alone is vindicated by this table: fixing
Welch t then would have left the arm calibrated and permanently underpowered.

### Change 1 — the pseudobulk arm's test is replaced (closes the question Amendment 1 left open; §3, §9 item 5)

**The selection rule, stated before its application and applied mechanically:**

> Among the tests **calibrated** — FP rate not significantly above α = 0.05 at *every* `sigma_het` level
> including 1.2, and λ ∈ [0.9, 1.1] — choose the one **maximising worst-case power across the grid**.
> Ties are broken by **runtime**.

*Applying it, with the arithmetic shown.*

**Step 1 — calibration filter.** All six donor-level tests pass (table (i) above). `wald` fails and is
eliminated. No test is eliminated by calibration, so the choice is made entirely on power.

**Step 2 — worst-case power. The rule as literally written is degenerate on this grid, and that is
reported rather than quietly repaired.** Every test — moderated and unshrunken alike — has a cell with
power ≈ 0.000 (the `n_donors = 4`, `sigma_donor = 0.5` corner, where nothing works: the best test in that
corner reaches 0.154). A criterion that is a 6-way tie at zero selects nothing, and would hand the whole
decision to the runtime tie-break, which would be an accident rather than a decision.

The rule is therefore evaluated in its **minimax-regret** form — worst-case power *shortfall from the
best test in the same cell*, which is the same intent (maximise worst-case performance) made
non-degenerate by measuring each test against what was achievable in that cell rather than against an
absolute scale that varies by two orders of magnitude across the grid. Over the 24 fully-populated cells:

| test | worst-case regret | mean regret | cells reaching power ≥ 0.60 |
|---|---|---|---|
| **ebayes_trend** | **0.0454** | 0.0038 | 6 |
| **ebayes** | **0.0470** | 0.0040 | 6 |
| voom | 0.0490 | 0.0058 | 6 |
| pooled_t | 0.2456 | 0.0534 | 5 |
| ttest | 0.2812 | 0.0756 | 4 |
| wilcoxon | 0.3512 | 0.1451 | 3 |

The three moderated tests separate cleanly from the three unshrunken ones by a factor of five, and are
within 0.004 of each other — which is a tie, not a ranking.

**Step 3 — the tie, made explicit.** Paired per-cell power differences over the 24 matched cells:

- `ebayes` − `ebayes_trend`: mean **−0.0002**, max |difference| **0.0061**, 10 wins / 13 losses.
- `ebayes` − `voom`: mean **+0.0018**, max |difference| 0.0159, 18 wins / 5 losses.

`ebayes` and `ebayes_trend` are indistinguishable (a mean difference of 2 parts in 10 000, sign-split
almost evenly); `voom` is very slightly behind both and never ahead on worst-case regret.

**Step 4 — runtime tie-break, with the numbers corrected.** Measured from the grid's own `seconds`
fields, paired cell-by-cell (`voom` ÷ `ebayes` on the identical configuration, n = 24 matched cells):
**median 9.85×, range 1.06×–11.98×**; grid totals `voom` 9933.0 s vs `ebayes` 3390.8 s (2.93× overall,
the ratio being compressed by the expensive tier-C cells where both are dominated by simulation rather
than by the test). *A figure of "≈11× slower, ~850 s vs ~75 s per cell" was in circulation while this
amendment was being drafted; it does not survive contact with `summary.json` and is corrected here. The
correct statement is that `voom` costs an order of magnitude more per cell in the cheap tiers and about
3× more over the grid as a whole.* `ebayes` vs `ebayes_trend`: 3390.8 s vs 3413.1 s over the grid —
`ebayes` is faster by 0.7%, which is within noise and is therefore **not** the deciding consideration
between them.

**Step 5 — the decision.** `voom` is eliminated: equal-or-worse power at ~10× the per-cell cost.
Between `ebayes` and `ebayes_trend`, power ties and runtime ties; the tie is broken toward `ebayes` as
the **strictly simpler estimator** — `ebayes_trend` adds a natural-cubic-spline trend in mean log-CPM
whose only justification would be a power gain, and the measured gain is −0.0002. An unjustified moving
part in the denominator arm of the whole study is a liability, not a feature.

**Selected: `ebayes`.** From this amendment forward the pseudobulk arm's test is:

> **Moderated t with empirical-Bayes shrunk variances (limma-style; Smyth 2004,
> *Stat. Appl. Genet. Mol. Biol.* 3, Article 3), on log2(CPM + 1) donor pseudobulk profiles.**
>
> Per gene, ordinary least squares for the two-group design `~ 1 + x` on the donor × gene matrix of
> log2(CPM + 1) — CPM taken against the row sums of the **universe-restricted** matrix, so no gene
> outside the frozen universe leaks into the normalisation — giving the contrast estimate `beta_g`, the
> residual variance `s2_g` on `d = n_donors − 2` degrees of freedom, and the unscaled variance
> `(X'X)^-1[1,1]`. A scaled-F prior (`s2_g | sigma2_g ~ sigma2_g · chisq_d/d`, `d0·s02/sigma2_g ~
> chisq_d0`, so marginally `s2_g ~ s02·F(d, d0)`) is fitted by Smyth's `fitFDist` method of moments on
> `log s2` via digamma/trigamma. The posterior variance is
> `s~2_g = (d0·s02 + d·s2_g)/(d0 + d)`, the statistic is `t_g = beta_g / sqrt(s~2_g · unscaled_var_g)`,
> referred to `t_{d + d0}`. `d0 = +inf` (complete pooling — the between-gene variance not separable from
> sampling noise) is handled in the limit `s~2 → s02`, `t_{inf} → N(0,1)`, and is reported, never
> silently clipped. **No mean-variance trend** (`covariate = None`): the trended variant was measured
> and rejected in step 5. Genes with non-finite `s2` yield `p = NaN` and are dropped by BH rather than
> given a p-value derived entirely from the prior.

This is a **pure-Python re-implementation of published estimators**, not a wrapper: there is no R, rpy2,
limma or edgeR on the machine and pbcheck targets the scanpy/scverse ecosystem, so an R dependency is
not available. The implementation being adopted is the *same code that produced this grid* — it is
lifted out of `scripts/pb_calibration_probe.py` into `src/pbcheck/methods/moderated.py` unchanged, with
a test pinning the extraction against the probe's own output, so the arm the study ships is numerically
the arm the selection evidence was measured on. It carries the probe's four numerical self-validations
(d0 → 0 reduces to the ordinary t; homoscedastic-null calibration by KS against `t_{d+d0}`; a power gain
over the ordinary t in the homoscedastic case; recovery of a known d0).

**DESeq2-Wald is retired from the arm but its implementation is retained.** `deseq_from_pdata` stays in
the tree, documented as superseded, so that every number in Amendment 1 and in `PILOT_FINDINGS.md`
remains reproducible from this repository. Retiring a method by deleting the code that produced the
retired results would make the record unfalsifiable.

Spec §3's `DeseqDataSet` / `DeseqStats` block, and C1 (`poscounts` size factors) and C4 (covariate df
rule) with it, no longer describe the arm's test. C1 and C4 are DESeq2-specific and lapse with it; the
log2(CPM + 1) transform normalises by library size directly, which is the moderated arm's counterpart to
a size factor. C2's requirement — no asymmetric NA between the arms — **survives and binds**, and is now
enforced by Change 2 rather than by a DESeq2 flag.

### Change 2 — fairness erratum: the paired BH was specified, implemented, tested, and never wired (§5, C2)

This is **not a spec change**. Spec §5 requires both arms to be BH-corrected over one identical tested
set, and item 4 requires that the two BH input vectors be asserted to be of the same length with
identical membership. `src/pbcheck/mtc.py` implements exactly that as `bh_both_arms` / `PairedBH`, it has
tests, and its own docstring (`mtc.py:73-76`) instructs *"Use `bh_both_arms` for anything whose output is
compared across arms."*

**No production caller ever did.** `permutation.run_null` (both arms, `permutation.py:126` and `:144`)
and `scripts/synthetic_gate.py` (`:35-36`) both call the per-arm `bh_over_universe`. Every gate number
published to date — including the ones in Amendment 1, `PILOT_FINDINGS.md`, and the README table — was
computed with **each arm BH-corrected over its own non-NaN subset**, not over one common set.

The direction of the error is known and is against the study's own thesis in one place and for it in
another: the naive Wilcoxon arm essentially never returns NaN, while the pseudobulk arm does (DESeq2 on
degenerate strata; the moderated arm on genes with non-finite `s2`). The pseudobulk arm was therefore
corrected over a **smaller** m — made per-gene systematically *more liberal* than the arm it is the
control for — and nothing reported it. Since the pseudobulk arm is the denominator of every inflation
number, this is not cosmetic.

Henceforth `run_null` and the gate use `bh_both_arms`, and its bookkeeping (`n_tested_common`,
`n_na_naive`, `n_na_pseudobulk`, `dropped_for_fairness`) is carried into the output rather than
discarded. **All gate numbers are recomputed under the paired BH and will shift.** Prior published
numbers are superseded by the recomputation, not by an argument that the shift is small — that has not
been measured yet at the time of writing.

### Change 3 — B5 restored, and its literal construction measured and rejected (§6, B5)

B5 pins that λ be computed from *empirical* permutation-null p-values — "rank of the observed statistic
within the donor-permutation null" — rather than from analytic p-values, "so lambda reflects clustering
miscalibration, not tie/discreteness artifacts". `metrics.py:39-45` cites `(spec B5)` while computing λ
from analytic p-values; the citation was an overclaim (what the function actually does is median λ
*across* permutations, which mitigates single-realization discreteness but is not B5's construction).

**B5's construction as literally written cannot serve the purpose §6 assigns to λ, and this was
measured, not argued.** Under the donor-permutation sharp null the real labeling is exchangeable with
every permuted one, so the rank of the real statistic within the permutation null is uniform *by
construction* — for any test, however miscalibrated. Measured on the null oracle (G = 400 surviving
genes, 8v8 donors, `sigma_donor` = 0.5, 120 cells/donor, dispersion 0.2, seed 1, 30 permutations):

| λ_naive computed as | value |
|---|---|
| analytic p-values, real labels | 22.93 |
| analytic p-values **under the permutation null** (what the gate currently uses) | **26.08** |
| **B5's literal empirical permutation p-values** | **0.93** |
| the same construction applied to a *held-out permuted* labeling | 0.93 (median over 30) |

The last row is the proof: B5's λ returns the identical value whether it is fed the real labeling or an
arbitrary permuted one, i.e. it carries **zero information** about the arm's calibration. It measures
the exchangeability of the permutation set, not the inflation of the test.

**Resolution.** The binding λ criterion (decision rule item 1, §8(a)) remains λ computed from each arm's
own p-values **under the permutation null** — which is, on reflection, exactly the clustering
miscalibration B5 wanted: the analytic p-value is (asymptotically) the *cell*-permutation p-value, so
its λ against the *donor*-permutation null measures precisely how much wider the true donor-level null
is than the cell-level null the naive test assumes. B5's tie/discreteness concern is separately
addressed by `tie_correct=True` (pinned by §2, wired in R0) and by taking the median across
permutations. The empirical construction is **implemented and reported as a named companion diagnostic**
— a validity check on the permutation machinery, which must sit near 1 and whose departure from 1 would
indicate a non-exchangeable permutation set — and is explicitly **not** a calibration criterion. The
`(spec B5)` citations in `metrics.py` are corrected to say what the code does.

### Change 4 — C5 restored, with a moderated analog of `fitType` (§1 item 5, §3, C5)

C5 requires a minimum frozen-universe gene count enforced before the pseudobulk test runs, and the
`fitType` actually used to be persisted per stratum. `frozen_universe(min_size=...)` exists and is
tested, but **no caller passes it** (`scripts/synthetic_gate.py:48` and the probe both call it bare), so
the gate was decorative in exactly the sense C5 was written to prevent.

1. All callers now pass `min_size` (default 200 genes, per inclusion-gate item 5).
2. `fitType` is DESeq2's dispersion-trend fitting mode and lapses with DESeq2. Its **moderated analog is
   the realised prior**, and that is what is now persisted in the gate output in its place: the fitted
   prior degrees of freedom `d0` (with the `prior_is_complete_pooling` boolean that disambiguates
   `d0 = +inf` from a missing value), the shrinkage factor `d0/(d0+d)`, the residual df `d`, whether a
   mean-variance trend was used (`s0_squared_is_trended`, `False` for the selected arm), `s0^2`, and the
   count of genes rendered untestable by a non-finite `s2`. These are the moderated arm's equivalent of
   "which model actually got fitted here", and a run whose `d0` collapsed or whose prior went to complete
   pooling is now visible in the artifact instead of being invisible.

### Change 5 — C3 is SUPERSEDED (retired with rationale) (C3)

C3 requires the pseudobulk DEG count to be reported **both ways** — DESeq2's native
independent-filtered `padj` and re-BH over the frozen universe — to show the choice changes neither
calibration nor the GO decision. Independent filtering is a DESeq2 feature. With DESeq2 retired from the
arm (Change 1) there is no native filtered `padj` to cross-check against, and the check has no referent.
C3 is therefore **retired, not deferred**: it is not work outstanding, it is a check whose subject no
longer exists.

Two things are put in its place, both cheap, and both serving C3's actual purpose (that the arm's
significance calls not be an artifact of one p-value-processing choice):

1. **Prior-strength disclosure** — the realised `d0` and shrinkage factor from Change 4. The moderated
   analog of "how much did the filtering choice matter" is "how much did the prior matter", and it is
   now reported rather than assumed negligible.
2. **A p-value uniformity check on the permutation null** — the moderated arm's null p-values are
   tested against Uniform(0,1), reported alongside λ. Where independent filtering would have distorted
   the null, this shows it directly.

Note that C2's motivation for `cooks_filter=False` — no asymmetric NA between the arms — is unaffected
and is now carried by Change 2's paired BH, which handles NAs from *any* source rather than from
DESeq2's filters specifically. That is a strictly stronger guarantee than the flag it replaces.

### Change 6 — A2 is DEFERRED to Phase 1, honestly (§4, A2)

A2 requires the permutation set to be **cell-count stratified** — restricted to (or reweighted toward)
label assignments whose per-group total cell count is close to the real split. What
`permutation.py:148-163` actually does is log per-group cell totals and check that the real split lies
inside the permutation range (`real_split_inside_perm_range`, `real_split_percentile_in_perms`). That is
a range check, not stratification, and the code comment labelling the block `---- A2:` reads as though
the correction were satisfied. **It is not, and the comment is corrected to say so.**

Full stratification is deferred to Phase 1, with the rationale stated so the deferral can be argued
with:

- **It is real statistical work, not wiring.** Restricting or reweighting a permutation set on a
  covariate changes the null distribution being sampled and requires deciding the matching tolerance,
  what to do when the restricted set is empty, and whether the resulting null is still exact. Doing it
  by feel inside an engineering pass is how a pre-registered study acquires an unrecorded protocol
  change.
- **It is infeasible at the donor counts the spec's own inclusion gate admits.** At the 3v3 minimum
  there are C(6,3) − 2 = 18 permutations in total; filtering them on cell-count proximity leaves single
  digits, and the floor's Monte-Carlo error would swamp the confound it is meant to remove.
- **The confound it targets is asymmetric between the arms, and the pseudobulk arm is the one that
  matters here.** Donor pseudobulk profiles are aggregated once and are label-invariant: a permutation
  changes only the donor→condition map, so per-group cell-count imbalance cannot move the pseudobulk
  arm's inputs at all. It is the naive per-cell arm whose statistic scales with cell count. Since the
  binding validity gate this amendment exists to repair is about the *pseudobulk* arm, A2 does not gate
  it. It does bear on the naive arm's floor, which is a Phase 1 headline quantity — hence Phase 1, not
  "never".
- The range check is retained and reported meanwhile, and is what it says it is: a check that the real
  split is not outside the permutation cell-count distribution. Where it fails, the floor comparison is
  size-confounded and must be read as such.

**A2 is recorded here as partially addressed. Any docstring or comment claiming otherwise is corrected
in the same change.**

### Change 7 — the thin-donor filter (`min_cells` / `min_counts`) is implemented manually (§1 item 2, §3)

Spec §1 inclusion-gate item 2 pins "≥ 10 cells per donor in that cell_type (donors below threshold are
**dropped, not merged**)" and §3 pins the aggregation call
`dc.pp.pseudobulk(..., min_cells=10, min_counts=1000)`. R0 documented the resulting hole: **decoupler
2.x removed those parameters**, so `build_pseudobulk` accepted `min_cells` / `min_counts` and silently
never applied them. The pre-registered filter has never run.

Implementing an equivalent filter ourselves is a protocol-affecting change, which is why R0 stopped at
documenting it and deferred it here. **Decision: implement it.** After aggregation, a donor × cell_type
pseudobulk profile is **dropped** — never merged into another donor, never back-filled — when it derives
from fewer than `min_cells` cells or when its summed counts are below `min_counts`. Defaults stay at the
pre-registered 10 and 1000. The number of profiles dropped and the reason is recorded, so a stratum
losing donors to the filter is visible rather than silently smaller; the existing "≥ 3 donors per group
post-aggregation, else SKIP" rule (§3) then applies to what survives.

This is a **re-implementation of the pre-registered intent in a library that no longer offers it**, not
a new criterion: the thresholds, the semantics ("dropped, not merged") and the placement (after
aggregation, before the universe is frozen) are all as §1/§3 already specify.

### What this does NOT settle

- **`sigma_donor` remains unanchored to real data.** Amendment 1's closing concern stands in full and is
  *not* relaxed by the test change. The power frontier above makes it sharper, not softer: the selected
  arm reaches power 0.60 at 8 donors per group when `sigma_donor` = 0.35, and **nowhere in the grid** at
  `sigma_donor` = 0.5 or 0.7. If real strata carry `sigma_donor` ≈ 0.5–0.7 then the pre-registered
  oracle (log2FC = 1.0, K = 200) is unreachable at realistic donor counts by *any* test in the grid, and
  §8(c) itself — not the choice of test — is what must be amended. `sigma_donor` must still be pinned to
  a real empirical mean-dispersion / donor-variance trend (§8(b)) before any minimum-donor
  stratum-inclusion rule is pre-registered.
- **A2 stratification is deferred, not solved** (Change 6). The naive arm's floor remains a
  cell-count-confounded quantity guarded only by a range check.
- **The real-data anchor is still untouched.** Oracle (d), Mathys 2019 (§8(d)), has not been run. Every
  number in this amendment — every FP rate, every λ, every power figure, and the entire selection
  argument — comes from our own generative model. A moderated test calibrated on synthetic NB data with
  a log-normal donor random effect is *not* thereby shown calibrated on real snRNA-seq, where the
  mean-variance relation, zero inflation and donor heterogeneity are all ours to get wrong. The binding
  check remains the real anchor, and it is still pending.
- **The selection rule was applied to data already seen** (opening section). Re-running the grid with
  different seeds would be a genuine, and so far unperformed, robustness check on the choice.
- **Whether the arm passes its gate under the paired BH is not known at the time of writing.** Change 2
  shifts every gate number; the grid evidence predicts the moderated arm is calibrated at these regimes,
  but the gate has not yet been run under the combined changes. If it fails, that is a result and it
  comes back here — no threshold in `PHASE0_SPEC.md` is to be touched to make it pass.

*Author attests: the synthetic evidence above is all that was available; no real data informed this
amendment. The selection data was seen in full before the rule was written, which is disclosed above and
is why the rule is arithmetic and re-runnable rather than discretionary.*

---

## Amendment 3 (2026-08-15) — §8(c)'s power criterion is re-scoped from a point to a declared operating envelope; the gate's permutation count is raised 40 → 200

Amendment 2 replaced the pseudobulk arm's test with moderated eBayes. The gate rerun that followed
**passed every calibration criterion and failed the power criterion** — sensitivity 0.35 against the
binding ≥ 0.60 — at the gate's operating point (`sigma_donor` = 0.5, 8 v 8 donors). Amendment 1 had
already named that outcome and its consequence, in its closing paragraph:

> If real strata carry `sigma_donor` ≈ 0.5–0.7, then log2FC = 1.0 / K = 200 is unachievable at realistic
> donor counts by **any** test, and §8(c) itself — not only the choice of test — must be amended.

The committed grid has since shown that no test in it reaches power 0.60 at `sigma_donor` = 0.5 at any
donor count the grid tested. This amendment is the entry Amendment 1 pointed at. It is written and
committed **before** any of the code that applies it, and before the rerun whose numbers it changes.

Spec sections touched: **§8(c)** — the *scope* over which the synthetic-positive criterion binds. Its
effect size (log2FC = 1.0, K = 200) is **UNCHANGED**, and Amendment 1 Change 2's prohibition on
substituting an easier oracle stands unmodified and is not weakened by anything below.
**Decision rule item 1 / §8(a)** — the *power* half of the pseudobulk validity gate; the calibration half
is untouched, stays binding, and is now explicitly pinned to the hardest regime available.
**§1 (inclusion gate)** — a new per-stratum requirement for the real sweep. Change 2 additionally alters
`gate_config`'s **not**-pre-registered instrument-sanity block, and is labelled as such there and here.

### Data visible at the time of this amendment (full disclosure)

Nothing new was measured for this amendment. It is written against evidence already committed, all of it
re-derived rather than quoted from memory:

1. **The 2026-08-15 gate run** under the Amendment 2 arm (`scripts/synthetic_gate.py` at defaults: seed 1,
   40 naive permutations of which 40 paired, 1500 genes, 8 v 8 donors, 250 cells/donor, dispersion 0.2,
   `sigma_donor` = 0.5), recorded in `docs/PILOT_FINDINGS.md`'s addendum (`1e26cec`), in the README status
   table (`5b6ad53`) and in `pilot/README.md` (`b53f044`). **No JSON artifact was written for that run** —
   `--out` was not passed — so the surviving record of it is prose. It was therefore **re-run verbatim
   while writing this entry**, on the same commit at the same defaults, and every published figure
   reproduced exactly: `lambda_naive` 54.77, `lambda_pseudobulk` 1.02, perm-null FP rate 0.05 (MC SE
   0.034), naive floor 1164 / 1500 (77.6 %), pseudobulk floor 0 (mean 0.05), power 0.35, verdict
   `INSTRUMENT NEEDS ATTENTION`, in 119.3 s against the recorded 117.9 s. The prose record is therefore
   sound and is treated as data here.
2. **The full 146-cell test-selection grid**, `pilot/testsel/summary.{csv,json}` at `72dec7b`, read by
   `scripts/analyze_test_selection.py`, which was re-run while writing this entry to re-derive the power
   frontier rather than copy it out of Amendment 2.
3. **Amendment 1's power frontier**: `n* = 4 / 8 / 13 / 23` donors per group at `sigma_donor` =
   0.2 / 0.35 / 0.5 / 0.7, "derived, then validated numerically, |error| < 0.033" — quoted verbatim from
   Amendment 1's closing section, checked against its text, and used below as a load-bearing input.
4. Everything that was visible for Amendments 1 and 2.
5. **No real data.** Oracle (d), the Mathys 2019 anchor (§8(d)), remains untouched. No CELLxGENE stratum
   has been run and no stratum list has been pre-registered. Every number below is synthetic, and every
   `sigma_donor` in it is a **free knob of our own simulator**, not a measurement of anything real. That
   fact is the entire reason Change 1 takes the form it does.

### What we found

**(i) The instrument is calibrated, and underpowered at exactly one point on one axis.** The 2026-08-15
run, at `sigma_donor` = 0.5 / 8 v 8:

| criterion (decision rule item 1) | reading | required | |
|---|---|---|---|
| `lambda_pseudobulk` | **1.02** | ∈ [0.9, 1.1] | PASS |
| pseudobulk perm-null FP rate | **0.05** (2 / 40 perms) | ≤ α = 0.05 | PASS, marginally — see (iv) |
| pseudobulk permutation floor | **0** median #DEG | ≈ 0 | PASS |
| `lambda_naive` (instrument sanity) | **54.77** | > 1.5 | PASS |
| naive false-positive floor | **1164 / 1500 genes (77.6 %)** | > 30 % | PASS |
| pseudobulk power at log2FC = 1.0, K = 200 | **0.35** | ≥ 0.60 | **FAIL** |

Five of six criteria pass and the sixth fails by a wide margin at one setting of one simulator knob.

**(ii) It is not a test-choice problem, and it is not an engineering gap.** From the committed grid
(`ebayes`, the selected arm, `sigma_het` = 0, generative arm `directnb`, power = mean sensitivity over the
probe's 25 positive-oracle replicates at the pre-registered oracle):

| `sigma_donor` | donors/group | `ebayes` power | calibrated? |
|---|---|---|---|
| 0.35 | 8 | **0.793** | yes (FP 0.043, λ 1.010) |
| 0.5 | 4 | 0.009 | yes (FP 0.029, λ 1.016) |
| 0.5 | 8 | 0.194 | yes (FP 0.031, λ 1.010) |
| 0.5 | 12 | 0.486 | yes (FP 0.034, λ 1.010) |
| 0.7 | 8 | 0.003 | yes (FP 0.015, λ 1.020) |

and the analyzer's frontier over *all* seven tests, recomputed on calibrated cells only:

| test | σ_donor = 0.35 | σ_donor = 0.5 | σ_donor = 0.7 |
|---|---|---|---|
| ttest / pooled_t / wilcoxon | > grid | > grid | > grid |
| ebayes / ebayes_trend / voom | **8** | > grid | > grid |

At `sigma_donor` = 0.5 **no test in the grid clears 0.60 at any donor count the grid tested** (the largest
was 12 v 12). The threshold the gate is failing is not reachable there by changing the test, by changing
the implementation, or by any means available inside the study except adding donors or changing the
regime. Amendment 1 said this in advance; the grid confirms it.

**(iii) Two independent derivations of the donor frontier agree.** Amendment 1's analytic frontier
(n\* = 4 / 8 / 13 / 23 at σ = 0.2 / 0.35 / 0.5 / 0.7) was derived from the power algebra and validated
numerically to |error| < 0.033. The grid, run later, on a differently-seeded and differently-generated
arm, brackets it from both sides where it has cells:

* at σ = 0.35 the grid's smallest tested count, 8, already reaches 0.793 — so the grid alone gives
  n\* ≤ 8, and Amendment 1's derivation puts it at exactly 8;
* at σ = 0.5 the grid's largest tested count, 12, reaches only 0.486 — so the grid alone gives n\* > 12,
  and Amendment 1's derivation puts it at 13;
* at σ = 0.7 the grid's 8 v 8 reaches 0.003, consistent with a frontier far above 8 (derived: 23).

Neither source is real data. What they establish jointly is that the frontier is a **steep, reproducible
function of a knob we have never measured** — which is the fact Change 1 acts on.

**(iv) The FP criterion passed at the edge of the run's resolution.** 0.05 is 2 rejecting permutations out
of 40. The Monte-Carlo standard error of a binomial rate at p = 0.05 over n = 40 is
√(0.05 · 0.95 / 40) = **0.034**. A reading of 0.05 at that resolution is not distinguishable from a true
rate of ~0.12. Calling the calibration half of the validity gate "met" on this evidence, while
simultaneously narrowing the instrument's claimed domain on power evidence, would be applying two
different standards of proof in the same amendment. Change 2 fixes the resolution instead.

### Change 1 — §8(c)'s power criterion is re-scoped from a point to a declared operating envelope (§8(c); decision rule item 1 / §8(a); §1)

**The problem, stated precisely.** §8(c) pins an effect size (log2FC = 1.0, K = 200) and a threshold
(sensitivity ≥ 0.60). It does **not** pin the regime at which that threshold is evaluated — the
`sigma_donor` and donor count of the oracle are nowhere in the frozen spec. The gate supplied them from
`gate_config.ORACLE_SIM`, whose own docstring already says `donor_sigma` "is a FREE KNOB of the simulator
and is **not** anchored to real data". So the instrument is currently declared invalid against a
threshold that is unreachable by any known test, at a knob setting chosen by us, which corresponds to no
measured property of any real stratum. That is not a measurement of the instrument; it conflates
**instrument validity** with **domain of applicability**.

The standard measurement-science treatment of exactly this situation is to state an operating range. A
thermometer is not invalid because it cannot read 2000 K; it is a thermometer with a stated range, and
using it outside that range is the user's error, not the instrument's. What is *not* acceptable is
leaving the range unstated, or discovering it after the fact from the data one wishes to include.

**The re-scoping.** From this amendment forward:

> The pseudobulk arm is declared **valid for strata whose (`sigma_donor`, donors-per-group) lies inside
> the operating envelope** — the region in which the selected test's power at the **unchanged**
> pre-registered oracle (log2FC = 1.0, K = 200) is ≥ 0.60. It is **not** declared valid outside it, and
> no result from a stratum outside it may be reported as a pbcheck measurement.

Concretely, the gate henceforth applies decision rule item 1 as three separable parts:

**(a) The calibration criteria are unchanged, binding, and evaluated at the hard regime.** λ ∈ [0.9, 1.1]
and perm-null FP ≤ α continue to be evaluated at `sigma_donor` = 0.5 — the *worst* donor variance in the
grid's fully-populated tier, and the conservative choice for a false-positive criterion, since the null
gets harder as donor variance rises. Nothing about calibration is relaxed, moved, or re-scoped. This is
the half of the gate the instrument passes, and it goes on being tested where it is hardest.

**(b) The power criterion becomes binding at the envelope boundary point `sigma_donor` = 0.35, 8 v 8.**
This is the *lowest-σ / smallest-n* point at which the selected arm is grid-shown to clear 0.60
(measured: 0.793 at 8 v 8, calibrated, §(ii) above) and it coincides with Amendment 1's analytic
n\*(0.35) = 8. Evaluating power there tests the claim actually being made — "this instrument delivers
≥ 0.60 sensitivity inside its stated envelope" — rather than a claim nobody is making.

**(c) The gate REPORTS the envelope, in its console output and in its JSON artifact.** The envelope is
not a footnote to be lost; it is a first-class output of every run, so that no reader can take
"INSTRUMENT VALID" as unconditional. The declared envelope, from Amendment 1's frontier with the grid's
corroboration attached:

| `sigma_donor` | minimum donors per group | source |
|---|---|---|
| 0.2 | **4** | Amendment 1 frontier (derived, validated \|err\| < 0.033); not in the grid |
| 0.35 | **8** | Amendment 1 frontier; grid: `ebayes` power 0.793 at 8 v 8 → n\* ≤ 8 |
| 0.5 | **13** | Amendment 1 frontier; grid: `ebayes` power 0.486 at 12 v 12 → n\* > 12 |
| 0.7 | **23** | Amendment 1 frontier; grid: `ebayes` power 0.003 at 8 v 8 |

**This NARROWS the instrument's claimed validity domain. It does not lower any bar.** The distinction is
the whole content of this change and is stated plainly so that it cannot be read the other way:

* The effect size is untouched: log2FC = 1.0, K = 200, exactly as §8(c) froze it, exactly as Amendment 1
  Change 2 forbade changing.
* The power threshold is untouched: ≥ 0.60.
* **Power 0.60 at `sigma_donor` = 0.5 with 8 donors per group remains UNMET and is NOT claimed.** The
  measured 0.35 stands on the record. Nothing here converts it into a pass.
* A stratum with `sigma_donor` ≈ 0.5 will require **≥ ~13 donors per group** or must be **excluded**.
  Both of those are stricter constraints on what pbcheck may report than existed before this amendment —
  previously the study had no donor-count rule at all beyond the inclusion gate's 3 v 3 minimum, which
  Amendment 1 already warned "would certify strata whose FP rate is 3–9× alpha" if used naively.
* Consequently, the honest one-line summary of the instrument's status after this amendment is
  *"valid within a stated envelope that most real strata may well fall outside"*, not *"valid"*.

**A disclosed extrapolation in (b).** The grid's (σ = 0.35, 8 v 8) cell is on the `directnb` generative
arm — counts simulated directly at the donor level. The gate runs the *cell-level log-normal* oracle
(`synthetic/oracles.py`), for which the grid has **no** σ = 0.35 cell (its lognormal cells exist only at
σ = 0.5). The two arms do not give the same power: at σ = 0.5 / 8 v 8 the grid measures `ebayes` power
0.194 on `directnb` and **0.4006** on `lognormal`, the latter being consistent with the gate's own
single-realisation 0.35. The lognormal arm is therefore the *more* favourable of the two at the one
regime where both were measured, which is why ≥ 0.60 at σ = 0.35 is expected to hold there — but it is an
**extrapolation across generative arms, not a measurement**, and it is recorded as such here, before the
rerun, so that a failure is a result rather than a surprise. If the rerun does not clear 0.60 at
σ = 0.35, the envelope boundary is wrong and comes back to this log; no threshold is to be moved to
accommodate it.

**Two further resolution caveats on the evidence cited above**, disclosed because they bound how hard the
envelope numbers can be pushed:

* The gate's power figure is a **single positive-oracle realisation**; the grid's are means over the
  probe's 25 replicate datasets. Single-realisation sensitivity at K = 200 carries real Monte-Carlo
  error, and the gate's 0.35-vs-grid-0.4006 agreement should be read with that in mind.
* `summary.json`'s `power_sd` and `n_power_reps` columns are **`None` in every row**, because
  `scripts/run_test_selection_grid.py` reads keys (`power_sd`, `n_power_reps`) that
  `scripts/pb_calibration_probe.py` does not write (it writes `power_sd_across_datasets` and
  `power_n_reps`). The power *means* are correct and are what is cited; their per-cell Monte-Carlo error
  is simply not recoverable from the committed summary. This is a reporting defect in the grid driver,
  noted here rather than silently worked around, and it is not fixed by this amendment because the grid
  is frozen evidence.

**Consequence for Phase 0 real-data stratum selection (part of this change, recorded now rather than
when it becomes convenient).** Stratum inclusion in the real sweep now requires, in addition to
everything §1 already pins:

1. a **per-stratum estimate of `sigma_donor`**, and
2. **membership in the operating envelope** at that stratum's own donors-per-group.

The **mechanism** for (1) exists and needs no new machinery: the moderated arm's own fit already produces,
per gene, the residual variance `s2_g` of log2(CPM + 1) **across donors within group** on `d = n − 2`
degrees of freedom, and the fitted prior location `s0^2` is its shrunken typical value — both already
persisted in the gate artifact as of Amendment 2 Change 4. Since the simulator's `donor_sigma` is the
standard deviation of a per-(gene, donor) log-normal random effect on the **natural-log** scale
(`synthetic/oracles.py`: `exp(N(−σ²/2, σ))`), a stratum's between-donor dispersion on the log2 scale
converts to that parameterisation by a factor of ln 2.

**That conversion is a mechanism, not an anchor, and the difference is load-bearing.** `s2_g` contains the
donor random effect **plus** residual NB sampling noise and the compression of the `+1` offset at low
expression, so `sqrt(s0^2) · ln 2` is an **upper bound** on `donor_sigma`, not an estimate of it. Deriving
and validating the correction — against the simulator, where the truth is known — is required work that
this amendment does **not** do and does **not** authorise skipping. **The `sigma_donor` anchoring demanded
by Amendment 1 therefore remains OPEN.** This amendment supplies the instrument's stated range and the
mechanism by which a stratum could be tested against it; it does not supply the empirical anchor, and no
stratum may be admitted to the real sweep on the strength of this entry alone.

### Change 2 — the gate's permutation count is raised 40 → 200 (instrument-sanity block; NOT pre-registered)

Declared here, **before** the rerun that will apply it, so that the resulting FP number is a
pre-committed measurement rather than a chosen one.

`gate_config.N_PERM_PB` goes 40 → 200. The Monte-Carlo standard error of the FP criterion at p = 0.05
falls from √(0.05 · 0.95 / 40) = **0.034** to √(0.05 · 0.95 / 200) = **0.015**, which separates 0.05 from
0.12 — the ambiguity §(iv) records. The marginal reading gets **resolved** instead of leaned on.

**`N_PERM` is raised 40 → 200 in the same change, and this is not incidental.**
`permutation.run_null` computes `n_paired = min(len(perms[:n_perm]), n_perm_pb)`, so raising `n_perm_pb`
alone would have changed **nothing at all** — the paired count, and with it the FP rate and its MC SE,
would have stayed pinned at 40 by `n_perm`. Raising only the named constant would have produced a run
that looked like it had 200 permutations and did not. Both constants move together, and the reason is
recorded here so that the pair cannot later be separated by someone reading only the headline.

**This is affordable only because of Amendment 2.** Under the retired DESeq2-Wald arm a fit cost ~2.8 s
against the moderated arm's ~5 ms (measured in the grid driver's own notes), so 200 paired permutations
would have cost hours and the resolution would have been unaffordable rather than merely unmeasured.
Amendment 2 was selected on calibration, power and runtime; this is the runtime dividend being spent on
resolution.

Both constants live in `gate_config`'s **`INSTRUMENT_SANITY`** block, which is **not** pre-registered and
is labelled so in the artifact's manifest. Spec §4 pins `n_perm = 1000` and `n_perm_pb ≥ 200` for the
**real** sweep; this change brings the synthetic gate's paired count up to that floor but claims no more
than that. It changes no threshold and no criterion — only the resolution at which an existing criterion
is measured.

**Stated in advance, so that it binds:** if the FP rate at 200 permutations resolves **above** α = 0.05,
the calibration half of the validity gate **fails**, the instrument is not valid, and that failure is the
information this change exists to surface. Nothing is to be tuned, no permutation count is to be walked
back, and the finding comes back to this log.

### What this does NOT settle

* **`sigma_donor` is still not anchored to real data.** This is the third consecutive amendment to close
  on it and it is *not* resolved here. Change 1 supplies the envelope and the estimation mechanism; it
  does not supply the estimate, the validated conversion from `s0^2` to `donor_sigma`, or the empirical
  mean-dispersion / donor-variance trend §8(b) asks for. Until that exists, the envelope is a statement
  about our simulator's coordinates, and **whether any real stratum falls inside it is unknown**.
* **Whether the real sweep is feasible at all.** If real strata cluster at `sigma_donor` ≈ 0.5–0.7, the
  envelope admits them only at ≥ 13–23 donors per group. Whether enough CELLxGENE strata clear that is an
  open empirical question, and a negative answer is a live outcome of this study, not a failure mode to
  be designed around.
* **A2 stratification remains deferred** (Amendment 2 Change 6). The naive arm's floor is still a
  cell-count-confounded quantity guarded only by a range check.
* **The real-data anchor is still untouched.** Oracle (d), Mathys 2019 (§8(d)), has not been run. Every
  number in this amendment is from our own generative model, whose mean-variance relation, zero inflation
  and donor heterogeneity are all ours to get wrong.
* **The GO/NO-GO decision is not taken**, and nothing here moves it. A gate that passes within a declared
  envelope licenses the *measurement*, not the conclusion.
* **The gate has not been re-run under these changes at the time of writing.** The envelope boundary
  point is an extrapolation across generative arms (Change 1) and the FP rate at 200 permutations is
  unmeasured (Change 2). Either may fail. If either does, it is a result and it comes back here — no
  threshold in `PHASE0_SPEC.md`, and no number in this entry, is to be touched to make it pass.

*Author attests: the synthetic evidence above is all that was available; no real data informed this
amendment. Every figure quoted was re-derived from committed artifacts or re-run while writing, not
copied from an earlier entry; where a number could not be verified — the grid's per-cell power
Monte-Carlo error — that is said rather than glossed.*

---

## Amendment 4, Part A (2026-08-16) — a per-stratum `sigma_donor` estimator and an operating-envelope membership rule: the derivation, the estimator, the gating rule, and PRE-DECLARED validation criteria; Amendment 3's "upper bound" is corrected

Amendment 3 Change 1 declared an operating envelope and made stratum inclusion in the real sweep
conditional on two new things: a **per-stratum estimate of `sigma_donor`**, and **membership in that
envelope** at the stratum's own donors-per-group. It supplied a *mechanism* for the first —
`sqrt(s0^2) · ln 2` from the moderated arm's own fit — declared the mechanism unvalidated, and closed
by recording that

> Deriving and validating the correction — against the simulator, where the truth is known — is
> required work that this amendment does **not** do and does **not** authorise skipping.

This entry is that work. It is split into two parts, and the split is the point.

### Why this entry is in two parts, and what binds when

**Part A — this section — is written and committed BEFORE the validation grid is run.** It fixes the
estimand, the derivation, the estimator, the aggregation rule, the gating rule, the validation grid,
and the numeric criteria **V1–V9** that decide whether the estimator is fit to gate stratum
admission. **Part B is a dated addendum written after the run**: the V1–V9 outcome table, the
mechanical application of the functional-selection rule with its arithmetic, the resulting constant
in `gate_config`, and the consequences of any failure.

The ordering is not a formality. Amendment 1 exists because `scripts/synthetic_gate.py` had quietly
substituted an easier oracle and dropped two binding criteria, and `gate_config`'s own docstring
names that episode as the reason its thresholds live in one file instead of scattered as literals. A
validation whose PASS criteria are written after its numbers are read is not a validation; it is a
description. Every threshold below is therefore stated now, with its arithmetic. **No threshold below
is to be moved after the run.** A criterion that fails is a result, is reported as a failure, and
comes back to this log — the formula Amendments 2 and 3 both closed on.

Two consequences of the split are recorded so they cannot be dropped quietly. First, the constants
below may enter `gate_config`'s `PRE_REGISTERED` block as a new `SIGMA_GATE` group only **after** this
Part A is committed, and never in the same commit as the run that uses them. Second, the functional
slot of Change 3 is left open **in code as well as in prose** — `envelope_membership(...)` takes
`functional` as a required argument with **no default** until Part B supplies one — so that no
implicit choice can be made by a default value while the choice is formally open. That single slot is
the only thing in this entry left to be decided later, and it is decided by a rule stated here, not by
judgement exercised there.

Spec sections touched: **§1 (inclusion gate)** — Amendment 3's two new per-stratum requirements
acquire a definition and a decision procedure. **Nothing else in the frozen spec is touched.** §8(c)'s
effect size (log2FC = 1.0, K = 200) and threshold (≥ 0.60), §8(a)'s calibration criteria, and
Amendment 3's envelope table are unchanged by this entry and are not weakened by it.

### Data visible at the time of this amendment (full disclosure)

1. **Everything visible for Amendments 1, 2 and 3**, unchanged.
2. **The frozen stratum list**, `docs/PREREGISTRATION_STRATUM_LIST.md` and
   `pilot/preregistration/stratum_list_2026-08-16.json`, committed **earlier the same day**: 251
   stratum-contrasts over 12 independent datasets, `admitted_to_sweep = False` on every row, and **no
   metric computed on any of them**. This amendment is the work its §9 item 2 names as the blocker,
   and its §6 tiers become decidable exactly when Part B lands. It was read while writing this entry,
   and its distributional figures are used below as *the range the estimator must cover*: group
   medians of 11.0 … 6671.5 cells per donor with 85 of 502 group medians in the `[10, 30)` bin, median
   counts per cell 284 … 56 841.5, and `min(n_A, n_B)` from 3 to 39.
3. **The simulator's development seed range, `seed ∈ [1, 999]`, is disclosed as seen.** The estimator
   will be developed and debugged against it, and the probe numbers in Correction 1 below were already
   measured on it while writing this entry. The confirmatory grid runs on a declared, non-overlapping
   seed range (Change 5) and has **not** been run.
4. **The confirmatory grid has not been run and no estimator has been validated.** Nothing in this
   entry claims that the estimator works. `src/pbcheck/sigma_donor.py` does not exist at the time of
   writing.
5. **No real data.** Oracle (d), the Mathys 2019 anchor (§8(d)), remains untouched; the freeze records
   that its Census path is closed and the Synapse path gated on a ROSMAP data-use agreement begun
   2026-08-16. No CELLxGENE stratum has been loaded. Every number below is synthetic, and
   `sigma_donor` remains a free knob of our own simulator.
6. **What is new here and was measured for it**: ten probe cells of Amendment 3's quantity on
   development seeds — six over 16 seeds, two over 8, two on a single seed — tabulated under Correction
   1. They are the evidence that Amendment 3's claim is false, and they make the expected outcome of
   V9b known in advance, which is disclosed where V9b is stated. They are **not** the validation, which
   is V1–V9 and is unrun.

### What Amendment 3 left open, and precisely which part of it this closes

Amendment 3's own closing item, quoted rather than paraphrased:

> **`sigma_donor` is still not anchored to real data.** This is the third consecutive amendment to
> close on it and it is *not* resolved here. Change 1 supplies the envelope and the estimation
> mechanism; it does not supply the estimate, the validated conversion from `s0^2` to `donor_sigma`,
> or the empirical mean-dispersion / donor-variance trend §8(b) asks for. Until that exists, the
> envelope is a statement about our simulator's coordinates, and **whether any real stratum falls
> inside it is unknown**.

**What Part A closes: the middle clause — the validated conversion.** Part A specifies the estimator
that replaces the mechanism, derives it from the simulator's generative model term by term, and fixes
the criteria under which it may be called validated. Part B either validates it or records that it
failed. On success, `sigma_donor_estimate` and `envelope_membership` stop being `PENDING` columns in
`census_select`'s manifest and become computed columns with a stated error direction.

**What Part A does not close, stated with the same precision.**

* **§8(b)'s real empirical mean–dispersion / donor-variance trend.** The estimator is derived in the
  simulator's coordinates and validated against the simulator's truth. That establishes that it
  recovers `donor_sigma` *from data the simulator generated*; it establishes nothing about whether
  real snRNA-seq between-donor variation is the log-normal random effect the simulator imposes. A
  perfectly validated estimator of the wrong model's parameter is still wrong. §8(d) remains the
  binding real check and remains unrun.
* **The pooling question.** `pooled` is `unresolved` on all 1197 candidates and all 251 frozen strata:
  the pinned Census exposes no library/pool identifier, so D3's answer is its "where pooling cannot be
  resolved" state. Donor pseudobulk is therefore a **lower bound** on the correct replication unit
  throughout, and every `sigma_donor` estimate inherits that caveat rather than resolving it. Part A
  carries the caveat into the manifest and does not pretend to fix it.
* **Whether any real stratum lands inside the envelope at all.** This amendment supplies the measuring
  instrument, not the measurement. The freeze's §6 scenario table stays arithmetic under a
  hypothetical σ until the estimator has run on real strata, and at σ ≈ 0.5 the surviving set is 5 of
  12 datasets, below §1's own 8–12 floor. A negative answer remains a live outcome.

---

### Correction 1 to Amendment 3 Change 1 — `sqrt(s0^2) · ln 2` is NOT an upper bound on `donor_sigma`, and the direction of its error is the dangerous one

This is placed first rather than folded into the derivation, because it reverses a statement this log
has already published and that two further committed documents have since repeated.

**The statement being corrected**, quoted from Amendment 3 Change 1:

> **That conversion is a mechanism, not an anchor, and the difference is load-bearing.** `s2_g`
> contains the donor random effect **plus** residual NB sampling noise and the compression of the `+1`
> offset at low expression, so `sqrt(s0^2) · ln 2` is an **upper bound** on `donor_sigma`, not an
> estimate of it.

**Why it is wrong.** The sentence names three ingredients and then treats all three as additive
contaminations of the donor effect. Two of them are. The third is not: the `+1` in `log2(CPM + 1)` is
a **multiplicative attenuation** applied to everything inside the logarithm, the donor random effect
included. Writing it out — the full derivation is Change 1, this is only its consequence — with
`a = CPM/(1 + CPM)`:

```
Y = log2(1 + C)      =>      dY/d ln C = (1/ln 2) · C/(1 + C) = a / ln 2
Var_within-group(Y)  ≈  (a² / ln 2²) · (σ² + v_tech)
```

so the arm's own fit sees `s2_g ≈ a_g² (σ² + v_g) / ln 2²`, and Amendment 3's quantity is

```
sqrt(s0²) · ln 2   ≈   ā · sqrt(σ² + v̄)
```

There are two distortions and they point in **opposite** directions: `+ v̄` inflates, `ā ≤ 1` deflates.
Amendment 3 named both and kept only the first. The quantity is an upper bound on σ if and only if

```
ā² (σ² + v̄)  ≥  σ²        ⟺        v̄  ≥  σ² (1 − ā²)/ā²  =  σ² (1/ā² − 1)
```

and it is an **under**statement whenever the technical variance is smaller than that. Nothing
guarantees that it is not.

**Worked arithmetic, at a point that is ordinary rather than adversarial.** Take `ā = 0.9` — a gene at
CPM = 9. In a 15 000-gene universe the mean CPM is 10⁶/15 000 = 66.7 and the distribution is heavily
right-skewed, so under the simulator's own gene-mean distribution CPM = 9 is the 14th percentile, and
in real data, whose low expression tail is far heavier, it is commoner still. Take `v̄ = 0.011` (a donor
of a few thousand cells; see Change 1.3) and a true `σ = 0.35`:

* required for the bound to hold: `v̄ ≥ 0.1225 × (1/0.81 − 1) = 0.1225 × 0.234568 = 0.028735`;
* actual: `ā² v̄ = 0.81 × 0.011 = 0.00891` against `σ²(1 − ā²) = 0.1225 × 0.19 = 0.023275` — the
  inequality fails by a factor of 2.6;
* value returned: `sqrt(0.81 × (0.1225 + 0.011)) = sqrt(0.108135) = 0.3288`.

**The mechanism returns 0.329 against a truth of 0.350 — a 6 % understatement — and Amendment 3 would
have a reader treat it as a ceiling.**

**Where the reversal lives, and why the intuition points the wrong way.** Substituting the derived
technical term `v ≈ 1/T + φ·r_d` (Change 1.3) with `T = C·L/10⁶`, `L` the donor's universe-restricted
library size, the per-gene condition `a² v ≥ σ²(1 − a²)` becomes, exactly,

```
10⁶ · C / L   +   φ · r_d · C²    ≥    σ² · (1 + 2C)
```

For genes with `C ≫ 1` and a donor deep enough that the dispersion term `φ r_d C²` is subdominant,
both sides grow linearly in `C`, the gene's own expression cancels, and what remains is a statement
about the **donor's library size alone**:

```
the bound survives only while     L  ≲  10⁶ / (2σ²)
        σ = 0.2  →  1.25e7                 σ = 0.35 →  4.08e6
        σ = 0.5  →  2.00e6                 σ = 0.7  →  1.02e6
```

This is the opposite of what "the compression of the `+1` offset **at low expression**" suggests, and
the inversion is worth stating plainly, because a reader checking the claim will otherwise look in the
wrong place:

* **Shallower data does not produce the failure — it prevents it.** Thinning a stratum raises `1/T`,
  raises `v̄`, and pushes the quantity *up*, away from the reversal. The failure needs technical noise
  to be *small*. In the noiseless limit the quantity converges to exactly `ā·σ`, which is below `σ` for
  every `ā < 1`. **Amendment 3's bound is least true on the cleanest data.**
* **What sets `ā` is the universe, not the depth.** CPM is a *within-universe relative* quantity:
  scaling every gene's mean by a constant scales `T` and `L` together and leaves CPM unchanged. The
  attenuation is therefore governed by how many genes share the million — median CPM ≈
  `(10⁶/G) · e^(−s²/2)` for log-normal gene means of log-sd `s` — and not by sequencing depth at all.
  At the simulator's `G = 1500` the median CPM is ≈ 300 and `ā ≈ 0.997`; at a realistic frozen universe
  of `G ≈ 15 000` it is ≈ 30 and `ā ≈ 0.97`, with a long low tail beneath.

**The sign of Amendment 3's error is therefore a property of the stratum, not a constant**, and the
frozen list straddles the threshold in both directions: its group medians run 11.0 … 6671.5 cells per
donor and 284 … 56 841.5 counts per cell, so donor libraries across the 251 span several orders of
magnitude around every `L` tabulated above. (Those two ranges are marginals from the freeze §4.2/§5
and are not paired per stratum; the point is the span, not a product.) That is precisely why a
*validated estimator* is required and a bound of unknown sign is no substitute for one.

**Why this is the dangerous direction and not a technicality.** Understating σ admits strata for which
the pseudobulk arm is not declared valid. That arm is the denominator of every inflation number in the
study; spec §10 ranks "pseudobulk denominator collapse" as risk 1, the single most dangerous failure,
and decision rule item 1 makes it **VOID → NO-GO**. A ceiling that is really a floor converts the
envelope from a protection into a rubber stamp — and it fails silently, because the reader has been
told the number errs the other way.

**What this correction does and does not change.**

* Amendment 3's **envelope is untouched**. The table of `(sigma_donor, min_donors_per_group)` is a
  statement about power at given `(σ, n)`; nothing here bears on it. Change 1 (a), (b) and (c) stand in
  full, as does Change 2.
* The **status of the mechanism** changes: `sqrt(s0²)·ln 2` is demoted from "upper bound on
  `donor_sigma`" to **an audit quantity of unknown error sign**. It goes on being computed and written
  beside every estimate (Change 4 item 6), because comparing it against `σ̂` is how a reader sees that
  the correction did work — but it may not gate anything, in either direction.
* **No measurement is retracted.** No stratum was ever admitted on the strength of it:
  `admitted_to_sweep` is `False` on all 2190 rows of the candidate manifest and all 251 frozen strata,
  and `sigma_donor_estimate` is `PENDING` on every one. The error has not propagated into any published
  number. That is a debt the freeze's refusal to admit anything happens to have covered, not evidence
  that the claim was harmless.
* **Three committed texts repeat the erroneous claim, and they are named here rather than silently
  edited**: `src/pbcheck/census_select.py`'s `PENDING_FIELDS["sigma_donor_estimate"]` ("states it is an
  UPPER BOUND, not an estimate"), and `docs/PREREGISTRATION_STRATUM_LIST.md` §6 ("an unvalidated upper
  bound") and §9 item 2 ("stated in terms that the quantity is an **upper bound**"). Both documents
  were correct to *cite* Amendment 3 and are wrong only by inheritance. They are corrected in the
  commit that next touches their owning code — Part B — and not before, so that the correction and its
  evidence land together.

**Criterion V9 exists to demonstrate this empirically rather than leave it resting on the algebra
above.** It is pre-declared in Change 5, with both clauses and both numbers.

#### The measurement behind Correction 1 (development seeds, disclosed)

Ten probe cells, run while writing this entry through the arm's own code path (`build_pseudobulk` →
`frozen_universe` → `moderated.log_cpm` → `wls_two_group` → `fit_f_dist`) at 8 v 8 donors, dispersion
0.2, simulator defaults otherwise. `L̃` is the median donor universe-restricted library size; `ā` the
median of `CPM/(1 + CPM)`. Rows are ordered by `L̃/L_crit`, with `L_crit = 10⁶/(2σ²)` from the closed
form above — a quantity computed from the design, not fitted to the results. Six cells were run over 16
development seeds, two over 8, and two on a single seed; the seed count is a column, and no row is read
for more than its seed count supports.

| universe / depth | σ | `L̃` | med. CPM | `ā` | seeds | mean `sqrt(s0²)·ln 2` | ratio to σ | below σ | `L̃/L_crit` |
|---|---|---|---|---|---|---|---|---|---|
| 1000 genes, 30 cells, `mean_log_mu = 0` | 0.35 | 6.12e4 | 461.7 | 0.9978 | 16 | 0.4216 | 1.205 | 0/16 | 0.015 |
| 1000 genes, 30 cells | 0.35 | 1.63e5 | 461.3 | 0.9978 | 16 | 0.3886 | 1.110 | 0/16 | 0.040 |
| 1000 genes, 300 cells, `mean_log_mu = 0` | 0.35 | 6.02e5 | 464.1 | 0.9978 | 16 | 0.3599 | 1.028 | 0/16 | 0.15 |
| 1000 genes, 300 cells | 0.35 | 1.63e6 | 463.6 | 0.9978 | 16 | 0.3536 | 1.010 | 1/16 | 0.40 |
| 1500 genes, 250 cells | 0.35 | 2.11e6 | 302.3 | 0.9967 | 16 | 0.3533 | 1.009 | 1/16 | 0.52 |
| **`ORACLE_SIM`: 1500 genes, 250 cells** | **0.50** | 2.12e6 | 282.3 | 0.9965 | 16 | **0.5006** | **1.001** | **7/16** | **1.06** |
| 15 000 genes, 300 cells | 0.20 | 2.62e7 | 31.5 | 0.9692 | 1 | **0.1955** | **0.977** | 1/1 | 2.10 |
| 15 000 genes, 100 cells | 0.35 | 8.91e6 | 30.2 | 0.9682 | 8 | **0.3415** | **0.976** | **8/8** | 2.18 |
| 15 000 genes, 300 cells | 0.35 | 2.62e7 | 30.3 | 0.9683 | 8 | **0.3318** | **0.948** | **8/8** | 6.42 |
| 15 000 genes, 300 cells | 0.50 | 2.60e7 | 28.4 | 0.9660 | 1 | **0.4669** | **0.934** | 1/1 | 13.0 |

Four things are read off this table and nothing further is claimed from it.

1. **The reversal is real, and it is ordered by exactly the quantity the derivation says orders it.**
   Reading down the table, `L̃/L_crit` climbs from 0.015 to 13 and the ratio falls monotonically from
   1.205 to 0.934, crossing 1 where the closed form says it should. The threshold is a heuristic — it is
   taken at the median gene and drops the `φ r C²` term, whereas `fit_f_dist`'s location is a log-scale
   central value over the whole gene distribution (its estimator is a mean of
   `e_g = log s2_g − digamma(d/2) + log(d/2)`), so the aggregate is nearer
   `GM_g(a_g) · sqrt(GM_g(σ² + v_g))` than the median-gene expression, and two approximation errors of
   opposite sign are partly cancelling. It is quoted for orientation only and **no code reads it**. What
   it does earn is this: the single cell it places on the boundary is the single cell that measures as
   being on the boundary.
2. **That crossing sits inside the gate's own operating point.** `gate_config.ORACLE_SIM` is 1500 genes,
   250 cells per donor, 8 v 8, `donor_sigma` = 0.5 — the calibration regime the gate has run at since
   Amendment 3 — and there the quantity averages **0.5006 against 0.5000** over 16 seeds (sd 0.0021,
   SE of the mean 0.0005), with **7 of 16 individual realisations below the truth**. A quantity that
   falls under the value it is supposed to bound on nearly half of single draws is not an upper bound;
   it is an estimate whose bias happens to be about +0.1 % at this one geometry. At σ = 0.35 on the same
   geometry the margin is +0.9 %, one realisation in 16 below. Neither is a margin an admission rule can
   rest on — and the fact that the claim is *nearly* true in the only regime that has ever been
   exercised is a large part of why it went unchallenged for an amendment.
3. **Shallow is the wrong place to look, exactly as derived.** The two `mean_log_mu = 0` cells — the
   intuitive "low expression" probe — move *away* from the reversal, to ratios 1.028 and 1.205 with 0 of
   16 realisations below σ in either, because thinning adds technical noise without touching CPM.
   Anyone hunting this failure by thinning the data will conclude the bound is safe, and the thinner
   they make it the safer it will look.
4. **Against a realistic universe the shortfall is not a rounding error, and it is not noise.** At
   15 000 genes — an order of magnitude closer to a real frozen universe than the simulator's 1500 — the
   quantity reads **2.3 % to 6.6 % below the truth at every σ probed**. The two cells measured over 8
   seeds are below σ on **8 of 8** realisations with a replicate sd of **0.0003**, because `s0²` over
   15 000 genes is a very stable quantity: at 300 cells the mean is 0.3318 against a truth of 0.35, a
   gap of 0.018 against a standard error of the mean of 0.0001. This is a systematic offset, not a
   sampling accident. The remaining two rows are single realisations far from the boundary and are read
   only for their sign.

These are development-seed probes, disclosed as such, and each row is read only as far as its seed count
allows. Together they establish that Amendment 3's claim is false — a counterexample suffices for that.
The pre-registered quantitative statement is V9's, on disjoint confirmatory seeds.

---

### Change 1 — the estimand, and the derivation term by term (§1)

#### 1.1 The estimand

> **A stratum's `sigma_donor`** is the value that, substituted into `synthetic/oracles.py::simulate`,
> reproduces that stratum's between-donor dispersion as the pseudobulk arm's own statistic sees it:
> the standard deviation, on the **natural-log** scale, of the per-(gene, donor) log-normal random
> effect `re[g,d] = exp(N(−σ²/2, σ²))`.

Three properties of the definition are load-bearing.

**It is pooled over both groups.** The arm's statistic is `wls_two_group`'s pooled residual variance on
`d = n_A + n_B − 2` degrees of freedom, and Amendment 3's envelope was computed for exactly that
statistic. Per-group estimates `σ̂_A` and `σ̂_B` are reported as diagnostics and gate nothing.

**Where the groups are unbalanced, the envelope's `n` is `min(n_A, n_B)`.** That is already this
repository's convention — `census_select.envelope_max_sigma_supported`, from which the frozen
artifact's column of the same name is built — it is conservative, and it is kept unchanged.

**The simulator has one σ for all genes; a real stratum has a distribution `{σ_g}`.** The estimand is
therefore a *functional* of that distribution, and which functional is the single question Part A
deliberately leaves open (Change 3). Everything else here is fixed.

#### 1.2 The generative model, quoted from the code

For cell `c` of donor `d`, gene `g`, group `k`, from `synthetic/oracles.py::simulate`:

```
counts[c,g] ~ NB(mean = λ, Var = λ + φ λ²),      λ = depth_c · μ_g · re[g,d] · fc[g,k]
depth_c = exp(N(0, s_dep²)),      s_dep = depth_log_sigma = 0.3
re[g,d] = exp(N(−σ²/2, σ²)),      σ = donor_sigma, NATURAL log
μ_g     = exp(N(mean_log_mu, mean_log_sigma²)),      defaults 1.0 and 1.2
```

Pseudobulk is decoupler's `mode='sum'`: `T[d,g] = Σ_c counts[c,g]`. The library is the row sum of the
**universe-restricted** matrix, `L_d = Σ_{g∈U} T[d,g]` (`moderated.log_cpm`) — universe restriction is
part of the arm and is preserved by the estimator. `C = 10⁶·T/L`, and the tested scale is
`Y = log2(C + 1)`. Write `b_gd = ln re[g,d] ~ N(−σ²/2, σ²)`, `S_d = Σ_c depth_c`, `Q_d = Σ_c depth_c²`.

#### 1.3 T1 — NB sampling noise of the pseudobulk sum

Conditional on `(re, depths)`:

```
E[T | ·]   = μ_g · re_gd · S_d
Var[T | ·] = μ_g · re_gd · S_d  +  φ · (μ_g · re_gd)² · Q_d
```

and the delta method on `ln T` gives

```
v_lnT[d,g]  ≈  1/E[T | ·]  +  φ · r_d ,          r_d = Q_d / S_d²
```

Both terms fall as `1/n_cells`: for iid depths `r_d ≈ (1 + cv_depth²)/n_c`, and with
`depth = exp(N(0, 0.3²))` the squared coefficient of variation is `e^0.09 − 1 = 0.0942`, so
`r_d ≈ 1.0942/n_c`.

Magnitudes, because this term's importance is entirely regime-dependent. At 250 cells per donor with
φ = 0.2, `φ r_d = 0.2 × 1.0942/250 = 8.75e-4`, negligible against `σ² = 0.1225` at σ = 0.35. At 30
cells and a moderate gene (`T ≈ 30`), `1/T + φ r_d = 0.0333 + 0.0073 = 0.0406`; at 10 cells,
`0.0333 + 0.0219 = 0.0552`. Both are **comparable to `σ² = 0.04` at σ = 0.2**. The correction carries
its whole weight at the thin end, which is where real strata live: 85 of the frozen list's 502 group
medians sit in the `[10, 30)` cells-per-donor bin and its smallest per-donor cell count is 10.

#### 1.4 T2 — library-composition noise

`ln C = ln T − ln L + const`. The common factor `ln S_d` enters both and **cancels exactly**: a donor's
overall depth scale does not appear in CPM. What survives is

```
ln L_d  ≈  ln S_d + ln( Σ_g w_g e^{b_gd} ) + (NB noise of L),      w_g = μ_g / Σμ
Var(ln L̃_d)  ≈  (e^{σ²} − 1) · Σ_g w_g²
```

Bounding it: for log-normal gene means with `mean_log_sigma = 1.2`, `Σ w² ≈ e^{1.44}/G = 4.221/G`, so
at `G = 1500` it is `2.81e-3`; at σ = 0.5 that gives `(e^{0.25} − 1) × 2.81e-3 = 0.284 × 2.81e-3 =
8.0e-4`, which is **0.32 % of σ²**. The NB noise of `L` itself is O(1/L) and vanishing. The covariance
`Cov(b_gd, ln L̃_d) ≈ w_g·(…)` is bounded by the same `Σ w²` and is material only for the handful of
top-expressed genes.

**T2 is derived as a bounded small addition and is NOT corrected for**; the bound is checked by
validation rather than assumed. In real data its analogue — donor-specific normalisation shifts — is
*genuine* between-donor variance that the arm actually experiences. It therefore belongs inside σ and
must not be subtracted; subtracting it would flatter the estimate in the dangerous direction.

#### 1.5 T3 — the `+1` compression

`Y = log2(1 + C)`, `dY/d ln C = (1/ln 2) · C/(1 + C) ≡ a/ln 2`, with `a_dg = C_dg/(1 + C_dg)`:

```
Var_within-group(Y[d,g])  ≈  (a_dg² / ln 2²) · [ σ² + v_lnT[d,g] + Var(ln L̃) − 2 Cov(b, ln L̃) ]
```

This is the term Correction 1 is about. It is a *multiplicative* attenuation of everything inside the
logarithm — the donor effect included — which is why the estimator must **de-attenuate**, dividing by
an attenuation factor, and not merely subtract a noise term.

#### 1.6 What the arm's own `s2_g` is an expectation of, written exactly

`wls_two_group` fits `~ 1 + x` with `x` a 0/1 group indicator, so the hat-matrix diagonal is `1/n_A`
for a donor in group A and `1/n_B` in group B, and

```
E[s2_g]  =  Σ_d w_d · Var(Y[d,g]),        w_d = (1 − 1/n_{group(d)}) / (n_A + n_B − 2)
```

with `Σ_d w_d = [(n_A − 1) + (n_B − 1)]/(n_A + n_B − 2) = 1`. Substituting T3:

```
E[s2_g] · ln 2²  =  σ² · ( Σ_d w_d a_dg² )  +  ln 2² · ( Σ_d w_d v_Y[d,g] )   (+ the bounded T2 term)
```

where `v_Y[d,g]` is the technical variance **already on the Y scale**, i.e. carrying its own
`a_dg²/ln 2²`. The estimator inverts this, per gene:

```
                     s2_g  −  Σ_d w_d · v̂_Y[d,g]
σ̂²_g   =   ln 2² · ───────────────────────────────  ,        ā²_g  ≡  Σ_d w_d a_dg²
                                ā²_g
```

Two details in that line are not cosmetic, and they are pinned here because getting either wrong is
silent. **`ā²_g` is a weighted mean of squares, not the square of a weighted mean** — using `(Σ w a)²`
inflates σ̂ by `Var_d(a_dg)`, a safe direction but a wrong number. And **the weights are the df weights
`w_d`, not `1/n`** — at `n_A = 3` the factor `1 − 1/n_A = 2/3` is a 33 % correction, and the two
coincide only when the groups are balanced. The frozen list holds designs as skewed as 100 versus 10
donors (#5) and 20 versus 7 (#3), so this is a live case and not a limit.

#### 1.7 T4–T6, the terms that are named and not corrected

* **T4 — between-cell depth variation.** Cancels in CPM in the mean (through `S_d`, step 1.4) and
  enters only via `r_d` in T1. Not separately corrected.
* **T5 — second-order delta-method error.** The expansion truncates at first order in `b`, which
  matters most at σ = 0.7 where the log-normal is visibly skewed; and `E[1/T] > 1/E[T]` (Jensen) at
  small `T` biases a plug-in `v̂` **upward** and hence `σ̂` **downward** — the dangerous direction. It is
  suppressed by the estimation stratum's count floor (`ESTIMATION_MIN_MEDIAN_COUNT`, Change 2) rather
  than corrected analytically, and its residue is measured as net bias by V1 and V2.
* **T6 — selection bias from choosing the estimation stratum.** Genes are selected on the
  **label-agnostic median CPM pooled over all donors** — never on a minimum, never on any quantity
  involving the condition labels, and never on observed variance, so the selection is not made against
  the quantity being estimated. The residual direction and magnitude are measured by V1/V2, not
  asserted.

#### 1.8 What a real stratum can and cannot yield

**Estimable.** `T_dg` and `L_d` directly; `a_dg` as a plug-in from the observed CPM; `r_d` **exactly**
from the per-cell library sizes `ℓ_c`, since within a donor `ℓ_c ∝ depth_c` and the composition
constant cancels in `Σℓ²/(Σℓ)²`; `φ` from the per-cell counts within donor by method of moments with
`ℓ_c` as the depth covariate, pooled over genes; and `v_Y[d,g]` **directly**, by resampling the donor's
own cells (E2, Change 2) — which is the reason E2 rather than E1 is primary.

**Not estimable.** Per-gene `σ_g` separately from the library-composition term at single-gene
resolution; it is bounded above by T2's derivation and reported, not identified. And the state of donor
pooling: D3 is `unresolved` on 251 of 251 frozen strata, so donor pseudobulk is a lower bound on the
correct replication unit and every `σ̂` inherits that caveat. It is written into the manifest row beside
the estimate; it is not resolved here and cannot be.

---

### Change 2 — the estimator: E2-bootstrap primary, split-half as the scaling cross-check, E1 as a misspecification detector, E3 reported

All three candidates share the skeleton of 1.6 and differ only in the source of `v̂_Y`. All three are
computed in one pass; the cost is set by E2.

**E2 — non-parametric within-donor resampling. PRIMARY.** Two implementations, both required, each the
other's cross-check.

* **Bootstrap — the primary source of `v̂`.** `B_BOOT = 50` times, rebuild the donor's pseudobulk profile
  by multinomial resampling of its own cells (one sparse mat-vec per replicate), recompute the
  universe-restricted `L*`, `C*` and `Y*`, and take `v̂_Y[d,g] = Var_B(Y*_dg)`. This estimates the
  technical variance **at full depth with no rescaling**, and it automatically contains T1, T3, the NB
  noise of the library, and any within-donor structure the cells themselves carry. It assumes no
  parametric count model whatsoever.
* **Split-half — the scaling cross-check.** `SPLITHALF_R = 20` random partitions of the donor's cells
  into halves `n1 + n2 = n_c`. Both technical terms scale as `1/cells` **exactly in expectation** — the
  Poisson term through `S_d`, the dispersion term through `r_d`, both by 1.3 — so the
  half-depth-to-full-depth conversion is
  ```
  v̂_full  =  E_splits[ (Y_h1 − Y_h2)² ]  /  ( n_c · (1/n1 + 1/n2) )        [ = …/4 when n1 = n2 ]
  ```
  The unequal-halves form is **mandatory, not an optimisation**: `n_c` is odd about half the time, and
  at the inclusion gate's floor of 10 cells the difference between `/4` and the exact divisor is not a
  rounding error. The two-line proof of the `1/cells` scaling is carried in the function's docstring,
  beside the code that depends on it.

E2's per-donor `v̂` does not depend on the number of donors, so it does not degrade at `n = 3`; its
residual error is O(1/n_c) *relative*. Its cost is `B_BOOT` sparse mat-vecs per donor, and per-donor
slices and the resampling mat-vec are computed without densifying, after `io_counts._iter_value_blocks`.

**E1 — `nb_plugin`, the misspecification detector.** `v̂_Y[d,g] = (a_dg²/ln 2²) · (1/max(T_dg, 1) +
φ̂ · r_d)`, with `φ̂` estimated from the per-cell counts inside each donor by moments
(`Var_c(x_cg) = m̄ + φ m̄² (1 + cv_ℓ²)`, solved for φ, pooled robustly over genes and donors). It is
nearly free, and it is **not** primary for a reason stated rather than discovered: under zero inflation
it underestimates `v` and therefore overestimates σ, which is safe; but under within-donor substructure
that NB cannot see it also underestimates `v`, and there the direction is **not** guaranteed. Its role
is that of a dissenting witness — `|σ̂_E1 − σ̂_E2| > 0.05` raises an `nb_misspecification` flag, always
reported, and blocking only if V7's M4 arm fails.

**E3 — `plateau`, reported and never primary.** Restrict to genes whose predicted technical share is
negligible (`v̂_Y ≤ ε · median(s2)` and `a ≥ 0.95`) and read `s2_g · ln 2² ≈ σ²` off directly, located by
the unbiased log-scale construction `fit_f_dist` already uses. It is the simplest of the three and it
errs upward — residual technical noise only adds — which is the safe direction. But on thin or small
strata no such genes exist and it is **obliged to refuse**, and that is exactly the zone of real-data
risk. It is reported as a consistency check and gates nothing.

**Constants.** These enter `gate_config` as a new `SIGMA_GATE` block in `PRE_REGISTERED` once this Part
A is committed. `scripts/sigma_probe.py` carries its own copy and does **not** import `gate_config`,
following the anti-circularity pattern that module's docstring sets out for
`scripts/pb_calibration_probe.py`; a test pins the two copies equal, as `tests/test_gate_config.py`
already does for the existing constants.

| constant | value | why this value |
|---|---|---|
| `B_BOOT` | 50 | the relative Monte-Carlo error of a variance over 50 replicates is `sqrt(2/49)` = 20 % per gene per donor, averaged down by the df-weighted pooling over donors and by the ≥ 100-gene estimation stratum |
| `SPLITHALF_R` | 20 | cross-check only; not the primary `v̂` |
| `ESTIMATION_MIN_MEDIAN_CPM` | 20.0 | `a = 20/21 = 0.952 ≥ 0.95`, which bounds the de-attenuation divisor so it cannot amplify error without limit. It does **not** make the naive bound safe — at `a = 0.952` and σ = 0.35 that bound still needs `v̄ ≥ 0.1225 × (1/0.9070 − 1) = 0.0126` — which is why the estimator de-attenuates instead of assuming `a ≈ 1` |
| `ESTIMATION_MIN_MEDIAN_COUNT` | 10 | suppresses T5's Jensen bias on `1/T`, which pushes σ̂ **down** |
| `MIN_ESTIMATION_GENES` | 100 | below this the per-gene aggregate is not a distribution |
| `UCB_LEVEL` | 0.90 | Change 4 |
| `ENVELOPE_LOOKUP` | `"step_up"` | Change 4 item 3 |
| `E1_E2_DIVERGENCE_FLAG` | 0.05 | in σ units; the flag above |
| `MAX_TECHNICAL_SHARE` | 0.8 | Change 4 item 4 |

**The pre-declared reasons for an `indeterminate` verdict.** A stratum is `indeterminate` — excluded
from the sweep and **counted, with its reason, in the D4 excluded-strata statistic** — if and only if
one of the following holds. The list is closed as of this commit; adding to it later is an amendment.

1. Fewer than `MIN_ESTIMATION_GENES` = 100 genes pass the estimation-stratum filters.
2. Degenerate jackknife: some leave-one-donor-out replicate fails to produce a finite functional, or
   `SE_jack` is not finite.
3. Median technical share `median_g( Σ_d w_d v̂_Y[d,g] / s2_g ) > MAX_TECHNICAL_SHARE = 0.8` — the
   estimate would be mostly its own correction.
4. A Tier-C misspecification signature that V7 converts into a blocking flag. The **set** of candidate
   signatures is fixed now — LODO influence (M3), E1/E2 divergence (M4), the coefficient of variation
   of cells per donor (M2), and the realised `d0` of the corrected per-gene values (M1) — and *which* of
   them block is decided mechanically by which Tier-C arms fail V7.
5. Median cells per donor below a floor, **if and only if V2 fails at 10 cells**, in which case the
   floor is 30 and is written into the rule automatically rather than invented afterwards.

---

### Change 3 — the aggregation functional, and the mechanical rule that will select it

The simulator has one σ, the envelope is indexed by a scalar, and a real stratum has `{σ_g}`. Something
has to collapse.

**The primary aggregate is a pooled moment and not a quantile, and the reason is arithmetic.** Per-gene
`σ̂²_g` at `d = n − 2 ≤ 14` degrees of freedom is enormously noisy: `s2_g` has relative standard
deviation `sqrt(2/d)`, which is 71 % at 3 v 3 (`d = 4`) and 38 % at 8 v 8 (`d = 14`). A quantile of
`{σ̂²_g}` is therefore a quantile of the *convolution* of the true `σ_g` distribution with that noise,
which is systematically wider than the truth — a per-gene q75 reads high, a per-gene q25 reads low, and
neither estimates the corresponding quantile of `σ_g`. What **is** honestly estimable is the **mean of
the variances**, because the noise averages out:

```
M  =  mean_{g ∈ estimation stratum} σ̂²_g      estimates      mean_g σ_g²  ,        σ_rms = sqrt(M)
```

Negative per-gene values are **not clipped before aggregation** — clipping a symmetric noise
distribution at zero is a systematic upward bias, and it is exactly the sort of quiet safety margin
this log exists to prevent. Clipping happens once, on the final scalar, and the fraction of negative
per-gene mass is reported.

**The four candidate gate functionals** are `rms`; `trimmed_rms_10`, a two-sided 10 % trim of the
per-gene values; `median_log`, the `fit_f_dist` location over corrected values, i.e. approximately a
geometric mean, which under heterogeneity is expected to read low and therefore expected to be
rejected; and `q75`.

**There is an analytic prior in favour of `rms`, and it is labelled a prior.** The arm's sensitivity
over K random genes is `E_g[power(σ_g²)]`; in the operating region power is convex and decreasing in
σ², so by Jensen `power(mean_g σ_g²) ≤ E_g[power(σ_g²)]` — gating on the RMS is conservative relative
to the sensitivity actually achieved. That is an argument, not a measurement, and this amendment does
not let it decide.

**The choice is made by a mechanical rule, declared now and applied to data not yet seen.** This is
Amendment 2 Change 1's device, used here in the honest order rather than the retrospective one:

> Among the four functionals, discard any whose wrong-admission rate exceeds V5's threshold on **any**
> heterogeneous cell of Tier C arm M1 — where a cell's ground truth is the **measured** sensitivity of
> the `ebayes` arm on the pre-registered oracle (log2FC = 1.0, K = 200), the cell counting as `inside`
> iff that sensitivity is ≥ 0.60. Among the survivors, choose the one **maximising pooled deep-inside
> admission**, V6's statistic. Ties are broken toward the **simplest** functional, in the declared
> complexity order `rms < trimmed_rms_10 < median_log < q75`.

If the rule eliminates all four, that is the result: no functional is fit for the job, no
`GATE_FUNCTIONAL` constant is written, and the failure comes back to this log. If the rule turns out
degenerate on the realised grid, as Amendment 2's worst-case-power rule did, the degeneracy is reported
and the non-degenerate reformulation is argued in Part B in the open, as Amendment 2 did — never
repaired quietly.

**Reported alongside, and never in place of, the scalar**: quantiles of the corrected per-gene values,
carrying the convolution caveat above; the realised `d0` from `fit_f_dist` over corrected values as a
heterogeneity measure; per-group `σ̂_A` and `σ̂_B`; and the tail mass `P(σ̂²_g > k·M)`.

**What the scalar loses, stated rather than discovered.** It cannot represent bimodal `σ_g`. A stratum
with a low RMS and a heavy tail is admitted at a σ for which its tail genes are not covered; the tail
metric makes that partly visible and does not fix it. The limitation is repeated in *What this does NOT
settle*.

---

### Change 4 — the operating-envelope membership rule (§1 inclusion gate)

**The error being excluded, and why it is asymmetric.** Understating σ admits strata for which the arm
is not valid, and the arm is the denominator of every inflation number (§10 risk 1; decision rule item
1's VOID → NO-GO). Overstating σ merely excludes valid strata, which costs sample size and is reported
honestly in the D4 excluded-strata bookkeeping. **The project is defended against understatement**, and
every choice below breaks that way.

1. **The gate reads an upper confidence bound, not the point estimate.** `σ_gate = UCB_90(σ_f)`, a
   one-sided 90 % upper bound on the selected functional.

2. **The interval is a leave-one-donor-out jackknife, clustered on the donor.** For each of the
   `n_d = n_A + n_B` donors, drop it and recompute **everything** — `s2`, the `v̂` rows, `ā²`, and the
   aggregate — then
   ```
   SE_jack  = sqrt( (n_d − 1)/n_d · Σ_i (θ_(i) − θ̄)² ) ,     θ = the functional, on the σ² scale
   UCB(σ²)  = θ̂ + t_{0.90, n_d − 1} · SE_jack ,              σ_gate = sqrt(max(UCB, 0))
   ```
   with `t_{0.90, n_d−1}` = 1.476 / 1.415 / 1.341 / 1.316 at `n_d − 1` = 5 / 7 / 15 / 25, i.e. 3 v 3 /
   4 v 4 / 8 v 8 / 13 v 13. The donor is the right cluster because between-gene correlation *within* a
   donor is large in real data and any per-gene interval would be fictitiously narrow. The alternative
   — a parametric donor bootstrap — is named and **not** chosen: at `n = 3` it is more degenerate than
   the jackknife and it reintroduces a model assumption the estimator has otherwise avoided. Jackknife
   coverage at `n = 3…4` is not assumed; it is measured by V4.

3. **Between the envelope's four tabulated points: step up, do not interpolate.** `σ_gate` is rounded
   **up** to the nearest tabulated `sigma_donor` in `gate_config.OPERATING_ENVELOPE` — {0.2, 0.35, 0.5,
   0.7} — and that row's `min_donors_per_group` is applied. `σ_gate ≤ 0.2` uses the 0.2 row;
   `σ_gate > 0.7` is outside the envelope and there is **no extrapolation**. Power falls monotonically
   in σ across the grid, so stepping up is strictly conservative with respect to the envelope as
   declared and requires no new pre-registered power model. The rule **reads**
   `gate_config.OPERATING_ENVELOPE` and does not restate its numbers.

   **The cost of the branch not taken is stated, because it is real.** Interpolating along Amendment
   1's analytic frontier would admit strata in the wedges — a stratum estimated at σ̂ = 0.4 with 10
   donors per group, say, which step-up rejects — but doing that honestly requires pre-registering a
   power model *between* the tabulated points and validating it, which is a second grid of work.
   Step-up loses those strata and requires nothing; the loss will be visible in the manifest, since
   `envelope_max_sigma_supported` already exposes the wedge. If real strata turn out to cluster there,
   interpolation is a future amendment with its own validation — not a decision to be taken at analysis
   time.

4. **Four verdicts, and only four.** `inside` | `outside_underpowered` (`min(n_A, n_B)` below the
   stepped-up row's requirement) | `outside_above_envelope` (`σ_gate > 0.7`) | `indeterminate` (the
   closed list in Change 2). `indeterminate` is an exclusion, never a pass, and is counted in D4.

5. **`n` for the envelope is `min(n_A, n_B)`**, per 1.1.

6. **The audit column.** `sqrt(s0²)·ln 2` is written into the manifest row beside every estimate,
   together with the realised `ā`, the median `L` and the technical share, so a reader can see what the
   correction did and can locate the stratum relative to Correction 1's threshold. Per Correction 1 it
   is an **audit quantity of unknown error sign** and gates nothing in either direction.

7. **The pooling caveat rides along.** Every emitted estimate carries the freeze's D3 state —
   `pooled = unresolved`, donor pseudobulk a lower bound on the correct replication unit. A σ̂ from a
   stratum whose donors may share a lane is a σ̂ of the wrong unit, and no verdict above repairs that.

8. **Admission remains a separate act.** `fill_sigma_columns` fills `sigma_donor_estimate`,
   `envelope_min_donors_per_group` and `envelope_membership`, and leaves `admitted_to_sweep = False`.
   Two of the freeze's four blockers — `integer_check` and `frozen_universe_size` — are untouched by
   this amendment and still stand on all 251 rows.

---

### Change 5 — the validation grid, and its PRE-DECLARED criteria V1–V9

The generator for well-specified cells is `synthetic/oracles.py::simulate`, used **verbatim**:
`oracles.py` is the frozen correctness specification of the engine and is not modified. The
misspecification generators of Tier C live in `scripts/sigma_probe.py`, following the precedent that
the Amendment 2 selection grid's generative arms lived in `scripts/pb_calibration_probe.py`.

**Seeds.** Development runs on the disclosed range `seed ∈ [1, 999]`. The confirmatory grid runs on a
range that cannot overlap it: replicate `r` of cell `i` uses

```
seed(i, r)  =  1000 · (20260816 + i)  +  r ,        0 ≤ i < 176 ,  0 ≤ r < 1000
```

so the smallest confirmatory seed is 20 260 816 000. Cells are indexed in the tabulated order below,
Tier A then Tier B then Tier C, from zero.

**Tier A — verdicts and coverage.** `σ ∈ {0, 0.1, 0.2, 0.275, 0.35, 0.425, 0.5, 0.6, 0.7, 0.85}` ×
donors per group `∈ {3, 4, 6, 8, 13, 24}`; fixed at 1000 genes, 300 cells per donor, φ = 0.2, simulator
defaults otherwise. **60 cells × 200 replicates = 12 000 simulations.** The σ grid deliberately
includes 0, three points *between* the envelope's tabulated values (0.275, 0.425, 0.6) where step-up
does its work, and one point above the envelope (0.85); the donor grid includes 3 and 4, where the
jackknife is weakest.

**Tier B — bias and RMSE against cells per donor, dispersion, and universe width.** 50 replicates per
cell.

* `σ ∈ {0, 0.2, 0.35, 0.5}` × `n ∈ {4, 8}` × cells `∈ {10, 30, 100, 1000, 3000}` at φ = 0.2 — **40
  cells.** This spans and exceeds the 10 … 1000 range the inclusion gate admits and reaches the frozen
  list's real ceiling of 6671.5 median cells per donor from below.
* `σ ∈ {0, 0.2, 0.35, 0.5}` × `n ∈ {4, 8}` × `φ ∈ {0.05, 0.8}` × cells `∈ {30, 300}` — **32 cells.**
* `σ ∈ {0, 0.2, 0.35, 0.5}` × `n ∈ {4, 8}` × **`n_genes = 15 000`** × cells `∈ {100, 300}` — **16
  cells.** This is the **compression block**, and its design follows Correction 1 rather than intuition.
  The instinctive way to probe the `+1` compression is to make the data *shallow*; Correction 1's table
  shows that moves away from the effect, because CPM is invariant to a global scale on gene means while
  `1/T` is not. What sets the attenuation is how many genes share the million. A 15 000-gene universe
  is also the realistic one: `gate_config.ORACLE_SIM`'s 1500 genes is roughly an order of magnitude
  narrower than a real frozen universe, which is the structural reason Amendment 3's claim looked safe.

**88 cells × 50 replicates = 4400 simulations.**

**Tier C — misspecification: regimes the simulator was not built for.** 100 replicates per cell.

* **M1, gene-varying σ.** `σ_g` log-normal with median `∈ {0.2, 0.35, 0.5}` × heterogeneity
  `∈ {0.4, 0.8}` × `n ∈ {4, 8, 13}` — 18 cells. Each replicate **additionally** measures the `ebayes`
  arm's sensitivity on the pre-registered oracle (log2FC = 1.0, K = 200); that measurement is the ground
  truth for Change 3's functional-selection rule and exists for no other purpose.
* **M2, unequal cells per donor.** `n_c ~ LogUniform(18, 3000)` — 2 cells.
* **M3, donor outlier.** One donor per group displaced by 3σ — 2 cells.
* **M4, zero inflation.** Extra dropout `π ∈ {0.15, 0.3}` — 4 cells.
* **M5, unbalanced groups.** 3 v 8 and 4 v 13 — 2 cells.

**28 cells × 100 replicates = 2800 simulations.** Expected directions are declared in advance, so a
surprise is recognisable as one: M4 — E1 underestimates `v` and therefore **overestimates** σ (safe),
E2 approximately unbiased; M3 — σ̂ up (safe), with the LODO influence diagnostic lighting up; M2 — E2
carries, E1 degrades.

**Total: 176 cells, 19 200 simulations.** The budget is deliberately not stated as a promise: 160 of the
176 cells cost about a second each, but the 16 compression cells carry fifteen times the genes of every
other cell and will dominate the wall clock by a wide and, at the time of writing, unmeasured margin.
That is precisely why the driver (`scripts/run_sigma_grid.py`) is **resumable and takes a time budget**,
after `scripts/run_test_selection_grid.py`, and why the grid runs on the PC as the 146-cell grid did
rather than in CI. The summary lands in `pilot/sigmaval/summary.{csv,json}` and is committed; per-cell
JSON is not, following `pilot/testsel/`. A run that exhausts its budget short of the full grid is
reported as an incomplete grid with the missing cells named — never as a grid whose scope was the cells
that finished.

**Metrics.** Bias, relative bias and RMSE of `σ̂_rms`; one-sided coverage of the UCB; and — the deciding
metric — **the frequency of a wrong verdict in each direction, separately**. In Tiers A and B the
verdict's ground truth is obtained by applying the *same* step-up procedure to the **true** σ, so
correctness is measured against the envelope as declared and not against an unstated power model. In
Tier C arm M1 the ground truth is the measured sensitivity of Change 3.

#### The criteria. Every number below is fixed by this commit and none is to be moved after the run

**V1 — bias, core.** Every Tier-A cell with `n ≥ 4`: `|mean σ̂_rms − σ| ≤ max(0.03, 0.10·σ)`.
*On failure*: the estimator is not fit to gate admission; no `GATE_FUNCTIONAL` is written and the
failing region is reported. It is **not** narrowed to the cells that happened to pass.

**V2 — bias, thin end.** Tier-B cells at 10 and 30 cells per donor: `|bias| ≤ max(0.05, 0.15·σ)`.
*On failure at 10 cells only*: the pre-declared consequence fires automatically — strata whose median
cells per donor is below 30 receive `indeterminate` (Change 2, reason 5). The threshold is not relaxed;
the estimator's domain is narrowed, and the narrowing is written into the gating rule rather than
invented after the fact. *On failure at 30 cells*: as V1 — the estimator fails, because `[30, 100)` is
the frozen list's modal cells-per-donor bin.

**V3 — RMSE.** Reported for every cell; **not binding**. Its decision-relevant projection is V5/V6, and
a second binding threshold on the same quantity would be double counting.

**V4 — UCB coverage.** Every Tier-A cell with `n ≥ 4`: `P(σ_UCB ≥ σ) ≥ 0.85` against the nominal 0.90;
pooled over the tier, `≥ 0.88`. The Monte-Carlo SE of a coverage estimate at 0.90 over 200 replicates is
`sqrt(0.9 × 0.1/200) = 0.021`, so 0.85 sits a little over two SE below nominal — tight enough to catch a
broken interval, loose enough not to fail on Monte-Carlo error. `n = 3` is **reported and not binding**:
the envelope's easiest row already demands 4 donors per group, so a 3 v 3 stratum is never admitted
whatever its interval does.
*On failure*: the jackknife is not delivering the coverage the gate assumes, the gate is not usable as
specified, and the alternative named in Change 4 item 2 goes on the table as a future amendment with its
own validation.

**V5 — false admission. The main criterion.** On deep-outside cells — those whose true σ is at least one
tabulated step above what the cell's `n` supports — `P(verdict = inside) ≤ 0.05` in **every** such cell
and `≤ 0.02` pooled.
*On failure*: the gate admits strata the arm is invalid for, which is spec §10 risk 1. The estimator does
not gate admission. No threshold moves.

**V6 — the price of over-exclusion.** On deep-inside cells — those where `n` suffices even for the next
tabulated σ upward — `P(inside) ≥ 0.80` in every such cell and `≥ 0.90` pooled. Borderline cells are
**reported and not binding**: an upper-confidence-bound gate is *supposed* to cut them, and binding there
would penalise the construction for working.
*On failure*: the instrument is honest but unusable, excluding nearly everything. That is reported as
such, and the sweep-feasibility question (freeze §6) is answered in the negative for reasons of
instrument conservatism rather than of real σ. It is not repaired by loosening the UCB.

**V7 — misspecification.** In each Tier-C arm, wrong-admission must be `≤ 2 ×` the corresponding
well-specified Tier-A cell.
*On failure*: **the threshold does not move, and the failing arm's pre-declared signature becomes a
blocking `indeterminate` flag** (Change 2, reason 4) — LODO influence for M3, E1/E2 divergence for M4,
cells-per-donor coefficient of variation for M2, realised `d0` of the corrected values for M1. The
estimator's domain shrinks to the strata that do not show the signature, and the shrinkage is reported as
a result. **This is one of the two failures whose consequence is a new blocking flag rather than a moved
threshold.**

**V8 — double use of the data.** Among Tier-A replicates whose verdict is `inside`, the `ebayes` arm's
fresh-null false-positive rate — the fraction of replicates with at least one BH rejection — must be
`≤ α + 2 · MC SE`, with `MC SE = sqrt(α(1 − α)/n_inside)` computed from the realised count of `inside`
replicates. A stratum selected for having drawn a *quiet* realisation of its donors must not thereby
become anti-conservative conditional on that selection.
*On failure*: the selection is not benign and admission must be conditioned differently. **This is the
second failure whose consequence is a new blocking flag rather than a moved threshold**: the estimator
does not gate admission at all until a selection-aware rule is pre-registered in a further amendment.

**V9 — the correction is not a no-op, and Amendment 3's bound is demonstrated to reverse.** Two clauses,
both binding.

* **V9a.** At the Tier-B cell (1000 genes, 30 cells per donor, φ = 0.2, σ = 0.35, `n = 8`): the bias of
  `sqrt(s0²)·ln 2` must be at least **3 ×** the absolute bias of `σ̂`. Over 16 development seeds that
  cell puts `sqrt(s0²)·ln 2` at 0.3886 (sd 0.0031), a bias of +0.0386, so the clause demands
  `|bias σ̂| ≲ 0.013`. **That is roughly four times tighter than V2 permits at the same cell, and the
  asymmetry is intentional** — it is stated here so a V9a failure cannot later be re-read as a V2 pass.
* **V9b.** At the Tier-B compression cell (`n_genes = 15 000`, 300 cells per donor, φ = 0.2, σ = 0.35,
  `n = 8`): the **mean over replicates** of `sqrt(s0²)·ln 2` must be `< σ`. This is the obligatory
  demonstration that Amendment 3's "upper bound" stops being an upper bound in the regime real universes
  occupy. The full 50-replicate distribution is reported, not merely its mean. **Its expected outcome is
  disclosed rather than dramatised**: the development probes already put this cell at 0.3318 with a
  replicate sd of 0.0003 and 8 of 8 seeds below σ, so V9b is expected to pass. It is pre-registered on
  disjoint confirmatory seeds all the same, because its job is to pin a reproducible demonstration into
  a committed artifact rather than to manufacture suspense — and because a criterion written after its
  own run would be worth nothing regardless of which way it went.

*On failure of V9a*: the correction is not earning its complexity at the thin end and the estimator's
advantage over the mechanism is not demonstrated there — reported, with the estimator's fate decided by
V1/V2/V5 as usual. *On failure of V9b*: Correction 1's empirical demonstration has failed at the cell
chosen for it. Correction 1's algebra stands on the derivation and on the counterexamples in its own
table regardless — but the failure is recorded, the cell's inadequacy is diagnosed, and the demonstration
is **not** quietly moved to a cell where it works.

**Change 3's functional-selection rule** is applied mechanically to Tier C arm M1 together with Tier A,
and its outcome is fixed in Part B. It is not a criterion and cannot pass or fail; it selects, or it
eliminates everything.

---

### What this does NOT settle

* **The real anchor §8(b) asks for.** This amendment validates an estimator *in the simulator's
  coordinates*. It does not pin `sigma_donor` to a real empirical mean–dispersion / donor-variance
  trend, and a validated estimator of the wrong model's parameter is still wrong. Oracle (d), Mathys
  2019 (§8(d)), remains **binding and unrun**, and is now additionally blocked on a ROSMAP data-use
  agreement begun 2026-08-16 and on a second, unwritten loader (freeze §8).
* **The pooling question.** `pooled = unresolved` on 251 of 251. Donor pseudobulk remains a lower bound
  on the correct replication unit, the D3 gold-standard claim cannot be made on any stratum in the list,
  and every σ̂ this machinery produces inherits that caveat. Nothing here improves it.
* **Whether any real stratum lands inside the envelope.** The freeze's §6 tiers stay a scenario analysis
  until the estimator runs on real data. At σ ≈ 0.5 the surviving set is 5 of 12 datasets, below §1's own
  8–12 floor, and a negative answer remains a live outcome of this study rather than a failure mode to be
  designed around.
* **Bimodal or heavy-tailed `σ_g`.** A single scalar cannot represent it (Change 3). A stratum with a low
  RMS and a heavy tail is admitted at a σ its tail genes are not covered by; the tail metric makes that
  partly visible and nothing here fixes it.
* **Interpolation inside the envelope's wedges** is not provided, and neither is any extension below
  σ = 0.2. Both are named as possible future amendments with their own validation grids, and neither may
  be improvised at analysis time.
* **Which functional gates.** That is Part B's, decided by Change 3's rule and by nothing else.
* **Whether the estimator passes at all.** Part A claims no result. If V1, V2, V4, V5, V6 or V7 fails,
  admission stays closed and the failure is reported here.
* **A2 stratification remains deferred** (Amendment 2 Change 6); the naive arm's floor is still
  cell-count confounded and guarded only by a range check.
* **The GO/NO-GO decision is not taken.** An estimator that gates admission licenses a measurement, not a
  conclusion.

*Author attests: every figure in this entry was derived while writing it. The ten probes under
Correction 1 were run against this repository's own arm code on disclosed development seeds; the algebra
of Correction 1 and of Change 1 was checked term by term against `synthetic/oracles.py` and
`src/pbcheck/methods/moderated.py` rather than quoted from an earlier entry; the frozen list's
distributional figures are taken from `docs/PREREGISTRATION_STRATUM_LIST.md` as committed. No real data
informed this amendment. The confirmatory grid has not been run and no criterion V1–V9 has an outcome.
Every probe row is read only as far as its seed count allows, and where a number is a single realisation
rather than a mean it is said so. The `L ≷ 10⁶/(2σ²)` threshold that orders the probe table is labelled
a heuristic, its two cancelling approximations are named, and it is read by no code. The correction to
Amendment 3 is stated as the reversal of a published claim of this log, not as a clarification of it.*
