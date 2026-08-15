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
