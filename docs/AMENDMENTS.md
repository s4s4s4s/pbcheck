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
and the numeric criteria **V1–V11** that decide whether the estimator is fit to gate stratum
admission. **Part B is a dated addendum written after the run**: the V1–V11 outcome table, the
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

### What an adversarial read of this entry changed before it was committed

An entry whose whole claim is *"these criteria were written before the numbers"* has to survive being
read by someone trying to break it. This one was, and it did not survive intact. The changes are listed
here rather than folded in silently, because a pre-registration that quietly improved between drafts is
indistinguishable from one that was written to fit — and because two of the defects below are of a kind
this log has charged its own predecessors with.

**One was blocking.** Reason 3's technical-share test fired **deterministically at σ = 0**, where
Change 1.6's own algebra makes the share 1 by construction, so V6 would have failed on three
deep-inside cells before any data existed and V6's declared consequence is a negative feasibility
verdict for the whole sweep. Two of this entry's definitions, meeting, would have pre-committed the
study to its own answer. Repaired in Change 2 with an exemption keyed to the **upper confidence bound**,
derived there, with its cost located in V5 and every other reason and criterion checked against it one
at a time.

**Seven were substantive.**

1. **The de-attenuation divisor — this entry's central novelty — was falsifiable by nothing.** No
   criterion bound tightly enough, at a wide enough universe, for an estimator that never divided by
   `ā²` to fail. **V10** is added, on the compression block, sized against the measured magnitude of
   the divisor there.
2. **The df-weight scheme was correct and unidentifiable by the grid**: the difference from `1/n`
   weights vanishes in expectation on every balanced-cells design, and every cell in the grid was one.
   Tier C arm **M6** is added — group imbalance *and* group-asymmetric cells per donor — under a
   binding bias criterion, **V11**.
3. **V7 was not computable as written**: four of its five arms had no declared coordinates, its M1 arm
   compared two incommensurable ground truths, its ratio had no floor and would have failed one run in
   five on a correct estimator, and three of the four blocking misspecification signatures had no
   declared threshold —
   which would have left Part B to invent the number that decides how many real strata get excluded.
   All four repaired, in Change 5 and in V7.
4. **Change 3's functional-selection rule rewarded liberality** — "maximise pooled deep-inside
   admission" selects whichever functional reads lowest — and carried an escape hatch letting Part B
   argue a reformulation after seeing the numbers. Screen widened, conservative tie-break and declared
   tie tolerances added, escape hatch deleted.
5. **deep-inside / deep-outside / borderline were prose and ambiguous in three places**, and ten
   vacuous `n` = 3 cells were carrying 30 % of V5's pooled denominator. All 60 Tier-A cells are now
   tabulated with their class; vacuous cells are excluded from pooled denominators in V4, V5 and V6.
6. **One quantitative justification was simply false** — the split-half divisor's, at the inclusion
   gate's floor of 10 cells. Corrected in Change 2, with the true magnitudes.
7. **Correction 1's evidence existed only as prose**, with no code, no artifact and no named seeds —
   precisely the defect Amendment 3 charged its predecessor with ("the surviving record of it is
   prose"). It is now `scripts/check_upper_bound_claim.py` and a committed artifact on named seeds, and
   the live code string that repeated the corrected claim is fixed **in this commit** rather than
   deferred to an undated Part B.

**Two more were found while checking those, and are not smaller.**

8. **E1's method-of-moments equation for `φ` was missing a term** — the pure depth-variation
   contribution `m̄² cv_ℓ²`, which is present even at φ = 0. Solving the form as drafted returns
   `φ + cv_ℓ²/(1 + cv_ℓ²)`, a **+43 %** overstatement at the simulator's own depth spread, enough to
   make the dissenting witness dissent from everything at σ = 0. Derived and corrected in Change 2,
   with the measurement beside it.
9. **56 of the 88 Tier-B cells were under no binding bias criterion at all**, because V1 was scoped to
   Tier A and V2 to 10 and 30 cells per donor, and nothing covered the gap between them — the whole
   compression block included. V2 now extends V1's band to every Tier-B cell at 100 cells per donor and
   above.

Nothing below was loosened to make anything pass. Where a criterion changed, it changed because its
*construction* was wrong; V2, V4, V5 and V6 are strictly harder than they were, V7's M1 arm becomes
binding where it was an unusable ratio, and two new binding criteria are added. The one place a rule now
permits something it did not — reason 4's exemption — is derived, its cost is named, and the risk it
moves lands on criteria that can measure it.

One live-code defect surfaced by the same work is fixed here rather than filed: `ebayes_from_pdata`'s
universe restriction was O(G²) and cost 113.7 s of a 115.5 s call at 15 000 genes. It is named where it
was found, under Correction 1's measurement.

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
6. **What is new here and was measured for it**: ten probe cells of Amendment 3's quantity, on the
   **named** development seeds 1–16, regenerated by `scripts/check_upper_bound_claim.py` into the
   committed artifact `pilot/upper_bound_check/`, and tabulated under Correction 1. They are the
   evidence that Amendment 3's claim is false, and they make the expected outcome of V9b known in
   advance, which is disclosed where V9b is stated. They are **not** the validation, which is V1–V11
   and is unrun.
7. **Three further development-seed measurements, made while revising this entry and disclosed here
   rather than only where they are used.** (a) The attenuation `ā` over the estimation stratum at 1000
   and 15 000 genes, with and without the CPM filter — this sizes **V10**. (b) The median technical
   share at `donor_sigma` = 0 — this corroborates the algebra behind reason 4's exemption, which does
   not depend on it. (c) The df-weight-versus-`1/n` gap on the group-asymmetric unbalanced designs —
   this establishes that Tier C arm **M6** can falsify the weighting choice, which no existing cell
   could. All three were made with a throwaway prototype of the estimator that is deliberately **not
   committed**, on the disclosed seed range, and none of them is validation: they size criteria, and
   the criteria are then judged on the disjoint confirmatory range. `src/pbcheck/sigma_donor.py` still
   does not exist and the confirmatory grid is still unrun.

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
* **Four committed texts repeat the erroneous claim, and they are named here rather than silently
  edited**: `src/pbcheck/census_select.py`'s `PENDING_FIELDS["sigma_donor_estimate"]` ("states it is an
  UPPER BOUND, not an estimate") and its module docstring ("that quantity is an *upper bound* on
  `donor_sigma`"); `docs/PREREGISTRATION_STRATUM_LIST.md` §6 ("an unvalidated upper bound") and §9
  item 2 ("stated in terms that the quantity is an **upper bound**"); and `docs/PILOT_FINDINGS.md`'s
  closing list ("the conversion from the fitted prior to the simulator's parameterisation is an upper
  bound that still has to be derived and validated"). All were correct to *cite* Amendment 3 and are
  wrong only by inheritance. **The live code and the dated documents are handled differently, and the
  difference is a rule rather than a convenience.**

  **The live code is corrected in this commit.** An earlier draft of this entry deferred it to Part B
  "so that the correction and its evidence land together" — but the evidence lands *here*, in
  `pilot/upper_bound_check/`, and Part B has no date. A string that a reader of the manifest sees
  today, telling them the number errs in the safe direction when this entry has just shown it does not,
  is the exact failure Correction 1 is about; leaving it in place for an undated future commit would
  reproduce the error while describing it. `census_select.py`'s pending-field text and module docstring
  now say *audit quantity of unknown error sign*, and `tests/test_census_select.py` pins that wording
  so it cannot regress — the assertion is written to reject any surviving sentence that asserts the
  bound while permitting the ones that name it as retracted.

  **The two dated `.md` snapshots are treated differently, and the difference is publication, not
  age.** The rule this repository corrects dated documents by — **live code is read as current truth
  and is fixed in place; a published dated snapshot is read as history and is corrected by the
  log** — turns on whether the document has already been published, because only then does editing it
  destroy a record someone could have read.

  `docs/PILOT_FINDINGS.md` is dated 2026-07-19 and has been public since. It **stays exactly as
  written**, and its repetition of the bound is retracted here, in the log it cites — exactly as
  Amendment 1 retracted its statistics error and left the document standing. A published dated
  document edited to match a later correction is no longer a record of what was believed, and is
  worth less than an uncorrected one.

  `docs/PREREGISTRATION_STRATUM_LIST.md` is dated the same day as this entry and **has not been
  pushed**: it carries no external timestamp, nobody has read it, and it goes out in the same push as
  this correction. There is no history to preserve, and publishing a pre-registration that asserts a
  bound its own amendment log retracts a thousand lines away would be a defect on the day of
  publication rather than a faithful record of anything. Its §6 and §9 item 2 therefore carry the
  retraction **in place**, at the point of the claim. From the push onward it is a published dated
  snapshot and the first rule governs it.

**Criterion V9 exists to demonstrate this empirically rather than leave it resting on the algebra
above.** It is pre-declared in Change 5, with both clauses and both numbers.

#### The measurement behind Correction 1 (development seeds, disclosed)

Ten probe cells, run through the arm's own code path (`build_pseudobulk` → `frozen_universe` →
`moderated.log_cpm` → `wls_two_group` → `fit_f_dist`) at 8 v 8 donors, dispersion 0.2, simulator
defaults otherwise. `L̃` is the median donor universe-restricted library size and `ā` the median of
`CPM/(1 + CPM)`; those two and median CPM are medians over seeds. Rows are ordered by `L̃/L_crit`, with
`L_crit = 10⁶/(2σ²)` from the closed form above — a quantity computed from the design, not fitted to
the results.

**These are a committed reproducer and an artifact, not prose, and the change from the draft is
deliberate.** An earlier draft of this entry tabulated these ten cells from an interactive session,
with seed counts given as "16 / 8 / 1" and the seeds themselves **unnamed**. Amendment 3 charged its own
predecessor with precisely that defect — *"the surviving record of it is prose"* — and an entry whose
opening act is to retract a published claim of this log cannot rest the retraction on numbers a reader
cannot regenerate. The cells are therefore regenerated by **`scripts/check_upper_bound_claim.py`** on
the **named development seeds 1–16, for all ten cells**, into
**`pilot/upper_bound_check/upper_bound_check_2026-08-16.json`** (with a `.csv` summary beside it), which
carries every per-seed replicate record, the argv, the commit, and the SHA-256 of each source file the
measurement passes through — so every aggregate below is recomputable from the artifact alone, and the
script's `--smoke` mode re-checks its one deviation from the arm's default call path on each run.

| universe / depth | σ | `L̃` | med. CPM | `ā` | seeds | mean `sqrt(s0²)·ln 2` (sd) | ratio to σ | below σ | `L̃/L_crit` |
|---|---|---|---|---|---|---|---|---|---|
| 1000 genes, 30 cells, `mean_log_mu = 0` | 0.35 | 6.28e4 | 464.7 | 0.9979 | 16 | 0.4216 (0.0044) | 1.205 | 0/16 | 0.015 |
| 1000 genes, 30 cells | 0.35 | 1.70e5 | 463.3 | 0.9978 | 16 | 0.3886 (0.0031) | 1.110 | 0/16 | 0.042 |
| 1000 genes, 300 cells, `mean_log_mu = 0` | 0.35 | 6.31e5 | 465.9 | 0.9979 | 16 | 0.3599 (0.0026) | 1.028 | 0/16 | 0.155 |
| 1000 genes, 300 cells | 0.35 | 1.71e6 | 465.0 | 0.9979 | 16 | 0.3536 (0.0026) | 1.010 | 1/16 | 0.420 |
| 1500 genes, 250 cells | 0.35 | 2.15e6 | 312.5 | 0.9968 | 16 | 0.3533 (0.0019) | 1.009 | 1/16 | 0.526 |
| **`ORACLE_SIM`: 1500 genes, 250 cells** | **0.50** | 2.15e6 | 291.9 | 0.9966 | 16 | **0.5006** (0.0021) | **1.001** | **7/16** | **1.076** |
| 15 000 genes, 300 cells | 0.20 | 2.64e7 | 31.96 | 0.9697 | 16 | **0.1955** (0.00021) | **0.977** | **16/16** | 2.109 |
| 15 000 genes, 100 cells | 0.35 | 8.85e6 | 30.63 | 0.9684 | 16 | **0.3415** (0.00032) | **0.976** | **16/16** | 2.169 |
| 15 000 genes, 300 cells | 0.35 | 2.63e7 | 30.65 | 0.9684 | 16 | **0.3319** (0.00033) | **0.948** | **16/16** | 6.446 |
| 15 000 genes, 300 cells | 0.50 | 2.64e7 | 28.77 | 0.9664 | 16 | **0.4676** (0.00044) | **0.935** | **16/16** | 13.187 |

**Reconciliation with the draft's table, because numbers moved.** The six cells the draft ran at 16
seeds reproduce to **every printed digit** — 0.4216 / 0.3886 / 0.3599 / 0.3536 / 0.3533 / 0.5006, and
0 / 0 / 0 / 1 / 1 / 7 realisations below σ — which also establishes what the draft's unnamed seeds were:
seeds 1…N. The four heavy cells were the under-sampled ones and they are the ones that changed. c09
goes 0.3318 (8 seeds) → **0.3319 (16)** and 8/8 → **16/16** below σ; c10 0.4669 (1 seed) → **0.4676
(16)**, 1/1 → **16/16**; c07 0.1955 (1) → 0.1955 (16), 1/1 → **16/16**; c08 0.3415 (8) → 0.3415 (16),
8/8 → **16/16**. The design columns moved 1–3 % because the draft reported them from **seed 1** where
they are now medians over seeds — `L̃/L_crit` 1.06 → 1.076 at `ORACLE_SIM`, 13.0 → 13.187 at the last
row — which changes the ordering nowhere and no conclusion anywhere. The draft's hedge that "no row is
read for more than its seed count supports" is retired: every row now carries 16.

**Why the draft under-sampled the four cells that mattered most, since the answer is a live defect.**
`ebayes_from_pdata`'s universe restriction rebuilt `set(pdata.var_names)` once per gene, making it
O(G²): profiled at G = 15 000, that one line took **113.7 s** of a 115.5 s replicate, 98 % of the cost,
against 0.02 s when the set is hoisted. It is fixed in `src/pbcheck/methods/moderated.py` in this
commit — a hoist that alters no formula, constant, branch or guard, verified bit-identical on a
15 000-gene replicate down to the p-values, the `log2fc` column, the gene order and every entry of the
`moderation` dict. The same construction survives at `src/pbcheck/methods/pseudobulk.py` (the
superseded DESeq2 arm) and `scripts/pb_calibration_probe.py` (the deliberately frozen reference
instrument); both are named here and neither is touched, since neither is on the path this amendment's
grid runs and the probe is pinned as frozen evidence.

Four things are read off this table and nothing further is claimed from it.

1. **The reversal is real, and it is ordered by exactly the quantity the derivation says orders it.**
   Reading down the table, `L̃/L_crit` climbs from 0.015 to 13.2 and the ratio falls monotonically from
   1.205 to 0.935, crossing 1 where the closed form says it should. The threshold is a heuristic — it is
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
   quantity reads **2.3 % to 6.5 % below the truth at every σ probed**, and **64 of 64 realisations
   across those four cells fall below σ**, with replicate sds of 0.0002 – 0.0004. `s0²` over 15 000
   genes is a very stable quantity: at 300 cells and σ = 0.35 the mean is 0.3319 against a truth of
   0.35, a gap of 0.0181 against a standard error of the mean of 0.00008 — **220 standard errors.** This
   is a systematic offset, not a sampling accident, and it is now measured at 16 seeds on all four rows
   rather than at 8, 8, 1 and 1.

These are development-seed probes, disclosed as such, on named seeds and with a committed artifact.
Together they establish that Amendment 3's claim is false — a counterexample suffices for that. The
pre-registered quantitative statement is V9's, on disjoint confirmatory seeds.

---

### Change 1 — the estimand, and the derivation term by term (§1)

#### 1.1 The estimand

> **A stratum's `sigma_donor`** is the value that, substituted into `synthetic/oracles.py::simulate`,
> reproduces that stratum's between-donor dispersion as the pseudobulk arm's own statistic sees it:
> the standard deviation, on the **natural-log** scale, of the per-(gene, donor) log-normal random
> effect `re[g,d] = exp(N(−σ²/2, σ²))`.

**A notation note, because this log has written the same object two ways.** `N(μ, τ²)` here is
*(mean, variance)*. Amendment 3 wrote this random effect as `exp(N(−σ²/2, σ))`, which is
*(mean, standard deviation)* — the `(loc, scale)` convention `numpy.random.Generator.normal` takes,
and therefore the one `oracles.py` is literally written in (`rng.normal(-0.5 * donor_sigma**2,
donor_sigma, ...)`, where the second argument is the **sd**). The two expressions denote the
identical distribution and no code changes with the choice. **This entry uses `N(mean, variance)`
throughout**; wherever a formula below is checked against a line of `oracles.py`, the code's second
argument is the square root of this entry's.

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

**But the difference is *identifiable* only in a narrower case than that sentence implies, and saying
so is what forced a grid cell to exist for it.** The two schemes differ by
`Σ_d (w_d − 1/n_d)·v_Y[d,g]` with `n_d = n_A + n_B`, and that difference is a reweighting **between
the two groups**: the df scheme gives group A a total weight of `(n_A − 1)/(n_d − 2)` where `1/n`
gives it `n_A/n_d` — 0.760 against 0.741 at 20 versus 7, and 0.9167 against 0.9091 at 100 versus 10.
The difference therefore **vanishes in expectation** whenever `E[v_Y]` is the same in both groups,
which it is on every balanced-cells design: every Tier-A cell, every Tier-B cell, M2, and both of
M5's unbalanced cells. It is non-zero only when the groups are unbalanced **and** the per-donor
technical variance differs systematically between them — i.e. when the two arms have different cells
per donor or different counts per cell, which is the ordinary situation in real data, where a small
disease arm is rarely matched cell for cell against a large control arm. Tier C arm **M6** (Change 5)
exists to create exactly that combination and **V11** binds on it. Without M6 the df weighting would
have been correct, derived, and untested by anything in the grid.

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
  The unequal-halves form is used because it is **exact and free**, and the quantitative justification
  an earlier draft of this entry gave for it was **false and is withdrawn here rather than quietly
  dropped**. That draft said the form was "mandatory, not an optimisation" because "at the inclusion
  gate's floor of 10 cells the difference between `/4` and the exact divisor is not a rounding error".
  Ten is even. At even `n_c` the exact divisor **is** 4 and the difference is exactly zero. Written
  out: with `n1 = n2` the divisor `n_c(1/n1 + 1/n2)` is 4 identically, and with `n_c` odd and
  `n1 = (n_c − 1)/2`, `n2 = (n_c + 1)/2` it is `4n_c²/(n_c² − 1)`, a relative excess of
  `1/(n_c² − 1)`: **+0.833 % at `n_c` = 11** — the smallest odd count above the `MIN_CELLS` = 10
  floor and the largest deviation reachable there — then +0.595 % at 13, +0.446 % at 15, +0.104 % at
  31, and 0 at every even `n_c`. Using `/4` at odd `n_c` would bias `v̂` **low** by that amount and
  therefore σ̂ **high**, which is the safe direction, and the split-half is the cross-check and never
  the primary `v̂`, so nothing in the gate would have broken. **The exact form is kept because it is
  right, not because the error would have mattered**; the claim that it would have was arithmetic
  written to make a design choice sound forced, which is the failure mode this log exists to catch,
  and it is corrected here for the same reason Correction 1 is. The two-line proof of the `1/cells`
  scaling is carried in the function's docstring, beside the code that depends on it.

E2's per-donor `v̂` does not depend on the number of donors, so it does not degrade at `n = 3`; its
residual error is O(1/n_c) *relative*. Its cost is `B_BOOT` sparse mat-vecs per donor, and per-donor
slices and the resampling mat-vec are computed without densifying, after `io_counts._iter_value_blocks`.

**E1 — `nb_plugin`, the misspecification detector.** `v̂_Y[d,g] = (a_dg²/ln 2²) · (1/max(T_dg, 1) +
φ̂ · r_d)`, with `φ̂` estimated from the per-cell counts inside each donor by moments, pooled robustly
over genes and donors, from

```
Var_c(x_cg)  =  m̄_g  +  m̄_g² · [ φ (1 + cv_ℓ²)  +  cv_ℓ² ]
        =>     φ̂  =  [ (Var_c(x_cg) − m̄_g)/m̄_g²  −  cv_ℓ² ] / (1 + cv_ℓ²)
```

**with `cv_ℓ²` the squared coefficient of variation of the donor's per-cell library sizes, and the
trailing `+ cv_ℓ²` inside the bracket is not optional.** An earlier draft of this entry wrote this
moment equation as `Var_c(x_cg) = m̄ + φ m̄²(1 + cv_ℓ²)`, dropping that term. It is the pure
depth-variation contribution, `Var(E[x | depth]) = m̄² cv_ℓ²`, which is present **even at φ = 0**: with
`x | depth ~ NB(mean = depth·m, var = depth·m + φ(depth·m)²)` the law of total variance gives
`Var(x) = E[Var(x|depth)] + Var(E[x|depth]) = m̄ + φ m̄²(1 + cv_ℓ²) + m̄² cv_ℓ²`. Solving the dropped
form returns not φ but `φ + cv_ℓ²/(1 + cv_ℓ²)`, and at the simulator's `depth_log_sigma` = 0.3 that
offset is `cv_ℓ² = e^{0.09} − 1 = 0.09417` over `1.09417` = **+0.0861**, predicting 0.286 in place of
φ = 0.2 — a **+43 %** overstatement. Measured on development seeds 1–4, median over (donor, gene) pairs
with `m̄ ≥ 1`: the dropped form returns **0.281–0.284** against a truth of 0.2 at every geometry probed
(+41 % to +42 %, the realised `cv_ℓ²` running a little under its nominal 0.0942), the corrected form
**0.195**. The error is not cosmetic: at the reference geometry (1000 genes, 300 cells) the `1/T` and
`φ r_d` terms are comparable, so a 43 % overstatement of φ inflates `v̂_E1` by about 16 %, and at
`donor_sigma` = 0 that is enough to drive the pre-clip aggregate `M` negative in **8 of 8** seeds and
make E1 return exactly 0 — a dissenting witness that dissents from everything. **A second limit of E1
is recorded at the same time**: at 10 cells per donor `Var_c` over 10 cells is too noisy for a
median-of-ratios and the corrected form reads 0.130 against 0.2, so E1's φ̂ is not usable at the
inclusion gate's own floor. That matters beyond E1's own accuracy, because E1/E2 divergence is a
candidate blocking signature (reason 5): a thin stratum could be flagged for E1's noise rather than for
misspecification. It is why the divergence flag blocks **only** if V7's M4 arm fails, and why M4's cells
are run at the reference geometry's 300 cells per donor, where φ̂ is sound. E1 is a detector, not an
estimator, and it is now declared not to be one at the thin end either.

It is nearly free, and it is **not** primary for a reason stated rather than discovered: under zero inflation
it underestimates `v` and therefore overestimates σ, which is safe; but under within-donor substructure
that NB cannot see it also underestimates `v`, and there the direction is **not** guaranteed. Its role
is that of a dissenting witness — `|σ̂_E1 − σ̂_E2| > 0.05` raises an `nb_misspecification` flag, always
reported, and blocking only if V7's M4 arm fails.

**E3 — `plateau`, reported and never primary.** Restrict to genes whose predicted technical share is
negligible (`Σ_d w_d v̂_Y[d,g] ≤ ε · median_g(s2_g)` with `ε = ESTIMATION_PLATEAU_EPS` = **0.05**,
declared in the constants table below rather than left as a free symbol, and `ā_g ≥ 0.95`) and read
`s2_g · ln 2² ≈ σ²` off directly, located by
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
| `E1_E2_DIVERGENCE_FLAG` | 0.05 | in σ units; the flag above, and — if V7's M4 arm fails — the blocking threshold of reason 5's E1/E2 signature. One number, two roles, and it does not move between them |
| `MAX_TECHNICAL_SHARE` | 0.8 | reason 4 below, with the exemption stated there |
| `ESTIMATION_PLATEAU_EPS` | 0.05 | E3's `ε`. A gene joins the plateau read-off only if its df-weighted technical variance is ≤ 5 % of the median gene's total variance. This is a *different kind of number* from `MAX_TECHNICAL_SHARE`, not a tighter setting of it: 0.8 is the point past which a **stratum's** corrected estimate is mostly correction, 0.05 is the point below which a **gene** needs no correction at all. E3 gates nothing, so this constant enters no verdict |
| `CELLS_PER_DONOR_CV_BLOCK` | 1.0 | reason 5's M2 signature; derived under reason 5 below |
| `IMBALANCE_RATIO_BLOCK` | 0.35 | reason 5's M6 signature; derived under reason 5 below |
| `SIGNATURE_QUANTILE` | 0.99 | the well-specified quantile that sets reason 5's M3 and M1 thresholds (the M1 one at `1 − 0.99` = 0.01, since low `d0` is the extreme end); reason 5 below |

**The pre-declared reasons for an `indeterminate` verdict.** A stratum is `indeterminate` — excluded
from the sweep and **counted, with its reason, in the D4 excluded-strata statistic** — if and only if
one of the following holds. **They are evaluated in the order listed and the first one that fires is
the recorded reason**, so that D4 carries exactly one reason per excluded stratum and so that a test
which needs a quantity an earlier reason would have found missing is never reached with that quantity
absent. (An earlier draft of this entry left the list unordered, which left the recorded reason
undefined whenever two of them held at once, and made reason 4's exemption below unevaluable because
it reads a quantity reason 3 exists to certify.) The list is closed as of this commit; adding to it
later is an amendment.

1. Fewer than `MIN_ESTIMATION_GENES` = 100 genes pass the estimation-stratum filters.
2. Median cells per donor below a floor, **if and only if V2 fails at 10 cells**, in which case the
   floor is 30 and is written into the rule automatically rather than invented afterwards.
3. Degenerate jackknife: some leave-one-donor-out replicate fails to produce a finite functional, or
   `SE_jack` is not finite.
4. **Median technical share above `MAX_TECHNICAL_SHARE`, except where the upper confidence bound has
   already placed the stratum in the envelope's easiest row.** Precisely, both of

   ```
   median_g( Σ_d w_d v̂_Y[d,g] / s2_g )  >  MAX_TECHNICAL_SHARE = 0.8
   σ_gate = UCB_90(σ_f)                 >  min( row["sigma_donor"] for row in OPERATING_ENVELOPE )
   ```

   The second line reads `gate_config.OPERATING_ENVELOPE` and does not restate its numbers, exactly as
   Change 4 item 3 does; as the envelope stands that lowest tabulated value is 0.2.
5. A Tier-C misspecification signature that V7 or V11 converts into a blocking flag. The **set** of
   candidate signatures and the **threshold rule for each** are both fixed now — see the table below —
   and *which* of them block is decided mechanically by which Tier-C arms fail.

**Why reason 4 carries an exemption, and why the exemption is not a loosened threshold.** As the
technical-share reason stood in the earlier draft of this entry — unexempted, at 0.8, and numbered 3
there — **it fired deterministically at `σ = 0`, for a reason that has nothing to do with the
estimator.** Set σ = 0 in 1.6: the σ² term vanishes identically and

```
E[s2_g]  =  Σ_d w_d · v_Y[d,g]
```

*exactly*. The technical share is therefore **1 in expectation at σ = 0 by construction**, for any `v̂`
that is consistent for `v_Y` — that is not a symptom of an unreliable estimate, it is the definition of
a stratum with no donor effect. Tier A contains σ = 0 at every donor count, and by the classification
table in Change 5 three of those cells — `n` = 8, 13 and 24 — are **deep-inside**, where V6 demands
`P(inside) ≥ 0.80`. The realised probability would have been 0 in all three, V6 would have failed, and
V6's declared consequence is that "the sweep-feasibility question (freeze §6) is answered in the
negative". **Two of this entry's own definitions, meeting, would have pre-committed the study to a
negative feasibility verdict before a single simulation ran.** That is a defect in the construction of
the rule, not a discovery about the estimator, and it is repaired before the entry binds.

**The algebra above is sufficient, and it was measured anyway.** On development seeds 1–8 with the E2
bootstrap `v̂`, at 1000 genes and 300 cells per donor, the median technical share at σ = 0 reads
**1.0455** at 8 v 8, **1.0292** at 13 v 13 and **1.0124** at 24 v 24, and **1.0153** at 8 v 8 with the
thin-donor floor's 10 cells — **8 of 8 seeds above 0.8 in every one of them**, the lowest single
realisation anywhere being 0.9626. The same cells at σ = 0.1 read 0.1774 and at σ = 0.2 read 0.0515,
with **0 of 8 above 0.8**. The threshold is not a close call in either direction: it separates σ = 0
from σ ≥ 0.1 by a factor of five with no overlap, and the failure is confined to σ = 0 exactly. The
excess above 1.0 is not estimator bias either — it is the median/mean gap of `s2_g ~ σ² χ²_d/d`, whose
`1/median(χ²_d/d)` is 1.0495 / 1.0284 / 1.0147 at `d` = 14 / 24 / 46 against the 1.0455 / 1.0292 /
1.0124 measured. And the **clipping** bias at σ = 0, which V1 binds on at its absolute leg of 0.03,
measures 0.0031 ± 0.0034 in σ units at 300 cells per donor, an order of magnitude inside it.

**The repair is not to move 0.8, and no threshold in this entry moves.** Reason 4 exists to catch *"the
estimate is mostly its own correction, therefore unreliable"*. That reading is right when the corrected
quantity has to carry a verdict. It is wrong when the **conservative** quantity has already placed the
stratum in the envelope's easiest row: if `σ_gate ≤ 0.2` then step-up (Change 4 item 3) selects the 0.2
row whatever else is true, the only question left is whether the stratum has 4 donors per group, and no
additional precision in σ̂ can change any verdict the rule can reach. Refusing it buys nothing. It also
costs the one control the whole simulator is built around: `oracles.py` calls `donor_sigma = 0` "the
single most important control", the falsification control that separates a pseudoreplication signal
from a power artifact, and an instrument that answers `indeterminate` on it is not being conservative,
it is broken. **A high technical share at genuinely low σ is the correct reading of the data, not an
unreliable one**, and the direction that actually endangers the study — an estimate that is too *low* —
is carried by the UCB gate of Change 4 item 1, not by this reason.

**The exemption reads the UCB and never the point estimate**, for the same reason the gate does: the
point estimate is not the conservative quantity, and an exemption keyed to it would be an exemption
keyed to optimism.

**What the exemption costs, stated rather than waved past.** It admits — to the 0.2 row and to no other
— strata whose `v̂` has **over**-subtracted: where the technical correction is too large, σ̂² is pushed
toward or below zero, and the UCB follows it down through the very inequality that grants the
exemption. That is a real failure mode and it is the dangerous direction. It is not left unmeasured. It
is exactly what **V5** counts: every deep-outside Tier-A cell with `n ≥ 4` is a cell where an
over-subtracting estimator reads low, clears the exemption and is admitted, and V5 bounds that at 0.05
per cell and 0.02 pooled over 19 cells and 3800 replicates. **V4** bounds the coverage of the interval
that the exemption now leans on more heavily than before. The exemption therefore moves risk off a rule
that could not measure it and onto two criteria that can, which is a different act from relaxing it.

**Interaction with every other reason and every criterion, checked one at a time rather than asserted.**

* **Reason 1** (gene count) is upstream, reads neither `v̂` nor the UCB, and is unaffected.
* **Reason 2** (thin donors) is a property of the design, fires before any estimate exists, unaffected.
* **Reason 3** (degenerate jackknife) is unaffected, and the exemption *depends* on it having passed —
  which is why the order above puts it first. A stratum with no finite `SE_jack` has no `σ_gate` to
  compare against 0.2 and is recorded under reason 3, never under reason 4.
* **Reason 5** (misspecification signatures) is unaffected: none of its five statistics reads the
  technical share or the exemption.
* **V1 / V2 / V3 / V9 / V10 / V11** are bias, RMSE and ablation criteria on σ̂ and on
  `sqrt(s0²)·ln 2`. None of them reads a verdict, so none of them changes.
* **V4** (coverage) reads the interval, not the verdict; unchanged. σ = 0 is separately vacuous there,
  for an unrelated reason given under V4.
* **V5** (false admission) is where the exemption's cost lands, above.
* **V6** (over-exclusion) is what the repair fixes; the σ = 0 cells enter it as ordinary deep-inside
  cells rather than as three guaranteed zeros.
* **V7** (misspecification ratios) is a ratio between a Tier-C cell and a well-specified Tier-A cell,
  and the exemption applies identically on both sides, so it cannot manufacture a Tier-A/Tier-C
  asymmetry. None of M2, M3 or M4 is run at σ = 0, so no V7 arm is exercised by the exemption at all.
* **V8** (double use of the data) counts `inside` replicates, and the exemption adds the σ = 0
  replicates to that pool. At σ = 0 the `ebayes` arm's fresh null is the easiest case it ever faces, so
  this **dilutes** V8 rather than stressing it. V8's threshold does not move; the dilution is recorded
  under *What this does NOT settle*, together with the rest of V8's denominator problem, and a
  borderline-only companion figure is reported beside it.
* **Change 3's functional-selection rule** reads V5's and V6's statistics, and the repair changes both
  — in the direction that matters. Without it, all four candidate functionals would have scored an
  identical 0 on the three σ = 0 deep-inside cells, the rule's objective would have been decided
  entirely by the cells the defect happened not to reach, and the tie structure the rule turns on would
  have been an artifact of the unexempted technical-share rule.

**How σ = 0 enters V5 and V6, stated explicitly because it is a boundary case in both** (and in V1 and
V4, which are given here too so the four are in one place).

* **V5 — not at all.** A σ = 0 cell is **never** deep-outside: at every `n ≥ 4` the envelope's easiest
  row already covers it, so it is deep-inside (`n ≥ 8`), borderline-inside (`n` = 4, 6) or, at `n` = 3,
  structurally vacuous. σ = 0 contributes nothing to V5 at any donor count.
* **V6 — binding, at three cells.** σ = 0 at `n` ∈ {8, 13, 24} is deep-inside and binding at
  `P(inside) ≥ 0.80`. This is the hardest thing the grid asks of `v̂` and the right place to ask it: at
  σ = 0 the whole of `s2` is technical, so `σ̂² = (s2 − Σ w v̂)/ā²` is an undiluted test of whether `v̂`
  is right, with no donor variance for an error in it to hide inside. σ = 0 at `n` ∈ {4, 6} is
  borderline-inside, reported and not binding; at `n` = 3 it is vacuous.
* **V1 — binding, at the absolute leg.** The tolerance at σ = 0 is `max(0.03, 0.10 × 0) = 0.03`. The
  final scalar is clipped once at zero (Change 3), so at a true σ of 0 the estimator's error is
  one-sided by construction and V1 is a one-sided criterion there; the pre-clip aggregate `M` is
  reported beside σ̂ so the clipping's own contribution is visible rather than folded into the bias.
* **V4 — vacuous, and therefore excluded.** `σ_gate = sqrt(max(UCB, 0)) ≥ 0 = σ` identically, so
  coverage is 1 at σ = 0 whatever the interval does. Five such cells sit at `n ≥ 4`; leaving them in
  would put five free 1.0s into a pooled coverage denominator. See V4.

**Reason 5's signatures and the threshold rule for each, declared now.** An earlier draft of this entry
fixed the *set* of signatures and left every threshold but one unstated, which would have forced Part B
— reading the outcome — to invent the number that decides how many real strata get excluded. That is
the ordering this entry exists to prevent, so each threshold is either a literal fixed here or a
quantile of a distribution this run itself commits to an artifact.

| signature | arm | the statistic, exactly | blocking threshold |
|---|---|---|---|
| LODO influence | M3 | `max_i \|θ_(i) − θ̂\| / SE_jack`, on the σ² scale, over the `n_d` leave-one-donor-out replicates Change 4 item 2 already computes | the `SIGNATURE_QUANTILE` = **0.99** quantile of the same statistic over the well-specified Tier-A replicates at the **same donor count**, with `n ≥ 4` |
| E1/E2 divergence | M4 | `\|σ̂_E1 − σ̂_E2\|`, in σ units | `E1_E2_DIVERGENCE_FLAG` = **0.05** — already declared as the reporting flag; failure makes the same number blocking and does not move it |
| cells-per-donor CV | M2 | `sd_d(n_c,d) / mean_d(n_c,d)` over the stratum's donors | `CELLS_PER_DONOR_CV_BLOCK` = **1.0**, an absolute. Tier A cannot supply a quantile here: its cells per donor are constant, so this statistic is identically 0 on every well-specified cell and no reference distribution exists. 1.0 is placed *below* M2's own design CV so the rule blocks the demonstrated failure with margin instead of only at it — `n_c ~ LogUniform(18, 3000)` has `E[n_c] = (3000 − 18)/ln(3000/18) = 582.9`, `E[n_c²] = (3000² − 18²)/(2 ln(3000/18)) = 8.7956e5`, hence sd 734.7 and **CV = 1.2605** |
| realised `d0` | M1 | `fit_f_dist`'s `d0` over the **corrected** per-gene values; low `d0` means heterogeneous | the `1 − SIGNATURE_QUANTILE` = **0.01** quantile of the same statistic over the well-specified Tier-A replicates at the same donor count, with `n ≥ 4` (the low tail, because heterogeneity drives `d0` down) |
| group-imbalance ratio | M6, via V11 | `min(n_A, n_B) / max(n_A, n_B)` | `IMBALANCE_RATIO_BLOCK` = **0.35**, the ratio of M6's *milder* design (7 / 20), so the rule blocks both demonstrated cells and everything more extreme. Its cost is counted rather than gestured at: **38 of the 251 frozen strata** sit at or below 0.35 — 14 in #5, 12 in #3, 10 in #4, 2 in #8 — and **no dataset loses all of its strata**. #3 does hold 20 v 7 (8 strata) and #5 does hold 100 v 10 (5 strata), the two designs M6 is built from |

The two quantile rules need a reference distribution that outlives the run, so **the driver writes both
quantiles, per Tier-A cell with `n ≥ 4`, into the committed `pilot/sigmaval/summary.json`**. The
threshold is then fixed by an artifact rather than recomputed by whoever needs it, which is the same
device `gate_config` uses for the gate's own constants.

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

> **Step 1 — the conservatism screen.** Discard any functional whose wrong-admission rate exceeds V5's
> threshold — `≤ 0.05` in every cell and `≤ 0.02` pooled — on **either** the 19 deep-outside cells of
> Tier A (ground truth: step-up applied to the **true** σ) **or** any heterogeneous cell of Tier C arm
> M1 (ground truth: the **measured** sensitivity of the `ebayes` arm on the pre-registered oracle,
> log2FC = 1.0, K = 200, the cell counting as `inside` iff that sensitivity is ≥ 0.60).
>
> **Step 2 — the objective, with a declared tie tolerance.** Among the survivors, compute the pooled
> deep-inside admission rate — V6's statistic, over 15 cells and 3000 replicates — and keep every
> survivor within **0.01** of the largest.
>
> **Step 3 — the conservative tie-break.** Among those, compute the pooled wrong-admission rate over
> Tier A's deep-outside cells and M1's, and keep every survivor within **0.005** of the smallest.
>
> **Step 4 — simplicity.** Among those, take the simplest, in the declared complexity order
> `rms < trimmed_rms_10 < median_log < q75`.

**Why steps 1 and 3 are wrapped around step 2, and why an earlier draft of this entry was wrong to omit
them.** That draft screened on M1 alone and then chose the survivor "maximising pooled deep-inside
admission". Maximising admission **rewards the functional that reads lowest**, because a lower σ̂ steps
up to an easier envelope row and admits more. `median_log` is a geometric-mean-like location that
Change 3's own text expects to read *low* under heterogeneity — so the rule as drafted would have
selected it **for the property that makes it wrong**, in flat contradiction of Change 4's declared
asymmetry, in every realisation where the M1-only screen did not happen to cut it first. Step 1 makes
conservatism a hard constraint rather than a quantity traded against admission, and adds Tier A's
deep-outside cells to the screen because those are where a downward-biased functional is most visible
and where the M1-only screen was blindest. Step 3 makes conservatism the tie-break as well. Step 2 then
decides only among functionals that are *already* conservative enough, and the thing admission rate is
allowed to settle is the one thing it should: which of them wastes the fewest valid strata.

**Both tolerances are declared here, before the grid runs, and both are derived rather than picked.**
0.01 is a little under two Monte-Carlo standard errors of step 2's pooled proportion at 0.90
(`sqrt(0.9 × 0.1/3000)` = 0.0055); 0.005 is a little over two of step 3's at 0.02
(`sqrt(0.02 × 0.98/3800)` = 0.0023). Both are computed **as if** the four functionals were evaluated
independently. They are not — all four are computed on the identical replicates, so their *paired*
differences carry less error than that — which makes both tolerances deliberately generous and hands
the decision to the conservative tie-break in every case where the four are close. That is the intended
direction. Amendment 2 had to declare "0.004 is a tie" *after* seeing its grid; these numbers are on
the record before this one runs.

**The rule is total, and the escape hatch is deleted.** Step 4 is a strict order over four distinct
items, so steps 1–4 return exactly one functional, or — at step 1 — none. There is no degenerate
outcome left for anyone to reformulate. The clause an earlier draft of this entry carried, allowing
Part B to report a degeneracy and argue "the non-degenerate reformulation ... in the open, as Amendment
2 did", is therefore **removed rather than bounded**: a rule that may be rewritten by the party who has
already seen the numbers is not a pre-registered rule, and the fact that Amendment 2 did it honestly
under duress is not a licence to plan on doing it again. That clause was the one place the Part A /
Part B split leaked, and it is the only thing in this entry that could have made Part B's judgement
load-bearing.

If step 1 eliminates all four, that is the result: no functional is fit for the job, no
`GATE_FUNCTIONAL` constant is written, and the failure comes back to this log.

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

   **The destination row that will carry almost every real stratum is the one row with no grid support
   at all, and the conservatism argument above does not repair that.** `OPERATING_ENVELOPE`'s σ = 0.2
   row records its own `grid_support` as *"not in the grid; Amendment 1 frontier only"*: its
   `min_donors_per_group` = 4 comes from Amendment 1's analytic power frontier, validated numerically
   to |error| < 0.033 against that same derivation, and from **no measured cell** — the 146-cell
   selection grid has nothing at σ = 0.2. It is also where nearly everything lands if the anchor is
   optimistic: freeze §6 puts **11 of 12 datasets and 227 of 251 strata** in the σ ≈ 0.2 tier, against
   7 / 150 at 0.35, 5 / 94 at 0.5 and 3 / 30 at 0.7. And step-up sends *every* stratum with
   `σ_gate ≤ 0.2` there, σ = 0 included. So the claim above should be read exactly as narrow as it is:
   **stepping up is conservative with respect to the envelope as declared, and inherits whatever the
   declared envelope is wrong about.** Where the destination is the 0.2 row, "conservative" means
   "consistent with a frontier that has never been measured". Fixing that is a *power* measurement at
   (σ = 0.2, n = 4) on the pre-registered oracle, not a `sigma_donor` estimation problem; it is not in
   this amendment's grid and this amendment does not claim it. It is named here, at the point where the
   conservatism argument is made, as the largest unvalidated dependency the membership rule has.

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

8. **Admission remains a separate act.** The function that writes these columns — `fill_sigma_columns`,
   which **does not exist at the time of writing** and is Part B's to add beside `census_select`'s
   manifest writer, alongside `src/pbcheck/sigma_donor.py`, which does not exist either — **will** fill
   `sigma_donor_estimate`, `envelope_min_donors_per_group` and `envelope_membership`, and **will leave**
   `admitted_to_sweep = False`.
   Two of the freeze's four blockers — `integer_check` and `frozen_universe_size` — are untouched by
   this amendment and still stand on all 251 rows.

---

### Change 5 — the validation grid, and its PRE-DECLARED criteria V1–V11

The generator for well-specified cells is `synthetic/oracles.py::simulate`, used **verbatim**:
`oracles.py` is the frozen correctness specification of the engine and is not modified. The
misspecification generators of Tier C live in `scripts/sigma_probe.py`, following the precedent that
the Amendment 2 selection grid's generative arms lived in `scripts/pb_calibration_probe.py`.

**Seeds.** Development runs on the disclosed range `seed ∈ [1, 999]`. The confirmatory grid runs on a
range that cannot overlap it: replicate `r` of cell `i` uses

```
seed(i, r)  =  1000 · (20260816 + i)  +  r ,        0 ≤ i < 180 ,  0 ≤ r < 1000
```

so the smallest confirmatory seed is 20 260 816 000. Cells are indexed in the tabulated order below,
Tier A then Tier B then Tier C, from zero.

**Tier A — verdicts and coverage.** `σ ∈ {0, 0.1, 0.2, 0.275, 0.35, 0.425, 0.5, 0.6, 0.7, 0.85}` ×
donors per group `∈ {3, 4, 6, 8, 13, 24}`; fixed at 1000 genes, 300 cells per donor, φ = 0.2, simulator
defaults otherwise. **60 cells × 200 replicates = 12 000 simulations.** The σ grid deliberately
includes 0, three points *between* the envelope's tabulated values (0.275, 0.425, 0.6) where step-up
does its work, and one point above the envelope (0.85); the donor grid includes 3 and 4, where the
jackknife is weakest. This geometry — 1000 genes, 300 cells per donor, φ = 0.2 — is the **reference
geometry**, and every Tier-C cell below is run at it too unless its own definition overrides a
coordinate, so that "the corresponding well-specified Tier-A cell" is always a cell that exists.

#### Every Tier-A cell, classified — because V5 and V6 turn on words that are ambiguous as prose

V5 binds on *deep-outside* cells, V6 on *deep-inside* cells, and both leave *borderline* unbound. Left
as prose those three words are ambiguous in at least three places, and the resolution decides which
cells carry the study's main criterion. All 60 are therefore classified here, by a rule stated first
and applied mechanically.

Write `step_up(σ)` for the smallest tabulated envelope σ that is `≥ σ` (undefined above 0.7),
`req(σ)` for that row's `min_donors_per_group`, and `sup(n)` for the largest tabulated σ whose row `n`
already satisfies (undefined at `n` = 3, which satisfies no row).

* **deep-inside** — `n` satisfies not only `req(step_up(σ))` but also the requirement of the **next
  tabulated σ above** `step_up(σ)`.
* **deep-outside** — the true σ is **at or above the next tabulated σ above `sup(n)`**; where `sup(n)`
  is already the top row, that means σ > 0.7.
* **vacuous** — `n` = 3. `inside` is unreachable at 3 donors per group, because the envelope's easiest
  row demands 4, so `P(inside) = 0` identically and the cell measures nothing about the estimator.
* **borderline** — everything else, split into *borderline-inside* and *borderline-outside* by the true
  verdict.

**Three consequences of that rule are exactly the three places the prose was ambiguous**, and each is
resolved in the direction that makes the binding criteria harder to pass on cells that carry
information and stops them being padded with cells that do not.

1. **σ ∈ (0.5, 0.7] has no "next tabulated σ upward"**, so `(n = 24, σ = 0.6)` and `(n = 24, σ = 0.7)`
   are **not** deep-inside. The margin the definition asks for does not exist: the only step up from the
   0.7 row is *out of the envelope*, and Change 4 item 3 forbids extrapolating there. Both are
   borderline-inside — reported, not binding under V6. The other reading would have demanded
   `P(inside) ≥ 0.80` at the envelope's own outer edge, which is precisely where an upper-confidence
   gate is supposed to be cutting.
2. **`(σ = 0.275, n = 4)` and `(σ = 0.275, n = 6)` are borderline-outside, not deep-outside.** "At least
   one tabulated step above what the cell's `n` supports" is read as *the true σ reaches the next
   tabulated value* — 0.35 — and not merely as *`step_up(σ)` reaches it*. At σ = 0.275 with 4 donors the
   estimator only has to read above 0.2 to get the verdict right; the margin is 0.075, not a tabulated
   step, and calling that "deep" would put a near-boundary cell under V5's hard 0.05 and misname what
   V5 measures. Under the other reading these two would join V5's binding set; here they are reported.
3. **All ten `n` = 3 cells are vacuous, in both directions, and are excluded from the per-cell
   requirement and from the pooled denominator of V5 and V6 alike.** Without the exclusion, eight of
   them would count as deep-outside — 8 of 27 cells, **30 % of V5's pooled weight** — as eight
   guaranteed zeros that no estimator can fail, and the remaining two (σ = 0 and 0.1) as
   borderline-outside. V4 already declines to bind at `n` = 3 for the same structural reason; this
   makes it explicit for V5 and V6 too. **A vacuous cell in a pooled denominator is a free pass, and
   free passes are what pooled denominators are for hiding.**

| true σ | step-up row | donors required | n = 3 | n = 4 | n = 6 | n = 8 | n = 13 | n = 24 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.2 | 4 | vac | bi | bi | **DI** | **DI** | **DI** |
| 0.1 | 0.2 | 4 | vac | bi | bi | **DI** | **DI** | **DI** |
| 0.2 | 0.2 | 4 | vac | bi | bi | **DI** | **DI** | **DI** |
| 0.275 | 0.35 | 8 | vac | bo | bo | bi | **DI** | **DI** |
| 0.35 | 0.35 | 8 | vac | **DO** | **DO** | bi | **DI** | **DI** |
| 0.425 | 0.5 | 13 | vac | **DO** | **DO** | bo | bi | **DI** |
| 0.5 | 0.5 | 13 | vac | **DO** | **DO** | **DO** | bi | **DI** |
| 0.6 | 0.7 | 23 | vac | **DO** | **DO** | **DO** | bo | bi |
| 0.7 | 0.7 | 23 | vac | **DO** | **DO** | **DO** | **DO** | bi |
| 0.85 | above envelope | — | vac | **DO** | **DO** | **DO** | **DO** | **DO** |

`**DI**` deep-inside (V6 binds) · `**DO**` deep-outside (V5 binds) · `bi` borderline-inside · `bo`
borderline-outside (both reported, neither binding) · `vac` structurally vacuous.

**Counts: 15 deep-inside, 19 deep-outside, 12 borderline-inside, 4 borderline-outside, 10 vacuous — 60
in all.** So V6 binds on 15 cells and 3000 replicates, whose pooled proportion at 0.90 carries a
Monte-Carlo SE of `sqrt(0.9 × 0.1/3000)` = **0.0055**; V5 binds on 19 cells and 3800 replicates, whose
pooled proportion at 0.02 carries `sqrt(0.02 × 0.98/3800)` = **0.0023**. Those two numbers are what size
Change 3's two tie tolerances, and they are computed here rather than there so the classification and
the arithmetic that depends on it sit in one place.

**Tier B — bias and RMSE against cells per donor, dispersion, and universe width.** 50 replicates per
cell.

* `σ ∈ {0, 0.2, 0.35, 0.5}` × `n ∈ {4, 8}` × cells `∈ {10, 30, 100, 1000, 3000}` at φ = 0.2 and
  **`n_genes` = 1000** — **40 cells.** This spans and exceeds the 10 … 1000 range the inclusion gate
  admits and reaches the frozen list's real ceiling of 6671.5 median cells per donor from below.
* `σ ∈ {0, 0.2, 0.35, 0.5}` × `n ∈ {4, 8}` × `φ ∈ {0.05, 0.8}` × cells `∈ {30, 300}`, at
  **`n_genes` = 1000** — **32 cells.**
* `σ ∈ {0, 0.2, 0.35, 0.5}` × `n ∈ {4, 8}` × **`n_genes = 15 000`** × cells `∈ {100, 300}` — **16
  cells.** This is the **compression block**, and its design follows Correction 1 rather than intuition.
  The instinctive way to probe the `+1` compression is to make the data *shallow*; Correction 1's table
  shows that moves away from the effect, because CPM is invariant to a global scale on gene means while
  `1/T` is not. What sets the attenuation is how many genes share the million. A 15 000-gene universe
  is also the realistic one: `gate_config.ORACLE_SIM`'s 1500 genes is roughly an order of magnitude
  narrower than a real frozen universe, which is the structural reason Amendment 3's claim looked safe.

**88 cells × 50 replicates = 4400 simulations.**

**Tier C — misspecification: regimes the simulator was not built for.** 100 replicates per cell.
**Every Tier-C cell runs at the Tier-A reference geometry — 1000 genes, 300 cells per donor, φ = 0.2,
`σ` = 0.35, simulator defaults — except where the arm's own definition overrides a coordinate**, and the
overrides are the only thing tabulated per arm below. An earlier draft of this entry gave M2 … M5 as
counts only, with no σ and no `n`, which left V7's phrase "the corresponding well-specified Tier-A
cell" undefined for four of its five arms; the coordinates are declared here so the comparison V7 makes
is a comparison between two cells that both exist.

* **M1, gene-varying σ.** `σ_g` log-normal with median `∈ {0.2, 0.35, 0.5}` × heterogeneity
  `∈ {0.4, 0.8}` × `n ∈ {4, 8, 13}` — 18 cells. Each replicate **additionally** measures the `ebayes`
  arm's sensitivity on the pre-registered oracle (log2FC = 1.0, K = 200); that measurement is the ground
  truth for Change 3's functional-selection rule and for M1's own criterion, and exists for no other
  purpose. **M1 is not under V7** — see V7 for why its truth is not commensurable with Tier A's.
* **M2, unequal cells per donor.** `n_c ~ LogUniform(18, 3000)` per donor, replacing the fixed 300;
  `σ` = 0.35, `n ∈ {4, 8}` — 2 cells. Reference cells: Tier A (0.35, 4) and (0.35, 8).
* **M3, donor outlier.** One donor per group displaced by 3σ; `σ` = 0.35, `n ∈ {4, 8}` — 2 cells.
  Reference cells: Tier A (0.35, 4) and (0.35, 8).
* **M4, zero inflation.** Extra dropout `π ∈ {0.15, 0.3}`; `σ` = 0.35, `n ∈ {4, 8}` — 4 cells.
  Reference cells: Tier A (0.35, 4) and (0.35, 8).
* **M5, unbalanced groups, equal cells per donor.** `(n_A, n_B)` ∈ {(8, 3), (13, 4)}, `σ` = 0.35 — 2
  cells. These are under **V11**, not V7: at `min(n) = 3` the 8 v 3 cell's `inside` verdict is
  unreachable, so a wrong-admission criterion on it passes vacuously for the same structural reason the
  ten `n` = 3 Tier-A cells do. Their job is the bias criterion and the weighting ablation.
* **M6, unbalanced groups with group-asymmetric cells per donor — NEW, and added because nothing else
  in the grid can falsify the df weighting of Change 1.6.** Two designs × `σ ∈ {0.2, 0.35}` — 4 cells:
  * **M6a**: `(n_A, n_B) = (20, 7)`, cells per donor `~ LogUniform(100, 900)` in the 20-donor group and
    `~ LogUniform(10, 90)` in the 7-donor group.
  * **M6b**: `(n_A, n_B) = (100, 10)`, the asymmetry **reversed** — `~ LogUniform(10, 90)` in the
    100-donor group and `~ LogUniform(100, 900)` in the 10-donor group — so that a sign error in the
    weighting cannot pass both.

  The two donor ratios, 7/20 and 10/100, are the frozen list's own two most skewed designs (#3 at 20 v 7
  and #5 at 100 v 10). Both σ values are run because the weighting error enters through
  `Σ_d w_d v̂_Y`, whose weight relative to σ² grows as σ falls, so σ = 0.2 is the more discriminating of
  the two and σ = 0.35 is the envelope boundary. Reference cells for reporting: Tier A at the same σ and
  at the largest tabulated donor count not exceeding `min(n_A, n_B)`, i.e. `n` = 6 for M6a and `n` = 8
  for M6b. M6 is under **V11**, not V7 — M6a's design is deep-outside at σ = 0.35 and deep-inside at
  σ = 0.2, M6b's is borderline-inside at 0.35 and deep-inside at 0.2, so the two designs do not share a
  single verdict-error direction to build a ratio on. Its wrong-admission and over-exclusion rates are
  reported per cell all the same.

**32 cells × 100 replicates = 3200 simulations.** Expected directions are declared in advance, so a
surprise is recognisable as one: M4 — E1 underestimates `v` and therefore **overestimates** σ (safe),
E2 approximately unbiased; M3 — σ̂ up (safe), with the LODO influence diagnostic lighting up; M2 — E2
carries, E1 degrades; M5 — both weightings agree, because the difference between them vanishes in
expectation when cells per donor do not differ by group (1.6); M6 — they do **not** agree, and the df
weighting is the one that matches `wls_two_group`.

**Total: 180 cells, 19 600 simulations.** The budget is deliberately not stated as a promise: 164 of the
180 cells cost about a second each, but the 16 compression cells carry fifteen times the genes of every
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

**V2 — bias, Tier B.** Tier-B cells at 10 and 30 cells per donor: `|bias| ≤ max(0.05, 0.15·σ)`.
**Tier-B cells at 100 cells per donor and above take V1's band instead**, `|bias| ≤ max(0.03, 0.10·σ)`.
That second clause closes a gap in the earlier draft of this entry, which scoped V1 to Tier A and V2 to
10 and 30 cells and thereby left **56 of the 88 Tier-B cells under no binding bias criterion at all** —
including the entire 16-cell compression block, which runs at 100 and 300. Those cells are no harder
than Tier A's own 300-cell geometry and there was no reason for them to be looser; V3 (RMSE) does not
bind and could not have covered them. The 12 compression cells with σ > 0 additionally take **V10**'s
much tighter band.
*On failure at 10 cells only*: the pre-declared consequence fires automatically — strata whose median
cells per donor is below 30 receive `indeterminate` (Change 2, reason 2). The threshold is not relaxed;
the estimator's domain is narrowed, and the narrowing is written into the gating rule rather than
invented after the fact. *On failure at 30 cells*: as V1 — the estimator fails, because `[30, 100)` is
the frozen list's modal cells-per-donor bin (161 of its 502 group medians). *On failure at 100 cells per
donor or above*: as V1, with no domain narrowing available — those are the **easy** end of the tier, no
harder than the geometry V1 itself binds on, and a bias failure there is a defect rather than a limit.

**V3 — RMSE.** Reported for every cell; **not binding**. Its decision-relevant projection is V5/V6, and
a second binding threshold on the same quantity would be double counting.

**V4 — UCB coverage.** Every Tier-A cell with `n ≥ 4` **and σ > 0**: `P(σ_UCB ≥ σ) ≥ 0.85` against the
nominal 0.90; pooled over those cells, `≥ 0.88`. The Monte-Carlo SE of a coverage estimate at 0.90 over
200 replicates is `sqrt(0.9 × 0.1/200) = 0.021`, so 0.85 sits a little over two SE below nominal — tight
enough to catch a broken interval, loose enough not to fail on Monte-Carlo error. `n = 3` is **reported
and not binding**: the envelope's easiest row already demands 4 donors per group, so a 3 v 3 stratum is
never admitted whatever its interval does. **σ = 0 is excluded because it is vacuous, not because it is
awkward**: `σ_gate = sqrt(max(UCB, 0)) ≥ 0 = σ` identically, so those cells report coverage 1 whatever
the jackknife does. That is **45 binding cells, not 50** — the exclusion removes five guaranteed 1.0s
from a pooled denominator and therefore *tightens* the pooled criterion rather than relaxing it. The
per-cell threshold's exposure to multiplicity over 45 cells is disclosed under *What this does NOT
settle* and no multiplicity adjustment is applied to it.
*On failure*: the jackknife is not delivering the coverage the gate assumes, the gate is not usable as
specified, and the alternative named in Change 4 item 2 goes on the table as a future amendment with its
own validation.

**V5 — false admission. The main criterion.** On deep-outside cells — **the 19 cells marked `DO` in
Change 5's classification table**, whose true σ is at or above the next tabulated step past what the
cell's `n` supports — `P(verdict = inside) ≤ 0.05` in **every** such cell and `≤ 0.02` pooled over the
19 of them (3800 replicates). The ten `n` = 3 cells are **excluded from both the per-cell requirement
and the pooled denominator** as structurally vacuous, per that table's third note: eight of them would
otherwise be deep-outside and would contribute 30 % of the pooled weight as guaranteed zeros. σ = 0
never enters V5 at any donor count, for the reason given under reason 4.
*On failure*: the gate admits strata the arm is invalid for, which is spec §10 risk 1. The estimator does
not gate admission. No threshold moves.

**V6 — the price of over-exclusion.** On deep-inside cells — **the 15 cells marked `DI` in Change 5's
classification table**, where `n` suffices even for the next tabulated σ upward — `P(inside) ≥ 0.80` in
every such cell and `≥ 0.90` pooled over the 15 of them (3000 replicates). Borderline cells are
**reported and not binding**: an upper-confidence-bound gate is *supposed* to cut them, and binding there
would penalise the construction for working. `n` = 3 is vacuous and excluded, as under V5. The three
σ = 0 deep-inside cells (`n` = 8, 13, 24) **are** binding, and they are the strictest thing the grid asks
of `v̂` — see reason 4, where the exemption that makes them reachable at all is derived.
*On failure*: the instrument is honest but unusable, excluding nearly everything. That is reported as
such, and the sweep-feasibility question (freeze §6) is answered in the negative for reasons of
instrument conservatism rather than of real σ. It is not repaired by loosening the UCB.

**V7 — misspecification.** Binds on Tier-C arms **M2, M3 and M4**. In each of their cells,

```
P(verdict = inside | true verdict = outside)   ≤   max( 2 × p_ref ,  0.05 )
```

where `p_ref` is the same rate at the arm's declared **corresponding well-specified Tier-A cell**
(Change 5 tabulates one per cell; they are (0.35, 4) and (0.35, 8) throughout).

Three things about that line differ from the earlier draft of this entry, and each is a repair rather
than a relaxation.

**(i) The arms are named, and their coordinates exist.** The draft said "in each Tier-C arm", but only
M1 had declared (σ, `n`); M2 … M5 were counts. "The corresponding well-specified Tier-A cell" was
therefore undefined for four of the five arms, and V7 was not computable as written. Change 5 now
declares every Tier-C coordinate.

**(ii) `M1` is removed from V7, because its ground truth is not the same quantity as Tier A's.** In
Tiers A and B a cell's truth comes from applying step-up to the **true** σ. In M1 the truth is the
**measured** sensitivity of the `ebayes` arm on the pre-registered oracle, `inside` iff that sensitivity
is ≥ 0.60. Those two truths do not agree even on well-specified data, and this log has already measured
the gap: Amendment 3 records the same arm at the same nominal `(σ = 0.5, 8 v 8)` measuring power
**0.194** on the grid's `directnb` generative arm and **0.4006** on its `lognormal` arm — the arm
`oracles.py` actually is — while the envelope's `min_donors_per_group` comes from Amendment 1's analytic
frontier and from neither. A ratio between a rate computed against one truth and a rate computed
against the other is not a ratio of like quantities, and `2 ×` of it bounds nothing. **M1 is instead
judged against its own truth, by V5's own threshold**: `P(verdict = inside | measured sensitivity <
0.60) ≤ 0.05` in every M1 cell and `≤ 0.02` pooled. That is exactly the screen Change 3's step 1 already
applies to it, made a criterion in its own right. It makes M1 **binding** rather than unbinding it, and
it removes an incommensurable comparison instead of papering one over.

**(iii) The ratio has an absolute floor, and the floor is not a new number.** `2 × p_ref` is **zero
whenever the reference cell reads 0 of 200** — which a working estimator is expected to do at most
deep-outside cells — and against a bound of zero a single false admission in 100 Tier-C replicates
fails the criterion. At a true rate of 0.002 that happens with probability `1 − 0.998¹⁰⁰` = **0.181**: a
criterion that fails one run in five on a *correct* estimator is measuring Monte-Carlo error, not
misspecification. The floor is **V5's own per-cell deep-outside threshold, 0.05**, reused rather than
invented, so an arm that meets the absolute standard the well-specified grid is held to cannot fail V7
merely because its reference cell was clean.

**M5 and M6 are not under V7 either**, and Change 5 says why for each: M5's 8 v 3 cell has
`min(n) = 3`, where `inside` is structurally unreachable and a wrong-admission criterion passes
vacuously; M6's four cells do not share one verdict-error direction. Both arms are under **V11**, whose
criterion is bias — a quantity that is defined and non-vacuous at every one of their cells.

*On failure of any V7 arm*: **the threshold does not move, and the failing arm's pre-declared signature
becomes a blocking `indeterminate` flag** (Change 2, reason 5) — LODO influence for M3, E1/E2 divergence
for M4, cells-per-donor coefficient of variation for M2 — at the threshold declared for it in reason 5's
table, which is a literal or a quantile of this run's own committed artifact and is in no case a number
Part B chooses. The estimator's domain shrinks to the strata that do not show the signature, and the
shrinkage is reported as a result. *On failure of M1's criterion*: as V5 — the estimator does not gate
admission, and the realised `d0` signature blocks at reason 5's declared quantile. **This is one of the
two failures whose consequence is a new blocking flag rather than a moved threshold.**

**V8 — double use of the data.** Among Tier-A replicates whose verdict is `inside`, the `ebayes` arm's
fresh-null false-positive rate — the fraction of replicates with at least one BH rejection — must be
`≤ α + 2 · MC SE`, with `MC SE = sqrt(α(1 − α)/n_inside)` computed from the realised count of `inside`
replicates. A stratum selected for having drawn a *quiet* realisation of its donors must not thereby
become anti-conservative conditional on that selection.

**A companion figure is reported beside it and is explicitly not binding**: the same false-positive rate
over the `inside` replicates of the **borderline** cells alone — the 16 cells where admission was
actually in doubt and the conditioning V8 is named for can exist at all. The pooled form stays binding
and its threshold does not move; the reason the companion cannot also be binding is that its denominator
is the realised count of `inside` replicates in cells whose admission rate is unknown in advance, and a
threshold sized on an unknown denominator is a threshold invented after the run. The dilution of the
pooled form is recorded under *What this does NOT settle*, with its arithmetic. If the two figures
disagree materially, that is a result and it comes back to this log.

*On failure*: the selection is not benign and admission must be conditioned differently. **This is the
second failure whose consequence is a new blocking flag rather than a moved threshold**: the estimator
does not gate admission at all until a selection-aware rule is pre-registered in a further amendment.

**V9 — the correction is not a no-op, and Amendment 3's bound is demonstrated to reverse.** Two clauses,
both binding.

* **V9a.** At the Tier-B cell (1000 genes, 30 cells per donor, φ = 0.2, σ = 0.35, `n = 8`): the bias of
  `sqrt(s0²)·ln 2` must be at least **3 ×** the absolute bias of `σ̂`. Over development seeds 1–16 that
  cell puts `sqrt(s0²)·ln 2` at 0.3886 (sd 0.0031; `pilot/upper_bound_check/`), a bias of +0.0386, so the clause demands
  `|bias σ̂| ≲ 0.013`. **That is roughly four times tighter than V2 permits at the same cell, and the
  asymmetry is intentional** — it is stated here so a V9a failure cannot later be re-read as a V2 pass.
* **V9b.** At the Tier-B compression cell (`n_genes = 15 000`, 300 cells per donor, φ = 0.2, σ = 0.35,
  `n = 8`): the **mean over replicates** of `sqrt(s0²)·ln 2` must be `< σ`. This is the obligatory
  demonstration that Amendment 3's "upper bound" stops being an upper bound in the regime real universes
  occupy. The full 50-replicate distribution is reported, not merely its mean. **Its expected outcome is
  disclosed rather than dramatised**: the committed development artifact already puts this cell at
  0.3319 with a replicate sd of 0.00033 and **16 of 16** seeds below σ, so V9b is expected to pass. It is pre-registered on
  disjoint confirmatory seeds all the same, because its job is to pin a reproducible demonstration into
  a committed artifact rather than to manufacture suspense — and because a criterion written after its
  own run would be worth nothing regardless of which way it went.

*On failure of V9a*: the correction is not earning its complexity at the thin end and the estimator's
advantage over the mechanism is not demonstrated there — reported, with the estimator's fate decided by
V1/V2/V5 as usual. *On failure of V9b*: Correction 1's empirical demonstration has failed at the cell
chosen for it. Correction 1's algebra stands on the derivation and on the counterexamples in its own
table regardless — but the failure is recorded, the cell's inadequacy is diagnosed, and the demonstration
is **not** quietly moved to a cell where it works.

**V10 — the de-attenuation divisor. NEW, because without it this entry's central novelty was
falsifiable by nothing.** The estimator's one structural departure from Amendment 3's mechanism is that
it **divides by `ā²_g`** (Change 1.6). Every criterion above can be passed by an estimator that never
does. Measured on development seeds 1–8 at 8 v 8, σ = 0.35, φ = 0.2, with the E2 bootstrap `v̂`
(± is the sd over the 8 seeds):

| geometry | `rms_g(ā_g)`, estimation stratum | `rms_g(ā_g)`, all universe genes | σ̂ **with** the divisor | σ̂ **without** | rel. bias, with | rel. bias, without |
|---|---|---|---|---|---|---|
| 1000 genes, 300 cells | 0.99572 ± 0.00024 | 0.99537 ± 0.00030 | 0.350434 ± 0.002674 | 0.348924 ± 0.002685 | **+0.124 %** | **−0.308 %** |
| 15 000 genes, 300 cells | 0.97865 ± 0.00013 | 0.94484 ± 0.00056 | 0.350201 ± 0.000546 | 0.342678 ± 0.000545 | **+0.058 %** | **−2.092 %** |
| 15 000 genes, 100 cells | 0.97860 ± 0.00015 | 0.94445 ± 0.00057 | 0.350687 ± 0.000669 | 0.343124 ± 0.000658 | **+0.196 %** | **−1.965 %** |

The ratio `σ̂_with / σ̂_without` equals `1/rms_g(ā_g)` over the estimation stratum to within 2.4e-4
absolute across all 24 seed × geometry runs, so the divisor's whole effect is that one number and it is
predictable from the design. **At the reference geometry it is worth +0.43 % in σ — a sixth of V1's
tolerance at σ = 0.35, and invisible to every criterion above. At 15 000 genes it is worth +2.20 %, and
that is measurable**: an estimator omitting it reads 2.09 % low with a Monte-Carlo SE of 0.055 %, a
38-sigma bias.

**V10 binds on the 12 compression cells with σ > 0** (`n_genes` = 15 000, cells ∈ {100, 300},
`n` ∈ {4, 8}, σ ∈ {0.2, 0.35, 0.5}):

```
| mean_reps σ̂_rms  −  σ |   ≤   0.010 · σ
```

**The band is half the measured magnitude of the effect it must detect**, which is how it is sized and
the only reason it is 0.010 rather than a round guess. An estimator omitting the divisor sits at about
−2.1 %, **more than twice the band outside it**, and fails. A correct one sits at +0.06 % to +0.20 % on
the cells measured, five to sixteen times inside it, against a Monte-Carlo SE over 50 replicates of
roughly 0.02 % of σ. The four σ = 0 compression cells are outside V10 because relative bias is
undefined at σ = 0 and the divisor has nothing there to act on; they are covered by V2's second clause.

**Why V10 cannot be satisfied by tightening the filter instead of dividing.** The `median CPM ≥ 20`
filter already does the larger share of the de-attenuation: over *all* universe genes at 15 000 the
attenuation is `rms(ā)` = 0.94484, so an unfiltered, undivided estimator reads 5.84 % low, and the
filter absorbs that to 2.20 %. So the criterion has to be shown to be about the divisor and not about
the filter. It is, because **`ESTIMATION_MIN_MEDIAN_CPM` = 20.0 is fixed by Change 2's constant table
and moving it is an amendment.** The filter is a constant of this pre-registration; the divisor is the
only free thing left that can close the remaining 2.20 %. Both attenuation figures — filtered and
unfiltered — are written per cell into the committed summary, so the two contributions stay separable
to a reader.

**A mandatory ablation is reported at all 12 cells**: the identical estimator with `ā²_g ≡ 1`, on the
identical replicates. *If the ablation is not outside V10's band at some cell*, that cell has no power
to falsify the divisor and is reported as such — a statement about the cell, **not** an estimator
failure, and it does not move V10.

*On failure*: as V1 — the estimator is not fit to gate admission at the universe width real strata
actually have, and the failing region is reported. A V10 failure in the low direction is specifically
the finding that **the de-attenuation is not working**, and it may not be re-read as "within V1's
tolerance", exactly as V9a may not be re-read as a V2 pass.

**V11 — the degrees-of-freedom weighting, on unbalanced designs. NEW.** Change 1.6 pins the pooling
weights to `w_d = (1 − 1/n_group(d))/(n_A + n_B − 2)` and states that `1/n` weights are wrong. Nothing
in the grid, as the earlier draft of this entry defined it, could tell the two apart: their difference
is `Σ_d (w_d − 1/n_d)·v_Y[d,g]`, which vanishes in expectation whenever `E[v_Y]` is the same in both
groups — and every Tier-A cell, every Tier-B cell, M2, M3, M4 and both M5 cells have equal cells per
donor and therefore satisfy that. Tier C arm **M6** exists to break it, and V11 is what binds on it.

**V11a — bias.** On M5's 2 cells and M6's 4: `|mean σ̂_rms − σ| ≤ max(0.03, 0.10·σ)` — V1's band,
applied to the unbalanced designs V1 itself never reaches.

**V11b — the weighting ablation, and it must be read *paired*.** On the same replicates, re-aggregate
with `w_d = 1/(n_A + n_B)` in place of the df weights and report `Δ = σ̂_df − σ̂_1/n` **per replicate**,
then its mean and its **paired** standard error. `Δ` must be non-zero at `|Δ| ≥ 5` paired SE on each M6
cell, and **must carry opposite signs on M6a and M6b** — positive where the small group is the shallow
one, negative where it is the deep one — so that a sign error in the weighting cannot pass both.

The pairing is not a stylistic preference, and this is measured rather than asserted. On development
seeds 1–8 at σ = 0.35, 1000 genes, φ = 0.2, holding the cell-count design fixed across seeds:

| design | realised median cells/donor, A / B | group mean `v̂_Y` ratio B/A | `Δ = σ̂_df − σ̂_1/n` | paired SE | **paired SE apart** | unpaired SE apart |
|---|---|---|---|---|---|---|
| M6a — 20 v 7 | 336 / 29 | 19.4 | **+1.925e-3** | 4.04e-5 | **47.7** | 0.97 |
| M6b — 100 v 10 | 33 / 378 | 0.047 | **−5.326e-4** | 1.02e-5 | **52.1** | 0.41 |
| control — 20 v 7, 300 cells everywhere | 300 / 300 | 1.00 | +4.3e-7 | 2.6e-7 | 1.6 | 0.00 |

Read **unpaired** — two means over replicates compared with each other — the same eight seeds give 0.97
and 0.41 SE: no power whatsoever, and of the order of 330 and 2000 replicates would be needed to reach
3 SE. Read paired they give 48 and 52. **A criterion written as a comparison of two means would have
been unfalsifiable at the grid's 100 replicates while looking like a test.** The grid draws each
donor's cell count once per replicate from M6's stated `LogUniform`, as M2 does, which adds design
variance to the pairing; measured under that construction the separations are 6.0 and 37.8 paired SE at
8 development replicates, scaling to roughly 21 and 134 at the grid's 100. The 5-SE bar is set well
below both.

The control row is the whole point of M6's design: **group imbalance alone is not enough.** At 20 v 7
with equal cells per donor the two weightings agree to 4e-7 in σ, which is why M5 cannot do this job
and why the earlier draft's grid could not falsify the weighting anywhere.

**What V11b may NOT be read as.** It is a criterion on the **sign and size of the paired difference**,
never on which weighting lands nearer the truth. On M6b the measured df estimate is 0.95 % below σ and
the `1/n` estimate 0.80 % below — the **rejected** scheme is nearer the truth on that cell, because
M6b's residual bias comes from its 100 shallow donors (median 33 cells, technical share 0.158) and not
from the weighting at all. A criterion phrased as "the df weighting is the more accurate one" would have
rejected the correct weighting on the cell built to confirm it. That reading is foreclosed here, before
the run, so it cannot be adopted after it.

*On failure of V11a at any M6 cell*: the estimator is not validated on group-asymmetric unbalanced
designs, and the group-imbalance ratio becomes a blocking signature at `IMBALANCE_RATIO_BLOCK` = 0.35
(Change 2, reason 5). **That would exclude 38 of the 251 frozen strata across four datasets, and no
dataset entirely** — counted from the freeze here rather than discovered in Part B. **What M6 does not
cover is also stated**: the freeze's two most extreme designs are #5's 54 v 3 (ratio 0.056) and #4's
3 v 28 (0.107), both more skewed than M6b, and both have `min(n_A, n_B) = 3`, so the envelope's own
4-donor floor already refuses them whatever the weighting does. M6's range therefore spans the skew
that can actually reach a verdict. *On failure of V11b*, i.e. the paired difference not
separating: M6 has no power to falsify the weighting choice, that is reported, and the df weighting
stands on Change 1.6's algebra and on `wls_two_group`'s hat matrix alone rather than on the grid.
**Neither failure moves a threshold.**

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
* **The envelope's σ = 0.2 row itself is unvalidated, and it is the row nearly everything will land
  on.** Its `min_donors_per_group` = 4 is Amendment 1's analytic frontier with `grid_support` recorded
  as *"not in the grid"*, and freeze §6 puts 11 of 12 datasets and 227 of 251 strata in that tier.
  Step-up sends every stratum with `σ_gate ≤ 0.2` there. Nothing in this amendment measures power at
  (σ = 0.2, n = 4); the membership rule is exactly as sound as that unmeasured row. Change 4 item 3
  states this at the point where the conservatism argument is made, and it is repeated here because it
  is the largest thing a validated estimator would still leave open.
* **V4's per-cell threshold is multiplicity-exposed, and no correction is applied to it.** At *exactly*
  nominal coverage, a single cell reads below 0.85 with probability `P(Bin(200, 0.9) ≤ 169)` =
  **0.00951**, so the probability that **at least one** of V4's 45 binding cells does is
  `1 − (1 − 0.00951)⁴⁵` = **0.349**. (Over the 50 cells V4 covered before the σ = 0 vacuity exclusion it
  is 0.380.) A lone per-cell failure is therefore close to a coin flip under a perfectly calibrated
  interval. The threshold is **not** moved and no adjustment is applied — a multiplicity correction
  would be a loosening dressed as rigour, and the pooled `≥ 0.88` leg, whose own exceedance probability
  at nominal coverage over 9000 replicates is below 1e-9, is not exposed at all. What is recorded here
  is the number itself, so that a single failing cell can neither be read as proof of a broken interval
  nor waved away later without this arithmetic on the table.
* **V8's pooled average is dominated by the replicates that least need checking, and it barely
  constrains the selection effect it is named for.** V8 asks whether selecting a stratum for having
  drawn a quiet realisation of its donors makes the `ebayes` arm anti-conservative *conditional on that
  selection*. That conditioning can only exist near the admission boundary: a replicate that would be
  admitted whatever it drew was not selected on anything. But V8 pools all Tier-A `inside` replicates,
  and by Change 5's classification 15 of the 60 cells are deep-inside, where V6's own demand is
  `P(inside) ≥ 0.80` and the realised rate is meant to be near 1 — including six cells at σ ∈ {0, 0.1}
  and `n ≥ 8`, where the arm's fresh null is the easiest case it ever faces. Only the 16 borderline
  cells carry real conditioning, and they are a minority of the pooled denominator by an amount that
  will not be known until the run. So a pooled V8 pass is weak evidence that the selection is benign,
  and V8 should be read as a **screen against gross anti-conservativity**, which is what its threshold
  is sized for, rather than as a measurement of the selection effect. The borderline-only companion
  figure V8 now reports is the closer look; it is not binding, for the reason given there.
* **The realistic hard corner — a wide universe *and* shallow donors — is not in the grid at all.** The
  compression block widens the universe to 15 000 genes at the simulator's default gene means, and the
  arithmetic of that is `E[Σ_g μ_g] = G · e^{mean_log_mu + mean_log_sigma²/2} = 15000 · e^{1.72}` =
  **83 768 counts per cell** — *above* the frozen list's own maximum of 56 841.5 median counts per cell,
  and 295 times its minimum of 284. (At the reference geometry's 1000 genes the same arithmetic gives
  5585, comfortably inside the freeze's range.) So the block isolates the attenuation `ā` on donors
  deeper than any real one: exactly the right place to look for the `+1` compression, per Correction 1,
  and exactly the wrong place to look for its **interaction** with technical noise. A real stratum with
  a 15 000-gene universe and 284 counts per cell has `ā` ≈ 0.97 **and** a large `1/T` simultaneously,
  and no cell in this grid has both. Reaching it needs `mean_log_mu` about 4.5 natural-log units below
  the simulator's default at `n_genes` = 15 000, which no cell in this grid does. It is named rather
  than added: at 15 000 genes these are the grid's most expensive cells by a wide margin, and a
  two-cell probe would not validate the estimator there — it would only tell us whether to build a grid
  for it. **The honest statement is that the estimator is not validated in the regime a real
  wide-universe, shallow-donor stratum occupies**, and that is the first candidate for the next grid.
* **Whether the estimation stratum is even reachable on real strata is unmeasured, and the count gate
  binds harder than the CPM gate.** The two filters interact: `T = CPM · L/10⁶`, so
  `ESTIMATION_MIN_MEDIAN_COUNT` = 10 is the condition `CPM ≥ 10⁷/L`, and the **effective CPM floor** is
  `max(20, 10⁷/L)` — the count gate is the binding one whenever the donor library `L` is below 5 × 10⁵.
  From the freeze's own committed marginals, estimating each group's median donor library as
  (median cells per donor) × (median counts per cell) and taking the smaller of a stratum's two groups:
  **154 of the 251 strata, and 286 of the 502 group-arms, sit below 5 × 10⁵** and are therefore
  count-gated rather than CPM-gated. The effective floor has median **42.1**, 75th percentile **128**,
  90th percentile **347** and maximum **1760** CPM, and exceeds 100 on **78 of the 251**. Every one of
  the 12 datasets holds at least one stratum in the count-gated regime, and one (#7, Binvignat 2024
  rheumatoid-arthritis blood) has *all* 15 of its strata there.

  Whether ≥ `MIN_ESTIMATION_GENES` = 100 genes clear that floor cannot be computed from the freeze,
  because the per-stratum expression distribution and `frozen_universe_size` are `pending` on all 251
  rows. As a **scenario** — the simulator's own log-normal gene-mean model, under which
  `ln CPM_g ~ N(ln(10⁶/G) − s²/2, s²)` with `s = mean_log_sigma` = 1.2, independent of `mean_log_mu`,
  verified numerically before being used — the floor at which fewer than 100 genes survive is
  `f*(G) = 1145 / 794 / 632 / 535` CPM at `G` = 5000 / 10 000 / 15 000 / 20 000, and the number of the
  251 strata whose effective floor exceeds it is **1 / 3 / 6 / 8**. No dataset loses all its strata at
  any tested `G`; the failures concentrate in #7 and in #3 (Melms 2021, COVID-19 lung).

  Three caveats bound that, all in the same direction. The library estimate is a product of two
  marginals rather than a measured library; it uses whole-transcriptome counts per cell, so it
  **over**-states the universe-restricted `L` and therefore **under**-states the floor — the count gate
  binds at least as hard as this says, and the failure counts are a lower bound. And the gene-mean model
  is our simulator's, not any real stratum's, whose low-expression tail is heavier. What this arithmetic
  establishes is not a number to plan on but that the estimation stratum is **not obviously reachable on
  a sixth of the frozen list**, that reason 1 (`indeterminate`, fewer than 100 genes) is a live exclusion
  path rather than a formality, and that the exclusions it causes will be concentrated in the two
  shallowest datasets rather than spread evenly. It is arithmetic Part A owed and had omitted.
* **Which functional gates.** That is Part B's, decided by Change 3's rule and by nothing else.
* **Whether the estimator passes at all.** Part A claims no result. If V1, V2, V4, V5, V6, V7, V10 or
  V11a fails, admission stays closed and the failure is reported here.
* **A2 stratification remains deferred** (Amendment 2 Change 6); the naive arm's floor is still
  cell-count confounded and guarded only by a range check.
* **The GO/NO-GO decision is not taken.** An estimator that gates admission licenses a measurement, not a
  conclusion.

*Author attests: every figure in this entry is derived, measured on named seeds, or quoted from a
committed artifact. The ten probes under Correction 1 were run against this repository's own arm code
on disclosed development seeds 1–16 by `scripts/check_upper_bound_claim.py`, and the artifact they
produced is committed, so every aggregate in that table is recomputable rather than asserted; where a
figure moved against the draft's table, the move is reconciled beside it. The sizing measurements for
V10 and V11 and the corroboration of reason 4's exemption were made on development seeds 1–8 with an
uncommitted prototype, disclosed under "Data visible" item 7 and labelled as sizing rather than
validation. The algebra of Correction 1 and of Change 1 was checked term by term against
`synthetic/oracles.py` and `src/pbcheck/methods/moderated.py` rather than quoted from an earlier entry —
which is how the missing depth term in E1's moment equation was found, and it is corrected in Change 2
rather than carried. The frozen list's distributional figures were recomputed from
`pilot/preregistration/stratum_list_2026-08-16.json` and agree with
`docs/PREREGISTRATION_STRATUM_LIST.md` as committed on every figure quoted. No real data informed this
amendment. The confirmatory grid has not been run and no criterion V1–V11 has an outcome. The
`L ≷ 10⁶/(2σ²)` threshold that orders the probe table is labelled a heuristic, its two cancelling
approximations are named, and it is read by no code. The correction to Amendment 3 is stated as the
reversal of a published claim of this log, not as a clarification of it. Where this entry differs from
the draft that an adversarial read broke, the differences are listed at the top rather than folded in
silently, and no threshold among them was loosened to make anything pass.*

---

## Amendment 5, Part A (2026-08-16) — the §6 reportability split: naive-arm null quantities become reportable at min(n_A, n_B) ≥ 8 on measured per-stratum calibration; Amendment 3 Change 1's blanket reporting prohibition is NARROWED; the powered claim stays inside the envelope

Amendment 3 Change 1 declared an operating envelope and closed it with a blanket prohibition:

> It is **not** declared valid outside it, and no result from a stratum outside it may be reported
> as a pbcheck measurement.

That sentence is one day old and it is wider than the argument that produced it. The envelope is a
**power** statement about the *pseudobulk* arm, derived from a power frontier in
(`sigma_donor`, donors-per-group). Two of §6's quantities — `lambda_naive` and the naive arm's
permutation false-positive floor — are properties of the *naive* arm under the donor-permutation
null and have no pseudobulk hit in their numerator or denominator. They are gated today by a
sentence whose reasoning does not reach them.

**This entry narrows that prohibition. Narrowing a prohibition is a relaxation, and it is called
one here rather than described as a clarification** — the substitution that made Amendment 1
necessary. What is relaxed is exactly one sentence of one amendment. No threshold moves.

Spec sections touched: **§6** — the inflation metrics acquire a declared reportability partition,
and the cross-dataset matching target that B2 has always demanded is supplied as a number.
**§4 / A2** —
cell-count stratification, deferred by Amendment 2 Change 6, is **un-deferred at ≥ 8 v 8** and
becomes a blocking condition of the reporting this entry permits. **§1 (inclusion gate)** — a second
admission rule, for the map only, sitting beside (not replacing) Amendment 4 Change 4's envelope
membership. **Decision rule items 3 and 4** — item 3's floor becomes reportable outside the
envelope under the conditions below; item 4's `signal_above_floor` is demoted out of
decision-relevant use outside the envelope.

**Unchanged, and not weakened by anything below:** α = 0.05; the λ band [0.9, 1.1]; `POWER_TARGET`
= 0.60; the oracle's log2FC = 1.0 and K = 200; Amendment 3's envelope table; §1's floor of 8
independent datasets; and **decision rule item 1 in its entirety**, including its
VOID → NO-GO on a broken denominator.

### Why this entry is written BEFORE the sigma anchor, and what binds when

Amendment 4 Part A stated the principle this entry inherits:

> A validation whose PASS criteria are written after its numbers are read is not a validation; it
> is a description.

The same holds one level up, for scope of publication rather than for criteria. A rule about *what
may be reported*, written after the anchor lands, is not a pre-registration — it is a description
of the set that survived. `sigma_donor` is measured nowhere: `sigma_donor_estimate` is `pending` on
all 251 frozen rows and on all 2190 rows of the source manifest, `admitted_to_sweep` is `False` on
every frozen row, and `pilot/results/` is empty. So the rule can still be written in the only order
that makes it binding.

**Three specific things make the order load-bearing, not ceremonial.**

1. **The numbers a later rule would be chosen from are already public.** The freeze's §6 tier table
   — 227 / 150 / 94 / 30 strata and 11 / 7 / 5 / 3 datasets at the envelope's four rows — was
   committed earlier today and is disclosed below as visible. After the anchor, *any* donor-count
   rule is a selection from a table whose answer is known, and the freeze forbids exactly that:
   "Nothing above may be selected from", and shrinkage is "a reported outcome, never a
   re-selection" (freeze §6, §9 item 7).
2. **It is not yet too late in the mechanical sense.** Nothing has been admitted, no stratum has
   been loaded, and no metric has been computed on any real row. A rule fixed now costs nothing to
   fix and cannot be fitted to a result that does not exist.
3. **Amendment 1 exists because a gate was quietly made easier.** An entry that widens what may be
   published is the exact shape of the failure this log was started to catch. It is therefore
   obliged to be loud: to name the relaxation as a relaxation, to enumerate what it does *not*
   move, and to pre-commit the observables that would make it the wrong decision (below), so that
   abandoning it later is a rule being applied rather than a judgement being made.

**What Part B is.** Part B is a dated addendum written after the anchor: which strata actually
cleared the rules below, what the σ estimator returned, what landed in Tier 3, and — if any
abandonment trigger fired — the withdrawal. Part B decides nothing that is left open here.

### Data visible at the time of this amendment (full disclosure)

1. **Everything visible for Amendments 1, 2, 3 and 4**, unchanged.
2. **The committed synthetic gate**, `pilot/gate/synthetic_gate_2026-08-15.json`, read in full while
   writing: `lambda_naive` 54.5718, `lambda_pseudobulk` 1.00656, pseudobulk perm-null FP rate 0.035
   (200 paired permutations), naive permutation floor median 1161.5 of G = 1500 (77.43 %, IQR
   22.25, Monte-Carlo SE 1.106), pseudobulk floor median 0 (max 1, MC SE 0.013), real-label #DEG
   naive 1130 against pseudobulk 1, `power_sensitivity` 0.86 at empirical FDR 0.0444, verdict
   `INSTRUMENT VALID WITHIN THE STATED OPERATING ENVELOPE`.
3. **The 146-cell test-selection grid**, `pilot/testsel/summary.json` at `72dec7b`. Re-read for this
   entry, and one property of it is load-bearing below: **the only `sigma_donor` values in the grid
   are 0.35, 0.5 and 0.7. There is no σ = 0.2 cell at all.**
4. **The frozen stratum list**, `docs/PREREGISTRATION_STRATUM_LIST.md` and
   `pilot/preregistration/stratum_list_2026-08-16.json`: 251 strata over 12 datasets, all four
   admission blockers standing on every row. **Its §6 tier table was read before the rule below was
   chosen, and is reproduced here so that the disclosure cannot be softened later** — re-derived
   from the JSON rather than copied: `min(n_A, n_B)` ≥ 4 / 8 / 13 / 23 gives **227 / 150 / 94 / 30**
   strata over **11 / 7 / 5 / 3** datasets, matching the freeze element for element. Also read:
   the D1 bins, the balance of the list (42 of 251 strata have `n_A = n_B`; 21 of the 150 at
   ≥ 8 v 8 do), the per-tier depth (median of the group medians 103.5 / 118.0 / 206.75 / 200.0 /
   310.0 at ≥ 3 v 3 / 4 v 4 / 8 v 8 / 13 v 13 / 23 v 23) and the A2 confound figures tabulated in
   Change 6.
5. **The pinned candidate manifest**,
   `pilot/preregistration/census_candidates_run31910799023_2026-08-15.json`: 1197 candidate rows
   over 68 datasets in 50 collections, with the donor-count ladder re-derived here —
   ≥ 4 v 4: **1017 strata / 62 datasets / 46 collections**; ≥ 7 v 7: **630 / 38 / 30**;
   ≥ 8 v 8: **554 / 33 / 25**; ≥ 13 v 13: **311 / 21 / 15**; ≥ 23 v 23: **134 / 12 / 10**. The
   dataset and collection columns agree with the freeze's own `MANIFEST_TIER_CENSUS` and
   `MANIFEST_TIER_COLLECTIONS` at every tier it declares.
6. **`sigma_donor` is measured nowhere**, `pooled` is `unresolved` on 251 of 251 and 1197 of 1197,
   and the manifest header's `pool_columns_detected` is `[]`. `integer_check` and
   `frozen_universe_size` are `pending` on all 251. **Every stratum count in this entry is therefore
   an upper bound that a real run can only reduce.**
7. **What was measured for this entry, and is new here.** (a) The leak coefficients of Change 4,
   computed exactly from the hypergeometric law over the permutation set
   `pbcheck.permutation.build_perms` actually produces. (b) The exact-binomial thresholds and
   false-exclusion arithmetic of Change 3. (c) The A2 confound distributions of Change 6, over the
   frozen 150 and the manifest's 554. (d) A fresh-null probe of the moderated arm at
   `sigma_donor` = 0.2 on the gate's own geometry, on the declared seed range 1000–1999, reported
   under "What this does NOT settle" — it is a check that **failed to confirm** a claim that was in
   circulation while this entry was drafted, and the failure is recorded rather than the claim.
8. **No real data.** Oracle (d), Mathys 2019 (§8(d)), remains binding and unrun. No CELLxGENE
   stratum has been loaded. Every synthetic number here comes from our own generative model.

### What an adversarial read of this entry changed before it was committed

A draft of this entry rested on four supports that did not survive being attacked. They are listed
rather than folded in, on the precedent of Amendment 4: a pre-registration that quietly improved
between drafts is indistinguishable from one written to fit.

1. **The claim that 8 v 8 is where the permutation null becomes resolvable is FALSE, and is
   deleted.** The draft argued that the admission threshold coincides with the point at which the
   permutation set reaches spec §4's `n_perm` = 1000. It does not: C(12, 6) = 924 < 1000 but
   **C(14, 7) = 3432 > 1000**, so resolution closes at **7 v 7**, one tier below. The manifest
   confirms it — 630 candidate strata over 38 datasets have `min(n_A, n_B)` ≥ 7, and the smallest
   `permutation_count` among them is exactly 3432. Two independent-sounding reasons for the same
   threshold were in fact one reason and one coincidence-that-isn't. **8 v 8 rests on correction A1
   alone**, and Change 2 says so.
2. **All three external-precedent citations are deleted, and none is cited anywhere in this
   entry.** (a) A two-arm *paired* figure from the published commentary was being read as a
   single-arm permutation demonstration; it is a paired construction in which the second panel is
   pseudobulk under the *same* permuted labels, not the naive arm standing alone. (b) The eLife
   review record was being read as endorsement of that permutation demonstration; the reviewers
   **constrained** it — explicitly declining to conclude that all of the re-analysed DEGs were
   wrong — and used the pseudobulk panel *against* pseudobulk at small cell counts. (c) A
   "size-then-power" principle was cited under a name that does not exist in the methodological
   literature, and the supporting quotation was truncated at the clause that reverses it: its next
   item re-couples test size to power, and its antecedent is the *conservative* test, which in this
   comparison is pseudobulk. **The case below rests on nothing external.** Where a reader wants
   precedent, the honest statement is that we have none and are not claiming any.
3. **The argument that "the grid shows calibration failing at n = 4, so our rule is stricter than
   the envelope where it matters" is deleted, on three independent grounds.** (i) **Wrong null** —
   the grid driver measures *fresh* nulls, independently simulated data with real labels; the map
   stands on the donor-permutation null, which is a different object. (ii) **Wrong σ** — every cell
   the argument cited sits at `sigma_donor` = 0.5, where the envelope demands 13 donors per group
   and not 4, and the grid has **no σ = 0.2 cell whatsoever** (verified: its σ values are 0.35,
   0.5, 0.7). (iii) **It bites only where there is no crisis** — the argument needs the anchor to
   land near 0.2, and at 0.2 the envelope keeps 11 of 12 datasets and 227 of 251 strata and the
   present entry is unnecessary. The rule below therefore does **not** claim to be stricter than
   the envelope at the likeliest destination; the four-tier table in the next section is where that
   comparison is made, and it shows the rule stricter at σ ≈ 0.2 and looser at 0.5 and 0.7.
4. **`signal_above_floor` is demoted out of decision-relevant use outside the envelope.** The draft
   had it in the reportable tier without qualification. Its numerator is a count of *real-label*
   naive hits, and a real-label hit count is the one place where the naive arm's map touches a
   quantity the pseudobulk arm is the control for. It stays printable — with its leak coefficient
   beside it — and it enters no GO/NO-GO logic outside the envelope. Measuring and deciding are
   different acts, and only the first is being permitted here.

### Correction to Amendment 3 Change 1 — the prohibition is NARROWED, and that is a relaxation

**The sentence being corrected**, quoted rather than paraphrased:

> It is **not** declared valid outside it, and no result from a stratum outside it may be reported
> as a pbcheck measurement.

**Henceforth:** the clause "no result … may be reported" is replaced by the partition of Change 1.
The declaration of *validity* is untouched — the pseudobulk arm remains declared valid only inside
the envelope, and nothing outside it is called valid.

**This is a relaxation of a prohibition.** After this entry, some results from strata outside the
envelope may be reported that could not be reported before. Describing that as a clarification, or
as "making explicit what was always intended", would be precisely the move that Amendment 1 exists
to punish, and it is refused here in the same words.

**What the relaxation costs and what it buys, per envelope tier** (re-derived from the frozen JSON;
"envelope gives" is the freeze's §6 table, "map rule gives" is `min(n_A, n_B)` ≥ 8 on the same 251):

| anchor lands at σ ≈ | envelope demands | envelope gives | map rule gives | direction |
|---|---|---|---|---|
| 0.2 | ≥ 4 v 4 | 227 strata / 11 datasets | 150 / 7 | **−77 strata, −4 datasets — STRICTER** |
| 0.35 | ≥ 8 v 8 | 150 / 7 | 150 / 7 | identical |
| 0.5 | ≥ 13 v 13 | 94 / 5 | 150 / 7 | +56 / +2 — looser |
| 0.7 | ≥ 23 v 23 | 30 / 3 | 150 / 7 | +120 / +4 — looser |

The rule is inert or stricter under the optimistic anchor and permissive under the pessimistic one.
That is the signature of a conditional pre-commitment rather than of a reaction to a result — but
it is a signature, not a proof, and the reader who distrusts it should read the abandonment
triggers below, which are the part that can be checked.

**And the rule does not repair coverage.** 7 datasets is below spec §1's own floor of 8, and no
arrangement of the frozen twelve fixes that. Where the map runs on the twelve alone it runs below
the spec's dataset floor, and that must be reported as a limitation on every occasion, exactly as
the freeze reports it.

### Change 1 — the §6 partition

§6's metrics are partitioned into three tiers. The partition is by **what a quantity's numerator
and denominator are made of**, not by convenience.

**Tier 1 — reportable at `min(n_A, n_B)` ≥ 8 outside the envelope**, subject to every condition in
Changes 2–8:

* `lambda_naive`, the genomic-inflation factor of the naive arm's own p-values under the
  donor-permutation null.
* The **naive-arm permutation false-positive floor**: median count, fraction of G, IQR, and
  Monte-Carlo SE, at that stratum's own cells-per-donor.
* The **floor-versus-cells-per-donor curve** (B2, D1), which is the map's only legitimate
  cross-dataset axis.
* Four per-cell controls, printed on every cell and never summarised away: `lambda_pseudobulk`,
  the pseudobulk permutation FP rate with its exact-binomial p-value, the pseudobulk permutation
  floor, and the B5 exchangeability diagnostic (which must sit near 1 and whose departure
  invalidates every other number from the same null).

Tier 1 is admissible because none of these quantities has a pseudobulk rejection anywhere in it.
`metrics.genomic_inflation` and `metrics.perm_floor` take no argument derived from the pseudobulk
arm's power, and the naive arm's per-permutation p-value vector is read from its own raw table
before any BH runs. The pseudobulk arm appears in Tier 1 **only as the negative control** that
establishes the null construction is sound.

**Tier 2 — reportable, never decision-relevant outside the envelope**: `signal_above_floor`, with
the leak coefficient of Change 4 printed on the same line. It is read from the naive arm's own
solo-BH floor (Change 3's plumbing requirement), and it enters no GO/NO-GO logic outside the
envelope.

**Tier 3 — inside the envelope only, with no change of any kind**: `real_label_ratio`,
`concordance`, and **any sentence asserting that a published finding is false**. These have a
pseudobulk hit count in the denominator or a claim about biology in the predicate, and decision
rule item 1 governs them exactly as Amendment 3 left it.

**Why Tier 3 is where it is, measured rather than argued.** In the committed gate artifact — a
synthetic null with donor structure and **truth = 0 DE**, real labels, `sigma_donor` = 0.5, 8 v 8,
G = 1500 — the naive arm calls **1130** genes and the moderated arm calls **1**
(`pilot/gate/synthetic_gate_2026-08-15.json`, `real_label_ndeg`). The real-label ratio there is
1130 on data containing no signal at all. A number that large, produced by the denominator's
behaviour on a stratum where the truth is known to be zero, cannot be reported as a count of false
discoveries in anyone's published analysis. That is decision rule item 1's point, and it is why
nothing in Tier 3 moves.

**The reportability contract.** These two sentences are printed under every map panel, verbatim,
with the angle-bracketed slots filled from that cell. They are the contract: a panel that cannot
fill every slot is not published.

> In \<dataset_id\> × \<cell_type\> (\<n_A\> v \<n_B\> donors, median \<M_A\>/\<M_B\> cells per
> donor, G = \<G\> frozen-universe genes), the naive per-cell Wilcoxon test is miscalibrated
> against the donor-permutation null by a genomic-inflation factor of λ_naive = \<L\>, and under
> that null — every cell keeping its real donor and its real counts, only the donor→condition map
> reshuffled at fixed group sizes, so the correct answer is no association — it calls a median of
> \<F\> genes (\<P\> % of G; IQR \<I\>; Monte-Carlo SE \<S\>) at BH-FDR < 0.05, against the ≈ 0
> rejections a calibrated test gives under the complete null, while the donor-level moderated arm
> run on the identical universe and the identical permutations rejects a median of 0 there
> (λ_pb = \<l\> ∈ [0.9, 1.1]; permutation false-positive rate \<r\>, one-sided exact-binomial
> p = \<p\> against α = 0.05).

> This is a sharp-null upper bound on the pseudoreplication false-positive floor at this stratum's
> own cells-per-donor (permutation-to-truth leak E|corr| = \<c\> at \<n_A\> v \<n_B\>; cell-count
> stratification per spec A2 applied with tolerance \<t\>), not a count of false discoveries in any
> published analysis; the donor-level arm appears here solely as the negative control establishing
> that the null construction is sound and is NOT declared powered at this stratum, so no
> naive-to-pseudobulk ratio, no concordance figure, and no statement that any specific published
> finding is false may be read from this panel; and because the pinned Census exposes no
> library/pool identifier, donor pseudobulk is a lower bound on the correct replication unit, which
> makes this floor a lower bound on the true pseudoreplication floor.

**The plumbing this requires, and the live defect it exposes.** The map reads the naive arm's
**solo-BH** floor explicitly, with a per-stratum assertion that the pseudobulk arm NA'd nothing
(`n_na_pseudobulk == 0`) wherever a paired count is also present. No spec change is needed for
that: §5 item 3 requires BH "over that same G-length set", and `mtc`'s own note records that the
single-arm path "remains correct for genuinely within-arm quantities (e.g. the naive arm's own
floor over permutations with no pseudobulk counterpart)".

**But a real defect must be fixed first, and it is invisible today.** At `2d9092b`,
`permutation.run_null` fills its `naive_ndeg` array from the **paired** BH while a paired pseudobulk
fit exists and from the naive arm's **own** BH above that index, and `scripts/synthetic_gate.py`
computes its headline floor over that whole hybrid array while `monte_carlo.naive_floor_median` is
computed over the paired prefix only. The two agree today for one reason and one only:
`gate_config.N_PERM == N_PERM_PB == 200`, so the array has no second half. **At spec §4's
pre-registered counts for the real sweep — `n_perm` = 1000, `n_perm_pb` ≥ 200 — four fifths of that
array come from the other BH convention and the headline floor becomes a median over two different
conventions at once.** No map cell may be computed until a #DEG series carries the convention it was
produced under. This is a latent defect, not a wrong published number, and it is named here rather
than filed.

### Change 2 — the map's admission rule: `min(n_A, n_B)` ≥ 8, on A1 alone

**The rule.** A stratum may carry Tier 1 quantities outside the envelope only if
`min(n_donors_A, n_donors_B)` ≥ 8.

**Its sole justification is correction A1**, frozen 2026-07-19, a month before `sigma_donor` became
a question at all:

> High-donor strata (≥ 8 vs 8, where balanced permutations are ~orthogonal to the true grouping)
> are weighted for the headline floor and are the ONLY strata where the signal-above-floor
> "kill-switch" ratio (metric 4) is treated as decision-relevant.

**The permutation-resolution argument is explicitly NOT part of this justification**, for the
reason given above: C(14, 7) = 3432 already exceeds spec §4's `n_perm` = 1000, so resolution closes
at 7 v 7 and cannot be what picks 8.

**What A1's threshold looks like when the leak is computed rather than asserted.** Define, for a
balanced donor-label permutation drawing `n_A` donors from `D = n_A + n_B`, the Pearson correlation
between the true and permuted condition indicators:

```
corr(k) = (k·D − n_A²) / (n_A · n_B),      k = |true test set ∩ permuted test set| ~ Hypergeometric
```

and take `E|corr|` over the permutation set `build_perms` actually returns — the identity always
removed, its exact complement removed as well when `n_A = n_B`. On balanced designs:

| donors per group | 3 v 3 | 4 v 4 | 5 v 5 | 6 v 6 | 7 v 7 | **8 v 8** | 9 v 9 | 10 v 10 |
|---|---|---|---|---|---|---|---|---|
| E\|corr\| | 0.3333 | 0.2353 | 0.2800 | 0.2148 | 0.2327 | **0.1902** | 0.2015 | 0.1719 |
| P(corr = 0) | 0 | 0.5294 | 0 | 0.4338 | 0 | **0.3808** | 0 | 0.3437 |

**The leak is a sawtooth, not a curve**, and 8 v 8 is a **strict local minimum** — below 7 v 7
(0.2327) and below 9 v 9 (0.2015), and the smallest value anywhere at or below 9 donors per group.
A1's threshold, chosen for a verbal reason, lands on the best point in its neighbourhood. That is
corroboration, not the argument.

8 is also strictly stricter than every donor-count floor the study already carries: the inclusion
gate's 3 per group, and the envelope's most permissive row, which demands 4.

### Change 3 — measured per-stratum calibration replaces extrapolated calibration

**No cell may be published on the strength of the synthetic gate's calibration.** Today the only
evidence that the moderated arm is calibrated is one synthetic point at
(`sigma_donor` = 0.5, 8 v 8, G = 1500). Extending publication outside the envelope on an
extrapolation from one simulated cell would be the same error the envelope was invented to stop.
Each cell of the map earns its own calibration, on its own data, under its own permutations.

**Both criteria, always both. The word "calibrated" may never appear on a cell that has not passed
both.**

1. `lambda_pseudobulk` ∈ [0.9, 1.1] — the pre-registered band, unchanged.
2. The pseudobulk permutation false-positive rate is **not significantly above α**, by an exact
   one-sided binomial test at a **per-cell level of 0.01** — not by a raw comparison against 0.05.

**Why the raw comparison is refused, with the arithmetic.** Let a stratum be *perfectly* calibrated,
so its true FP rate is exactly α = 0.05. The rule "exclude if the observed rate exceeds 0.05"
excludes it with probability

```
n_perm = 200   :  P(X > 10)  = P(X ≥ 11) = 0.4169
n_perm = 1000  :  P(X > 50)  = P(X ≥ 51) = 0.4625
n_perm = 10000 :  P(X > 500) = P(X ≥ 501) = 0.4881          (X ~ Binomial(n_perm, 0.05))
```

**The false-exclusion rate rises with computation and tends to 1/2**, because the estimate
concentrates on α and the rule discards every cell that lands on the wrong side of it. A criterion
that gets worse the harder you work is not a criterion.

**The thresholds, and why they are the smallest that hold the level.** At a per-cell level of 0.01,
the cell fails when

```
n_perm = 200   :  k ≥ 19   (19/200 = 0.0950)     P(X ≥ 19 | 200) = 0.00582 ;  P(X ≥ 18 | 200) = 0.01209
n_perm = 1000  :  k ≥ 68   (68/1000 = 0.0680)    P(X ≥ 68 | 1000) = 0.00741 ;  P(X ≥ 67 | 1000) = 0.01059
```

One step lower breaches 0.01 in both cases, so these are the strictest thresholds the declared level
admits; neither is rounded and neither is chosen.

**Multiplicity is declared, not corrected.** At the nominal per-cell level of 0.01 over the 150
strata of the ≥ 8 v 8 map, P(at least one false exclusion) = 1 − 0.99¹⁵⁰ = **0.7785** and the
expected count is **1.5**. At the *realised* discrete thresholds the same arithmetic gives
**0.5836** and **0.874** at `n_perm` = 200, and **0.6722** and **1.111** at `n_perm` = 1000. Both
are stated because the first is the level we declared and the second is what the integer threshold
actually delivers. **No correction is applied**, in Amendment 4's own words: *a multiplicity
correction here would be a loosening dressed as rigour.* The direction of the error is exclusion,
which costs sample size and is reported in D4's bookkeeping — the direction Amendment 4 Change 4
already fixed the project's asymmetry on: **"The project is defended against understatement."**

**A weakness in the test itself, and it is an order of magnitude larger than a disclosure of its
direction would suggest.** The exact binomial treats the `n_perm` permutations of a stratum as
independent Bernoulli trials. They are not: every permutation reuses the same donors and the same
counts, so the rejection indicators are positively dependent, the count's variance exceeds the
binomial variance, and the true tail is heavier than tabulated. That much follows from the
construction. **The magnitude does not, so it was measured rather than left as a direction.**

On 100 simulated strata at the gate's own geometry (G = 1500, 8 v 8, 250 cells per donor,
dispersion 0.2), each run under 200 balanced donor-label permutations of its own pre-aggregated
donor profiles — **the exact object this criterion tests, on the permutation null and not on a
fresh one** — the realised per-cell false-exclusion rate of `k ≥ 19 of 200` on **perfectly
calibrated** strata is:

| `donor_sigma` | strata failing `k ≥ 19` | realised rate | declared rate | ratio | P(that many or more at the declared rate) |
|---|---|---|---|---|---|
| 0.2 | **15 of 100** | 0.150 | 0.00582 | **25.8×** | 4.7 × 10⁻¹⁷ |
| 0.5 | **5 of 100** | 0.050 | 0.00582 | **8.6×** | 3.2 × 10⁻⁴ |

Mean per-stratum permutation-null FP rate 0.0647 at σ = 0.2 and 0.0476 at σ = 0.5, over 20 000
permutations per σ; seeds 20 261 003 000–099 and 20 261 004 000–099. The permutation construction
was **transcribed from `build_perms` at `2d9092b` rather than imported** — balanced sets of size
`n_A`, identity and exact complement excluded, rejection-sampled because C(16, 8) = 12 870 > 200 —
because that module is under concurrent edit in this tree. A future reader who diffs the two
constructions should find this note rather than an unexplained second implementation; if the two
ever disagreed, the measurement would be of something else. Paired and solo BH coincide at this
geometry (neither arm NaNs), so the choice of path does not enter.

**So the 0.01 is nominal and the realised level is one to two orders of magnitude above it.** On a
150-cell map that projects to roughly **22 VOID cells at σ = 0.2 and 8 at σ = 0.5 caused by the
criterion's own dependence structure rather than by any stratum's miscalibration.** Three
consequences, pre-declared here so that none of them can be decided after a run:

1. **No threshold moves.** `k ≥ 19 of 200` and `k ≥ 68 of 1000` stand exactly as derived above. The
   arithmetic that produced them is correct arithmetic under a stated assumption; the assumption is
   what fails. The response is to measure it and publish the cost, not to re-tune the level until
   the exclusion count looks comfortable — that would be the loosening this log exists to catch.
2. **The direction is exclusion, which is the safe side — but safe is not free.** Losing a
   calibrated stratum costs sample size and is counted in D4; admitting an invalid one breaks the
   denominator, and Amendment 4 Change 4 already fixed the project's asymmetry that way. At 15 %,
   however, the criterion discards about one cell in seven for the instrument's reason, and that
   loss is an instrument property and must be reported as one, never as a property of the strata.
3. **The grid must measure this before Change 3 runs on real data**, at every envelope σ including
   0.2, and the measured rate — not the binomial nominal — is what the abandonment triggers below
   are evaluated against. Until that measurement exists, the criterion may not be applied.

**A second weakness, about the arm rather than the test.** Change 3 admits a cell on a measurement
of the *arm's* calibration at that cell. Where the shipped arm is itself mildly anticonservative, a
cell fails for a reason belonging to the instrument and not to the stratum. That is not
hypothetical: the arm's fresh-null rejection rate at σ = 0.2 sits above nominal (tabulated under
"What this does NOT settle"), so part of the 15 % above belongs to the arm and part to the
dependence, and **this measurement does not separate them.** Separating them is the grid's job, and
it is named as such rather than assumed away.

**No number crosses between the two nulls.** Everything in this section is the **donor-permutation**
null — donor pseudobulk profiles held fixed, only the labels moving. The figures under "What this
does NOT settle" are **fresh** nulls, independently simulated data with real labels. The two have
different dependence structures, they need not agree, and here they do not: at σ = 0.5 the fresh
null reads nominal while the permutation null still excludes at 8.6× the declared rate. A reader
must not carry a value from one to the other, and this entry never does.

### Change 4 — the per-cell leak coefficient, keyed by `(n_A, n_B)` and never by `min`

Every cell prints `E|corr|` computed at its **own** `(n_A, n_B)`, by the formula of Change 2, over
the permutation set that cell actually used.

**Keying it by `min(n_A, n_B)` would mis-order the map, and by how much is measurable.** Within the
frozen list, the label "min = 8" covers four realised shapes — 8 v 9, 8 v 10, 8 v 24, 8 v 27 — whose
leak coefficients run **0.1344 (8 v 24) to 0.2015 (8 v 9 and 8 v 10)**, a spread of 50 % inside one
label. Over the pinned manifest's 27 distinct min-8 shapes the same label spans **0.0730 (8 v 138)
to 0.2015**, a factor of 2.8. Across the whole frozen ≥ 8 v 8 set — 76 distinct shapes over 150
strata — it runs 0.0737 to 0.2015. And the map is overwhelmingly unbalanced: only **21 of the 150**
frozen ≥ 8 v 8 strata have `n_A = n_B` (42 of all 251), so on **129 of the 150** a `min` key would
hand the reader a leak figure that is not the cell's own.

**Exact orthogonality has an exact condition, and it is not "n is even".** `corr = 0` requires
`k = n_A²/D`, i.e. **`D | n_A²`**. Two consequences that a `min` key hides:

* 8 v 24 satisfies it — D = 32 divides 64 — and **P(corr = 0) = 0.3583**: better than a third of
  its permutations are exactly orthogonal to the truth.
* 9 v 9 does not — D = 18 does not divide 81 — and **P(corr = 0) = 0 exactly**: no permutation of a
  9 v 9 design is orthogonal to the true grouping, despite 9 being larger than 8.

Only **20 of the frozen 150** ≥ 8 v 8 strata admit an exactly orthogonal assignment at all. The
coefficient is a per-cell property and is printed as one.

### Change 5 — the per-cell void rule

Decision rule item 1's VOID is a *study-level* verdict and there is no per-stratum analogue anywhere
in the spec or in Amendments 1–4. The quantities exist — `run_null` already emits `pb_fp_rate`, the
pseudobulk floor, both λs and the Monte-Carlo errors — and no rule consumes them. That gap is closed
here, before any cell is computed.

**A cell is VOID when its own `lambda_pseudobulk` falls outside [0.9, 1.1], or its own pseudobulk
permutation FP rate fails Change 3's binomial test, or its A2 stratification cannot be constructed
under Change 6.** A VOID cell publishes exactly four things and nothing else: its identifiers, its
`lambda_pseudobulk`, its FP count `k` out of `n_perm` with the exact-binomial p-value, and the word
VOID. It contributes to **no** curve, **no** aggregate, **no** per-dataset summary and **no**
figure, and it is counted in D4's excluded-strata bookkeeping with its reason.

A void cell is a result. The count of void cells is reported beside the map, and if that count is
large the abandonment trigger below fires rather than the map being trimmed to its survivors.

### Change 6 — A2 is UN-deferred at ≥ 8 v 8

**This is the price of the relaxation and it is a blocking condition, not an aspiration. If A2
stratification is not implemented and passing at ≥ 8 v 8, the map is not published — inside the
envelope or outside it.**

**Amendment 2 Change 6 deferred A2 on two arguments. One was about 3 v 3 and inverts here; the other
is correct and is exactly why A2 now binds.**

The infeasibility argument, quoted:

> It is infeasible at the donor counts the spec's own inclusion gate admits. At the 3v3 minimum
> there are C(6,3) − 2 = 18 permutations in total; filtering them on cell-count proximity leaves
> single digits, and the floor's Monte-Carlo error would swamp the confound it is meant to remove.

That is a true statement about 3 v 3 and it does not survive the move to 8 v 8. **The smallest
member of the ≥ 8 v 8 tier holds C(16, 8) = 12 870 label assignments** — verified as the minimum
`permutation_count` over the manifest's 554 candidates at that tier — and 13 v 13 holds
C(26, 13) = 10 400 600. There is room to restrict on cell-count proximity and still sample
`n_perm` = 1000.

The asymmetry argument, quoted:

> The confound is asymmetric between the arms, and the pseudobulk arm is the one that matters here.
> … It is the naive per-cell arm whose statistic scales with cell count. Since the binding validity
> gate this amendment exists to repair is about the *pseudobulk* arm, A2 does not gate it.

Correct, and it is the reason A2 binds now. **The map is made of the naive arm's floor — precisely
the quantity A2 protects.** Amendment 2 said so in the same breath: "It does bear on the naive
arm's floor, which is a Phase 1 headline quantity."

**The range check that stands in for A2 today has never fired where it could.** On the synthetic
oracle every donor has the same number of cells by construction, so the check is degenerate and its
passing carries no information about real strata.

**On real strata the confound is large, and here it is measured** (re-derived from the committed
JSONs; within-group figures are over the two group-arms of each stratum):

| population | within-group cells/donor max÷min | between-group total-cell ratio |
|---|---|---|
| frozen 150 at ≥ 8 v 8 (300 group-arms) | median **18.0**, **70.7 %** above 10, max **524.5** | median **1.43**, **21.3 %** above 3, max **40.4** |
| manifest 554 at ≥ 8 v 8 (1108 group-arms) | median **30.4**, **81.8 %** above 10, max **3279.6** | median **2.32**, **40.1 %** above 3, max **130.6** |

Two in five of the strata the mechanical extension would admit have one group carrying more than
three times the other's cells. A floor compared against permutations of a systematically different
size is a size measurement, not a replication-unit measurement.

**What must be delivered before any cell is published**, and it is statistical work rather than
wiring: the matching tolerance, stated as a number before the run; the behaviour when the restricted
permutation set is thin, stated as a rule and not decided per stratum; a statement of whether the
restricted null remains exact or becomes approximate; and the declared consequence — a stratum whose
restricted set cannot be built under the tolerance is VOID by Change 5, never silently unstratified.

### Change 7 — the common cells-per-donor target T, pre-registered as a number

B2 makes cross-dataset floor comparison legal only "after bootstrapping cells to a common per-donor
target", and B3 restricts λ to a binary/ordinal flag. So the map's cross-dataset object is the
**floor-versus-cells-per-donor curve**, never a single floor, and the matched comparison needs a
target. **No amendment and no section of the spec has ever supplied one**, which leaves T a free
knob of the analyst that moves both the size of the map and every floor in it.

**T = 30 cells per donor.**

Two reasons, in order. First, 30 is the **smallest pre-registered D1 bin edge strictly above the
inclusion gate's own floor of 10**, and a *lower* target is the stricter direction against this
study's own thesis: the floor and `lambda_naive` grow with cells per donor (spec §10 risk 2), so
subsampling deeper makes decision rule item 2's persistence requirement harder, not easier. This is
the same reasoning the freeze used to put the bin floor at 10 rather than higher. Second, 10 itself
is unusable as a matching target: a donor sitting exactly at the inclusion gate's floor cannot be
subsampled at all without being dropped.

**Its reachability is disclosed as a consequence, not used as the reason.** Of the frozen 150 at
≥ 8 v 8, **94.3 %** of the 300 group-arms have a per-group median at or above 30 (against 64.0 % at
100); of the manifest's 554, **88.9 %** of 1108 do. The extension of Change 8 therefore does not
move T, and T is fixed once.

**One honest qualification.** Reaching T requires *per-donor* counts at or above T, and the frozen
artifact carries only per-group medians: only **45.0 %** of the frozen ≥ 8 v 8 group-arms have a
*minimum* per-donor count at or above 30. How many donors a cell loses at matching is a load-time
fact, and it is reported per cell rather than assumed away.

### Change 8 — the selection rule if the map goes past the frozen twelve

Declared now, as a number, because after the anchor it could not be.

**If the map is extended beyond the frozen twelve, it is extended by exactly one mechanical rule:
every candidate stratum in the pinned manifest with `min(n_A, n_B)` ≥ 8. That is 554 strata over 33
datasets in 25 collections.** Never a list of datasets, never a subset chosen for coverage, never a
substitution for a stratum that failed a gate. The freeze's prohibition binds this extension exactly
as it binds the twelve: shrinkage is a reported outcome, never a re-selection.

Collections, not datasets, are the honest unit — two datasets of one collection are not two
independent choices under D2 — which is why the collection count is declared beside the dataset
count and must be reported beside it.

**Whether and when to extend is a separate decision and is NOT taken here.** Both admission maps are
pre-declared so that the choice between them cannot be made after the anchor is read: the frozen
twelve's map is 150 strata over 7 datasets, and the mechanical extension is 554 over 33 datasets in
25 collections. Declaring both costs nothing now and removes a degree of freedom later. Execution
order — twelve first, manifest later, or one pass — is an engineering decision that Part B or a
later entry records with its reasons.

### Pre-declared counts, so the surviving set cannot be chosen later

Every count below is re-derived in this entry from the committed JSON artifacts, and every one is
an **upper bound**: `integer_check` and `frozen_universe_size` are `pending` on all 251 rows and can
only remove strata.

| declared object | strata | datasets | collections |
|---|---|---|---|
| the map on the frozen twelve, `min(n_A, n_B)` ≥ 8 | **150** | **7** (below §1's floor of 8) | — |
| the mechanical extension, pinned manifest, `min(n_A, n_B)` ≥ 8 | **554** | **33** | **25** |
| envelope tiers on the frozen twelve, ≥ 4 / 8 / 13 / 23 v same | 227 / 150 / 94 / 30 | 11 / 7 / 5 / 3 | — |
| manifest ladder, ≥ 4 / 7 / 8 / 13 / 23 v same | 1017 / 630 / 554 / 311 / 134 | 62 / 38 / 33 / 21 / 12 | 46 / 30 / 25 / 15 / 10 |

### Pre-declared abandonment triggers

Each is an observable with a number attached and a declared action. They are here, before the
anchor, because a trigger invented after a disappointing measurement is not a trigger.

1. **The anchor lands at σ ≤ 0.2 with a tight interval.** The envelope then holds 227 strata over 11
   datasets, spec §1's floor of 8 is met, and this entry's rule is pure narrowing — it would cost 77
   strata and 4 datasets for nothing. **Action: the rule is WITHDRAWN and everything is published
   under the envelope.** This is the cleanest of the triggers and it is directly observable.
2. **The A2-stratified floor disagrees with the unstratified floor.** If, on **more than 20 %** of
   the ≥ 8 v 8 strata, the stratified median floor differs from the unstratified median floor by
   **more than 2 × the Monte-Carlo SE**, the floor is an artifact of cell counts rather than of
   pseudoreplication, and a map of it is not a map of pseudoreplication. **Action: nothing is
   published, inside the envelope or outside it.** This trigger kills the map outright and it is the
   one the whole relaxation is bought against.
3. **Calibration does not survive contact with real data.** If the number of ≥ 8 v 8 strata failing
   Change 3's binomial leg exceeds **5 × the expected count** on a perfectly calibrated set,
   "calibrated" is not a property that transfers to real data and the map has no reference arm.

   **The baseline for that expectation is NOT the binomial nominal, and the reason is measured.**
   A draft of this trigger set E from the exact-binomial tail — E = 0.874 of 150 at `n_perm` = 200,
   E = 1.111 at `n_perm` = 1000 — which would have fired at 5 and 6 respectively. Change 3's
   measurement shows that baseline is wrong by 8.6× to 25.8× on *calibrated synthetic* strata,
   because the permutations of one stratum are dependent. Against it, a perfectly calibrated map
   would fire this trigger immediately and for the instrument's reason rather than the data's.
   **E is therefore the grid's measured per-cell false-exclusion rate at the σ of the map's own
   tier, and the trigger cannot be evaluated until that measurement exists.** Provisionally, from
   this entry's 100-stratum probe, E would be ≈ 22 of 150 at σ = 0.2 and ≈ 8 at σ = 0.5 — quoted as
   the order of magnitude to expect, never as the baseline itself.

   **Two independent probes of that rate disagree by a factor of two, which is itself the argument
   for waiting on the grid.** At σ = 0.2, 8 v 8, `n_perm` = 200: this entry's 100-stratum probe
   reads 15 of 100 (0.150, 25.8 × the nominal), and a separate 40-stratum probe on a disjoint seed
   block (777000–777039) reads 3 of 40 (0.075, 12.9 × the nominal, exact-binomial p = 1.7 × 10⁻³
   against the nominal), with mean per-stratum permutation FP 0.0601 against the first probe's
   0.0647. The two are consistent under sampling — 3/40 against 15/100 does not separate — and the
   qualitative finding is the same in both: the declared 0.582 % is wrong by an order of magnitude
   and the criterion voids calibrated strata at a rate between roughly 8 % and 15 %. But the
   *number* is not settled to better than a factor of two by either, and Trigger 8's 10 % threshold
   sits between them. That is deliberate: a threshold the present evidence straddles is one the
   grid can genuinely fail, which is the only kind worth pre-registering.

   **The direction of this correction is stated because it is the uncomfortable one.** Raising E
   makes the trigger *harder* to fire, which is the loosening direction, and it is done only
   because the number being replaced is measured-wrong rather than merely inconvenient. The
   compensating control is that it is now **blocking**: with no measured E the criterion may not be
   applied to real data at all, so the trigger cannot be quietly skipped by leaving E unmeasured.
   **Action: the map stops, and the reference arm becomes the subject of a new amendment.**
4. **`lambda_naive` at matched cells-per-donor does not exceed `LAMBDA_NAIVE_GO` = 2.0 in a majority
   of independent datasets.** Decision rule item 2 has then failed and there is nothing to publish
   at any scope. **Action: no map, narrow or wide; the finding is the failure.**
5. **The frozen universe comes in below `MIN_UNIVERSE_SIZE` = 200 on a large share of ≥ 8 v 8
   strata.** G then varies so much that even fraction-of-G stops being comparable between cells.
   No arithmetic on this gate exists yet — `frozen_universe_size` is `pending` on all 251 rows and
   only loading X decides it. The nearest thing on the record is about a **different** gate and is
   cited as such: Amendment 4's scenario for the σ-estimation stratum's own 100-gene minimum puts
   **1 / 3 / 6 / 8** of the 251 strata below it at G = 5000 / 10 000 / 15 000 / 20 000, with the
   failures concentrated in datasets #7 and #3. It bounds nothing here; it names where to look.
   **Action: cells below `MIN_UNIVERSE_SIZE` are SKIPs with their measured size reported, never
   rounded up; if more than a quarter of the ≥ 8 v 8 cells SKIP on it, the cross-cell comparison is
   withdrawn and only within-stratum numbers stand.**
6. **§8(d) runs and the permutation floor does NOT account for most naive calls at ≥ 8 v 8.** The
   binding real anchor then contradicts the thesis, and what is in question is the whole study, not
   the scope of publication. **Action: back to this log, with the GO/NO-GO rule, not the map, as the
   subject.**
7. **Engineering.** If the property tests of any accelerated permutation engine do not show exact
   agreement with the shipped `scanpy` statistic, the cost model behind Change 8 is void. **Action:
   the mechanical extension is deferred and only the frozen twelve's 150 are attempted.** No
   core-hour figure is quoted here: the speed-up claim in circulation while this entry was drafted
   has not been checked adversarially, and an unverified number has no place in a pre-registration.
8. **Change 3's criterion voids cells for the instrument's reason rather than the stratum's.** If
   the grid's measured per-cell false-exclusion rate of `k ≥ 19 of 200` — on **calibrated
   synthetic** strata, at the σ of the map's own tier, on the **permutation** null — exceeds
   **10 %**, then more than one map cell in ten is discarded by the instrument, the admission rule
   is not a calibration test of the stratum, and the map's coverage figure means nothing.
   **This trigger is live rather than hypothetical**: this entry's own 100-stratum probe already
   reads **15 % at σ = 0.2** and **5 % at σ = 0.5**, so the grid may well fire it, and it is
   written at a threshold the current evidence straddles rather than at one chosen to clear.
   Scope: the trigger is on the **criterion and the arm**, not on the σ = 0.2 row — the σ = 0.2
   confinement holds on the fresh null (χ² = 9.90, p = 0.0071) but *not* on the permutation null,
   where σ = 0.5 still runs at 8.6 × the declared rate. **Action: the map is not published; the
   arm's null behaviour becomes the subject of its own amendment; and no threshold in Change 3 is
   adjusted to bring the rate down.**

### What this NARROWS and what it RELAXES

**RELAXED — one thing, named without euphemism.** Amendment 3 Change 1's sentence "no result from a
stratum outside it may be reported as a pbcheck measurement" no longer holds for Tier 1 quantities
at `min(n_A, n_B)` ≥ 8. Some results from outside the envelope become reportable that were not.

**NARROWED — five things that did not exist before this entry.**

1. A donor-count floor of 8 for the map, against the inclusion gate's 3 and the envelope's most
   permissive row's 4.
2. Measured per-stratum calibration on both criteria, replacing an extrapolation from one synthetic
   cell — and a per-cell VOID rule that discards a failing cell entirely.
3. A leak coefficient printed per cell and keyed by `(n_A, n_B)`, where nothing was printed before.
4. A2 stratification, deferred since Amendment 2, made a blocking condition of publication at
   ≥ 8 v 8.
5. A pre-registered matching target T and a pre-registered extension rule, both of which were free
   analyst knobs until this entry.

**MOVED — nothing.** α, the λ band, `POWER_TARGET`, the oracle's log2FC and K, the envelope table,
§1's 8-dataset floor and decision rule item 1 all stand exactly as they stood at `2d9092b`.

### What this does NOT settle

* **The powered claim is untouched and stays inside the envelope.** `real_label_ratio`,
  `concordance` and any sentence asserting that a published finding is false remain Tier 3 under
  decision rule item 1 as Amendment 3 left it. This entry licenses a *measurement* of the naive
  arm's own null behaviour; it licenses no conclusion about anyone's published result.
* **The shipped arm's fresh-null rejection rate at `sigma_donor` = 0.2 sits ABOVE nominal, in three
  independent seed blocks, and this entry does not resolve it.** A draft reported the excess; a
  first check here failed to reject it and was written up as a failure to *reproduce* it. **That
  was a category error and the sentence is withdrawn** — a p of 0.21 is not evidence of absence.
  The quantity is P(≥ 1 BH rejection at FDR 0.05) on **fresh** nulls with real labels, at the
  gate's own geometry (G = 1500, 8 v 8, 250 cells per donor, dispersion 0.2):

  | seed block | n | P(≥ 1 rejection) | one-sided exact binomial vs α |
  |---|---|---|---|
  | 1000–1999 | 1000 | 0.0560 | 0.210 |
  | 5000–5299 | 300 | 0.0900 | 0.0026 |
  | 20 261 000 000–20 261 000 999 | 1000 | 0.0700 | 0.0035 |
  | **pooled** | **2300** | **0.0665** | **0.00029** |

  **Homogeneity first: the three blocks do not separate** — χ² = 4.64, df = 2, **p = 0.098**. The
  apparent conflict between the first two (Fisher p = 0.043) does not survive the third block, so
  the correct description is one rate near 0.066, not three rates in disagreement. Pooled, the
  excess is real (p = 0.00029) and mild — about 1.3 × nominal.

  **Across σ, again homogeneity before interpretation.** At n = 1000 per σ on disjoint blocks:
  σ = 0.2 → 0.0700 (p = 0.0035); σ = 0.35 → 0.0420 (p = 0.89); σ = 0.5 → 0.0440 (p = 0.83).
  **χ² = 9.90, df = 2, p = 0.0071**, and σ = 0.2 against the other two pooled gives Fisher
  p = 0.0022. The blocks **do** separate: on the fresh null the excess is confined to the low-σ end
  and is not uniform across the envelope.

  **That is the convenient reading, and it is the one to handle most carefully.** σ = 0.2 is the
  envelope row with **no grid support at all** — `OPERATING_ENVELOPE` records its own `grid_support`
  as "not in the grid", Amendment 4 calls it "the largest unvalidated dependency the membership rule
  has", and freeze §6 puts 11 of 12 datasets and 227 of 251 strata in that tier. An arm-side excess
  sitting exactly on the unmeasured row that nearly every real stratum will land on earns **more**
  disclosure, not less, however reassuring σ = 0.35 and 0.5 look.

  **And on the null that actually gates, it is not confined to σ = 0.2.** Change 3's measurement of
  the permutation-null counterpart excludes calibrated strata at 8.6 × the declared rate even at
  σ = 0.5, where the fresh null reads nominal. Fresh-null and permutation-null figures are different
  objects with different dependence structure; neither may be read as the other.

  **Three things this is not.** Not a broken arm: under the complete null BH's P(any rejection) is
  exactly the quantity FDR control bounds, so 0.066 against 0.05 is mild anticonservatism,
  plausibly a dependence effect of the shared donor random effect and the library-size
  normalisation. Not a licence to touch anything: no threshold moves for it, the resulting
  exclusions run in the safe direction, and it returns to this log as a result. And not settled —
  the row must be measured by the grid, on both nulls, before Change 3's criterion is applied to
  real data, and trigger 8 puts a number on it.

  **A note on the seed base, so the constraint is visible rather than asserted.** The third block
  uses 20 261 000 000 + i and not the natural 20 260 816 000 + i, because the latter is exactly
  `seed(0, r)` of Amendment 4 Change 5's declared confirmatory grid
  (`1000·(20260816 + i) + r`, i.e. 20 260 816 000–20 260 995 999). Spending it here would have made
  the first cell of an unrun pre-registered grid seen for a different estimand. The base was forced
  by that collision, not preferred; all blocks used here are disjoint from it and from each other.

  The map does not depend on the σ = 0.2 row, which is its one genuine advantage over the envelope
  path — and the advantage does not extend to Tier 3, which still needs the envelope and therefore
  still needs the row.
* **`pooled` is `unresolved` on 251 of 251 and 1197 of 1197**, and `pool_columns_detected` is empty
  in the manifest header. Donor pseudobulk is a lower bound on the correct replication unit, so
  every floor in the map is a lower bound on the true pseudoreplication floor. The disclaimer above
  carries that on every cell rather than once in the verdict, because a caveat stated once and
  attached nowhere is not a caveat. No choice of datasets fixes it; it is a property of the pin.
* **§8(d), Mathys 2019, remains binding and unrun**, additionally blocked on a ROSMAP data-use
  agreement begun 2026-08-16 and on a second, unwritten loader. Nothing here reduces its standing as
  the binding real check.
* **`integer_check` and `frozen_universe_size` are `pending` on all 251 rows and all 1197
  candidates.** Every stratum count in this entry is an upper bound that a real run can only reduce.
* **D5 is not harmonised, and harmonising it would re-cut the frozen strata.** The freeze records
  that collapsing 124 cell-type labels to a common ontology depth merges strata and therefore
  changes what the 251 are. A cross-dataset map is hard to read without it and this entry does not
  supply it.
* **`lambda_naive` is not a hidden σ coordinate, and the reason it is not is measured rather than
  assumed.** At fixed geometry it rises monotonically with `sigma_donor`, which is why the worry
  arises; but at fixed σ it moves by orders of magnitude with cells per donor, so it does not
  identify σ and cannot be read as a back door into the envelope's coordinate. Change 7's matched
  target and the per-cell σ̂, once it exists, are how the confusion is kept out of the map.
* **A2's tolerance is not chosen here.** Change 6 makes A2 blocking and states what must be
  delivered; it does not pick the number, and the number must be pre-registered before the run in
  its own entry or in Part B.
* **The GO/NO-GO decision is not taken, and nothing here moves it.** Widening what may be measured
  licenses a measurement, not a conclusion.

*Author attests: every figure in this entry was re-derived while writing it, from the committed
artifacts at `2d9092b`, using this repository's own virtual environment — the tier tables and
donor-count ladders from `pilot/preregistration/stratum_list_2026-08-16.json` and
`census_candidates_run31910799023_2026-08-15.json`, the gate readings from
`pilot/gate/synthetic_gate_2026-08-15.json`, the grid's σ coverage from `pilot/testsel/summary.json`
— or computed from first principles and shown here with its arithmetic. The leak coefficients are
exact hypergeometric expectations over the permutation set `build_perms` constructs, not
simulations. The binomial thresholds and the false-exclusion rates were computed here and the one
step below each threshold is shown so the choice can be checked. Where this entry differs from the
draft an adversarial read broke, the differences are listed at the top rather than folded in
silently: a false resolution argument, three external precedents, one grid argument and one
demotion — and no threshold among them was loosened to make anything pass. Where a figure in
circulation could not be reproduced, the failure to reproduce is reported in place of the figure
and the probe that failed to find it is described with its seeds; where a claim rests on unverified
work, no number from it is quoted. The relaxation this entry performs is stated as a relaxation of
a prohibition, not as a clarification of one. No real data informed this amendment, no stratum has
been admitted, and `sigma_donor` remains unmeasured on every row.*
