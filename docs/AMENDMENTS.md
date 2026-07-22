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
