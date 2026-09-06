# pbcheck audit report

## 1. Header

pbcheck audit of '(no file)'

pbcheck 0.1.0, generated 2026-09-06T07:18:48Z.

Status: design_only (at least one donor was measured under both conditions, which makes this a paired design; a paired design needs a paired or mixed model, which pbcheck does not implement in this release, and the donor-permutation null would treat the two halves of one donor as independent)

> This is a pbcheck v0.1.0 audit of one stratum of one file. It measures how a naive per-cell test and a donor-pseudobulk test behave on this data under a donor-permutation null. It is run outside pbcheck's pre-registered Phase 0 protocol, is not a Phase 0 measurement, and makes no claim about any publication or about whether any reported result is true.

## 2. Read-out

No detection arms were run at this status: at least one donor was measured under both conditions, which makes this a paired design; a paired design needs a paired or mixed model, which pbcheck does not implement in this release, and the donor-permutation null would treat the two halves of one donor as independent.

> Amendment 3 declares an operating envelope for the pseudobulk arm on synthetic oracles. For each donor-variance point it states the donor count per group at which the power target is reached, and each row says what the committed grid supports there: a point measured on the grid, or a count derived or extrapolated from it.  
> sigma_donor 0.2: at least 4 donors per group (power at least 0.6 at log2FC 1.0 in 200 genes; grid support 'not in the grid; Amendment 1 frontier only')  
> sigma_donor 0.35: at least 8 donors per group (power at least 0.6 at log2FC 1.0 in 200 genes; grid support 'ebayes power 0.793 at 8v8 (calibrated) -> n* <= 8')  
> sigma_donor 0.5: at least 13 donors per group (power at least 0.6 at log2FC 1.0 in 200 genes; grid support 'ebayes power 0.486 at 12v12, the largest n tested -> n* > 12')  
> sigma_donor 0.7: at least 23 donors per group (power at least 0.6 at log2FC 1.0 in 200 genes; grid support 'ebayes power 0.003 at 8v8 -> n* far above 8')  
> The arm's calibration was evaluated at one hard regime (sigma_donor 0.5, 8 against 8 donors) and nowhere else. pbcheck does not estimate sigma_donor for real data, so whether this stratum lies inside that envelope is not determined here.

> Everything in this report is a diagnostic of this file: the naive arm's inflation factor and permutation floor, the pseudobulk arm's inflation factor, floor and false-positive rate as its negative control, and the bookkeeping of the shared gene universe. None of it is a Phase 0 result, and pbcheck does not rate this file against the Phase 0 decision rule.

## 3. What these words mean

lambda: the genomic inflation factor: the ratio of observed test statistics to the null expectation, with 1.0 meaning no inflation.

permutation floor: the number of genes a test calls under donor permutation, when condition labels carry no real signal; it is the test's own false-positive baseline on this data.

donor-permutation null: the null distribution built by reassigning condition labels between donors, keeping every donor's cells together, and rerunning the test on each reassignment.

gene universe: the fixed set of genes both arms are tested and corrected over, frozen before either arm sees the real labels.

replication unit: the unit whose independent draws the statistics assume; for donor data that unit is the donor, not the cell.

thin-donor filter: the rule that drops a donor's pseudobulk profile when it is built from too few cells or too few counts, rather than keeping a noisy profile.

operating envelope: the region of donor count and donor-to-donor variability that Amendment 3 declares for the pseudobulk arm on synthetic data; each of its points states the donor count per group at which the power target is reached, on the committed grid or by the derivation the grid support of that point names.

sigma_donor: a knob of pbcheck's synthetic simulator for how much donors differ from each other; it cannot be measured on your data, which is why the envelope question is left open.

Monte-Carlo SE: the standard error of a quantity estimated from a finite number of permutations; it shrinks as more permutations are drawn.

BH convention (solo vs paired): solo BH corrects an arm's p-values over the whole gene universe on its own; paired BH corrects both arms together over the genes common to both, so their real-label counts are directly comparable.

## 4. Design audit

**Donors per group**

| group | donors (pre-filter) | pseudobulk profiles (post-filter) |
| --- | --- | --- |
| ctrl | 8 | n/a |
| stim | 8 | n/a |

**Cells per donor (top 8 of 8; 0 more not shown)**

| donor | cells |
| --- | --- |
| patient_1488 | 2745 |
| patient_1256 | 2042 |
| patient_1244 | 2000 |
| patient_1015 | 1786 |
| patient_1016 | 963 |
| patient_101 | 820 |
| patient_1039 | 503 |
| patient_107 | 379 |

No batch columns were provided.

**Design flags**

| Field | Value |
| --- | --- |
| donor nests in condition | False |
| imbalance ratio | 1.021 |
| usable for pseudobulk | False |
| flags | donor spans multiple conditions — not a standard case/control design |

**Dropped cells**

| Field | Value |
| --- | --- |
| dropped, other condition | 0 |
| dropped, other cell type | 13435 |
| dropped, missing donor | 0 |
| dropped, missing condition | 0 |
| dropped, missing cell type | 0 |

## 5. Counts check

Counts check: not run, because at least one donor was measured under both conditions, which makes this a paired design; a paired design needs a paired or mixed model, which pbcheck does not implement in this release, and the donor-permutation null would treat the two halves of one donor as independent.

## 6. Naive per-cell arm

Naive per-cell arm: not run, because at least one donor was measured under both conditions, which makes this a paired design; a paired design needs a paired or mixed model, which pbcheck does not implement in this release, and the donor-permutation null would treat the two halves of one donor as independent.

## 7. Donor-pseudobulk arm

Donor-pseudobulk arm: not run, because at least one donor was measured under both conditions, which makes this a paired design; a paired design needs a paired or mixed model, which pbcheck does not implement in this release, and the donor-permutation null would treat the two halves of one donor as independent.

## 8. Shared universe bookkeeping

**Frozen gene universe**

| Field | Value |
| --- | --- |
| size | 0 |
| minimum size | 0 |
| builder | n/a |

Builder rule: 

Thin-donor filter: not run.

## 9. Settings and provenance

**Tool settings**

| setting | value |
| --- | --- |
| n_perm_requested | 1000 |
| n_perm_pb_requested | 200 |
| seed | 0 |
| top_n | 25 |
| design_only | no |
| naive_method | wilcoxon |
| naive_engine | fast |
| pseudobulk_method | moderated_ebayes |
| trend | no |
| min_donors_per_group | 3 |
| universe_min_total_count | 15 |
| universe_min_prop | 0.5 |
| fallback_universe_min_prop | 0.5 |
| fallback_universe_min_size | 200 |
| min_profiles_per_group_after_thin_filter | 3 |

**Protocol constants**

| protocol constant | value |
| --- | --- |
| alpha | 0.05 |
| lambda_band | [0.9, 1.1] |
| min_universe_size | 200 |
| min_cells | 10 |
| min_counts | 1000 |

**Operating envelope**

| sigma_donor | min donors per group | grid support |
| --- | --- | --- |
| 0.2 | 4 | not in the grid; Amendment 1 frontier only |
| 0.35 | 8 | ebayes power 0.793 at 8v8 (calibrated) -> n* <= 8 |
| 0.5 | 13 | ebayes power 0.486 at 12v12, the largest n tested -> n* > 12 |
| 0.7 | 23 | ebayes power 0.003 at 8v8 -> n* far above 8 |

**Package versions**

| package | version |
| --- | --- |
| numpy | 2.5.2 |
| scipy | 1.18.0 |
| pandas | 3.0.3 |
| anndata | 0.13.2 |
| scanpy | 1.12.2 |
| statsmodels | 0.14.6 |
| decoupler | 2.2.0 |
| pydeseq2 | 0.5.4 |

**Platform**

| Field | Value |
| --- | --- |
| platform | Windows-11-10.0.26200-SP0 |
| python | 3.12.10 |

**Runtime by stage**

| stage | seconds |
| --- | --- |
| load | n/a |
| prepare | 0.3069 |
| design | 0.01704 |
| counts | n/a |
| pseudobulk_build | n/a |
| universe | n/a |
| real_label | n/a |
| permutation_null | n/a |
| render | n/a |

> The engine was measured on one synthetic oracle point (sigma_donor 0.5, 8 against 8 donors, 1500 genes), recorded in pilot/gate/synthetic_gate_2026-08-15.json; that measurement is not repeated on this file. Of the settings below, only the ones listed as protocol constants are pre-registered values taken from pbcheck.gate_config (alpha, lambda_band, min_universe_size, min_cells, min_counts); the permutation counts, the universe filter parameters, the fallback universe rule and the display settings are the tool's own and are not protocol values.

## 10. Footer

Generated by pbcheck 0.1.0. Protocol: docs/PHASE0_SPEC.md and docs/AMENDMENTS.md in the pbcheck repository.
