# pbcheck audit report

## 1. Header

pbcheck audit - scaled_liver_atlas.h5ad

pbcheck 0.1.0, generated 2026-09-01T12:00:00Z.

Status: naive_only (no counts matrix passed the raw-count check)

> This is a pbcheck v0.1.0 audit of one stratum of one file. It measures how a naive per-cell test and a donor-pseudobulk test behave on this data under a donor-permutation null. It is run outside pbcheck's pre-registered Phase 0 protocol, is not a Phase 0 measurement, and makes no claim about any publication or about whether any reported result is true.

## 2. Read-out

Under donor permutation, with no true signal to find, the naive per-cell test calls a median of 41 genes (0.5% of the 6120-gene universe; Monte-Carlo SE 0.63) over 1000 permutations; on the real labels it calls 63.

The naive arm's inflation factor lambda is 3.21 (IQR 0.44): inflated against the band 0.9-1.1.

The donor-pseudobulk arm was not run: no counts matrix passed the raw-count check.

11 donors in lupus_nephritis, 10 in healthy; 352716 distinct donor splits exist.

Real-label calls over the permutation floor: naive 1.54 times the floor.

> The pseudobulk arm's calibration and power were established on synthetic oracles only inside the operating envelope declared in Amendment 3 (minimum donors per group 4 / 8 / 13 / 23 at sigma_donor 0.2 / 0.35 / 0.5 / 0.7). pbcheck does not estimate sigma_donor for real data, so whether this stratum lies inside that envelope is not determined here.

> Everything in this report is a diagnostic of this file: the naive arm's inflation factor and permutation floor, the pseudobulk arm's inflation factor, floor and false-positive rate as its negative control, and the bookkeeping of the shared gene universe. None of it is a Phase 0 result, and pbcheck does not rate this file against the Phase 0 decision rule.

> The count matrix failed the raw-count check (values are not integers). The pseudobulk arm was dropped, never rounded, so this run has no negative control. The naive arm was run on the matrix as found, with the naive pipeline's own normalisation and log transform applied on top of it; its numbers describe that pipeline on this matrix and are not comparable to a run on raw counts.

## 3. What these words mean

lambda: the genomic inflation factor: the ratio of observed test statistics to the null expectation, with 1.0 meaning no inflation.

permutation floor: the number of genes a test calls under donor permutation, when condition labels carry no real signal; it is the test's own false-positive baseline on this data.

donor-permutation null: the null distribution built by reassigning condition labels between donors, keeping every donor's cells together, and rerunning the test on each reassignment.

gene universe: the fixed set of genes both arms are tested and corrected over, frozen before either arm sees the real labels.

replication unit: the unit whose independent draws the statistics assume; for donor data that unit is the donor, not the cell.

thin-donor filter: the rule that drops a donor's pseudobulk profile when it is built from too few cells or too few counts, rather than keeping a noisy profile.

operating envelope: the region of donor count and donor-to-donor variability where the pseudobulk arm's calibration and power were established on synthetic data.

sigma_donor: a knob of pbcheck's synthetic simulator for how much donors differ from each other; it cannot be measured on your data, which is why the envelope question is left open.

Monte-Carlo SE: the standard error of a quantity estimated from a finite number of permutations; it shrinks as more permutations are drawn.

BH convention (solo vs paired): solo BH corrects an arm's p-values over the whole gene universe on its own; paired BH corrects both arms together over the genes common to both, so their real-label counts are directly comparable.

## 4. Design audit

**Donors per group**

| group | donors (pre-filter) | pseudobulk profiles (post-filter) |
| --- | --- | --- |
| lupus_nephritis | 11 | n/a |
| healthy | 10 | n/a |

**Cells per donor (top 20 of 21; 1 more not shown)**

| donor | cells |
| --- | --- |
| donor_00 | 620 |
| donor_01 | 605 |
| donor_02 | 590 |
| donor_03 | 575 |
| donor_04 | 560 |
| donor_05 | 545 |
| donor_06 | 530 |
| donor_07 | 515 |
| donor_08 | 500 |
| donor_09 | 485 |
| donor_10 | 470 |
| donor_11 | 455 |
| donor_12 | 440 |
| donor_13 | 425 |
| donor_14 | 410 |
| donor_15 | 395 |
| donor_16 | 380 |
| donor_17 | 365 |
| donor_18 | 350 |
| donor_19 | 335 |

No batch columns were provided.

**Design flags**

| Field | Value |
| --- | --- |
| donor nests in condition | True |
| imbalance ratio | 1.07 |
| usable for pseudobulk | True |
| flags | none |

**Dropped cells**

| Field | Value |
| --- | --- |
| dropped, other condition | 1500 |
| dropped, other cell type | 600 |
| dropped, missing donor | 90 |
| dropped, missing condition | 40 |
| dropped, missing cell type | 20 |

## 5. Counts check

**Counts check**

| Field | Value |
| --- | --- |
| source | n/a |
| passed | False |
| reason | 3812 of 4,000,000 scanned values are not integers |
| dtype | float32 |
| values scanned | 4000000 |

Examples of the values checked: 0.0, 1.34, 2.71, 0.0.

## 6. Naive per-cell arm

**Real-label counts**

| Field | Value |
| --- | --- |
| real-label calls, solo BH | 63 |

**Top genes, naive arm**

| gene | pval | padj | log2fc | pct group | pct reference |
| --- | --- | --- | --- | --- | --- |
| GOLGA8A | 2.1e-06 | 0.00014 | 1.83 | 71.2 | 22.5 |
| <b> | 3.4e-05 | 0.00089 | -1.21 | 12 | 44.3 |
| TXN2 | 9e-05 | 0.0019 | 0.97 | 55.5 | 30.1 |

**Naive inflation factor**

| Field | Value |
| --- | --- |
| lambda | 3.21 |
| lambda IQR | 0.44 |
| class | inflated |

**Solo permutation floor**

| Field | Value |
| --- | --- |
| median count | 41 |
| median fraction | 0.005 |
| IQR count | 12 |
| bh mode | solo |
| Monte-Carlo SE | 0.63 |

The paired floor is not shown: the pseudobulk arm left 6120 genes without a value, so the paired series is not comparable; the solo floor above stands alone.

**Permutations, naive arm**

| Field | Value |
| --- | --- |
| requested | 1000 |
| achieved | 1000 |

> B5 machinery check, not a calibration criterion: empirical-permutation-p lambda 1.02.

## 7. Donor-pseudobulk arm

Donor-pseudobulk arm: not run: no counts matrix passed the raw-count check.

## 8. Shared universe bookkeeping

**Frozen gene universe**

| Field | Value |
| --- | --- |
| size | 6120 |
| minimum size | 200 |
| builder | naive_detection_fallback |

Builder rule: a gene is kept when it has a value greater than 0 in at least one cell of at least half the donors, computed on the matrix as found; no count-sum rule applies to a non-integer matrix

Thin-donor filter: not run.

**Shared universe, real-label bookkeeping**

| Field | Value |
| --- | --- |
| common tested | 0 |
| dropped for fairness | 0 |
| no value in naive arm | 0 |
| no value in pseudobulk arm | 6120 |
| pseudobulk arm complete for every gene | False |

**Real split against the permutation range**

| Field | Value |
| --- | --- |
| real split inside permutation range | n/a |
| real split percentile in permutations | n/a |

## 9. Settings and provenance

**Tool settings**

| setting | value |
| --- | --- |
| n_perm_requested | 1000 |
| n_perm_pb_requested | 200 |
| seed | 0 |
| top_n | 20 |
| design_only | no |
| naive_method | wilcoxon |
| naive_engine | fast |
| pseudobulk_method | moderated_ebayes |
| trend | no |
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
| numpy | 1.26.4 |
| scipy | 1.13.1 |
| pandas | 2.2.2 |
| anndata | 0.10.7 |
| scanpy | 1.10.1 |
| statsmodels | 0.14.2 |
| decoupler | 1.7.0 |
| pydeseq2 | 0.4.10 |

**Platform**

| Field | Value |
| --- | --- |
| platform | Windows-10-10.0.26200 |
| python | 3.12.4 |

**Runtime by stage**

| stage | seconds |
| --- | --- |
| load | 3.1 |
| prepare | 0.4 |
| design | 0.2 |
| counts | 0.6 |
| pseudobulk_build | 1.8 |
| universe | 0.3 |
| real_label | 5.2 |
| permutation_null | 30.9 |
| render | 0.2 |

> This instrument was calibrated on synthetic oracles (pilot/gate/synthetic_gate_2026-08-15.json). Of the settings below, only those listed under 'protocol constants' are pre-registered values taken from pbcheck.gate_config (alpha, lambda_band, min_universe_size, min_cells, min_counts); the permutation counts, the universe filter parameters, the fallback universe rule and the display settings are the tool's own and are not protocol values.

## 10. Footer

Generated by pbcheck 0.1.0. Protocol: docs/PHASE0_SPEC.md and docs/AMENDMENTS.md in the pbcheck repository.
