# pbcheck audit report

## 1. Header

pbcheck audit - rare_cohort_biopsy.h5ad

pbcheck 0.1.0, generated 2026-09-01T12:00:00Z.

Status: complete

> This is a pbcheck v0.1.0 audit of one stratum of one file. It measures how a naive per-cell test and a donor-pseudobulk test behave on this data under a donor-permutation null. It is run outside pbcheck's pre-registered Phase 0 protocol, is not a Phase 0 measurement, and makes no claim about any publication or about whether any reported result is true.

## 2. Read-out

Under donor permutation, with no true signal to find, the naive per-cell test calls a median of 41 genes (0.5% of the 4310-gene universe; Monte-Carlo SE 0.412) over 20 permutations, coarse because few distinct donor splits exist; on the real labels it calls 63.

The naive arm's inflation factor lambda is 3.21 (IQR 0.44): inflated against the band 0.9-1.1.

The donor-pseudobulk arm's lambda is 1.02 (calibrated); its permutation false-positive rate is 0.041 (MC SE 0.014) and its floor a median of 1 genes; on the real labels it calls 5 genes.

3 donors in lupus_nephritis, 4 in healthy; 20 distinct donor splits exist.

After the thin-donor filter (fewer than 10 cells or 1000 counts), 3 and 4 pseudobulk profiles remain.

Real-label calls over the permutation floor: naive 1.54 times the floor; pseudobulk 5 times the floor.

> On this file, with condition labels shuffled between donors and therefore no real signal to find, the per-cell test still calls a median of 41.0 of 8213 genes at FDR 0.05 (0.5%); on the real labels it calls 63. A gene list produced by a per-cell test on this data cannot be separated from that floor. The donor is the replication unit this design supports; the donor-pseudobulk test called 5 genes on the real labels against a permutation median of 1.0.

> The pseudobulk arm's calibration and power were established on synthetic oracles only inside the operating envelope declared in Amendment 3 (minimum donors per group 4 / 8 / 13 / 23 at sigma_donor 0.2 / 0.35 / 0.5 / 0.7). pbcheck does not estimate sigma_donor for real data, so whether this stratum lies inside that envelope is not determined here.

> Everything in this report is a diagnostic of this file: the naive arm's inflation factor and permutation floor, the pseudobulk arm's inflation factor, floor and false-positive rate as its negative control, and the bookkeeping of the shared gene universe. None of it is a Phase 0 result, and pbcheck does not rate this file against the Phase 0 decision rule.

> Fewer than 8 donors in at least one group: the permutation null has few distinct donor splits, floors are coarse, and pbcheck's own protocol treats floor-based comparisons below 8 donors per group as contaminated by the per-cell leak. Do not compare these numbers with a run on another file.

> Only 20 distinct donor splits exist for this design, so the permutation null has 20 draws, not the 1000 requested; the floor's Monte-Carlo standard error (0.412) is correspondingly large.

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
| lupus_nephritis | 3 | 3 |
| healthy | 4 | 4 |

**Cells per donor (top 7 of 7; 0 more not shown)**

| donor | cells |
| --- | --- |
| donor_00 | 200 |
| donor_01 | 190 |
| donor_02 | 180 |
| donor_03 | 170 |
| donor_04 | 160 |
| donor_05 | 150 |
| donor_06 | 140 |

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
| source | X |
| passed | True |
| reason | n/a |
| dtype | float32 |
| values scanned | 283654210 |

Examples of the values checked: 0.0, 1.0, 3.0, 12.0.

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
| Monte-Carlo SE | 0.412 |

The paired floor is not shown: the pseudobulk arm left 812 genes without a value, so the paired series is not comparable; the solo floor above stands alone.

**Permutations, naive arm**

| Field | Value |
| --- | --- |
| requested | 1000 |
| achieved | 20 |

> B5 machinery check, not a calibration criterion: empirical-permutation-p lambda 1.02.

## 7. Donor-pseudobulk arm

**Real-label counts**

| Field | Value |
| --- | --- |
| real-label calls, paired BH | 5 |

**Top genes, donor-pseudobulk arm**

| gene | pval | padj | log2fc |
| --- | --- | --- | --- |
| GOLGA8A | 0.0004 | 0.0061 | 1.75 |
| TXN2 | 0.0012 | 0.011 | 0.88 |

**Pseudobulk inflation and false-positive rate**

| Field | Value |
| --- | --- |
| lambda | 1.02 |
| lambda IQR | 0.09 |
| class | calibrated |
| false-positive rate | 0.041 |
| false-positive rate Monte-Carlo SE | 0.014 |

**Pseudobulk permutation floor**

| Field | Value |
| --- | --- |
| median count | 1 |
| median fraction | 0.00012 |
| bh mode | paired |
| Monte-Carlo SE | 0.09 |

**Permutations, pseudobulk arm**

| Field | Value |
| --- | --- |
| requested | 200 |
| achieved | 20 |

> B5 machinery check, not a calibration criterion: empirical-permutation-p lambda 1.01.

**Moderated eBayes technical detail**

| quantity | value | meaning |
| --- | --- | --- |
| d0 | 3.21 | moderated eBayes prior degrees of freedom |
| shrinkage factor | 0.62 | how strongly a gene's own variance is pulled toward the prior |
| complete pooling | no | whether every gene's variance was replaced by the prior outright |
| residual df | 19 | residual degrees of freedom feeding the moderated test |

## 8. Shared universe bookkeeping

**Frozen gene universe**

| Field | Value |
| --- | --- |
| size | 4310 |
| minimum size | 200 |
| builder | pseudobulk_frozen |

Builder rule: a gene is kept when its total count is at least 15 in at least half the pseudobulk profiles of each group, after the thin-donor filter

**Thin-donor filter**

| Field | Value |
| --- | --- |
| min_cells | 10 |
| min_counts | 1000 |
| donors_dropped | 1 |

**Shared universe, real-label bookkeeping**

| Field | Value |
| --- | --- |
| common tested | 8102 |
| dropped for fairness | 111 |
| no value in naive arm | 0 |
| no value in pseudobulk arm | 812 |
| pseudobulk arm complete for every gene | False |

**Real split against the permutation range**

| Field | Value |
| --- | --- |
| real split inside permutation range | yes |
| real split percentile in permutations | 90 |

**Monte-Carlo detail**

| Field | Value |
| --- | --- |
| floor_gap_over_mc_se | 4.8 |
| n_perm_naive | 1000 |
| n_perm_pb | 200 |

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
