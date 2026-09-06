# pbcheck audit report

## 1. Header

pbcheck audit of '(no file)'

pbcheck 0.1.0, generated 2026-09-06T07:37:08Z.

Status: complete

> This is a pbcheck v0.1.0 audit of one stratum of one file. It measures how a naive per-cell test and a donor-pseudobulk test behave on this data under a donor-permutation null. It is run outside pbcheck's pre-registered Phase 0 protocol, is not a Phase 0 measurement, and makes no claim about any publication or about whether any reported result is true.

## 2. Read-out

Under donor permutation, with no true signal to find, the naive per-cell test calls a median of 6915 genes (86.44% of the 8000-gene universe; Monte-Carlo SE 1.018) over 1000 permutations; on the real labels it calls 6892. Both counts are corrected over the whole universe on their own (solo BH).

The naive arm's inflation factor lambda is 144.7 (IQR 4.847): above the band 0.9 to 1.1, which is the donor-pseudobulk arm's band, shown here to describe the naive number and not as the naive arm's own criterion.

The donor-pseudobulk arm's lambda is 1 (inside the band 0.9 to 1.1); its permutation false-positive rate is 0.05 (MC SE 0.01541) and its floor a median of 0 genes; on the real labels it calls 0 genes, corrected across both arms together (paired BH).

8 donors in 'disease', 8 in 'ctrl'; 12868 distinct donor splits exist.

After the thin-donor filter (fewer than 10 cells or 1000 counts), 8 and 8 pseudobulk profiles remain.

Real-label calls over the permutation floor: the per-cell arm calls 0.9967 times its own floor, both counts corrected over the whole universe on their own (solo BH). The donor-pseudobulk arm calls 0 times its floor, both counts corrected across the two arms together (paired BH).

> On this file, with condition labels shuffled between donors and therefore no real signal to find, the per-cell test still calls a median of 6915 of 8000 genes at FDR 0.05 (86.4%), corrected over the whole gene universe on its own (solo BH); on the real labels, corrected the same way, it calls 6892. A gene list produced by a per-cell test on this data cannot be separated from that floor. The donor is the replication unit this design supports. The donor-pseudobulk test called 0 genes on the real labels against a permutation median of 0, both corrected across the two arms together over the genes they have in common (paired BH).

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
| ctrl | 8 | 8 |
| disease | 8 | 8 |

**Cells per donor (top 16 of 16; 0 more not shown)**

| donor | cells |
| --- | --- |
| d00 | 625 |
| d01 | 625 |
| d02 | 625 |
| d03 | 625 |
| d04 | 625 |
| d05 | 625 |
| d06 | 625 |
| d07 | 625 |
| d08 | 625 |
| d09 | 625 |
| d10 | 625 |
| d11 | 625 |
| d12 | 625 |
| d13 | 625 |
| d14 | 625 |
| d15 | 625 |

No batch columns were provided.

**Design flags**

| Field | Value |
| --- | --- |
| donor nests in condition | True |
| imbalance ratio | 1 |
| usable for pseudobulk | True |
| flags | none |

**Dropped cells**

| Field | Value |
| --- | --- |
| dropped, other condition | 0 |
| dropped, other cell type | 0 |
| dropped, missing donor | 0 |
| dropped, missing condition | 0 |
| dropped, missing cell type | 0 |

## 5. Counts check

**Counts check**

| Field | Value |
| --- | --- |
| source | X |
| passed | True |
| reason | n/a |
| dtype | int32 |
| values scanned | 80000000 |

Examples of the values checked: none.

## 6. Naive per-cell arm

**Real-label counts**

| Field | Value |
| --- | --- |
| real-label calls, solo BH | 6892 |
| real-label calls, paired BH | 6892 |

**Top genes, naive arm**

| gene | pval | padj | log2fc | pct group | pct reference |
| --- | --- | --- | --- | --- | --- |
| g0232 | 0 | 0 | 0.8765 | 0.9918 | 0.953 |
| g0500 | 0 | 0 | -0.9131 | 0.968 | 0.9938 |
| g0636 | 0 | 0 | 0.7821 | 0.9996 | 0.9932 |
| g0752 | 0 | 0 | -0.7551 | 0.9828 | 0.9954 |
| g1287 | 0 | 0 | -1.071 | 0.8268 | 0.948 |
| g1391 | 0 | 0 | 1.118 | 0.9988 | 0.9884 |
| g1629 | 0 | 0 | -1.189 | 0.8302 | 0.959 |
| g1976 | 0 | 0 | 0.9745 | 0.993 | 0.961 |
| g2398 | 0 | 0 | -1.049 | 0.888 | 0.9556 |
| g2753 | 0 | 0 | 0.7929 | 0.9998 | 0.9994 |
| g3053 | 0 | 0 | 1.116 | 0.997 | 0.9586 |
| g3452 | 0 | 0 | -0.9193 | 0.9392 | 0.9798 |
| g3643 | 0 | 0 | 0.7922 | 0.9994 | 0.9984 |
| g3779 | 0 | 0 | 1.3 | 0.8092 | 0.5564 |
| g3932 | 0 | 0 | 0.9441 | 0.9996 | 0.9962 |
| g4249 | 0 | 0 | 1.003 | 0.9986 | 0.9936 |
| g4458 | 0 | 0 | 1.436 | 0.8492 | 0.5708 |
| g4521 | 0 | 0 | -1.033 | 0.8462 | 0.9438 |
| g4581 | 0 | 0 | 1.163 | 0.9124 | 0.801 |
| g4820 | 0 | 0 | 0.964 | 0.9996 | 0.996 |
| g5247 | 0 | 0 | 1.019 | 0.96 | 0.8844 |
| g5317 | 0 | 0 | 1.082 | 0.9436 | 0.799 |
| g5405 | 0 | 0 | -1.101 | 0.9508 | 0.9944 |
| g5844 | 0 | 0 | -0.8106 | 0.9986 | 1 |
| g6185 | 0 | 0 | 0.9639 | 0.9952 | 0.9742 |

**Naive inflation factor**

| Field | Value |
| --- | --- |
| lambda | 144.7 |
| lambda IQR | 4.847 |
| class | above_band |

**Solo permutation floor**

| Field | Value |
| --- | --- |
| median count | 6915 |
| median fraction | 0.8644 |
| IQR count | 44 |
| bh mode | solo |
| Monte-Carlo SE | 1.018 |

**Paired permutation floor**

| Field | Value |
| --- | --- |
| median count | 6912 |
| median fraction | 0.864 |
| IQR count | 37.5 |
| bh mode | paired |
| Monte-Carlo SE | 2.3 |

**Permutations, naive arm**

| Field | Value |
| --- | --- |
| requested | 1000 |
| achieved | 1000 |

> Machinery check of the permutation engine, not a criterion of any kind: the inflation factor of the empirical permutation p-values is 0.9884.

## 7. Donor-pseudobulk arm

**Real-label counts**

| Field | Value |
| --- | --- |
| real-label calls, paired BH | 0 |

**Top genes, donor-pseudobulk arm**

| gene | pval | padj | log2fc |
| --- | --- | --- | --- |
| g4458 | 9.85e-05 | 0.6984 | 1.373 |
| g6837 | 0.0004187 | 0.6984 | -1.254 |
| g3878 | 0.0006038 | 0.6984 | -1.218 |
| g6862 | 0.0006301 | 0.6984 | 1.237 |
| g3779 | 0.0006696 | 0.6984 | 1.199 |
| g2522 | 0.0006998 | 0.6984 | 1.219 |
| g1629 | 0.0007639 | 0.6984 | -1.177 |
| g6632 | 0.000964 | 0.6984 | 1.184 |
| g7623 | 0.0009681 | 0.6984 | -1.171 |
| g3053 | 0.001182 | 0.6984 | 1.15 |
| g5003 | 0.001205 | 0.6984 | -1.153 |
| g5405 | 0.001209 | 0.6984 | -1.142 |
| g3343 | 0.001257 | 0.6984 | -1.133 |
| g1287 | 0.001361 | 0.6984 | -1.141 |
| g5317 | 0.001455 | 0.6984 | 1.14 |
| g1277 | 0.00147 | 0.6984 | 1.131 |
| g1118 | 0.001484 | 0.6984 | -1.142 |
| g0485 | 0.001598 | 0.71 | 1.121 |
| g1391 | 0.001701 | 0.7163 | 1.125 |
| g3425 | 0.0022 | 0.8801 | 1.094 |
| g1717 | 0.002339 | 0.891 | 1.08 |
| g6102 | 0.002559 | 0.9302 | 1.077 |
| g3973 | 0.002682 | 0.9302 | 1.067 |
| g6808 | 0.00279 | 0.9302 | 1.057 |
| g4581 | 0.003256 | 0.9951 | 1.042 |

**Pseudobulk inflation and false-positive rate**

| Field | Value |
| --- | --- |
| lambda | 1 |
| lambda IQR | 0.03131 |
| class | in_band |
| false-positive rate | 0.05 |
| false-positive rate Monte-Carlo SE | 0.01541 |

**Pseudobulk permutation floor**

| Field | Value |
| --- | --- |
| median count | 0 |
| median fraction | 0 |
| bh mode | paired |
| Monte-Carlo SE | 0.01765 |

**Permutations, pseudobulk arm**

| Field | Value |
| --- | --- |
| requested | 200 |
| achieved | 200 |

> Machinery check of the permutation engine, not a criterion of any kind: the inflation factor of the empirical permutation p-values is 0.9884.

**Moderated eBayes technical detail**

| quantity | value | meaning |
| --- | --- | --- |
| d0 | 320.8 | moderated eBayes prior degrees of freedom |
| shrinkage factor | n/a | how strongly a gene's own variance is pulled toward the prior |
| complete pooling | n/a | whether every gene's variance was replaced by the prior outright |
| residual df | n/a | residual degrees of freedom feeding the moderated test |

## 8. Shared universe bookkeeping

**Frozen gene universe**

| Field | Value |
| --- | --- |
| size | 8000 |
| minimum size | 200 |
| builder | pseudobulk_frozen |

Builder rule: A gene enters the frozen universe when its total count across the surviving donor pseudobulk profiles is at least 15 and it is detected (count greater than 0) in at least ceil(0.5 * n_profiles) of them. The rule ignores the condition labels, so the tested gene set cannot shift when the labels are permuted (pbcheck.gene_universe.frozen_universe, applied after the thin-donor filter).

**Thin-donor filter**

| Field | Value |
| --- | --- |
| min_cells | 10 |
| min_counts | 1000 |
| n_profiles_before | 16 |
| n_profiles_after | 16 |
| n_dropped | 0 |
| dropped_profiles | [] |
| per_threshold | {'min_cells': 0, 'min_counts': 0} |
| semantics | dropped, not merged (spec section 1 item 2; Amendment 2 Change 7) |

**Shared universe, real-label bookkeeping**

| Field | Value |
| --- | --- |
| common tested | 8000 |
| dropped for fairness | 0 |
| no value in naive arm | 0 |
| no value in pseudobulk arm | 0 |
| pseudobulk arm complete for every gene | True |

**Real split against the permutation range**

| Field | Value |
| --- | --- |
| real split inside permutation range | yes |
| real split percentile in permutations | 1 |

**Monte-Carlo detail**

| Field | Value |
| --- | --- |
| naive_floor_median | 6912 |
| naive_floor_mc_se | 2.3 |
| naive_floor_median_solo | 6915 |
| pb_floor_median | 0 |
| pb_floor_mc_se | 0.01765 |
| pb_fp_rate | 0.05 |
| pb_fp_rate_mc_se | 0.01541 |
| bh_mode | paired (mtc.bh_both_arms) over one common tested set — spec section 5 |
| naive_floor_median_solo_bh_mode | solo (mtc.bh_over_universe) — naive arm only |
| floor_gap | 6912 |
| floor_gap_over_mc_se | 3005 |

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
| prepare | 0.2831 |
| design | 0.01699 |
| counts | 0.2682 |
| pseudobulk_build | 2.384 |
| universe | 0.02562 |
| real_label | 23.35 |
| permutation_null | 62.79 |
| render | n/a |

> The engine was measured on one synthetic oracle point (sigma_donor 0.5, 8 against 8 donors, 1500 genes), recorded in pilot/gate/synthetic_gate_2026-08-15.json; that measurement is not repeated on this file. Of the settings below, only the ones listed as protocol constants are pre-registered values taken from pbcheck.gate_config (alpha, lambda_band, min_universe_size, min_cells, min_counts); the permutation counts, the universe filter parameters, the fallback universe rule and the display settings are the tool's own and are not protocol values.

## 10. Footer

Generated by pbcheck 0.1.0. Protocol: docs/PHASE0_SPEC.md and docs/AMENDMENTS.md in the pbcheck repository.
