# pbcheck — Phase 0 Pre-registration of the Stratum List (2026-08-16)

Spec [`PHASE0_SPEC.md`](PHASE0_SPEC.md) §1 closes with a sentence that is an instruction and not a
description:

> **Pre-register the stratum list before computing any metric.**

This document is that act. It is **not** an amendment: it changes no threshold, supersedes no
section and relaxes nothing. It executes a step the frozen spec demands and had left undone, and it
closes one numerical gap in §1 that the spec references three times and never states (the
cells-per-donor bins, §5 below).

Its machine-readable half is [`../pilot/preregistration/stratum_list_2026-08-16.json`](../pilot/preregistration/stratum_list_2026-08-16.json)
(and its CSV twin), emitted by [`../scripts/freeze_stratum_list.py`](../scripts/freeze_stratum_list.py)
from the committed candidate manifest. Neither half is authoritative alone: the script re-derives
the load-bearing figures below — the frozen set, the per-dataset counts and ceilings, the bin
occupancy, the Layer B subsets — and refuses to run if any of them disagrees with what is written
here. §10 lists what else was re-derived by hand while writing, and what was not.

---

## 1. What this act is, and what it binds

| | |
|---|---|
| Date of the freeze | **2026-08-16** |
| Census | **`2025-01-30`**, the §1 pin — `open_census()` rejects `latest` and `stable` by name |
| Independent datasets | **12** |
| Frozen strata (stratum-contrasts) in the analysis set | **251** |
| Named siblings, frozen as within-collection controls and excluded from the D2 denominator | **2** datasets, **27** strata |
| Metrics computed on any of these strata at the time of writing | **none** |
| Rows admitted to the sweep | **0** — `admitted_to_sweep = False` on all 2190 rows of the source manifest, and on all 251 here |

**What it binds.** From this commit forward, the analysis set of the Phase 0 real-data sweep is the
251 strata named here, and the "independent datasets" denominator of decision rule item 2 (D2) is
the 12 dataset ids named here. Any addition, removal, re-selection or re-interpretation — including
dropping a dataset because its numbers turn out unhelpful, and including adding one because the
survivors look thin — requires a dated, numbered entry in [`AMENDMENTS.md`](AMENDMENTS.md), written
**before** the change is applied, in the form the log has used three times already.

**What it deliberately does not bind.** It does not admit anything to the sweep (§9 below), it does
not claim that any of these strata are inside the operating envelope, and it does not predict any
result. It fixes *what will be measured*, so that *what is found* cannot retroactively decide it.

**Why now, and why this is worth doing at all.** [Amendment 2](AMENDMENTS.md) had to open with the
disclosure that its selection data "has been seen in full, before the selection rule was written
down… the opposite of the pre-registration ideal". That disclosure was honest and it was also
expensive: the defence of that choice cost a page of argument and still leaves a reader entitled to
distrust it. This freeze is the cheap version of the same guarantee — the list is fixed while the
metrics are, verifiably, not yet computed — and the cost of getting it wrong is that Phase 0's
headline becomes unfalsifiable.

---

## 2. Provenance

Everything below is derived from one artifact, committed here in full rather than referenced:

| | |
|---|---|
| File | [`../pilot/preregistration/census_candidates_run31910799023_2026-08-15.json`](../pilot/preregistration/census_candidates_run31910799023_2026-08-15.json) |
| sha256 | `33f8a800229dccc5f58f311e7d0c493655068d43563b31ff53fdaebb3b44e4b4` |
| Size | 6 630 446 bytes |
| CSV twin | [`../pilot/preregistration/census_candidates_run31910799023_2026-08-15.csv`](../pilot/preregistration/census_candidates_run31910799023_2026-08-15.csv) |
| CSV sha256 | `09eb110dd308155f64e10b2b05beff36854f7125b1434699935707d3551f12d6` |
| CSV size | 4 513 660 bytes |
| Header `generated_utc` | `2026-08-15T22:18:37+00:00` |
| Header `census_version` | `2025-01-30` |
| Produced by | `scripts/census_candidates.py` via `.github/workflows/census-candidates.yml`, **GitHub Actions run `31910799023`** (`workflow_dispatch`, `dry_run: false`) |
| Shape | 2190 stratum-contrasts over 73 datasets; 1197 `candidate`, 981 `excluded_inclusion_gate`, 12 `excluded_confound` |

The manifest is normally a **CI artifact and never a commit** — `pilot/results/` is gitignored
precisely so that "a candidate list that reached git by way of a CI job would have pre-registered
itself by accident". That reasoning is why it is committed *here*, under `pilot/preregistration/`,
by a deliberate human act with this document attached: a pre-registration whose evidence lives in an
expiring CI artifact is not auditable, and an evidence file whose bytes cannot be checked is not
evidence. `.gitattributes` marks the directory `-text` so no end-of-line conversion can move the
hashes.

**The proposal document is committed too, and its status is narrower than the manifest's.**

| | |
|---|---|
| File | [`../pilot/preregistration/stratum_list_proposal_2026-08-16.md`](../pilot/preregistration/stratum_list_proposal_2026-08-16.md) |
| sha256 | `50872414b0727c129a824b0c65ed179674ac5d6c9ecaac53327568b3eae6fb48` |
| Size | 92 589 bytes |
| Status | **the reasoning behind the choice of twelve, committed for auditability — not part of the binding act** |

The binding content of this pre-registration is the list of §3.2 and the rule of §3.1, and nothing
else. The proposal is the working document those were chosen from: it carries the per-dataset
rationale, the datasets considered and rejected, the **five named reserves**, the third 5′ candidate
and the method of the Mathys search. It is referenced by §3.2, §4.2 and §8, and until this commit it
existed only in a scratch directory — so the one act of discretion in the whole freeze was justified
by a file no reader could open. The reserves are the sharpest case: §9.7 forbids substituting a
replacement for a stratum that fails the counts gate, and a list of five candidate replacements must
not live outside the record where nobody can check whether one was quietly used.

It is committed **byte-for-byte as it was circulated**, including its Russian prose and the two
figures this document corrects (§3.2's Rexach ceiling, and §6's claim about what the manifest could
support). Translating or patching it would have produced a hash of a rendering rather than of the
thing that was acted on, which is the same mistake as not committing it. Where it disagrees with
this document, **this document governs** and the disagreement is recorded in §3.2, §6 and §10.

`scripts/freeze_stratum_list.py` checks the size and then the sha256 **before it parses anything**,
and then the header's `generated_utc`, `census_version` and row count, aborting on the first
mismatch. The hash guards against the file changing; the header stamps guard against the pinned
constant being edited to match some other file. Either check alone can be walked around by one edit.

**Two things are not derivable from this artifact and are labelled wherever they appear.** Assay,
suspension type, tissue and publication DOI are **not manifest columns** — `census_select` screens
assay and suspension for confounding but does not emit their levels — so §4's coverage table takes
them from the CELLxGENE Discover curation API index (`GET /curation/v1/datasets`, 2216 datasets,
read 2026-08-16). That index is not committed; every value is checkable against the public Discover
record for the dataset id. And the strong/subtle axis of §1 (i)/(ii) is a literature judgement, not
a measurement, marked as such in §4.

**Regenerating.** `python scripts/freeze_stratum_list.py` rewrites the two artifacts;
`--check` verifies **both** committed halves — JSON and CSV — against a fresh run without writing.
The output carries no generation timestamp, no package versions and no environment, so regeneration
is byte-identical on any platform, and `tests/test_stratum_list_freeze.py` enforces that against the
committed bytes.

---

## 3. The selection rule, and the frozen set

### 3.1 The rule

> **Every row of the source manifest whose `gate_status == "candidate"` and whose `dataset_id` is
> one of the twelve datasets in §3.2 is in the analysis set; every other row is out.**

That is the whole rule, and it is deterministic. The twelve dataset ids are the only judgement
exercised anywhere in this freeze. Given them, the 251 strata follow by arithmetic — no stratum,
cell type, disease term or donor-count tier is chosen by hand, and no stratum of a listed dataset is
dropped for being small, skewed, noisy or awkward. `select_strata()` in the freeze script *is* the
rule; there is no second place where membership is decided.

The reason to bind ourselves this tightly is specific and is named in the spec's own risk list
(§10 risk 13, *stratum cherry-picking*). Stratum-level discretion is the cheapest way to manufacture
a GO: the manifest offers 1197 candidates across 68 datasets, cells-per-donor spans three orders of
magnitude, and λ_naive is expected to grow with cells per donor (§10 risk 2). A rule that admitted
per-stratum choice would let the study pick its own answer while every individual choice looked
defensible.

`gate_status == "candidate"` is `census_select`'s own verdict: the stratum cleared the obs-decidable
half of the §1 inclusion gate (≥ 3 donors per group after the thin-donor drop, donor present,
non-constant and nested within condition) and was not excluded by the §1 confound pre-screen. **It
is not an admission.** See §9.

**The "recommended strata" tables of the 2026-08-16 proposal document are not this set.** That
document listed four to six strata per dataset as reading aids for a human deciding which datasets
to take. They are a strict subset of the 251, they were never a selection, and they carry no status
of any kind. Stated plainly here because the failure mode is obvious: a later reader finds the
smaller, tidier table, mistakes it for the pre-registration, and reports a study whose stratum list
was in fact chosen after the datasets were.

### 3.2 The twelve datasets

Ceiling = the largest `min(n_donors_A, n_donors_B)` over the dataset's candidate strata, i.e. the
best-powered design it contains. `≥ 8v8` and `3v3` count strata, not datasets. Every figure in this
table is computed from the manifest by the freeze script; the two it also declares in advance —
strata and ceiling — abort the freeze on a disagreement rather than being quietly re-derived.

| # | `dataset_id` | Short | Strata | Ceiling | ≥ 8v8 | exactly 3v3 | Disease terms vs `normal` |
|---|---|---|---|---|---|---|---|
| 1 | `6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3` | Gabitto 2024, SEA-AD DLPFC | 18 | **39** | 18 | 0 | dementia |
| 2 | `ac0c6561-7a48-4185-af6f-af799f699172` | Rexach 2024 Cell, cross-dementia | 27 | **10** | 23 | 2 | Alzheimer disease; Pick disease; progressive supranuclear palsy |
| 3 | `d8da613f-e681-4c69-b463-e94f5e66847f` | Melms 2021 Nature, lethal COVID-19 lung | 28 | 7 | 0 | 0 | COVID-19 |
| 4 | `2a498ace-872a-4935-984b-1afa70fd9886` | Yoshida 2022 Nature, PBMC | 47 | 20 | 34 | 0 | COVID-19; post-COVID-19 disorder |
| 5 | `ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3` | Ahern 2022 Cell, COMBAT blood atlas | 25 | 10 | 21 | 0 | COVID-19; influenza |
| 6 | `f1606894-59df-4794-a37f-baa7c6fb6de1` | Linna-Kuosmanen 2024, PERIHEART right atrium | 11 | 25 | 11 | 0 | atrial fibrillation |
| 7 | `d18736c3-6292-4379-919a-d6d973204c87` | Binvignat 2024, rheumatoid arthritis blood | 15 | 18 | 14 | 0 | rheumatoid arthritis |
| 8 | `a12ccb9b-4fbe-457d-8590-ac78053259ef` | KPMP adult human kidney snRNA-seq v1.5 | 37 | 24 | 29 | 0 | acute kidney failure; chronic kidney disease |
| 9 | `19e46756-9100-4e01-8b0e-23b557558a4c` | Heimlich 2024, clonal haematopoiesis PBMC | 7 | 7 | 0 | 0 | clonal hematopoiesis |
| 10 | `c893ddc3-f25b-45e2-8c9e-155918b4261c` | Phan 2024, opioid use disorder striatum | 10 | 6 | 0 | 0 | opiate dependence |
| 11 | `8e47ed12-c658-4252-b126-381df8d52a3d` | Elmentaite 2020, paediatric gut (Crohn) | 18 | 7 | 0 | 0 | Crohn disease |
| 12 | `4b6af54a-4a21-46e0-bc8d-673c0561a836` | Wang 2023, emphysema non-immune | 8 | **3** | 0 | 8 | pulmonary emphysema |
| | | **Total** | **251** | | **150** | **10** | 15 distinct terms |

Aggregates over the 251, all re-derived: 4 609 595 cells; 182 distinct `(dataset_id × cell_type)`
strata carrying 251 binary contrasts (a dataset with two disease terms contributes two contrasts
against the same `normal` group); 124 distinct cell-type labels; `residual_df` 4 … 108;
`permutation_count` 20 … 7.28 × 10²³, of which **45 strata fall below 1000 and therefore require
full enumeration** (§4, "Small D → enumerate all").

**A correction, recorded rather than absorbed.** The proposal document circulated on 2026-08-16
states Rexach's envelope ceiling as `min(A,B) = 11`. **It is 10.** The manifest's best-powered Rexach
strata are A = 11 versus B = 10 (the progressive supranuclear palsy arm, six cell types), and no
Rexach control group anywhere in the dataset exceeds 10 donors. The 11 is `max(n_donors_A)` read as
if it were the ceiling; the two coincide for balanced designs and diverge for every skewed one. The
pre-registration carries **10**. Nothing downstream changes — Rexach was already outside the σ = 0.5
tier at either value — but a pre-registration that silently improved a number it had published
would be worth less than one that did not.

**Where the 12 came from.** The 68 candidate-bearing datasets were read for coverage of the axes
§1 (iii) names, not for expected inflation. The rationale per dataset, the datasets considered and
rejected, and five named reserves are in the 2026-08-16 proposal document, committed at
[`../pilot/preregistration/stratum_list_proposal_2026-08-16.md`](../pilot/preregistration/stratum_list_proposal_2026-08-16.md)
and hash-pinned in §2. It is the reasoning, not the act: what binds is the list above and the rule
in §3.1, and **none of its five reserves may be substituted for anything** (§9.7). What matters for
auditing is stated in §4 — the coverage claim is checkable against the frozen set itself.

### 3.3 The two siblings — frozen as controls, and why they are not datasets 13 and 14

| `dataset_id` | What it is | Sibling of | Frozen strata | Role |
|---|---|---|---|---|
| `c2876b1b-06d8-4d96-a56b-5304f815b99a` | SEA-AD, middle temporal gyrus (`min(A,B)` 27 … 42) | #1 | 18 | `within_collection_control` |
| `1e5bd3b8-6a0e-4959-8d69-cafed30fe814` | Emphysema Cell Atlas, immune cells (all 3v3) | #12 | 9 | `within_collection_control` |

D2 clusters evidence by dataset because same-dataset strata share donors, batch and assay. Two
datasets from one collection share the **cohort and the laboratory** as well: SEA-AD MTG is the same
cohort and laboratory as SEA-AD DLPFC in a different cortical region, and the emphysema immune split
is the same three-versus-three cohort's other compartment. Counting either toward "majority of
independent datasets" would inflate the effective n exactly the way D2 exists to prevent.

*Not* "the same donors": the analysed data does not support that. DLPFC's best design is 39 v 44 and
MTG's is 42 v 46, and MTG's `min(A,B)` runs 27 … 42 against DLPFC's 33 … 39 — the cohorts overlap,
the donor sets in the strata do not coincide, and the artifact's own wording (*same cohort and
laboratory*) is the one that is true.

**They are frozen too, by the same rule, and emitted with an explicit role.** Every candidate row of
theirs — 18 + 9 = **27 strata** — is in
[`stratum_list_2026-08-16.json`](../pilot/preregistration/stratum_list_2026-08-16.json) under
`within_collection_control_rows`, carrying `role = "within_collection_control"`, where the 251 carry
`role = "analysis_set"`. Naming them without freezing them would have left 27 runnable, unlisted
strata outside the pre-registration — a set someone could reach for after seeing the results, which
is the freedom §3.1 exists to remove. The CSV twin carries both blocks, 278 data rows, told apart by
that same column.

Three rules attach to the control set and none of them is discretionary:

1. A result from a control stratum is **reported as a within-collection control**, labelled as such
   wherever it appears.
2. It **never enters the D2 denominator** and never counts toward a majority, whatever it shows.
3. **Promoting one to an independent dataset is an amendment** — dated, numbered, written before the
   change. An unnamed sibling can be promoted later by someone who did not know it was one; a named
   and frozen one cannot be promoted quietly at all.

They remain admissible as within-collection reproducibility controls, which is a useful thing to
have: a result that fails to reproduce between two regions of one cohort is informative, and it is
informative in a way that says nothing about independent replication.

---

## 4. Coverage against §1 (iii), re-derived

§1: *"First pass = 8–12 datasets chosen to SPAN the outcome space (not cherry-pick wins): (i) 2–3
with a biologically strong expected effect (pseudobulk shown non-null), (ii) 2–3 subtle/low-effect,
(iii) deliberate variation in assay (10x 3′ vs 5′), tissue, donor count (some exactly 3v3, some
≥ 8v8), and cells-per-donor spanning the pre-registered bins (D1)."*

### 4.1 (i) and (ii) — the effect-size axis is a literature judgement, and is labelled as one

No effect in this table has been measured by us and none can be before the sweep runs. The axis
exists for one purpose: to show that the list is not picked for wins. It is **not** a prediction of
λ_naive and must never be scored against the results as though it were a hypothesis. It is recorded
per `(dataset, disease term)` rather than per dataset, because four datasets carry more than one
disease arm against **one shared control group**, and in two of them (#4, #5) those arms differ in
expected strength — the cheapest effect-size control available to us, and one a per-dataset label
destroys. The labels are frozen in the machine-readable artifact (`expected_effect`, one value per
row) so they cannot be reassigned after the fact.

**§1 (i)'s parenthetical is not satisfied by anything checked here, and that is a gap, not a
formality.** §1 asks for datasets "with a biologically strong expected effect **(pseudobulk shown
non-null)**". The Basis column below gives biological expectation — what kind of disease this is,
in what tissue, against what control — and **not** evidence that donor-level pseudobulk was shown
non-null in those publications. No such check was performed for this freeze: the literature was read
for study design, not for the presence of a donor-aggregated differential-expression result. So the
*(pseudobulk shown non-null)* clause is **unmet**, and the strong / moderate / subtle labels rest on
biological expectation alone. It is recorded here rather than repaired, because repairing it means a
per-publication evidence review that this document did not do and must not pretend to have done.
The consequence is bounded and worth stating: the labels are a *coverage* claim about the list's
span, and a coverage claim that rests on expectation is weaker than one that rests on published
pseudobulk results — but it is not a claim about any result of ours, so nothing downstream inherits
an unearned number.

| Expected | Datasets | Arms | Basis — biological expectation from each dataset's own publication, **not** a pseudobulk result |
|---|---|---|---|
| **strong** | 5 (#1, #2, #3, #4, #5) | dementia; AD, Pick, PSP; lethal COVID-19 lung; acute COVID-19 PBMC ×2 | Cortical neurodegeneration and tauopathy with glial response (Gabitto 2024 `10.1038/s41593-024-01774-5`; Rexach 2024 `10.1016/j.cell.2024.08.019`); diffuse alveolar damage in autopsy lung (Melms 2021 `10.1038/s41586-021-03569-1`); acute severe viral infection in blood (Yoshida 2022 `10.1038/s41586-021-04345-x`; Ahern 2022 `10.1016/j.cell.2022.01.012`) |
| **moderate** | 3 (#5, #8, #12) | influenza; acute kidney failure; chronic kidney disease; pulmonary emphysema | Milder viral arm of the same atlas; injured proximal-tubule states, acute then chronic (KPMP v1.5, **DOI not established** — Discover returns `collection_doi = null`); pulmonary emphysema against an alveolar control (Wang 2023 `10.1016/j.immuni.2023.01.032`) — but see the note on #12 below, whose label is the weakest thing in this table |
| **subtle** | 6 (#4, #6, #7, #9, #10, #11) | post-COVID-19 disorder; atrial fibrillation; rheumatoid arthritis; clonal haematopoiesis; opiate dependence; Crohn disease | Clinically heterogeneous post-viral state; atrial remodelling against a cardiac-surgery control that is not a healthy heart (Linna-Kuosmanen 2024 `10.1016/j.xcrm.2024.101556`); subpopulation-level blood signatures (Binvignat 2024 `10.1172/jci.insight.178499`); a pre-clinical clonal state with no phenotype (Heimlich 2024 `10.1182/bloodadvances.2023011445`); cell-type-specific programmes rather than degeneration (Phan 2024 `10.1038/s41467-024-45165-7`); paediatric ileal mucosa (Elmentaite 2020 `10.1016/j.devcel.2020.11.010`) |

**"moderate" is ours, not §1's.** §1 names two categories — (i) strong and (ii) subtle/low-effect —
and a three-valued vocabulary is a labelling convenience of this freeze, not a spec category. It is
frozen as `expected_effect_vocabulary` in the artifact so it cannot grow later. What §1 actually
demands is counted over its own two categories, and the middle tier is counted toward **neither**:
5 datasets carry a strong arm and 6 carry a subtle arm, both above §1's 2–3, without any moderate
arm being borrowed to reach the floor.

**#12's row is a design justification wearing an effect label.** Wang 2023 is in the list to satisfy
§1 (iii)'s "some exactly 3v3" — all 8 of its strata are 3 v 3 — and its inclusion was decided on
that ground, not on an expected effect size. Its `moderate` label is therefore the one entry in the
table with the least behind it, and the honest reading is that #12 spans a *donor-count* axis and
was given the middle label because the vocabulary requires one, not because the emphysema literature
places it there. It is called out rather than quietly assigned so that no later reader treats #12's
result as bearing on the effect-size axis at all.

§1 asks for 2–3 of each; the list carries 5 and 6. That is a floor being exceeded, not a bar being
moved, and the excess is not free padding: it comes from the within-dataset arms (#4 strong + subtle,
#5 strong + moderate, #8 acute + chronic, #2 three tauopathies), which cost no additional dataset
and hold assay, tissue, laboratory and control group fixed across the effect-size contrast.

One consequence is worth pre-registering as a stated expectation rather than discovered later:
**#9 (clonal haematopoiesis) is expected to behave as a no-effect stratum set.** §4/A1 says the
clean truth-zero guarantee "lives in oracle (b) and in no-effect strata"; if that role is claimed
after seeing the result it is worthless, so it is claimed now.

### 4.2 (iii) — assay, tissue, suspension, donor count, cells per donor

Assay, suspension and tissue are from the Discover index (§2), not from the manifest.

| Axis | §1 requirement | Frozen set | Verdict |
|---|---|---|---|
| Assay 3′ | variation | 10 datasets: #1, #2, #3, #6, #7, #8, #9, #10, #11, #12 | met |
| Assay 5′ | variation | **2 datasets: #4 and #5, both 10x 5′ v1** | met, and thinly — see below |
| Assay chemistry within 3′ | not required | v3 in 9 datasets; v2 in #11 (pure) and #2 (mixed v2 + v3); Discover additionally lists `10x multiome` for #1 — see below | — |
| Tissue | variation | **6 organ systems**: brain (#1, #2, #10), lung (#3, #12), blood (#4, #5, #7, #9), heart (#6), kidney (#8), gut (#11) | met |
| Suspension | not named in §1; D3-adjacent | `nucleus` 6 (#1, #2, #3, #6, #8, #10) / `cell` 6 (#4, #5, #7, #9, #11, #12) | balanced |
| Donor count — "some exactly 3v3" | required | **10 strata**, in #12 (all 8) and #2 (2 T-cell strata) | met |
| Donor count — "some ≥ 8v8" | required | **150 strata** across 7 datasets (#1, #2, #4, #5, #6, #7, #8) | met |
| Donor count — range | — | `min(A,B)` from 3 to 39; 24 strata at `min(A,B) = 3`, of which 10 are exactly 3v3; largest design 39 v 44 | — |
| Cells per donor | "spanning the pre-registered bins" | group medians **11.0 … 6671.5**, all six bins occupied — §5 | met, against bins defined in §5 |
| Counts per cell | not in §1; D1-adjacent | group medians **284 … 56 841.5**, over two orders of magnitude; the extremes are #7 (RA blood) and #10 (striatal medium spiny neurons), both 10x 3′ v3 | — |

**The 5′ arm is the thinnest axis in the list, and it is thin in a way that matters.** Both 5′
datasets are PBMC, both are 10x 5′ v1, and both carry a COVID-19 arm. If the 3′-versus-5′ comparison
produces a difference, this list cannot separate "5′ chemistry" from "blood", "PBMC" or "acute viral
infection". Widening it was possible — the proposal names a third 5′ candidate — and was not done,
because that candidate uses *in-vitro-stimulated* PBMCs, which inserts a layer between disease and
transcriptome that seemed worse than the confound it would relieve. That trade is recorded here so
it can be disagreed with, and it is a limitation of the frozen list, not of the data.

**A second limitation, of the same kind.** The 3′/5′ split coincides almost exactly with the
suspension split among the donor-rich datasets, and blood is the only tissue represented by both
5′ datasets. Any assay-level statement from this list is therefore partially confounded with tissue
by construction. §7 (depth-matching, permutation null) does not fix this; only more datasets would,
and there are not more inside the 12-dataset ceiling §1 sets.

**The multiome entry is a Discover property, not a property of the frozen strata, and the two
disagree.** Discover lists `10x multiome` alongside `10x 3' v3` for SEA-AD DLPFC at the *dataset*
level. The manifest says something different about the *strata*: `confound_cramers_v["assay"]` is
**`null` on all 18 SEA-AD strata**, which is what `census_select` writes when assay is constant
within a stratum — there is no association to measure. Rexach (#2), whose declared v2 + v3 mix is
real, shows a **non-zero V on 24 of its 27 strata** (0.029 … 0.085), so the column does register a
within-stratum assay mix when there is one. Two readings survive: either the multiome cells are
absent from the 18 frozen strata (dropped upstream, or confined to cell types that did not clear the
gate), or they are present and the §1 confound pre-screen missed them. **`obs` cannot decide
between them**, so nothing here is asserted about the frozen set; the multiome value is recorded as
Discover's dataset-level metadata and is carried as such in the artifact's `assay` field. §9.7 lists
this among the things only an X load can settle.

---

## 5. The cells-per-donor bins (D1), pre-registered here

`PHASE0_SPEC.md` references "the pre-registered bins" three times — decision rule item 2 ("as
cells-per-donor is subsampled down toward the floor of the pre-registered bins"), §1 (iii)
("cells-per-donor spanning the pre-registered bins (D1)") and §7 item 3 ("pre-register the
cells-per-donor bins") — and **nowhere states them numerically**. Neither does any amendment.
Strictly, §1 (iii) is currently unsatisfiable: there is no object to span.

**This is a gap in §1 and it is closed here, in the same act that freezes the list. It is not a
change to §1**: no threshold moves and nothing §1 defines is redefined. A number the spec demanded
and never supplied is supplied.

It is **not** true, however, that nothing downstream is affected, and the earlier draft of this
paragraph claimed it. Decision rule item 2 requires the naive effect to persist "as cells-per-donor
is subsampled down toward **the floor of the pre-registered bins**", so the floor is an argument to
the GO/NO-GO rule and setting it at 10 parameterises that rule. The direction is the unhelpful one:
a **higher** floor would have made item 2 easier, because the effect would only have to persist down
to a shallower subsampling depth, where λ_naive is expected to be larger (§10 risk 2). Choosing 10 —
the lowest value §1's own inclusion gate permits — is therefore the strictest available choice and
not a neutral one, and it is fixed here before any curve is computed. Anyone who thinks the floor
should be higher is asking for a weaker item 2 and must say so in an amendment.

> **The pre-registered cells-per-donor bins are the half-decade log bins**
> `[10, 30) [30, 100) [100, 300) [300, 1000) [1000, 3000) [3000, ∞)`.

The lower edge is 10 because that is the §1 inclusion gate's own floor (`gate_config.MIN_CELLS`;
item 2 drops donors below it), so no admitted donor can fall outside the bins from below. Each step
is a factor of about √10, giving five closed bins per two decades. The top bin is open because the
Census's upper tail is unbounded and a closed top bin would be a threshold with nothing behind it.

**The construction is deliberately independent of our data**, so that no boundary can be suspected
of having been fitted to make an occupancy table look full: half-decade log spacing anchored at a
threshold the spec already fixed is decidable without opening the manifest, and it was. The
occupancy below is a *consequence*, reported after the fact and not used to choose the edges. If the
result had been three empty bins, that would have been the reported finding.

The binned quantity is the **per-group median cells per donor** of a stratum
(`cells_per_donor_by_group[A|B]["median"]`), so each of the 251 strata contributes two values.

| Bin | Group medians | | Bin | Group medians |
|---|---|---|---|---|
| `[10, 30)` | 85 | | `[300, 1000)` | 106 |
| `[30, 100)` | 161 | | `[1000, 3000)` | 38 |
| `[100, 300)` | 99 | | `[3000, ∞)` | 13 |

**All six bins are occupied** over 502 group medians spanning **11.0 … 6671.5** — a factor of 600.
Every one of the twelve datasets occupies at least four bins on its own, and three — #1 (SEA-AD
DLPFC), #2 (Rexach) and #5 (COMBAT) — occupy **all six**, which means the D1 curve can be traced
within a single dataset at fixed assay, tissue, cohort and laboratory before it is ever traced
across datasets. That is the form in which D1 is least confounded, and it is available here by
construction rather than by luck.

Two honest qualifications. The distribution is not uniform — the two extreme bins hold 85 and 13
values against 161 in the mode — so cross-dataset comparison at matched cells-per-donor will be
better resolved in the middle of the range than at either end. And the *median* is what is binned:
within a stratum the per-donor spread is wide (the smallest per-donor cell count anywhere in the
frozen set is 10, the largest 16 383), so a bin label describes the stratum's centre and not its
donors.

---

## 6. Layer B — the pre-declared truncation

[Amendment 3](AMENDMENTS.md) declared an operating envelope: the pseudobulk arm is valid only where
`(sigma_donor, donors-per-group)` puts the selected `ebayes` test at power ≥ 0.60 on the unchanged
pre-registered oracle. It also stated, three amendments running, that `sigma_donor` **is not
anchored to any real data** and that "whether any real stratum falls inside it is unknown".

When the anchor lands it will decide, mechanically, which of the 12 datasets survive. **That subset
is declared now, before the anchor exists**, so that it cannot be chosen later from among the
subsets that happen to be convenient. A dataset is in a tier if and only if it holds at least one
frozen stratum with at least the envelope's donors-per-group in **both** groups; the thresholds are
read from `gate_config.OPERATING_ENVELOPE`, not restated as literals, and the freeze script verifies
the declared membership against the manifest. **All four** envelope tiers are declared, not only the
two whose failure was obvious — an earlier draft declared 0.5 and 0.7 alone, and the two undeclared
rows of this table were consequently published without ever being checked against the manifest.

The last column is not typed. It is `below_spec_dataset_floor` from
[`stratum_list_2026-08-16.json`](../pilot/preregistration/stratum_list_2026-08-16.json), which the
freeze script computes as `n_datasets < SPEC_DATASET_FLOOR` with `SPEC_DATASET_FLOOR = 8`, and
`tests/test_stratum_list_freeze.py` parses this table and fails if any cell of it disagrees.

| σ_donor anchor lands at | Envelope demands | Surviving datasets | Surviving strata | Against §1's 8–12 floor |
|---|---|---|---|---|
| ≈ 0.2 | ≥ 4 v 4 | 11 of 12 (all but #12) | 227 | met |
| **≈ 0.35** | **≥ 8 v 8** | **7 of 12 — #1, #2, #4, #5, #6, #7, #8** | **150** | **BELOW** |
| **≈ 0.5** | **≥ 13 v 13** | **5 of 12 — #1, #4, #6, #7, #8** | **94** | **BELOW** |
| ≈ 0.7 | ≥ 23 v 23 | **3 of 12 — #1, #6, #8** | 30 | **BELOW** |

**Three of the four rows are the ones this section exists for. The list falls below §1's own
"8–12 datasets" floor at every tier except the most optimistic.** Seven is below eight, and an
earlier draft of this table called the σ ≈ 0.35 row *met, with no margin*, which was simply wrong
against the rule the freeze script has always encoded.

**The σ ≈ 0.35 failure is the one that matters, because 0.35 is not a pessimistic scenario.** It is
`gate_config.POWER_EVAL_SIGMA` — the envelope boundary that [Amendment 3](AMENDMENTS.md) Change 1(b)
makes binding, the point at which the shipped `ebayes` arm is grid-shown to clear the pre-registered
power target (0.793 at 8 v 8, calibrated), and therefore the instrument's own nominal operating
point. If the anchor lands there, the study is already below §1's dataset floor: **7 of 12, and it
gets worse from there — 5 at σ ≈ 0.5 and 3 at σ ≈ 0.7.** Only σ ≈ 0.2, the most optimistic tier in
the envelope, keeps 11. That is not a tail risk being noted; it is the expected case failing, and
softening it to "no margin" would have been the single most consequential piece of self-flattery
available in this document.

**And the shortfall is ours, not the data's.** The claim that no other list from this manifest could
have done better is false, and the artifact committed beside this document refutes it. Re-derived
from `census_candidates_run31910799023_2026-08-15.json`, over its 68 candidate-bearing datasets:

| Envelope demand | Datasets in the manifest clearing it anywhere | The frozen twelve | Most a §1 (iii)-compliant twelve could have kept |
|---|---|---|---|
| ≥ 4 v 4 | 62 | 11 | 12 |
| ≥ 8 v 8 | **33** | 7 | 12 |
| ≥ 13 v 13 | **21** | 5 | 11 |
| ≥ 23 v 23 | **12** | 3 | 11 |

Twenty-one is not a handful, and **21 > 12**: twelve datasets that all clear 13 v 13 could have been
picked out of this very manifest, so a list surviving σ ≈ 0.5 well above §1's floor was assemblable
from it. The last column is what a list optimised for donor counts would
have retained, and it is constructive rather than an upper bound nobody could reach. §1 (iii)
requires "some exactly 3v3", and **only 3 of the 68 candidate-bearing datasets hold an exactly-3v3
stratum at all** — Rexach (#2, 2 of its 27 strata), Wang (#12, all 8) and the Emphysema Cell Atlas
immune split (all 9), the last two being the two halves of one collection, so in truth only **two
independent collections** can supply the anchor. **None of the three clears 13 v 13**, so the anchor
costs a slot at the two hard tiers and eleven is the ceiling there. Rexach does clear 8 v 8, so at
the two easy tiers twelve is attainable. The witness is explicit: **Rexach as the 3v3 anchor plus
eleven of the twelve datasets that clear 23 v 23** is a twelve-dataset, §1 (iii)-compliant list
retaining **12 / 12 / 11 / 11** datasets across the four tiers — above §1's floor of eight at
*every* tier, where the frozen list is above it at one.

So the correct statement, which is both truer and worse for us than the one it replaces:

> **The truncation is a consequence of selecting for §1 (iii)'s coverage axes — tissue, assay
> chemistry, suspension, effect-size range and the mandatory 3v3 anchor — rather than for donor
> counts. A list optimised for donor counts would have survived the truncation and would have
> spanned less.**

Two things follow. First, the feasibility risk in §6 is a **cost of the coverage requirement**, and
if the anchor lands at 0.35 or worse the honest report is "§1's span requirement and §1's dataset
floor are not jointly satisfiable at this Census pin under this envelope" — a conflict inside the
spec, not a shortage of public data. Second, and this is why it is written here rather than
discovered later: we know this **before the anchor exists**. Saying it now is precisely what makes
it impossible to reassemble the list on donor counts afterwards and present that as the original
plan. The counterfactual is a statement of what was given up, not a reserve list; nothing in it may
be selected from (§9.7), and the three figures above are recomputed by the freeze script on every
run, which is what stops this paragraph from drifting back into the comfortable version.

Amendment 3 named this outcome in advance and the wording is quoted rather than paraphrased
(**emphasis added**), because it is the commitment being honoured:

> **Whether the real sweep is feasible at all.** If real strata cluster at `sigma_donor` ≈ 0.5–0.7,
> the envelope admits them only at ≥ 13–23 donors per group. Whether enough CELLxGENE strata clear
> that is an open empirical question, **and a negative answer is a live outcome of this study, not a
> failure mode to be designed around.**

Two things this table is **not**. It is not a measurement: `sigma_donor` is unmeasured for every
stratum here, so these are donor-count arithmetic under a hypothetical σ, i.e. a scenario analysis.
And it is not an admission at any tier — Amendment 3's conversion from the moderated fit's `s0²` to
`donor_sigma` is an unvalidated upper bound, and "no stratum may be admitted to the real sweep on
the strength of this entry alone".

Writing this down is the entire mechanism. Without it, an anchor at 0.35 would leave us free to
report "the seven datasets that were inside the envelope" as though seven had always been the plan.
With it, seven is visibly a *truncation of twelve*, it is visibly below §1's floor, and it is
reported as a limitation of the instrument's applicability domain.

---

## 7. Gating already fixed by the spec, restated where it will be applied

These rules are the spec's, restated because they bind at exactly the points this list is unusual,
and a rule that lives only in a document nobody re-reads at analysis time is not a control. One
thing here **is** new and is flagged rather than smuggled in: §4 says only "small D → enumerate all,
report the exact null distribution, flag coarse resolution" and never says what *small* is. The
threshold of **fewer than 1000 distinct label assignments** below is this document's
operationalisation of that looser text. It is pre-registered here, before any null is computed, and
it is chosen to coincide with §4's own `n_perm = 1000` for the sampled case — below it, sampling
would draw more permutations than exist — rather than being fitted to how many strata it catches.

**A1 — the kill-switch is decision-relevant only at ≥ 8 v 8, and is leak-contaminated at 3 v 3.**
The signal-above-floor ratio (decision rule item 4, §6) is treated as decision-relevant **only** in
strata with at least 8 donors per group, "where balanced permutations are ~orthogonal to the true
grouping"; on 3 v 3 strata it is "flagged leak-contaminated and excluded from GO logic". In this
list that partitions the 251 concretely: **150 strata across 7 datasets carry the kill-switch and
101 do not.** Of the 101, **24 have `min(A,B) = 3`** — among them the 10 that are exactly 3 v 3, the
case A1 names as leak-contaminated — and **77 sit between 4 v 4 and 7 v 7**, below A1's threshold
and therefore also outside the kill-switch, which is the reading A1's text forces even though it
names only the 3 v 3 case explicitly. Dataset #12 exists to satisfy §1's "some exactly 3v3" and
contributes nothing to item 4 by construction; that is its declared role, not a disappointment to
be discovered later.

The same boundary bites on permutation resolution. **45 of the 251 strata have fewer than 1000
distinct label assignments** (minimum 20, at 3 v 3, which §4 reduces to 10 after removing the
identity and its complement), so for those §4's "enumerate all, report the exact null distribution,
flag coarse resolution" applies and the null's resolution must be printed next to every number
derived from it.

**D2 — evidence clusters by dataset, not by stratum.** The majority in decision rule item 2 is over
**independent datasets**, and the denominator is the 12 of §3.2 — not 251, not 14. Same-dataset
strata are not independent replicates: #4 alone would contribute 47 of the 251, and #1, #4 and #8
together contribute 102. Per-stratum numbers are reported alongside, never in place of, per-dataset
aggregates, and the two siblings of §3.3 stay outside the denominator whatever they show.

**No pooled headline across incompatible granularity (D5).** Restated in §9 because it is a
consequence of this particular list rather than a general rule.

---

## 8. §8(d) — the Mathys 2019 oracle. Not a member of this list, and not weakened by it

**Mathys 2019 is not in the stratum list and must not be.** §8(d) makes it the *binding* real anchor
against the possibility that our simulators are optimistic: it is an oracle with a known target, not
a stratum whose inflation we are measuring. Putting it in the analysis set would convert the check
into one of the things being checked.

**The Census path is closed, and this is now checked against the pinned release itself rather than
inferred from Discover.** Mathys 2019 is absent from CELLxGENE Discover, and it is absent from
Census `2025-01-30`. Those are two claims and the second does not follow from the first, which the
earlier draft of this paragraph assumed it did.

*The Discover half.* The full curation index (`GET /curation/v1/datasets`, **2216 datasets**, read
2026-08-16) was searched case-insensitively for `mathys`, `rosmap`, `religious order`, `memory and
aging` and the DOI stem `s41586-019-1195` / `10.1038/s41586-019-1195-2` — first over
`title + collection_name + collection_doi + collection_doi_label + citation` (the method of Appendix
A.3 of the proposal document, §2), then over each record's *entire* JSON in case the provenance sat
in a field the original query did not read. **Zero hits, all six needles, both passes.**

*Why that was not enough.* The index was read in 2026; the Census is pinned to `2025-01-30`; Discover
can withdraw a dataset after a release is built. The control the earlier draft offered — all 73
manifest dataset ids resolve in Discover — cannot detect that by construction: those 73 ids are the
ones the Census reader already returned, so they are guaranteed to be in both. The claim needs
Discover to be a **superset of the pinned release**, and that is a different check.

*The pinned release, enumerated directly.* `pip install -e ".[census]"` is **not feasible on this
machine**: `cellxgene-census` requires `tiledbsoma`, which has never published a Windows wheel (PyPI
JSON API, **zero Windows artifacts across all 49 releases**) and whose sdist build fails without the
prebuilt native library, so `census["census_info"]["datasets"]` was unreachable. The release's own
contents were enumerated instead, from the public release bucket over plain HTTPS: every
`{dataset_id}.h5ad` published under
`s3://cellxgene-census-public-us-west-2/cell-census/2025-01-30/h5ads/`, which is the pinned build's
own bytes and not an index of them. **1573 datasets.** All 73 manifest ids appear in that listing,
which is the control that the listing is the right object.

*The result, which is not the tidy one.* **1567 of the 1573 resolve in Discover; 6 do not.** So the
inference "absent from Discover ⇒ absent from Census `2025-01-30`" is **false in general** — there
is a six-dataset blind spot, and the earlier sentence walked straight through it. The blind spot is
small enough to close by enumeration, and it was: each of the six was identified by reading
`uns/title` and `uns/citation` out of its own h5ad in the release bucket over HTTP range requests
(~1 MB transferred per file, no matrix read). Four are the lamina-propria / submucosa sorts of
collection `0c3f148e-02ff-4c81-8946-29beaaf5fa59` (`10.1101/2021.03.28.437379`) and two are the
cross-species pancreatic alpha- and beta-cell maps of collection
`0a77d4c0-d5d0-40f0-aa1a-5e1429bcbd7e` (`10.1016/j.molmet.2022.101595`). **None is Mathys, ROSMAP or
prefrontal cortex.** The conclusion therefore stands on the release's contents rather than on an
inference from an index: Mathys 2019 is not in Census `2025-01-30`.

§8(d)'s own parenthesis — "via Census if present, else Synapse syn18485175" — anticipated the
absence; what follows from it did not get anticipated, and is recorded here.

**What the remaining path costs, stated before it is walked.**

1. **A data-use agreement, on the critical path.** Synapse `syn18485175` is ROSMAP data and requires
   a ROSMAP DUA. This is a legal dependency, not a technical one, and it gates the binding check.
   The application was begun on **2026-08-16**, the date of this freeze; **access had not been
   granted at the time of writing**, and no data from the deposit has been seen. Elapsed time to
   access is a schedule risk with no engineering mitigation, and it is recorded here so that a long
   wait is a known cost rather than a discovery.
2. **A second loader, which is work not in §9's plan.** `io_counts.load_stratum(dataset_id,
   cell_type)` is Census-shaped throughout: it wraps `cellxgene_census.get_anndata` on the pinned
   version's raw layer and keys on `dataset_id`. A Synapse-hosted matrix shares none of that. A
   second loading path is therefore **required work that §9's implementation plan does not contain**,
   and it must carry the same integrality gate (item 4), the same thin-donor filter and the same
   frozen-universe check as the Census path, or the anchor is not running the same pipeline and
   proves nothing. Per public descriptions of the deposit — recorded as such, not verified here —
   the resource is ~80 660 nuclei from prefrontal cortex across 48 ROSMAP donors, 24 with AD
   pathology and 24 without, with a filtered matrix at `syn18681734`.
3. **The target is fixed and external.** §8(d) requires qualitative reproduction of Murphy & Skene
   2023 (eLife 90214): naive per-cell DE grossly exceeding pseudobulk, with the permutation-null
   floor accounting for most naive calls. That paper reports on the order of a 549-fold reduction in
   DEG count at FDR 0.05 — and it attributes that figure to its **whole corrected re-analysis**,
   corrected quality control *and* pseudobulk aggregation together, not to the change of replication
   unit alone. Attributing it to pseudobulk by itself, as an earlier draft did, overstates what the
   source claims and would set our anchor against a number nobody measured. It is **quoted from that
   publication and not re-derived in this repository**, and it is context rather than a threshold,
   since §8(d) asks for a qualitative pattern.

**This oracle does not block the freeze, and the freeze does not weaken the oracle.** Making the
stratum list wait on a DUA would stall the pre-registration behind a legal process while the
candidate manifest sat un-frozen and readable — the precise condition under which a list stops being
pre-registered. Conversely, nothing here downgrades §8(d): it remains **BINDING**, it remains
**unrun**, and no result from the 251 strata may be reported as validated by anything except the
oracle §8(d) names. If the DUA is refused, that is an amendment, and the amendment will have to
state what the study's binding real-data check is instead — not quietly drop the requirement.

---

## 9. What this freeze does NOT settle

Written in the amendments' own genre, and not as a formality: every item below is a live hole and
several of them can still shrink the list.

1. **Admission is still closed, and this document does not open it.** All four blockers stand on all
   2190 rows of the source manifest and on all 251 rows here: `integer_check`,
   `frozen_universe_size`, `sigma_donor_estimate`, `envelope_membership`. `admitted_to_sweep` is
   `False` everywhere. Freezing *what will be measured* is upstream of deciding *what may be
   reported*, and only the first has happened.

2. **`sigma_donor` is unanchored, so envelope membership is unknown for every stratum.** This is the
   fourth consecutive pre-registration document to close on it. Amendment 3 supplied the envelope
   and a *mechanism* (`sqrt(s0²) · ln 2` from the moderated fit) and stated in terms that the
   quantity is an **upper bound**, the conversion is unvalidated against the simulator where truth
   is known, and "the `sigma_donor` anchoring demanded by Amendment 1 therefore remains OPEN". §6's
   tiers are arithmetic under a hypothetical σ, not measurements, and no stratum here is inside the
   envelope until that work is done.

3. **`pooled` is `unresolved` on 1197 of 1197 candidates, and on all 251 frozen strata.** The pinned
   Census's `obs` exposes no library/pool/multiplex identifier at all, so D3's answer is not "not
   pooled" — it is the "where pooling cannot be resolved" state. Consequently **the donor-pseudobulk
   gold-standard claim cannot be made on any stratum in this list**: donor pseudobulk is a *lower
   bound* on the correct replication unit throughout, and every calibration statement inherits that
   caveat. This is a property of the Census pin, not of the selection; no choice of datasets fixes
   it, and choosing differently would only have hidden it behind different ids.

4. **The D5 granularity conflict is inside the list, by choice.** SEA-AD (#1) annotates at
   `L2/3-6 intratelencephalic projecting glutamatergic neuron` and `chandelier pvalb`; Rexach (#2)
   annotates at `astrocyte` and `glutamatergic neuron`. D5 forbids pooling a headline across strata
   of differing granularity, so **no pooled headline may be reported across #1 and #2**, and the
   same caution applies wherever the 124 cell-type labels differ in depth. `cell_type_ontology_depth`
   is `pending` on every row (it is a property of the CL graph, which `census_select` deliberately
   does not fetch), so the harmonised level D5 prefers is **not yet chosen**. Both datasets are in
   the list because dropping one to make pooling easy would trade evidence for tidiness; the cost is
   that the headline is per-dataset here and the aggregate is D2's, not a pool. **Choosing that
   harmonised level is an amendment**: collapsing 124 labels to a common depth merges strata, and
   merged strata are not the 251 frozen here — an open question that can silently change the frozen
   set is not an open question, it is an unlocked door.

5. **Tissue region within `dataset_id` is an unmodelled batch.** Rexach (#2) spans Brodmann area 4,
   insular cortex and primary visual cortex; KPMP (#8) spans kidney cortex, medulla and papilla;
   Phan (#10) spans caudate and putamen. `tissue_general` collapses all of these to one level, so
   the §1 pre-screen did not catch them — not because the confound is absent but because it is not
   in the column being screened. Whether region enters the *definition* of a stratum is **not
   decided here**, and it needs deciding before those datasets are analysed. **Deciding it is an
   amendment, because it changes what the 251 are**: splitting a dataset by region splits its
   strata, and it moves Layer B — KPMP at 24 donors' ceiling would not hold a 23 v 23 stratum inside
   a single region, so the σ ≈ 0.7 tier would lose one of its three survivors (and only 6 KPMP
   strata reach 23 v 23 pooled across regions at all). A decision that can shrink the pre-declared
   truncation is not a detail to settle at analysis time.

   **"Or the covariate set" is struck from this item, because the shipped arm has no covariate
   slot.** Item 8 below states it plainly: C4 lapsed with DESeq2 (Amendment 2 Change 1) and the
   moderated arm fits `~ 1 + x` via `wls_two_group`. Offering region-as-covariate as a live option
   in the same section that records the absence of any covariate slot was a contradiction inside one
   §, and the option is **unavailable** until some amendment supplies a design that has one.

6. **Sex, developmental stage and self-reported ethnicity are outside §1's confound screen.** §1
   screens exactly five covariates — assay, suspension type, `tissue_general`, a sequencing-depth
   bin and a pool id — and these three are materialised into `obs` but never tested. Two datasets
   make the omission concrete: Yoshida (#4) mixes children and adults in one cohort, and Elmentaite
   (#11) is entirely paediatric at 4–14 years, where developmental stage is an obvious candidate
   confounder with disease. **This is a hole in §1, not in the manifest**, and closing it is a spec
   change that would need an amendment; it is named here so that it is a known omission rather than
   an accident of the frozen list.

7. **`integer_check` and `frozen_universe_size` are unverified for all 12 datasets, so the list may
   still shrink.** Both are computed at X load and neither can be decided from `obs`
   (§1 items 4 and 5; §10 risk 6). A dataset whose raw layer is not integral is **dropped with its
   reason, never rounded**; a stratum whose frozen universe falls below 200 genes is a SKIP whose
   measured size is still recorded. The most likely casualty is identifiable in advance and is named
   now rather than after the fact: **#7 (rheumatoid arthritis blood) runs at 284–803 median counts
   per cell** against the frozen set's median of 3081 — roughly an order of magnitude thinner — and
   its smaller strata may not assemble a 200-gene universe. That is a stated expectation, not a
   reason to pre-emptively drop it: a dataset excluded now on a suspicion would be exactly the
   discretion §3.1 exists to remove. **Any such shrinkage is a reported outcome, never a
   re-selection**:
   the dropped stratum is reported with its reason and its measured value, and no replacement is
   substituted for it. Replacing a stratum that failed a counts gate with one that did not is
   selection on the data, arriving through the back door — and that prohibition binds the five
   reserves named in the proposal document (§2) exactly as it binds anything else. They are
   committed so that a substitution could be *detected*, not so that one becomes available.

   Two further questions can only be settled by loading X, and are listed here rather than left
   implicit. **(a) Whether SEA-AD's multiome cells are inside the frozen strata** (§4.2): Discover
   lists `10x multiome` for #1 while the manifest's `confound_cramers_v["assay"]` is null on all 18
   of its strata, and `obs` cannot say whether that means the multiome cells are absent from those
   strata or present and unscreened. Only reading `assay` per cell at load decides it. If they are
   present, #1's strata carry an unscreened assay mix and must be tagged as such — a reported
   property of a frozen stratum, never a reason to drop it. **(b) Whether the raw layer is integral
   for the multiome-derived cells specifically**, which is the same question as item 4 of the
   inclusion gate but with a named candidate for failing it.

8. **The confound tags on the frozen set are tags, and nothing downstream acts on them.** Of the 251
   strata, 231 carry a partial `sequencing_depth_bin` confound tag, 8 a near-confound tag
   (V ≥ 0.8), and 24 a partial `assay` tag; none carries a §1 exclusion, and none of the 12 datasets
   lost a stratum to the confound pre-screen (all 55 non-candidate contrasts among them failed the
   *inclusion gate*, mostly on donor counts). §1 would carry a tagged covariate into the design "only
   if C4's df rule allows", but C4 lapsed with DESeq2 (Amendment 2 Change 1) and the shipped
   moderated arm fits `~ 1 + x` with **no covariate slot**. A partially confounded stratum here is
   therefore neutralised **only** by the permutation null (§7.1), which is what the tag means and all
   it means. Depth in particular is not modelled anywhere.

9. **D4's excluded fractions are not whole-Census fractions.** The manifest's `counts` block reports
   an excluded fraction of 0.4534 over its 2190 contrasts, but those contrasts come only from the
   datasets pass 2 actually read; datasets the coarse pass-1 filter judged incapable of holding a
   gate-clearing stratum were never queried and appear in no denominator. Within the 12, the split is
   251 candidates against 55 inclusion-gate failures out of 306 contrasts. Any D4 statement must
   carry that scope.

10. **The freeze fixes the list, not the analysis.** `controls.py`, `decision.py` and `report.py`
    (§9 items 9, 12, 13) do not exist yet; the cells-per-donor sweep that D1 makes the primary
    conditioning axis is unwritten, and so is the code that applies the GO/NO-GO rule. A stratum list
    frozen before the analysis code is written is the right order — it is the only order in which the
    list cannot be shaped by what the code turns out to make easy — but it means nothing here has
    been executed end to end on real data.

11. **A2 stratification remains deferred** (Amendment 2 Change 6), so the naive arm's permutation
    floor is still a cell-count-confounded quantity guarded by a range check only. Nothing about this
    list improves that, and its skewed designs make it more visible: #5 runs 100 versus 10 donors and
    #3 runs 20 versus 7, where per-group cell totals differ systematically between the arms.

12. **`permutation_count` flatters skewed designs and must not be read as resolution.** It is
    C(D, n_A) over the *total* donor count: #5's COVID-19 strata reach 4.69 × 10¹³ permutations on a
    control group of at most 10 donors, and #1 reaches 7.28 × 10²³. The informative quantity is
    `min(n_donors_A, n_donors_B)`, which is why it is emitted as its own column in the frozen
    artifact. Any report that prints `permutation_count` must print `min_donors_per_group` beside it.

---

## 10. Author attestation

**Re-derived while writing this document**, from the committed artifact and not quoted from the
2026-08-16 proposal it acts on: the sha256 sums and byte counts of all three committed source files;
the manifest's
`generated_utc`, `census_version`, row count and gate-status split (2190 = 1197 + 981 + 12); the
membership of all twelve datasets and both siblings by full uuid; the 251 frozen strata and their
per-dataset counts (18, 27, 28, 47, 25, 11, 15, 37, 7, 10, 18, 8); the 27 within-collection control
strata and their split (18 + 9) with the siblings' ceilings (42 and 3); every donor ceiling in §3.2;
the ≥ 8v8 count of 150 and the 101 below it with its 24 / 77 split; the exactly-3v3 count of 10; the
six-bin occupancy over 502 group medians and its 11.0 … 6671.5 range; **all four** Layer B subsets —
≥ 4v4 (11 datasets, 227 strata), ≥ 8v8 (7, 150), ≥ 13v13 (5, 94), ≥ 23v23 (3, 30) — and the
`below_spec_dataset_floor` verdict of each; `admitted_to_sweep = False` on all 2190
rows; `pooled_flag = unresolved` on all 1197 candidates; `integer_check` and `frozen_universe_size`
`pending` on all 1197; the confound-tag counts 231 / 8 / 24 and the 55 inclusion-gate failures inside
the twelve; the 45 strata below 1000 permutations; the counts-per-cell range 284 … 56 841.5;
SEA-AD's null `confound_cramers_v["assay"]` on all 18 strata against Rexach's non-zero V on 24 of 27;
and the sibling donor spans (DLPFC 33 … 39, MTG 27 … 42) that make "the same donors" false.

**§6's counterfactual, re-derived and now recomputed on every run.** Of the manifest's 68
candidate-bearing datasets, **62 / 33 / 21 / 12** hold at least one stratum at ≥ 4v4 / ≥ 8v8 /
≥ 13v13 / ≥ 23v23; **3** hold an exactly-3v3 stratum (Rexach, Wang, the Emphysema immune sibling —
the last two being one collection), and **none of those three clears 13v13**, so §1 (iii)'s 3v3
requirement costs a dataset slot at the two hard tiers. A §1 (iii)-compliant twelve could therefore
have retained at most **12 / 12 / 11 / 11**, attained by Rexach plus eleven of the twelve ≥ 23v23
datasets. The earlier claim that "only a handful clear 13 v 13 anywhere, so the shortfall is a
property of the public data" was false against this same artifact and is retracted in §6.

**The pinned Census release, enumerated rather than inferred (§8).** 1573 datasets published under
`cell-census/2025-01-30/h5ads/`; 2216 in the Discover curation index read 2026-08-16; **zero hits**
for all six Mathys / ROSMAP needles over every field of every Discover record; all 73 manifest
dataset ids present in both. **1567 of the 1573 resolve in Discover and 6 do not**, so "absent from
Discover, therefore absent from Census" does not follow in general; the six were identified
individually from their own h5ads (two collections, `0c3f148e-…` and `0a77d4c0-…`) and none is
Mathys. `cellxgene-census` could not be installed to read `census_info/datasets` directly:
`tiledbsoma` publishes no Windows wheel at any version.

`scripts/freeze_stratum_list.py` recomputes the load-bearing
subset of these on every run and aborts on disagreement, so they cannot drift silently from the
document; `tests/test_stratum_list_freeze.py` re-derives the rest from the raw manifest rather than
through that script, so a bug in it cannot be reproduced by the test meant to catch it, and it
parses §6's table out of this file so a hand-typed verdict cannot disagree with the artifact.

**Two figures disagreed with what had been circulated, and are corrected rather than absorbed.**
Rexach's envelope ceiling is `min(A, B) = 10`, not the 11 stated in the 2026-08-16 proposal
document; §3.2 records the error, its magnitude and its cause. And the proposal's reading of the
manifest's donor-count distribution understated what a differently-chosen twelve could have
supported; §6 carries the corrected figures. The proposal is committed unmodified with both errors
in it (§2), because a hash of a corrected copy would not be a hash of the document that was acted
on.

**Judgement, not measurement.** The choice of these twelve datasets from 68 candidate-bearing ones;
the strong / moderate / subtle labels of §4.1 and the pre-declaration of #9 as an expected no-effect
anchor; the decision to accept a 5′ arm that is confounded with tissue rather than take an
in-vitro-stimulated alternative; the decision to keep #1 and #2 despite their D5 incompatibility; and
the choice of half-decade log bins in §5. The bins were fixed before their occupancy was computed;
the effect-size labels were fixed before any metric was computed; neither is a prediction of any
result, and the freeze script pins both so that neither can be revised after one.

**Could not be established.** The publication DOI of dataset #8 (KPMP kidney v1.5): the Discover
curation API returns `collection_doi = null` for that collection, so the provenance is recorded as
*not established* rather than inferred. Assay, suspension type and tissue for all twelve are taken
from the Discover index and are **not** re-derivable from the committed manifest, which does not
carry those columns. Whether SEA-AD's `10x multiome` cells are inside the 18 frozen strata is
undecidable from `obs` and is deferred to X load (§4.2, §9.7). Murphy & Skene's 549-fold figure is
quoted from the publication, which attributes it to the whole corrected re-analysis rather than to
pseudobulk alone, and was not re-derived here. The public description of the Mathys deposit
(~80 660 nuclei, 48 donors, 24 versus 24) is recorded as a public description and was not verified
against Synapse, since access requires the DUA that is still pending and has not been granted.

**Asserted on expectation, not evidence.** §1 (i)'s parenthetical *(pseudobulk shown non-null)* is
**not satisfied** by anything checked for this freeze: the strong / moderate / subtle labels rest on
the biology of each study's design and not on published donor-level pseudobulk results, no such
per-publication review was performed, and §4.1 records this as a known gap in what §1 asked for.

**Not claimed.** That any stratum in this list is inside the operating envelope; that any of them
will survive the counts gate; that donor pseudobulk is calibrated on any of them (D3 forbids the
claim while pooling is unresolved); or that the strong / subtle labels will match the measured
inflation. No metric has been computed on any of these 251 strata.
