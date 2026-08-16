# pbcheck — Phase 0 Pre-registration of the Stratum List (2026-08-16)

The frozen specification, [`PHASE0_SPEC.md`](PHASE0_SPEC.md), closes spec §1 with a sentence that is
an instruction and not a description:

> **Pre-register the stratum list before computing any metric.**

This document is that act. It is **not** an amendment: it changes no threshold, supersedes no
section and relaxes nothing. It executes a step the frozen spec demands and had left undone, and it
closes one numerical gap in spec §1 that the spec references three times and never states (the
cells-per-donor bins, §5 below).

**Cross-references, because they used to switch document silently.** `spec §N` is
[`PHASE0_SPEC.md`](PHASE0_SPEC.md); a bare `§N` or `§N.M` is a section of **this** document;
`Amendment N` is [`AMENDMENTS.md`](AMENDMENTS.md). `A1`, `A2`, `C4` and `D1`–`D5` are the spec's own
correction and decision codes. The convention is machine-enforced: `freeze_stratum_list.py` resolves
every `§` in this file against the headings of the document it names, and aborts if one points at
nothing. That guard exists because an earlier draft used a dotted subsection reference for
spec §7 item 1 — a notation naming nothing in either file, since the spec has no subsections and
this document's §7 has no numbered items — and because eight bare references changed document
relative to their neighbours.

Its machine-readable half is [`../pilot/preregistration/stratum_list_2026-08-16.json`](../pilot/preregistration/stratum_list_2026-08-16.json)
(and its CSV twin), emitted by [`../scripts/freeze_stratum_list.py`](../scripts/freeze_stratum_list.py)
from three hash-pinned sources. Neither half is authoritative alone, and neither is trusted: the
script re-derives every falsifiable figure below and then **parses this document's tables back out
of the file and compares them cell by cell**, refusing to run on a disagreement in either direction.

**The rule this round enforces, without exception.** Every falsifiable claim in this document is
recomputed by the freeze script from a source committed in the repository and pinned by sha256, or
it is labelled at the point it is made as a judgement, a quotation or an external fact. A claim that
is neither is deleted. §10 lists what was deleted and what is labelled. This is not a stylistic
preference: three adversarial reviews of this freeze found the same defect three times — a claim in
the prose that the artifact committed beside it refutes — and every one of them sat in the part of
the document that nothing recomputed.

---

## 1. What this act is, and what it binds

| | |
|---|---|
| Date of the freeze | **2026-08-16** |
| Census | **`2025-01-30`**, the spec §1 pin — `open_census()` rejects `latest` and `stable` by name |
| Independent datasets | **12**, in **12 distinct collections** |
| Frozen strata (stratum-contrasts) in the analysis set | **251** |
| Collection-siblings, frozen as within-collection controls and excluded from the D2 denominator | **5** datasets, **106** strata |
| Metrics computed on any of these strata at the time of writing | **none** |
| Rows admitted to the sweep | **0** — `admitted_to_sweep = False` on all 2190 rows of the source manifest, and on all 357 here |

**What it binds.** From this commit forward, the analysis set of the Phase 0 real-data sweep is the
251 strata named here, and the "independent datasets" denominator of decision rule item 2 (D2) is
the 12 dataset ids named here. Any addition, removal, re-selection or re-interpretation — including
dropping a dataset because its numbers turn out unhelpful, and including adding one because the
survivors look thin — requires a dated, numbered entry in [`AMENDMENTS.md`](AMENDMENTS.md), written
**before** the change is applied, in the form that log already uses.

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

### 2.1 The candidate manifest

Everything about strata, donors and cells is derived from one artifact, committed here in full
rather than referenced:

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
| Header `value_filter` | `is_primary_data == True and disease != 'na'` |
| Produced by | `scripts/census_candidates.py` via `.github/workflows/census-candidates.yml`, **GitHub Actions run `31910799023`** (`workflow_dispatch`, `dry_run: false`) |
| Shape | 2190 stratum-contrasts over 73 datasets; 1197 `candidate`, 981 `excluded_inclusion_gate`, 12 `excluded_confound` |

The manifest is normally a **CI artifact and never a commit** — `pilot/results/` is gitignored
precisely so that "a candidate list that reached git by way of a CI job would have pre-registered
itself by accident". That reasoning is why it is committed *here*, under `pilot/preregistration/`,
by a deliberate human act with this document attached: a pre-registration whose evidence lives in an
expiring CI artifact is not auditable, and an evidence file whose bytes cannot be checked is not
evidence. `.gitattributes` marks the directory `-text` so no end-of-line conversion can move the
hashes.

### 2.2 The two external indexes, pinned — and what that does not fix

The manifest carries **no collection, assay, suspension, tissue or DOI column**. Until this commit
those properties were read live from the CELLxGENE Discover API and never committed, so no reader
could check a single claim that rested on them — and both blocking defects of the third review were
exactly that: claims about collections. Two snapshots are now committed and hash-pinned, and every
such claim is recomputed from their bytes.

| | Discover index | Pinned release table |
|---|---|---|
| File | [`../pilot/preregistration/discover_index_2026-08-16.json`](../pilot/preregistration/discover_index_2026-08-16.json) | [`../pilot/preregistration/census_release_datasets_2025-01-30.json`](../pilot/preregistration/census_release_datasets_2025-01-30.json) |
| sha256 | `afc74c1c1ea8f22c9c86a7cd6a2e4eb8087b7db58c37b447f8383e76f4eaf416` | `b60e8e1920de09ef3e1a6de595d574bddca622ea64adf86498ae14bcfa26da0e` |
| Size | 2 343 258 bytes | 1 367 664 bytes |
| Source | `GET https://api.cellxgene.cziscience.com/curation/v1/datasets`, read **2026-08-16** | the release's own `census_info/datasets` array under `cell-census/2025-01-30/soma/`, read over anonymous S3 with TileDB-Py on **2026-08-16** |
| Records | **2216**, the whole index | **1573**, the whole release |
| Reduced to | `dataset_id`, `dataset_version_id`, `collection_id`, `collection_name`, `collection_doi`, `collection_doi_label`, `citation`, `title`, `assay`, `suspension_type`, `tissue`, `cell_count`, `is_primary_data` | `dataset_id`, `dataset_version_id`, `collection_id`, `collection_name`, `collection_doi`, `collection_doi_label`, `dataset_title`, `dataset_h5ad_path`, `citation` |

Both are written by [`../scripts/fetch_preregistration_evidence.py`](../scripts/fetch_preregistration_evidence.py),
which also records in each header the read date, the endpoint or URI, the full upstream size and the
sha256 of the raw upstream payload the reduction was made from. A file cannot carry its own sha256,
so each reduced file's hash is recorded here and pinned in the freeze script instead.

**The limitation this does not remove, stated plainly.** The Discover snapshot was read on
2026-08-16 and describes **Discover today**. The Census is pinned to `2025-01-30`. These are
measurably different objects, and the measurements are recomputed on every freeze:

| Comparison, over the 1573 datasets of the pinned release | |
|---|---|
| resolve in Discover by `dataset_id` | **1567** |
| do not resolve in Discover at all | **6** |
| `dataset_version_id` also present in Discover | **0 of 1573** |
| `collection_doi` differs between the two | **61** |
| `collection_id` differs between the two | **0** |

So the Discover snapshot is evidence **about Discover**, not about the release, and this document
says which of the two each claim rests on. Collection membership is the one place they can be
checked against each other, and they agree exactly — `collection_id` is identical for all 1567
datasets the two share, including all 73 of the manifest's — so §3.3's sibling set is derived from
the **release** table, which is the pin, and the agreement is asserted by the freeze script rather
than assumed.

**A recorded discrepancy, not argued away.** The manifest's own Census query filtered on
`is_primary_data == True`. Two of the seventeen datasets frozen here carry
`is_primary_data = [False]` in the 2026-08-16 Discover snapshot: **#8, the KPMP single-nucleus
kidney atlas, and its within-collection control, the KPMP single-cell arm** — both KPMP, and no
others. Discover today says neither holds a primary cell; the pinned release's `obs` returned cells
for both under a filter demanding one. Nothing here decides which is right, and no stratum is kept
or dropped on the strength of it. It is recorded because a pre-registration that noticed a
disagreement between its sources and did not write it down would be worth less than one that did.

### 2.3 The proposal document, and its redaction

| | |
|---|---|
| File, as committed | [`../pilot/preregistration/stratum_list_proposal_2026-08-16.redacted.md`](../pilot/preregistration/stratum_list_proposal_2026-08-16.redacted.md) |
| sha256, as committed | `5588ba845cc144b56ea27a25ca4599f3fbc33d69c5eeb8833da7772a9459f07d` |
| Size, as committed | 93 873 bytes |
| sha256, as circulated | `50872414b0727c129a824b0c65ed179674ac5d6c9ecaac53327568b3eae6fb48` |
| Size, as circulated | 92 589 bytes |
| Status | **superseded working document** — the reasoning behind the choice of twelve, committed for auditability; **not part of the binding act** |

The binding content of this pre-registration is the list of §3.2 and the rule of §3.1, and nothing
else. The proposal is the working document those were chosen from: it carries the per-dataset
rationale, the datasets considered and rejected, the **five named reserves**, the third 5′ candidate
and the method of the Mathys search. It is referenced by §3.2 and §8, and until this commit it
existed only in a scratch directory — so the one act of discretion in the whole freeze was justified
by a file no reader could open. The reserves are the sharpest case: §9 item 7 forbids substituting a
replacement for a stratum that fails the counts gate, and a list of five candidate replacements must
not live outside the record where nobody can check whether one was quietly used.

**What was redacted, so the redaction is itself auditable.** The circulated copy contained an
absolute filesystem path from the author's Windows account in **six** places — one in the backslash
form and five in the forward-slash form, all of them the same repository root. Each is replaced by
the literal `<REPO>`. **Nothing else is altered**: no prose, no number, no code block, no line
break. Both hashes are recorded above, so anyone holding the circulated copy can reproduce the
substitution byte for byte and confirm that only those six strings moved. The circulated copy is
deliberately **not** committed, because the point of the redaction is that the account path does not
enter the repository. The freeze script checks the redacted copy for any surviving drive-letter path
and refuses to run if it finds one.

The redacted copy opens with a banner saying it is superseded, that it has known errors, that this
document governs and that the freeze has already acted on it — so a reader who opens the file
directly on GitHub cannot mistake a decided question for a live one. **Its errors are enumerated in
§10**, with the correct value for each.

### 2.4 What the freeze script refuses to run without

`scripts/freeze_stratum_list.py` checks size then sha256 of all four pinned files **before it parses
anything**, then the manifest's `generated_utc`, `census_version`, `value_filter` and row count, then
the two snapshots' own headers, aborting on the first mismatch. The hash guards against a file
changing; the header stamps guard against a pinned constant being edited to match some other file.
Either check alone can be walked around by one edit.

**Regenerating.** `python scripts/freeze_stratum_list.py` rewrites the two artifacts;
`--check` verifies **both** committed halves — JSON and CSV — against a fresh run without writing,
and in both modes it verifies this document's tables against what it has just derived. The output
carries no generation timestamp, no package versions and no environment, so regeneration is
byte-identical on any platform, and `tests/test_stratum_list_freeze.py` enforces that against the
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
(spec §10 risk 13, *stratum cherry-picking*). Stratum-level discretion is the cheapest way to
manufacture a GO: the manifest offers 1197 candidates across 68 datasets, cells-per-donor spans
three orders of magnitude, and λ_naive is expected to grow with cells per donor (spec §10 risk 2). A
rule that admitted per-stratum choice would let the study pick its own answer while every individual
choice looked defensible.

`gate_status == "candidate"` is `census_select`'s own verdict: the stratum cleared the obs-decidable
half of the spec §1 inclusion gate (≥ 3 donors per group after the thin-donor drop, donor present,
non-constant and nested within condition) and was not excluded by the spec §1 confound pre-screen.
**It is not an admission.** See §9.

**The "recommended strata" tables of the 2026-08-16 proposal document are not this set.** That
document listed four to seven stratum-contrasts per dataset as reading aids for a human deciding
which datasets to take. They are a strict subset of the 251, they were never a selection, and they
carry no status of any kind. Stated plainly here because the failure mode is obvious: a later reader
finds the smaller, tidier table, mistakes it for the pre-registration, and reports a study whose
stratum list was in fact chosen after the datasets were.

### 3.2 The twelve datasets

Ceiling = the largest `min(n_donors_A, n_donors_B)` over the dataset's candidate strata, i.e. the
best-powered design it contains. `≥ 8v8` and `3v3` count strata, not datasets. Every figure in this
table is computed from the manifest by the freeze script, **and the script parses the table back out
of this file and compares it cell by cell** — the ranks, the ids, all four counts and the disease
terms — so a hand-typed cell cannot disagree with the artifact in either direction.

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

**The twelve occupy twelve distinct collections**, verified against the pinned release table. That
is what makes the D2 denominator clean: no two entries of it share a cohort and a laboratory before
a single sibling has been considered.

Aggregates over the 251, all recomputed by the freeze script and pinned as `attested_figures` in the
artifact: 4 609 595 cells; 182 distinct `(dataset_id × cell_type)` strata carrying 251 binary
contrasts (a dataset with two disease terms contributes two contrasts against the same `normal`
group); 124 distinct cell-type labels; `residual_df` 4 … 108; `permutation_count` 20 …
7.28 × 10²³, of which **45 strata fall below 1000 and therefore require full enumeration**
(spec §4, "Small D → enumerate all").

**A correction, recorded rather than absorbed.** The proposal document circulated on 2026-08-16
states Rexach's envelope ceiling as `min(A,B) = 11`. **It is 10.** The manifest's best-powered Rexach
strata are A = 11 versus B = 10 (the progressive supranuclear palsy arm, six cell types), and no
Rexach control group anywhere in the dataset exceeds 10 donors. The 11 is `max(n_donors_A)` read as
if it were the ceiling; the two coincide for balanced designs and diverge for every skewed one. The
pre-registration carries **10**. Nothing downstream changes — Rexach was already outside the σ = 0.5
tier at either value — but a pre-registration that silently improved a number it had published
would be worth less than one that did not. It is the first of the discrepancies §10 lists.

**Where the 12 came from.** The 68 candidate-bearing datasets were read for coverage of the axes
spec §1 (iii) names, not for expected inflation. The rationale per dataset, the datasets considered
and rejected, and five named reserves are in the 2026-08-16 proposal document (§2.3), and **none of
its five reserves may be substituted for anything** (§9 item 7). What matters for auditing is stated
in §4 — the coverage claim is checkable against the frozen set itself.

### 3.3 The five collection-siblings — frozen as controls, and computed rather than named

**This is where the third review found a blocker, and the fix is structural.** The previous version
of this section named **two** siblings. There are **five**. The other three were missed because
collection membership was not derivable from any committed source: the manifest has no collection
column, and Discover was read live. Three datasets holding **79 further runnable strata** were
therefore excluded from the D2 denominator's neighbourhood without being frozen — the exact set a
later reader could have reached for after seeing the results — and two of them clear every tier of
the operating envelope. The sibling set is now **derived** on every run from the pinned release
table, and the declaration is compared against that derivation in both directions.

> **A sibling is any candidate-bearing dataset of the manifest that shares a `collection_id` with
> one of the twelve.** That is the whole definition; there is no list to maintain by hand.

| `dataset_id` | What it is | Sibling of | Frozen strata | Ceiling | Role |
|---|---|---|---|---|---|
| `c2876b1b-06d8-4d96-a56b-5304f815b99a` | SEA-AD, middle temporal gyrus | `6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3` (#1) | 18 | 42 | `within_collection_control` |
| `edc8d3fe-153c-4e3d-8be0-2108d30f8d70` | Yoshida 2022, airway (bronchus, nasal cavity, trachea) | `2a498ace-872a-4935-984b-1afa70fd9886` (#4) | 25 | 30 | `within_collection_control` |
| `8f4f8502-9170-4ac2-9707-3b6985ebfe5f` | CAREBANK right atrium | `f1606894-59df-4794-a37f-baa7c6fb6de1` (#6) | 11 | 15 | `within_collection_control` |
| `dea717d4-7bc0-4e46-950f-fd7e1cc8df7d` | KPMP adult human kidney **scRNA**-seq v1.5 | `a12ccb9b-4fbe-457d-8590-ac78053259ef` (#8) | 43 | 26 | `within_collection_control` |
| `1e5bd3b8-6a0e-4959-8d69-cafed30fe814` | Emphysema Cell Atlas, immune cells | `4b6af54a-4a21-46e0-bc8d-673c0561a836` (#12) | 9 | 3 | `within_collection_control` |

D2 clusters evidence by dataset because same-dataset strata share donors, batch and assay. Two
datasets from one collection share the **cohort and the laboratory** as well: SEA-AD MTG is the same
cohort and laboratory as SEA-AD DLPFC in a different cortical region; the Yoshida airway split is
the other compartment of one paediatric SARS-CoV-2 cohort; CAREBANK and PERIHEART are two cohorts of
one right-atrium collection at the same tissue, assay and suspension; the KPMP pair is the
single-cell and single-nucleus arm of one kidney atlas; and the emphysema immune split is the same
three-versus-three cohort's other compartment. Counting any of them toward "majority of independent
datasets" would inflate the effective n exactly the way D2 exists to prevent.

*Not* "the same donors", for the pair where that was checked: DLPFC's best design is 39 v 44 and
MTG's is 42 v 46, and MTG's `min(A,B)` runs 27 … 42 against DLPFC's 33 … 39 — the cohorts overlap,
the donor sets in the strata do not coincide, and *same cohort and laboratory* is the wording that
is true. **The other four pairs have not been checked for donor overlap and no claim is made about
them**; `obs` would decide it and the freeze does not load `obs`.

**They are frozen, by the same rule, and emitted with an explicit role.** Every candidate row of
theirs — 18 + 25 + 11 + 43 + 9 = **106 strata** — is in
[`stratum_list_2026-08-16.json`](../pilot/preregistration/stratum_list_2026-08-16.json) under
`within_collection_control_rows`, carrying `role = "within_collection_control"`, where the 251 carry
`role = "analysis_set"`. The CSV twin carries both blocks, **357 data rows**, told apart by that same
column.

**The scope of what freezing them buys, stated exactly.** It means that no candidate stratum of any
collection represented in the analysis set is left unlisted and selectable after the fact. It is
**not** a claim about the whole Census pin: 1197 candidate strata over 68 datasets exist at this
pin, and only these 357 are frozen. The earlier artifact note said the broader thing and was wrong.

Three rules attach to the control set and none of them is discretionary:

1. A result from a control stratum is **reported as a within-collection control**, labelled as such
   wherever it appears.
2. It **never enters the D2 denominator** and never counts toward a majority, whatever it shows.
3. **Promoting one to an independent dataset is an amendment** — dated, numbered, written before the
   change. An unnamed sibling can be promoted later by someone who did not know it was one; a named
   and frozen one cannot be promoted quietly at all. That protection is now worth four times what it
   was: three of the five are donor-rich, and two of those three clear all four envelope tiers.

They remain admissible as within-collection reproducibility controls, which is a useful thing to
have: a result that fails to reproduce between two regions of one cohort is informative, and it is
informative in a way that says nothing about independent replication.

---

## 4. Coverage against spec §1 (iii), re-derived

Spec §1: *"First pass = 8–12 datasets chosen to SPAN the outcome space (not cherry-pick wins): (i)
2–3 with a biologically strong expected effect (pseudobulk shown non-null), (ii) 2–3
subtle/low-effect, (iii) deliberate variation in assay (10x 3′ vs 5′), tissue, donor count (some
exactly 3v3, some ≥ 8v8), and cells-per-donor spanning the pre-registered bins (D1)."*

### 4.1 (i) and (ii) — the effect-size axis is a literature judgement, and is labelled as one

No effect in this table has been measured by us and none can be before the sweep runs. The axis
exists for one purpose: to show that the list is not picked for wins. It is **not** a prediction of
λ_naive and must never be scored against the results as though it were a hypothesis. It is recorded
per `(dataset, disease term)` rather than per dataset, because four datasets carry more than one
disease arm against **one shared control group**, and in two of them (#4, #5) those arms differ in
expected strength — the cheapest effect-size control available to us, and one a per-dataset label
destroys. The labels are frozen in the machine-readable artifact (`expected_effect`, one value per
row) so they cannot be reassigned after the fact, and the freeze script parses the Datasets and Arms
columns of the table below back out of this file and compares them against that artifact.

**Spec §1 (i)'s parenthetical is not satisfied by anything checked here, and that is a gap, not a
formality.** Spec §1 asks for datasets "with a biologically strong expected effect **(pseudobulk
shown non-null)**". The Basis column below gives biological expectation — what kind of disease this
is, in what tissue, against what control — and **not** evidence that donor-level pseudobulk was
shown non-null in those publications. No such check was performed for this freeze: the literature
was read for study design, not for the presence of a donor-aggregated differential-expression
result. So the *(pseudobulk shown non-null)* clause is **unmet**, and the strong / moderate / subtle
labels rest on biological expectation alone. It is recorded here rather than repaired, because
repairing it means a per-publication evidence review that this document did not do and must not
pretend to have done. The consequence is bounded and worth stating: the labels are a *coverage*
claim about the list's span, and a coverage claim that rests on expectation is weaker than one that
rests on published pseudobulk results — but it is not a claim about any result of ours, so nothing
downstream inherits an unearned number.

| Expected | Datasets | Arms | Basis — biological expectation from each dataset's own publication, **not** a pseudobulk result. **Judgement, not measurement.** |
|---|---|---|---|
| **strong** | 5 (#1, #2, #3, #4, #5) | #1 dementia; #2 Alzheimer disease; #2 Pick disease; #2 progressive supranuclear palsy; #3 COVID-19; #4 COVID-19; #5 COVID-19 | Cortical neurodegeneration and tauopathy with glial response (Gabitto 2024 `10.1038/s41593-024-01774-5`; Rexach 2024 `10.1016/j.cell.2024.08.019`); diffuse alveolar damage in autopsy lung (Melms 2021 `10.1038/s41586-021-03569-1`); acute severe viral infection in blood (Yoshida 2022 `10.1038/s41586-021-04345-x`; Ahern 2022 `10.1016/j.cell.2022.01.012`) |
| **moderate** | 3 (#5, #8, #12) | #5 influenza; #8 acute kidney failure; #8 chronic kidney disease; #12 pulmonary emphysema | Milder viral arm of the same atlas; injured proximal-tubule states, acute then chronic (KPMP v1.5, **DOI not established** — Discover returns `collection_doi = null`); pulmonary emphysema against an alveolar control (Wang 2023 `10.1016/j.immuni.2023.01.032`) — but see the note on #12 below, whose label is the weakest thing in this table |
| **subtle** | 6 (#4, #6, #7, #9, #10, #11) | #4 post-COVID-19 disorder; #6 atrial fibrillation; #7 rheumatoid arthritis; #9 clonal hematopoiesis; #10 opiate dependence; #11 Crohn disease | Clinically heterogeneous post-viral state; atrial remodelling against a cardiac-surgery control that is not a healthy heart (Linna-Kuosmanen 2024 `10.1016/j.xcrm.2024.101556`); subpopulation-level blood signatures (Binvignat 2024 `10.1172/jci.insight.178499`); a pre-clinical clonal state with no phenotype (Heimlich 2024 `10.1182/bloodadvances.2023011445`); cell-type-specific programmes rather than degeneration (Phan 2024 `10.1038/s41467-024-45165-7`); paediatric ileal mucosa (Elmentaite 2020 `10.1016/j.devcel.2020.11.010`) |

The DOIs in the Basis column are the `collection_doi` values of the pinned Discover snapshot,
recomputed and compared against the declaration on every freeze; the prose beside them is the
judgement.

**"moderate" is ours, not the spec's.** Spec §1 names two categories — (i) strong and (ii)
subtle/low-effect — and a three-valued vocabulary is a labelling convenience of this freeze, not a
spec category. It is frozen as `expected_effect_vocabulary` in the artifact so it cannot grow later.
What spec §1 actually demands is counted over its own two categories, and the middle tier is counted
toward **neither**: 5 datasets carry a strong arm and 6 carry a subtle arm, both above the spec's
2–3, without any moderate arm being borrowed to reach the floor.

**#12's row is a design justification wearing an effect label.** Wang 2023 is in the list to satisfy
spec §1 (iii)'s "some exactly 3v3" — all 8 of its strata are 3 v 3 — and its inclusion was decided
on that ground, not on an expected effect size. Its `moderate` label is therefore the one entry in
the table with the least behind it, and the honest reading is that #12 spans a *donor-count* axis
and was given the middle label because the vocabulary requires one, not because the emphysema
literature places it there. It is called out rather than quietly assigned so that no later reader
treats #12's result as bearing on the effect-size axis at all.

Spec §1 asks for 2–3 of each; the list carries 5 and 6. That is a floor being exceeded, not a bar
being moved, and the excess is not free padding: it comes from the within-dataset arms (#4 strong +
subtle, #5 strong + moderate, #8 acute + chronic, #2 three tauopathies), which cost no additional
dataset and hold assay, tissue, laboratory and control group fixed across the effect-size contrast.

One consequence is worth pre-registering as a stated expectation rather than discovered later:
**#9 (clonal haematopoiesis) is expected to behave as a no-effect stratum set.** Spec §4 / A1 says
the clean truth-zero guarantee "lives in oracle (b) and in no-effect strata"; if that role is
claimed after seeing the result it is worthless, so it is claimed now. **This is a stated
expectation, not a measurement.**

### 4.2 (iii) — assay, tissue, suspension, donor count, cells per donor

Assay, suspension and tissue are recomputed from the pinned Discover snapshot (§2.2), not from the
manifest, and the freeze aborts if any of them disagrees with the declaration. Everything else in
this table is recomputed from the manifest.

| Axis | Spec §1 requirement | Frozen set | Verdict |
|---|---|---|---|
| Assay 3′ | variation | 10 datasets: #1, #2, #3, #6, #7, #8, #9, #10, #11, #12 | met |
| Assay 5′ | variation | **2 datasets: #4 and #5, both 10x 5′ v1** | met, and thinly — see below |
| Assay chemistry within 3′ | not required | v3 in 9 datasets; v2 in #11 (pure) and #2 (mixed v2 + v3); Discover additionally lists `10x multiome` for #1 — see below | — |
| Tissue | variation | **6 organ systems**: brain (#1, #2, #10), lung (#3, #12), blood (#4, #5, #7, #9), heart (#6), kidney (#8), gut (#11) | met |
| Suspension | not named in spec §1; D3-adjacent | `nucleus` 6 (#1, #2, #3, #6, #8, #10) / `cell` 6 (#4, #5, #7, #9, #11, #12) | balanced |
| Donor count — "some exactly 3v3" | required | **10 strata**, in #12 (all 8) and #2 (2 T-cell strata) | met |
| Donor count — "some ≥ 8v8" | required | **150 strata** across 7 datasets (#1, #2, #4, #5, #6, #7, #8) | met |
| Donor count — range | — | `min(A,B)` from 3 to 39; 24 strata at `min(A,B) = 3`, of which 10 are exactly 3v3; largest design 39 v 44 | — |
| Cells per donor | "spanning the pre-registered bins" | group medians **11.0 … 6671.5**, all six bins occupied — §5 | met, against bins defined in §5 |
| Counts per cell | not in spec §1; D1-adjacent | group medians **284 … 56 841.5**, a factor of 200; the extremes are #7 (RA blood) and #10 (striatal medium spiny neurons), both 10x 3′ v3 | — |

**"Organ system" is our grouping, not Discover's.** The tissue *labels* are Discover's and are
machine-checked; collapsing "dorsolateral prefrontal cortex", "Brodmann (1909) area 4", "insular
cortex", "primary visual cortex", "caudate nucleus" and "putamen" into one row called "brain" is a
labelling convenience of this document. Nothing downstream reads it.

**The 5′ arm is the thinnest axis in the list, and it is thin in a way that matters.** Both 5′
datasets are PBMC, both are 10x 5′ v1, and both carry a COVID-19 arm. If the 3′-versus-5′ comparison
produces a difference, this list cannot separate "5′ chemistry" from "blood", "PBMC" or "acute viral
infection". Widening it was possible — the proposal names a third 5′ candidate — and was not done,
because that candidate uses *in-vitro-stimulated* PBMCs, which inserts a layer between disease and
transcriptome that seemed worse than the confound it would relieve. **That trade is a judgement**,
recorded here so it can be disagreed with, and it is a limitation of the frozen list, not of the
data.

**A second limitation, of the same kind.** The 3′/5′ split coincides almost exactly with the
suspension split among the donor-rich datasets, and blood is the only tissue represented by both
5′ datasets. Any assay-level statement from this list is therefore partially confounded with tissue
by construction. Spec §7 item 2 (depth-matching) and the permutation null do not fix this; only more
datasets would, and there are not more inside the 12-dataset ceiling spec §1 sets.

**The multiome entry is a Discover property, not a property of the frozen strata, and the two
disagree.** Discover lists `10x multiome` alongside `10x 3' v3` for SEA-AD DLPFC at the *dataset*
level. The manifest says something different about the *strata*: `confound_cramers_v["assay"]` is
**`null` on all 18 SEA-AD strata**, which is what `census_select` writes when assay is constant
within a stratum — there is no association to measure. Rexach (#2), whose declared v2 + v3 mix is
real, shows a **non-zero V on 24 of its 27 strata** (0.029 … 0.085), so the column does register a
within-stratum assay mix when there is one. Two readings survive: either the multiome cells are
absent from the 18 frozen strata (dropped upstream, or confined to cell types that did not clear the
gate), or they are present and the spec §1 confound pre-screen missed them. **`obs` cannot decide
between them**, so nothing here is asserted about the frozen set; the multiome value is recorded as
Discover's dataset-level metadata and is carried as such in the artifact's `assay` field. §9 item 7
lists this among the things only an X load can settle.

---

## 5. The cells-per-donor bins (D1), pre-registered here

`PHASE0_SPEC.md` references "the pre-registered bins" three times — decision rule item 2 ("as
cells-per-donor is subsampled down toward the floor of the pre-registered bins"), spec §1 (iii)
("cells-per-donor spanning the pre-registered bins (D1)") and spec §7 item 3 ("pre-register the
cells-per-donor bins") — and **nowhere states them numerically**. Neither does any amendment.
Strictly, spec §1 (iii) is currently unsatisfiable: there is no object to span.

**This is a gap in spec §1 and it is closed here, in the same act that freezes the list. It is not a
change to spec §1**: no threshold moves and nothing spec §1 defines is redefined. A number the spec
demanded and never supplied is supplied.

It is **not** true, however, that nothing downstream is affected, and an earlier draft of this
paragraph claimed it. Decision rule item 2 requires the naive effect to persist "as cells-per-donor
is subsampled down toward **the floor of the pre-registered bins**", so the floor is an argument to
the GO/NO-GO rule and setting it at 10 parameterises that rule. The direction is the unhelpful one:
a **higher** floor would have made item 2 easier, because the effect would only have to persist down
to a shallower subsampling depth, where λ_naive is expected to be larger (spec §10 risk 2). Choosing
10 — the lowest value the spec §1 inclusion gate permits — is therefore the strictest available
choice and not a neutral one, and it is fixed here before any curve is computed. Anyone who thinks
the floor should be higher is asking for a weaker item 2 and must say so in an amendment.

> **The pre-registered cells-per-donor bins are the half-decade log bins**
> `[10, 30) [30, 100) [100, 300) [300, 1000) [1000, 3000) [3000, ∞)`.

The lower edge is 10 because that is the spec §1 inclusion gate's own floor (`gate_config.MIN_CELLS`;
its item 2 drops donors below it), so no admitted donor can fall outside the bins from below. Each
step is a factor of about √10 — **two closed bins per decade** — giving five closed bins over the
two and a half decades from 10 to 3000, plus an open top bin. The top bin is open because the
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

**All six bins are occupied** over 502 group medians spanning **11.0 … 6671.5** — a factor of more
than 600.
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
subsets that happen to be convenient. The thresholds are read from `gate_config.OPERATING_ENVELOPE`,
not restated as literals, and the freeze script verifies the declared membership against the
manifest. **All four** envelope tiers are declared, not only the two whose failure was obvious — an
earlier draft declared 0.5 and 0.7 alone, and the two undeclared rows of this table were consequently
published without ever being checked against the manifest.

**Both columns of the table below, defined, because only one of them used to be.**

* **Surviving datasets** — of the twelve, how many hold **at least one** frozen stratum with at
  least the tier's donors-per-group in **both** groups. A dataset-level any-stratum test.
* **Surviving strata** — of the 251, how many **individually** have
  `min(n_donors_A, n_donors_B)` at least the tier's threshold. This is *not* "the strata of the
  surviving datasets": at ≥ 4 v 4 the eleven surviving datasets hold 243 strata between them and the
  column reads 227, because sixteen of those strata are themselves below the threshold.
* **Against the floor** is not typed. It is `below_spec_dataset_floor` from
  [`stratum_list_2026-08-16.json`](../pilot/preregistration/stratum_list_2026-08-16.json), which the
  freeze script computes as `n_datasets < SPEC_DATASET_FLOOR` with `SPEC_DATASET_FLOOR = 8`.

| σ_donor anchor lands at | Envelope demands | Surviving datasets | Surviving strata | Against spec §1's 8–12 floor |
|---|---|---|---|---|
| ≈ 0.2 | ≥ 4 v 4 | 11 of 12 (all but #12) | 227 | met |
| **≈ 0.35** | **≥ 8 v 8** | **7 of 12 — #1, #2, #4, #5, #6, #7, #8** | **150** | **BELOW** |
| **≈ 0.5** | **≥ 13 v 13** | **5 of 12 — #1, #4, #6, #7, #8** | **94** | **BELOW** |
| ≈ 0.7 | ≥ 23 v 23 | **3 of 12 — #1, #6, #8** | 30 | **BELOW** |

**Three of the four rows are the ones this section exists for. The list falls below the spec's own
"8–12 datasets" floor at every tier except the most optimistic.** Seven is below eight, and an
earlier draft of this table called the σ ≈ 0.35 row *met, with no margin*, which was simply wrong
against the rule the freeze script has always encoded.

**The σ ≈ 0.35 failure is the one that matters, because 0.35 is not a pessimistic scenario.** It is
`gate_config.POWER_EVAL_SIGMA` — the envelope boundary that [Amendment 3](AMENDMENTS.md) Change 1(b)
makes binding, the point at which the shipped `ebayes` arm is grid-shown to clear the pre-registered
power target (0.793 at 8 v 8, calibrated), and therefore the instrument's own nominal operating
point. If the anchor lands there, the study is already below the spec's dataset floor: **7 of 12, and
it gets worse from there — 5 at σ ≈ 0.5 and 3 at σ ≈ 0.7.** Only σ ≈ 0.2, the most optimistic tier in
the envelope, keeps 11. That is not a tail risk being noted; it is the expected case failing, and
softening it to "no margin" would have been the single most consequential piece of self-flattery
available in this document.

**And the shortfall is ours, not the data's.** The claim that no other list from this manifest could
have done better is false, and the artifact committed beside this document refutes it. Re-derived
from `census_candidates_run31910799023_2026-08-15.json` over its **68 candidate-bearing datasets in
50 distinct collections**, with collection membership taken from the pinned release table:

| Envelope demand | Candidate-bearing datasets clearing it | In how many collections | Retained by the frozen twelve |
|---|---|---|---|
| ≥ 4 v 4 | 62 | 46 | 11 |
| ≥ 8 v 8 | **33** | **25** | 7 |
| ≥ 13 v 13 | **21** | **15** | 5 |
| ≥ 23 v 23 | **12** | **10** | 3 |

The last column is the "Surviving datasets" series of the table above, repeated here so the two can
be read against each other; the collections column is the honest unit, because two datasets of one
collection are not two independent choices and a dataset count silently offers slots D2 would refuse.

Twenty-one is not a handful, and **15 collections clear 13 v 13** against the frozen list's five
datasets. A list selected on donor counts would therefore have survived the truncation far better
than one selected for spec §1 (iii)'s coverage — and there is no need to construct a witness to say
so, because the counts say it. That trade was made deliberately, and it is recorded now so it cannot
be re-made later:

> **The truncation is a consequence of selecting for spec §1 (iii)'s coverage axes — tissue, assay
> chemistry, suspension, effect-size range and the mandatory 3v3 anchor — rather than for donor
> counts. A list optimised for donor counts would have survived the truncation and would have
> spanned less.**

Part of the cost is nameable: spec §1 (iii) requires "some exactly 3v3", and **only 3 of the 68
candidate-bearing datasets, in 2 collections, hold an exactly-3v3 stratum at all** — Rexach (#2, 2 of
its 27 strata), Wang (#12, all 8) and the Emphysema immune sibling (all 9), the last two being the
two halves of one collection. **None of the three clears 13 v 13**, so the anchor costs a slot at
both hard tiers. Rexach does clear 8 v 8, so it costs nothing at the two easy ones.

**An earlier version of this section built an explicit twelve-dataset witness, and it is deleted
rather than repaired.** It was meant to show the counterfactual attainable. It contained both SEA-AD
datasets and both KPMP datasets, so it double-counted two collections; its headline ≥ 23v23 figure
was 10 and it claimed 11; and it leaned on the Human Lung Cell Atlas, which owns all 12
`excluded_confound` rows in the manifest and whose donor-rich strata run assay Cramér's V up to 0.96
— a dataset the shipped covariate-less arm cannot analyse. The two tables above say everything the
section needs, are recomputed on every run, and contain nothing anyone can select from.

Two things follow. First, the feasibility risk recorded here is a **cost of the coverage
requirement**, and if the anchor lands at 0.35 or worse the honest report is "the spec's span
requirement and the spec's dataset floor are not jointly satisfiable at this Census pin under this
envelope" — a conflict inside the spec, not a shortage of public data. Second, and this is why it is
written here rather than discovered later: we know this **before the anchor exists**. Saying it now
is precisely what makes it impossible to reassemble the list on donor counts afterwards and present
that as the original plan. Nothing above may be selected from (§9 item 7), and **every cell of both
tables — 20 in the first, 16 in the second — is recomputed by the freeze script on every run and
compared against this file**, which is what stops this section from drifting back into the
comfortable version.

Amendment 3 named this outcome in advance and the wording is quoted rather than paraphrased
(**emphasis added**), because it is the commitment being honoured:

> **Whether the real sweep is feasible at all.** If real strata cluster at `sigma_donor` ≈ 0.5–0.7,
> the envelope admits them only at ≥ 13–23 donors per group. Whether enough CELLxGENE strata clear
> that is an open empirical question, **and a negative answer is a live outcome of this study, not a
> failure mode to be designed around.**

Two things this table is **not**. It is not a measurement: `sigma_donor` is unmeasured for every
stratum here, so these are donor-count arithmetic under a hypothetical σ, i.e. a scenario analysis.
And it is not an admission at any tier — Amendment 3's conversion from the moderated fit's `s0²` to
`donor_sigma` is unvalidated, and "no stratum may be admitted to the real sweep on the strength of
this entry alone".

> **Retracted, and marked here rather than left to be discovered a thousand lines away.** Amendment
> 3 called that conversion an **upper bound**, and [Amendment 4](AMENDMENTS.md) Part A, Correction 1
> withdraws the claim: `sqrt(s0²)·ln 2` is an audit quantity of **unknown error sign**, and the
> direction of the error is the dangerous one — it can fall *below* the σ it was supposed to bound,
> so a stratum can look inside the envelope when it is not. Nothing in this document was admitted on
> the strength of it, because nothing here is admitted at all.

Writing this down is the entire mechanism. Without it, an anchor at 0.35 would leave us free to
report "the seven datasets that were inside the envelope" as though seven had always been the plan.
With it, seven is visibly a *truncation of twelve*, it is visibly below the spec's floor, and it is
reported as a limitation of the instrument's applicability domain.

---

## 7. Gating already fixed by the spec, restated where it will be applied

These rules are the spec's, restated because they bind at exactly the points this list is unusual,
and a rule that lives only in a document nobody re-reads at analysis time is not a control. One
thing here **is** new and is flagged rather than smuggled in: spec §4 says only "small D → enumerate
all, report the exact null distribution, flag coarse resolution" and never says what *small* is. The
threshold of **fewer than 1000 distinct label assignments** below is this document's
operationalisation of that looser text. It is pre-registered here, before any null is computed, and
it is chosen to coincide with spec §4's own `n_perm = 1000` for the sampled case — below it, sampling
would draw more permutations than exist — rather than being fitted to how many strata it catches.

**A1 — the kill-switch is decision-relevant only at ≥ 8 v 8, and is leak-contaminated at 3 v 3.**
The signal-above-floor ratio (decision rule item 4, spec §6) is treated as decision-relevant **only**
in strata with at least 8 donors per group, "where balanced permutations are ~orthogonal to the true
grouping"; on 3 v 3 strata it is "flagged leak-contaminated and excluded from GO logic". In this
list that partitions the 251 concretely: **150 strata across 7 datasets carry the kill-switch and
101 do not.** Of the 101, **24 have `min(A,B) = 3`** — among them the 10 that are exactly 3 v 3, the
case A1 names as leak-contaminated — and **77 sit between 4 v 4 and 7 v 7**, below A1's threshold
and therefore also outside the kill-switch, which is the reading A1's text forces even though it
names only the 3 v 3 case explicitly. Dataset #12 exists to satisfy spec §1's "some exactly 3v3" and
contributes nothing to item 4 by construction; that is its declared role, not a disappointment to
be discovered later.

The same boundary bites on permutation resolution. **45 of the 251 strata have fewer than 1000
distinct label assignments** (minimum 20, at 3 v 3, which spec §4 reduces to 10 after removing the
identity and its complement), so for those spec §4's "enumerate all, report the exact null
distribution, flag coarse resolution" applies and the null's resolution must be printed next to
every number derived from it.

**D2 — evidence clusters by dataset, not by stratum.** The majority in decision rule item 2 is over
**independent datasets**, and the denominator is the 12 of §3.2 — not 251, and not 17. Same-dataset
strata are not independent replicates: #4 alone would contribute 47 of the 251, and #1, #4 and #8
together contribute 102. Per-stratum numbers are reported alongside, never in place of, per-dataset
aggregates, and the five siblings of §3.3 stay outside the denominator whatever they show.

**No pooled headline across incompatible granularity (D5).** Restated in §9 because it is a
consequence of this particular list rather than a general rule.

---

## 8. Spec §8(d) — the Mathys 2019 oracle. Not a member of this list, and not weakened by it

**Mathys 2019 is not in the stratum list and must not be.** Spec §8(d) makes it the *binding* real
anchor against the possibility that our simulators are optimistic: it is an oracle with a known
target, not a stratum whose inflation we are measuring. Putting it in the analysis set would convert
the check into one of the things being checked.

**The Census path is closed, and this is now read out of the pinned release itself.** Mathys 2019 is
absent from CELLxGENE Discover, and it is absent from Census `2025-01-30`. Those are two claims and
the second does not follow from the first, which an earlier draft of this paragraph assumed it did.

*The release, enumerated directly.* The release's own `census_info/datasets` dataframe — the object
`cellxgene_census` exposes as `census["census_info"]["datasets"]` — is an ordinary public TileDB
array under `cell-census/2025-01-30/soma/`, and it is committed here in reduced form (§2.2). It
holds **1573 datasets**, and all 73 of the manifest's dataset ids are among them, which is the
control that it is the right object. Searched over every retained field of every row for `mathys`,
`rosmap`, `religious order`, `memory and aging` and the DOI stems `s41586-019-1195` and
`10.1038/s41586-019-1195-2`: **zero hits, all six needles.** The freeze script re-runs that search on
every freeze and aborts on a hit.

*Why the earlier route did not need to be taken.* An earlier draft reported that
`census["census_info"]["datasets"]` was **unreachable** on this machine and enumerated the release
bucket's `h5ad` basenames instead. That sentence is withdrawn: it is false. `cellxgene-census`
requires `tiledbsoma`, which has published **no Windows artifact across any of its 52 PyPI
releases** — but `tiledb` (TileDB-Py) publishes Windows wheels at every recent version, and the
array is public. The dataset table was read directly with it. The h5ad enumeration gave the same
1573 and is no longer relied on.

*The Discover half, and why it is not enough on its own.* The full curation index
(`GET /curation/v1/datasets`, **2216 datasets**, read 2026-08-16, committed and pinned in §2.2) was
searched for the same six needles over `title`, `collection_name`, `collection_doi`,
`collection_doi_label` and `citation` — the method of Appendix A.3 of the proposal document —
yielding **zero hits**, and that search is recomputable from the committed bytes and is recomputed
on every freeze. A second pass over each record's *entire* JSON, before reduction, also returned zero
hits; **that pass is a read-time observation** recorded in the snapshot's header against the sha256
of the raw upstream payload, and it is not recomputable from the reduced file, because the reduction
dropped the fields it would have to search.

*The gap that makes Discover insufficient, measured.* **1567 of the release's 1573 datasets resolve
in Discover by `dataset_id`; 6 do not.** So the inference "absent from Discover ⇒ absent from Census
`2025-01-30`" is **false in general**, and the earlier sentence walked straight through it. The
control that draft offered — all 73 manifest dataset ids resolve in Discover — cannot detect the gap
by construction: those 73 ids are the ones the Census reader already returned, so they are guaranteed
to be in both. The blind spot is six datasets in two collections: four lamina-propria / submucosa
sorts of collection `0c3f148e-02ff-4c81-8946-29beaaf5fa59` (`10.1101/2021.03.28.437379`) and the
cross-species pancreatic alpha- and beta-cell maps of collection
`0a77d4c0-d5d0-40f0-aa1a-5e1429bcbd7e` (`10.1016/j.molmet.2022.101595`). Their ids and their two
collections are declared in the freeze script and recomputed from the two snapshots, and **none is
Mathys, ROSMAP or prefrontal cortex**. The conclusion now rests on the release's own contents rather
than on an inference from an index.

*The operational definition, stated because the two indexes are not the same object.* "Resolves in
Discover" means: the release row's `dataset_id` is a `dataset_id` in the Discover snapshot. Matching
on `dataset_version_id` instead gives **0 of 1573** — every dataset in the release has been
re-versioned since — and `collection_doi` differs for 61 of the 1567 that do resolve. Both reads are
dated 2026-08-16; the release they describe is dated 2025-01-30.

Spec §8(d)'s own parenthesis — "via Census if present, else Synapse syn18485175" — anticipated the
absence; what follows from it did not get anticipated, and is recorded here.

**What the remaining path costs, stated before it is walked.**

1. **A data-use agreement, on the critical path.** Synapse `syn18485175` is ROSMAP data and requires
   a ROSMAP DUA. This is a legal dependency, not a technical one, and it gates the binding check.
   The application was begun on **2026-08-16**, the date of this freeze; **access had not been
   granted at the time of writing**, and no data from the deposit has been seen. Elapsed time to
   access is a schedule risk with no engineering mitigation, and it is recorded here so that a long
   wait is a known cost rather than a discovery.
2. **A second loader, which is work not in the spec §9 plan.** `io_counts.load_stratum(dataset_id,
   cell_type)` is Census-shaped throughout: it wraps `cellxgene_census.get_anndata` on the pinned
   version's raw layer and keys on `dataset_id`. A Synapse-hosted matrix shares none of that. A
   second loading path is therefore **required work that the spec's implementation plan does not
   contain**, and it must carry the same integrality gate (inclusion-gate item 4), the same
   thin-donor filter and the same frozen-universe check as the Census path, or the anchor is not
   running the same pipeline and proves nothing. **Recorded as a public description and not
   verified here**: the resource is ~80 660 nuclei from prefrontal cortex across 48 ROSMAP donors,
   24 with AD pathology and 24 without, with a filtered matrix at `syn18681734`.
3. **The target is fixed and external.** Spec §8(d) requires qualitative reproduction of what it
   calls "Murphy & Skene 2023": naive per-cell DE grossly exceeding pseudobulk, with the
   permutation-null floor accounting for most naive calls. **The citation, verified against the
   article and recorded as an external fact:** Murphy AE, Fancy NN, Skene NG (2023), *Avoiding false
   discoveries in single-cell RNA-seq by revisiting the first Alzheimer's disease dataset*, **eLife
   12:RP90214**, `10.7554/eLife.90214.3`. It has three authors, not two; the spec's shorter name is
   left as the spec wrote it and is not a different paper. That article reports **549 times fewer
   DEGs at FDR 0.05** — 14 274 cell-level calls against 26 pseudobulk calls, exactly 549.0 — and it
   attributes the reduction to its **whole corrected re-analysis**, corrected quality control *and*
   pseudobulk aggregation together, not to the change of replication unit alone. Attributing it to
   pseudobulk by itself, as an earlier draft did, overstates what the source claims and would set our
   anchor against a number nobody measured. It is **quoted from that publication and not re-derived
   in this repository**, and it is context rather than a threshold, since spec §8(d) asks for a
   qualitative pattern. One further limit, recorded so the anchor is not over-sold: the article's own
   permutation experiment is reported qualitatively in its text ("consistently found high numbers of
   DEGs"), so the *quantitative* half of spec §8(d)'s criterion — that the null floor accounts for
   most naive calls — is the spec's own bar and not a figure this publication states.

**This oracle does not block the freeze, and the freeze does not weaken the oracle.** Making the
stratum list wait on a DUA would stall the pre-registration behind a legal process while the
candidate manifest sat un-frozen and readable — the precise condition under which a list stops being
pre-registered. Conversely, nothing here downgrades spec §8(d): it remains **BINDING**, it remains
**unrun**, and no result from the 251 strata may be reported as validated by anything except the
oracle spec §8(d) names. If the DUA is refused, that is an amendment, and the amendment will have to
state what the study's binding real-data check is instead — not quietly drop the requirement.

---

## 9. What this freeze does NOT settle

Written in the amendments' own genre, and not as a formality: every item below is a live hole and
several of them can still shrink the list.

1. **Admission is still closed, and this document does not open it.** All four blockers stand on all
   2190 rows of the source manifest and on all 357 rows here: `integer_check`,
   `frozen_universe_size`, `sigma_donor_estimate`, `envelope_membership`. `admitted_to_sweep` is
   `False` everywhere. Freezing *what will be measured* is upstream of deciding *what may be
   reported*, and only the first has happened.

2. **`sigma_donor` is unanchored, so envelope membership is unknown for every stratum.** Amendment 3
   supplied the envelope and a *mechanism* (`sqrt(s0²) · ln 2` from the moderated fit), stated that
   the conversion is unvalidated against the simulator where truth is known, and recorded that "the
   `sigma_donor` anchoring demanded by Amendment 1 therefore remains OPEN". §6's tiers are arithmetic
   under a hypothetical σ, not measurements, and no stratum here is inside the envelope until that
   work is done.

   > **Retracted.** Amendment 3 also called that quantity an **upper bound**, and
   > [Amendment 4](AMENDMENTS.md) Part A, Correction 1 withdraws the claim: it is an audit quantity
   > of **unknown error sign**, and the direction of the error is the dangerous one — it can sit
   > below the σ it was supposed to bound, so a stratum can appear inside the envelope when it is
   > not.

3. **`pooled` is `unresolved` on 1197 of 1197 candidates, and on all 357 frozen rows.** The pinned
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
   the spec §1 pre-screen did not catch them — not because the confound is absent but because it is
   not in the column being screened. Whether region enters the *definition* of a stratum is **not
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
   section, and the option is **unavailable** until some amendment supplies a design that has one.

6. **Sex, developmental stage and self-reported ethnicity are outside the spec §1 confound screen.**
   Spec §1 screens exactly five covariates — assay, suspension type, `tissue_general`, a
   sequencing-depth bin and a pool id — and these three are materialised into `obs` but never tested.
   Two datasets make the omission concrete: Yoshida (#4) mixes children and adults in one cohort, and
   Elmentaite (#11) is entirely paediatric at 4–14 years, where developmental stage is an obvious
   candidate confounder with disease. **This is a hole in the spec, not in the manifest**, and
   closing it is a spec change that would need an amendment; it is named here so that it is a known
   omission rather than an accident of the frozen list.

7. **`integer_check` and `frozen_universe_size` are unverified for all 17 datasets, so the list may
   still shrink.** Both are computed at X load and neither can be decided from `obs`
   (spec §1 items 4 and 5; spec §10 risk 6). A dataset whose raw layer is not integral is **dropped
   with its reason, never rounded**; a stratum whose frozen universe falls below 200 genes is a SKIP
   whose measured size is still recorded. The most likely casualty is identifiable in advance and is
   named now rather than after the fact: **#7 (rheumatoid arthritis blood) runs at 284–803 median
   counts per cell** against the frozen set's median of 3081 — roughly an order of magnitude thinner
   — and its smaller strata may not assemble a 200-gene universe. (Both figures are per-group
   medians; 3081 is the median of the 502 per-group medians of the analysis set, two per stratum, and
   the population has to be named because the per-stratum mean of the two groups gives 3126.75
   instead.) That is a **stated expectation**, not a reason to pre-emptively drop it: a dataset
   excluded now on a suspicion would be exactly the discretion §3.1 exists to remove. **Any such
   shrinkage is a reported outcome, never a re-selection**: the dropped stratum is reported with its
   reason and its measured value, and no replacement is substituted for it. Replacing a stratum that
   failed a counts gate with one that did not is selection on the data, arriving through the back
   door — and that prohibition binds the five reserves named in the proposal document (§2.3) exactly
   as it binds anything else. They are committed so that a substitution could be *detected*, not so
   that one becomes available.

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
   (V ≥ 0.8), and 24 a partial `assay` tag; none carries a spec §1 exclusion, and none of the 12
   datasets lost a stratum to the confound pre-screen (all 55 non-candidate contrasts among their 306
   failed the *inclusion gate*, mostly on donor counts). Spec §1 would carry a tagged covariate into
   the design "only if C4's df rule allows", but C4 lapsed with DESeq2 (Amendment 2 Change 1) and the
   shipped moderated arm fits `~ 1 + x` with **no covariate slot**. A partially confounded stratum
   here is therefore neutralised **only** by the permutation null (spec §7 item 1, which is where
   "neutralized by the permutation null" is written), which is what the tag means and all it means.
   Depth in particular is not modelled anywhere.

9. **D4's excluded fractions are not whole-Census fractions.** The manifest's `counts` block reports
   an excluded fraction of 0.4534 over its 2190 contrasts, but those contrasts come only from the
   datasets pass 2 actually read; datasets the coarse pass-1 filter judged incapable of holding a
   gate-clearing stratum were never queried and appear in no denominator. Within the 12, the split is
   251 candidates against 55 inclusion-gate failures out of 306 contrasts. Any D4 statement must
   carry that scope.

10. **The freeze fixes the list, not the analysis.** `controls.py`, `decision.py` and `report.py`
    (spec §9 items 9, 12, 13) do not exist yet; the cells-per-donor sweep that D1 makes the primary
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

13. **The Discover snapshot dates from after the release it describes.** §2.2 gives the measurement:
    0 of 1573 `dataset_version_id` match, 61 `collection_doi` differ, 6 datasets have no Discover
    record. Every assay, suspension, tissue and DOI value in §4 is therefore a 2026 Discover fact
    about a dataset id, not a property read out of the 2025-01-30 release, and the two KPMP
    `is_primary_data = [False]` records (§2.2) are an unresolved disagreement between them. Only
    reading `obs` at the pin would settle any of it, and this freeze does not read `obs`.

---

## 10. Author attestation

**What is machine-checked.** `scripts/freeze_stratum_list.py` recomputes, from four hash-pinned
files, and aborts on any disagreement: the sha256 sums and byte counts of all four; the manifest's
`generated_utc`, `census_version`, `value_filter`, row count and gate-status split
(2190 = 1197 + 981 + 12); the membership of all twelve datasets and all five siblings by full uuid;
the twelve distinct collections; the sibling set itself, derived from the pinned release table rather
than declared; the 251 frozen strata and their per-dataset counts and ceilings; the 106
within-collection control strata and their per-sibling counts and ceilings; every cell of §3.2's
table, §3.3's table, §4.1's Datasets and Arms columns, §5's occupancy table and both of §6's tables,
parsed back out of this file; the assay, suspension, tissue and DOI of all twelve against the
Discover snapshot; all four Layer B subsets and the `below_spec_dataset_floor` verdict of each; the
tier census in datasets and in collections; the exactly-3v3 census in datasets and in collections;
`admitted_to_sweep = False` on all 2190 rows; and the Mathys needle search over both snapshots.

**The eleven figures §10 used only to attest are recomputed too.** They were correct when a reviewer
checked them by hand, which is precisely the property that could not be relied on again: they are
now `attested_figures` in the artifact and the freeze aborts if any of them moves. They are
4 609 595 cells; 182 distinct `(dataset_id × cell_type)` strata; 124 cell-type labels; `residual_df`
4 … 108; `permutation_count` 20 … 7.28 × 10²³; per-donor cell counts 10 … 16 383; the 3081 median
counts per cell (over the 502 per-group medians, as §9 item 7 states); Rexach's assay V non-zero on
24 of 27 strata at 0.029 … 0.085; 306 contrasts inside the twelve, 251 candidates and 55
inclusion-gate failures; KPMP's 6 strata at ≥ 23v23; and "every dataset occupies at least four bins,
three occupy all six". Alongside them: the counts-per-cell range 284 … 56 841.5, the six-bin
occupancy over 502 group medians spanning 11.0 … 6671.5, the confound-tag counts 231 / 8 / 24, the
45 strata below 1000 permutations, the 24 / 77 split below 8 v 8, the 10 exactly-3v3 strata,
SEA-AD's null `confound_cramers_v["assay"]` on all 18 strata, and the sibling donor spans
(DLPFC 33 … 39, MTG 27 … 42) that make "the same donors" false for that pair.

`tests/test_stratum_list_freeze.py` re-derives the load-bearing figures from the raw manifest rather
than through the freeze script, so a bug in it cannot be reproduced by the test meant to catch it.

**Judgement, not measurement.** The choice of these twelve datasets from 68 candidate-bearing ones;
the strong / moderate / subtle labels of §4.1 and the pre-declaration of #9 as an expected no-effect
anchor; the grouping of Discover's tissue labels into six "organ systems" in §4.2; the decision to
accept a 5′ arm that is confounded with tissue rather than take an in-vitro-stimulated alternative;
the decision to keep #1 and #2 despite their D5 incompatibility; the choice of half-decade log bins
in §5; and the `role_note` prose describing each sibling relationship in §3.3. The bins were fixed
before their occupancy was computed; the effect-size labels were fixed before any metric was
computed; neither is a prediction of any result, and the freeze script pins both so that neither can
be revised after one.

**Quotation, not measurement.** Amendment 3's "a negative answer is a live outcome of this study"
(§6, emphasis added); Amendment 2's disclosure in §1; spec §1's first-pass sentence in §4; A1's
kill-switch wording and spec §4's "small D → enumerate all" in §7; spec §8(d)'s parenthesis in §8.

**External fact, not measurement.** Murphy, Fancy & Skene 2023, eLife 12:RP90214,
`10.7554/eLife.90214.3`, and its 549-fold figure (14 274 → 26 DEGs at FDR 0.05), verified against the
published article and attributed by it to corrected quality control *and* pseudobulk together; the
public description of the Mathys deposit (~80 660 nuclei, 48 donors, 24 versus 24), which was **not**
verified against Synapse, since access requires the DUA that is still pending; that `tiledbsoma` has
published no Windows artifact across its 52 PyPI releases while `tiledb` publishes Windows wheels;
and every value taken from the Discover snapshot, which is a 2026 read and not a property of the
2025-01-30 release (§2.2, §9 item 13).

**Read-time observation, not recomputable from the committed bytes.** The Mathys needle search over
*every* field of *every* Discover record before reduction (zero hits): the reduced snapshot no longer
carries the fields that pass would search. It is recorded in the snapshot's header against the sha256
of the raw upstream payload. The search over the retained fields is recomputable and is recomputed.

**Could not be established.** The publication DOI of dataset #8 (KPMP kidney v1.5): Discover returns
`collection_doi = null` for that collection, so the provenance is recorded as *not established*
rather than inferred. Whether SEA-AD's `10x multiome` cells are inside the 18 frozen strata is
undecidable from `obs` and is deferred to X load (§4.2, §9 item 7). Whether the four sibling pairs
other than SEA-AD share donors (§3.3). Which of Discover and the pinned release is right about the
two KPMP `is_primary_data` records (§2.2).

**Asserted on expectation, not evidence.** Spec §1 (i)'s parenthetical *(pseudobulk shown non-null)*
is **not satisfied** by anything checked for this freeze: the strong / moderate / subtle labels rest
on the biology of each study's design and not on published donor-level pseudobulk results, no such
per-publication review was performed, and §4.1 records this as a known gap in what spec §1 asked for.
The expectation that #7 may fail the frozen-universe gate (§9 item 7) is of the same kind.

**Not claimed.** That any stratum in this list is inside the operating envelope; that any of them
will survive the counts gate; that donor pseudobulk is calibrated on any of them (D3 forbids the
claim while pooling is unresolved); or that the strong / subtle labels will match the measured
inflation. No metric has been computed on any of these 251 strata.

### 10.1 The discrepancies established in the proposal document

§2.3 says the proposal has known errors. This is the list, and it is not a count taken on trust: each
entry was established by recomputation against the pinned manifest, the pinned snapshots or the
published article. The **wrong** value in each row is a quotation from a file pinned by sha256; the
**right** value is, wherever it is a number over a pinned source, one the freeze script recomputes.
An earlier draft of §2 said the proposal contained "two figures this document corrects". It contains
at least these twenty-seven, and four of the first seven are defects of reasoning rather than
arithmetic.

**Reasoning, not arithmetic.**

| # | What the proposal asserts | What is true |
|---|---|---|
| R1 | Mathys is absent from Discover **therefore** absent from Census `2025-01-30`; and the control that "all 73 manifest ids resolve in Discover" shows nothing was lost | The inference is invalid — 6 of the release's 1573 datasets have no Discover record (§8). The control is circular: those 73 ids are the ones the Census reader returned. The *conclusion* survives on other evidence |
| R2 | The 549-fold reduction is what re-analysing Mathys **by pseudobulk** gave | The article attributes it to the whole corrected re-analysis — corrected quality control *and* pseudobulk together (§8 item 3) |
| R3 | 7 of 12 datasets inside the envelope at σ ≈ 0.35 is "workable, but with no margin" | 7 < 8, the spec §1 floor. The governing verdict is **BELOW** (§6). The same sentence re-bases the D2 majority onto the 7 survivors; the denominator is the 12 (§7) |
| R4 | SEA-AD's absent `assay` confound flag shows "assay does not separate the conditions — confirmed, not assumed" | The value is `null`, which is what `census_select` writes when assay is *constant* within the stratum. Nothing was measured; `obs` cannot decide it (§4.2) |
| R5 | Yoshida's recommended COVID strata clear "0 of them under a strict ≥ 13v13", while stating their `min` is 13 | The envelope demands `≥ 13`, so all three clear it. The sentence contradicts itself |
| R6 | Depth or region can be carried into the design as a covariate if C4's `df ≥ 3` rule allows (five separate places) | C4 lapsed with DESeq2 (Amendment 2 Change 1) and the shipped arm fits `~ 1 + x` with no covariate slot at all (§9 item 8). Degrees of freedom are not the binding constraint; the absence of the slot is |
| R7 | The manifest is `census_candidates_full.json`, and "all commands are reproducible from scratch" | No such file exists in the repository; the artifact is `census_candidates_run31910799023_2026-08-15.json`. The proposal gives no hash or byte count for any evidence file |
| R8 | Three passages are quoted with bold emphasis (the manifest's `operating_envelope_source`, spec §8(d), Amendment 3) | None of the three sources carries that emphasis, and none of the three quotations says "emphasis added" |
| R9 | Two same-collection siblings, "+2 сиблинга" in the header | Five. The proposal's own body names two of the three it omits from its header (CAREBANK, KPMP scRNA) without drawing the D2 consequence — which is the defect §3.3 now fixes |
| R10 | The σ_donor anchor has been left OPEN "for the third amendment running"; two appendix expectations about `AMENDMENTS.md` | True when written, stale now: Amendment 4 Part A (2026-08-16) supplies an estimator and a membership rule |

**Arithmetic.**

| # | What the proposal states | What the pinned sources give |
|---|---|---|
| N1 | Rexach's envelope ceiling `min(A,B) = 11` (twice) | **10** — `max(n_donors_A)` read as the ceiling (§3.2) |
| N2 | Cells per donor "18 … 6071" (twice) | **11.0 … 6671.5** over the 502 per-group medians. The proposal's range is the group-**A** medians of its own recommended-strata tables, with the larger group-B value ignored beside it |
| N3 | Rexach has "three 3v3 T-cell strata" | **Two.** The third is 4 v 3 — visible in the same sentence's own "perm = 20/35". Its coverage table separately omits Rexach's two 3v3 strata altogether |
| N4 | "No other list from this manifest could satisfy spec §1 at σ = 0.5 — only 5 donor-rich datasets beyond those chosen" | **16** candidate-bearing datasets beyond the twelve clear 13 v 13; 21 in all, in 15 collections (§6) |
| N5 | An identical set of four admission blockers on all 2190 rows | 993 rows carry a fifth — 981 `excluded_inclusion_gate` and 12 `excluded_confound` |
| N6 | 279 rows fail donor nesting | **284** |
| N7 | Candidate **rows** flagged: depth 1179, assay 547, tissue 299 | Those are flag *occurrences*. The row counts are **1172 / 443 / 269** (suspension's 173 is right) |
| N8 | Median cells per donor over candidates, 96.5 | **96.25** — the two central values are 96.0 and 96.5 |
| N9 | SEA-AD is "the only dataset in the list alive in the σ = 0.7 tier" | **Three** are: #1, #6, #8 — as the proposal's own tier table says |
| N10 | SEA-AD DLPFC's group-B donors run 34 … 46 | **40 … 44.** 34 … 46 is the MTG sibling's |
| N11 | Yoshida's depth Cramér's V is 0.26 … 0.57 | **0.175 … 0.816** over its 47 strata (0.263 … 0.518 over the recommended seven) |
| N12 | KPMP reaches ≥ 23v23 in three strata, named | **Six** — as the proposal's own tier table says, and as §9 item 5 says |
| N13 | "300 counts/cell on 57 cells gives ≈ 24k counts per pseudosample" | 300 × 57 = **17 100**. The 24k uses the stratum's real 426 counts/cell |
| N14 | #10's counts per cell are "150× higher" than #7's | **200×** — 56 841.5 against 284.0 |
| N15 | #11 is "the smallest dataset in the list" at 22 502 cells | **#12 is smaller**, at 18 386 cells |
| N16 | Reserve `9f222629` is "the only source of Drop-seq / Seq-Well" | Two candidate-bearing datasets carry each, one of them another of the proposal's own reserves |
| N17 | The effect-size counts: 3 strong, 5 subtle | This document counts 5 strong and 6 subtle (§4.1). Both are judgements, but the records differ, and only this document discloses that spec §1 (i)'s parenthetical is unmet |

Five further statements describe a per-stratum quantity as if it were fixed for a dataset (control
group sizes given as "≈ 9–10 donors", "7 donors", "35 donors", "10 donors"; a near-confound said to
affect "half" of #12's strata when it affects 3 of 8). They are recorded for completeness and change
nothing.
