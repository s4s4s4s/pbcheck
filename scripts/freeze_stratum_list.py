"""Emit the frozen Phase 0 stratum list from the committed candidate manifest (spec §1).

Spec §1 ends with "**Pre-register the stratum list before computing any metric.**" This script is
the machine-readable half of that act; the human half is ``docs/PREREGISTRATION_STRATUM_LIST.md``,
which this file must be read next to. Neither is an amendment — the freeze applies §1 rather than
changing it — and once both are committed the only route by which any of it may change is a dated
entry in ``docs/AMENDMENTS.md``.

**Nothing here selects anything.** The one judgement in the freeze is the list of twelve
``dataset_id``s in :data:`FROZEN_DATASETS`, argued in the document. Given those twelve, the strata
follow mechanically: every row of the source manifest with ``gate_status == "candidate"`` belonging
to one of them is in the analysis set, and no row of any of them is out. There is no stratum-level
discretion to audit, because there is none to exercise — which is the point, and the reason this
runs as a script instead of being a table someone typed.

The two same-collection siblings are frozen by the same rule and emitted under
``role = within_collection_control``, separately from the 251. Naming them without freezing them
would have left 27 runnable strata unlisted and therefore selectable after the fact.

**Nothing here admits anything either.** The source manifest carries ``admitted_to_sweep = False``
on all 2190 of its rows and this script propagates that value unchanged; ``integer_check`` and
``frozen_universe_size`` are still ``pending`` (they are computed at X load, not from ``obs``), and
``sigma_donor`` is still unanchored, so envelope membership is unknown for every stratum here.
Freezing the list is upstream of admission and does not perform it.

**The output is byte-reproducible and that is load-bearing.** No timestamp of its own generation, no
package versions, no environment: everything in the header is either a pinned constant, a field of
the source manifest, or arithmetic over the source manifest's rows. Running this script twice, or
on another platform, must produce the same bytes, and ``tests/test_stratum_list_freeze.py`` compares
the committed artifact against a fresh run to enforce it. Line endings are pinned explicitly (LF for
the JSON, CRLF for the CSV, matching :func:`pbcheck.census_select.emit_manifest`'s writer) and
``.gitattributes`` marks ``pilot/preregistration/`` as non-text so git never rewrites them.

**Declared constants are verified, not trusted.** Every figure this script states in advance — the
per-dataset stratum counts and donor ceilings, the total, the Layer B subsets, the bin occupancy —
is recomputed from the manifest and compared against the declaration. A disagreement is a hard
failure with the two values printed, never a silent correction: a pre-registration that quietly
re-derives its own numbers when they stop matching is not a pre-registration.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from pbcheck import census_select as cs  # noqa: E402
from pbcheck import gate_config as gc  # noqa: E402

# ---------------------------------------------------------------------------
# The source artifact. Pinned by hash, not by path: the path is where we keep it, the hash is
# what makes it the same bytes the numbers below were derived from.
# ---------------------------------------------------------------------------

PREREG_DIR = REPO / "pilot" / "preregistration"

#: The GitHub Actions run of `.github/workflows/census-candidates.yml` that produced the manifest.
CI_RUN_ID = "31910799023"
CI_WORKFLOW = ".github/workflows/census-candidates.yml (workflow_dispatch, dry_run=false)"

SOURCE_JSON = PREREG_DIR / f"census_candidates_run{CI_RUN_ID}_2026-08-15.json"
SOURCE_CSV = PREREG_DIR / f"census_candidates_run{CI_RUN_ID}_2026-08-15.csv"

SOURCE_SHA256 = "33f8a800229dccc5f58f311e7d0c493655068d43563b31ff53fdaebb3b44e4b4"
SOURCE_CSV_SHA256 = "09eb110dd308155f64e10b2b05beff36854f7125b1434699935707d3551f12d6"
SOURCE_BYTES = 6_630_446
SOURCE_CSV_BYTES = 4_513_660

#: The proposal document the choice of twelve was made from, committed beside the manifest and
#: pinned here for the same reason: it is referenced by §3.2, §4.2 and §8 of the document — for the
#: per-dataset rationale, the rejected datasets, the FIVE NAMED RESERVES, the third 5' candidate and
#: the Mathys search method — and a justification nobody can open is not a justification.
#:
#: It is NOT part of the binding act and nothing here derives from it: what binds is
#: :data:`FROZEN_DATASETS` and :data:`SELECTION_RULE`. It is pinned so that a substitution from its
#: reserve list would be detectable, never so that one becomes available (§9 item 7). It is
#: committed byte-for-byte as circulated, errors included; where it disagrees with the document,
#: the document governs.
PROPOSAL_MD = PREREG_DIR / "stratum_list_proposal_2026-08-16.md"
PROPOSAL_SHA256 = "50872414b0727c129a824b0c65ed179674ac5d6c9ecaac53327568b3eae6fb48"
PROPOSAL_BYTES = 92_589
PROPOSAL_STATUS = (
    "The reasoning behind the choice of twelve, committed for auditability. NOT part of the "
    "binding act: the binding content is the frozen list and the selection rule. Committed "
    "unmodified, including the two figures the document corrects (Rexach's ceiling, and what the "
    "manifest could have supported at the envelope tiers); where the two disagree, the document "
    "governs."
)

#: The manifest's own header stamp. Checked, so that a file with the right hash but a header we did
#: not expect (i.e. a hash constant edited to match a different artifact) still fails.
SOURCE_GENERATED_UTC = "2026-08-15T22:18:37+00:00"
SOURCE_N_ROWS = 2190
SOURCE_N_CANDIDATES = 1197

#: The date of the freeze itself — the act's date, not the run's, and not "now". A generation
#: timestamp here would make the artifact non-reproducible for no gain: the run it derives from is
#: already stamped, and the freeze date belongs to the commit.
FREEZE_DATE = "2026-08-16"
OUT_STEM = f"stratum_list_{FREEZE_DATE}"

DOCUMENT = "docs/PREREGISTRATION_STRATUM_LIST.md"

# ---------------------------------------------------------------------------
# The rule and the frozen set.
# ---------------------------------------------------------------------------

SELECTION_RULE = (
    "Every row of the source manifest whose gate_status == 'candidate' and whose dataset_id is one "
    "of the twelve frozen_datasets below is in the analysis set; every other row is out. No "
    "stratum, cell type, disease term or donor-count tier is chosen by hand. The twelve dataset "
    "ids are the whole of the judgement; the 251 strata are its arithmetic consequence."
)

SELECTION_RULE_NOTES = (
    "The per-dataset 'recommended strata' tables that circulated in the 2026-08-16 proposal "
    "document were illustrative reading aids and are NOT this frozen set. They are a strict "
    "subset of it and carry no status. If a later reader finds a smaller list of strata attributed "
    "to this study, it is not the pre-registration.",
    "'candidate' is census_select's own gate_status: the stratum cleared the obs-decidable half of "
    "the §1 inclusion gate and was not excluded by the §1 confound pre-screen. It is not an "
    "admission — see the admission block below.",
    "The same rule is applied to the two same-collection siblings, whose rows are emitted under "
    "role = within_collection_control. They are frozen so that no runnable stratum of this Census "
    "pin is left unlisted, and they are NOT part of the analysis set: they never enter the D2 "
    "denominator and promoting one to an independent dataset is an amendment.",
)

#: The twelve independent datasets. ``n_strata`` and ``ceiling_min_donors_per_group`` are declared
#: here and re-derived from the manifest at run time; a mismatch aborts.
#:
#: ``assay`` / ``suspension`` / ``tissue`` / ``doi`` come from the CELLxGENE Discover curation API
#: index (``GET /curation/v1/datasets``, 2216 datasets, read 2026-08-16). They are NOT columns of
#: the manifest and cannot be re-derived from it — ``census_select`` screens assay and suspension
#: for confounding but does not emit their levels — so they are recorded here as declared metadata
#: with a named external source, checkable against the public Discover record for each dataset_id.
#: ``expected_effect`` is a LITERATURE JUDGEMENT per (dataset, disease term); see EXPECTED_EFFECT.
FROZEN_DATASETS = (
    {
        "rank": 1,
        "dataset_id": "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3",
        "short": "Gabitto 2024, SEA-AD DLPFC",
        "doi": "10.1038/s41593-024-01774-5",
        "assay": ("10x 3' v3", "10x multiome"),
        "suspension": ("nucleus",),
        "tissue": ("dorsolateral prefrontal cortex",),
        "n_strata": 18,
        "ceiling_min_donors_per_group": 39,
    },
    {
        "rank": 2,
        "dataset_id": "ac0c6561-7a48-4185-af6f-af799f699172",
        "short": "Rexach 2024 Cell, cross-dementia",
        "doi": "10.1016/j.cell.2024.08.019",
        "assay": ("10x 3' v2", "10x 3' v3"),
        "suspension": ("nucleus",),
        "tissue": ("Brodmann (1909) area 4", "insular cortex", "primary visual cortex"),
        "n_strata": 27,
        "ceiling_min_donors_per_group": 10,
    },
    {
        "rank": 3,
        "dataset_id": "d8da613f-e681-4c69-b463-e94f5e66847f",
        "short": "Melms 2021 Nature, lethal COVID-19 lung",
        "doi": "10.1038/s41586-021-03569-1",
        "assay": ("10x 3' v3",),
        "suspension": ("nucleus",),
        "tissue": ("lung",),
        "n_strata": 28,
        "ceiling_min_donors_per_group": 7,
    },
    {
        "rank": 4,
        "dataset_id": "2a498ace-872a-4935-984b-1afa70fd9886",
        "short": "Yoshida 2022 Nature, PBMC",
        "doi": "10.1038/s41586-021-04345-x",
        "assay": ("10x 5' v1",),
        "suspension": ("cell",),
        "tissue": ("blood",),
        "n_strata": 47,
        "ceiling_min_donors_per_group": 20,
    },
    {
        "rank": 5,
        "dataset_id": "ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3",
        "short": "Ahern 2022 Cell, COMBAT blood atlas",
        "doi": "10.1016/j.cell.2022.01.012",
        "assay": ("10x 5' v1",),
        "suspension": ("cell",),
        "tissue": ("blood",),
        "n_strata": 25,
        "ceiling_min_donors_per_group": 10,
    },
    {
        "rank": 6,
        "dataset_id": "f1606894-59df-4794-a37f-baa7c6fb6de1",
        "short": "Linna-Kuosmanen 2024, PERIHEART right atrium",
        "doi": "10.1016/j.xcrm.2024.101556",
        "assay": ("10x 3' v3",),
        "suspension": ("nucleus",),
        "tissue": ("right atrium auricular region",),
        "n_strata": 11,
        "ceiling_min_donors_per_group": 25,
    },
    {
        "rank": 7,
        "dataset_id": "d18736c3-6292-4379-919a-d6d973204c87",
        "short": "Binvignat 2024, rheumatoid arthritis blood",
        "doi": "10.1172/jci.insight.178499",
        "assay": ("10x 3' v3",),
        "suspension": ("cell",),
        "tissue": ("blood",),
        "n_strata": 15,
        "ceiling_min_donors_per_group": 18,
    },
    {
        "rank": 8,
        "dataset_id": "a12ccb9b-4fbe-457d-8590-ac78053259ef",
        "short": "KPMP adult human kidney snRNA-seq v1.5",
        # Discover returns collection_doi = null for this collection: the publication provenance is
        # NOT ESTABLISHED from the API and is recorded as such rather than guessed.
        "doi": None,
        "assay": ("10x 3' v3",),
        "suspension": ("nucleus",),
        "tissue": ("cortex of kidney", "kidney", "renal medulla", "renal papilla"),
        "n_strata": 37,
        "ceiling_min_donors_per_group": 24,
    },
    {
        "rank": 9,
        "dataset_id": "19e46756-9100-4e01-8b0e-23b557558a4c",
        "short": "Heimlich 2024, clonal haematopoiesis PBMC",
        "doi": "10.1182/bloodadvances.2023011445",
        "assay": ("10x 3' v3",),
        "suspension": ("cell",),
        "tissue": ("blood",),
        "n_strata": 7,
        "ceiling_min_donors_per_group": 7,
    },
    {
        "rank": 10,
        "dataset_id": "c893ddc3-f25b-45e2-8c9e-155918b4261c",
        "short": "Phan 2024, opioid use disorder striatum",
        "doi": "10.1038/s41467-024-45165-7",
        "assay": ("10x 3' v3",),
        "suspension": ("nucleus",),
        "tissue": ("caudate nucleus", "putamen"),
        "n_strata": 10,
        "ceiling_min_donors_per_group": 6,
    },
    {
        "rank": 11,
        "dataset_id": "8e47ed12-c658-4252-b126-381df8d52a3d",
        "short": "Elmentaite 2020, paediatric gut (Crohn)",
        "doi": "10.1016/j.devcel.2020.11.010",
        "assay": ("10x 3' v2",),
        "suspension": ("cell",),
        "tissue": ("ileal mucosa",),
        "n_strata": 18,
        "ceiling_min_donors_per_group": 7,
    },
    {
        "rank": 12,
        "dataset_id": "4b6af54a-4a21-46e0-bc8d-673c0561a836",
        "short": "Wang 2023, emphysema non-immune (3v3 anchor)",
        "doi": "10.1016/j.immuni.2023.01.032",
        "assay": ("10x 3' v3",),
        "suspension": ("cell",),
        "tissue": ("alveolus of lung",),
        "n_strata": 8,
        "ceiling_min_donors_per_group": 3,
    },
)

#: The two role values a frozen row can carry. Every emitted row states one of them, so a row
#: copied out of the artifact still says which set it belongs to.
ROLE_ANALYSIS_SET = "analysis_set"
ROLE_WITHIN_COLLECTION_CONTROL = "within_collection_control"

#: Same-collection datasets, named so that their exclusion is a recorded decision rather than an
#: oversight. D2 clusters evidence by dataset because same-dataset strata share donors, batch and
#: assay; two datasets of one collection share the cohort and the laboratory as well, so counting
#: both toward "majority of independent datasets" would inflate the effective n.
#:
#: They are **frozen too**, and by the same rule: every candidate row of theirs is emitted, under
#: ``role = within_collection_control``. Naming them without freezing them would have left 27
#: runnable strata outside the pre-registration, selectable after the fact — the exact freedom this
#: document exists to remove. ``n_strata`` and ``ceiling_min_donors_per_group`` are declared here
#: and re-derived from the manifest at run time, as for the twelve.
SIBLING_DATASETS = (
    {
        "dataset_id": "c2876b1b-06d8-4d96-a56b-5304f815b99a",
        "short": "SEA-AD MTG",
        "sibling_of": "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3",
        "role": ROLE_WITHIN_COLLECTION_CONTROL,
        "role_note": "different cortical region, same cohort and laboratory",
        "n_strata": 18,
        "ceiling_min_donors_per_group": 42,
    },
    {
        "dataset_id": "1e5bd3b8-6a0e-4959-8d69-cafed30fe814",
        "short": "Emphysema Cell Atlas, immune cells",
        "sibling_of": "4b6af54a-4a21-46e0-bc8d-673c0561a836",
        "role": ROLE_WITHIN_COLLECTION_CONTROL,
        "role_note": "immune compartment of the same 3v3 cohort",
        "n_strata": 9,
        "ceiling_min_donors_per_group": 3,
    },
)

#: What a control row's ``expected_effect`` says. The §1 (i)/(ii) coverage axis is a claim about the
#: twelve; giving a control row a strong/subtle label would smuggle it into that claim.
EXPECTED_EFFECT_NOT_APPLICABLE = "not_applicable"

#: §1 (i)/(ii): "2-3 with a biologically strong expected effect", "2-3 subtle/low-effect". The axis
#: is a LITERATURE JUDGEMENT about the design, made before any metric is computed and recorded here
#: so that it cannot be reassigned afterwards. It is not a prediction of lambda_naive and must never
#: be scored against the results as though it were a hypothesis. It is keyed per (dataset, disease
#: term) because four datasets carry a strong and a subtle arm against one shared control group,
#: which is the cheapest effect-size axis available to us and is lost by a per-dataset label.
EXPECTED_EFFECT_VOCABULARY = ("strong", "moderate", "subtle")
EXPECTED_EFFECT_SOURCE = (
    "literature judgement from each dataset's own publication (DOI above), recorded before any "
    "metric is computed; a coverage criterion for §1 (i)/(ii), never a prediction of the result"
)
EXPECTED_EFFECT = {
    ("6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3", "dementia"): "strong",
    ("ac0c6561-7a48-4185-af6f-af799f699172", "Alzheimer disease"): "strong",
    ("ac0c6561-7a48-4185-af6f-af799f699172", "Pick disease"): "strong",
    ("ac0c6561-7a48-4185-af6f-af799f699172", "progressive supranuclear palsy"): "strong",
    ("d8da613f-e681-4c69-b463-e94f5e66847f", "COVID-19"): "strong",
    ("2a498ace-872a-4935-984b-1afa70fd9886", "COVID-19"): "strong",
    ("2a498ace-872a-4935-984b-1afa70fd9886", "post-COVID-19 disorder"): "subtle",
    ("ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3", "COVID-19"): "strong",
    ("ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3", "influenza"): "moderate",
    ("f1606894-59df-4794-a37f-baa7c6fb6de1", "atrial fibrillation"): "subtle",
    ("d18736c3-6292-4379-919a-d6d973204c87", "rheumatoid arthritis"): "subtle",
    ("a12ccb9b-4fbe-457d-8590-ac78053259ef", "acute kidney failure"): "moderate",
    ("a12ccb9b-4fbe-457d-8590-ac78053259ef", "chronic kidney disease"): "moderate",
    ("19e46756-9100-4e01-8b0e-23b557558a4c", "clonal hematopoiesis"): "subtle",
    ("c893ddc3-f25b-45e2-8c9e-155918b4261c", "opiate dependence"): "subtle",
    ("8e47ed12-c658-4252-b126-381df8d52a3d", "Crohn disease"): "subtle",
    ("4b6af54a-4a21-46e0-bc8d-673c0561a836", "pulmonary emphysema"): "moderate",
}

# ---------------------------------------------------------------------------
# D1 — the cells-per-donor bins, pre-registered here (§1, §7 item 3, decision rule item 2).
# ---------------------------------------------------------------------------

#: Half-decade log bins over cells per donor. **Chosen independently of these data**, so that no
#: boundary can have been fitted to make an occupancy table look full: the sequence is
#: 10 -> 30 -> 100 -> 300 -> 1000 -> 3000, i.e. successive multiplications by ~sqrt(10), started at
#: the §1 inclusion gate's own floor of 10 cells per donor and left open at the top. The last bin is
#: open because the Census's upper tail is unbounded and a closed top bin would be a threshold with
#: nothing behind it.
#:
#: The spec references "the pre-registered bins" three times (decision rule item 2, §1 (iii), §7
#: item 3) and never states them numerically. That omission is a GAP IN §1 and it is closed here,
#: by pre-registration in the same act that freezes the list. It is not a change to §1: no threshold
#: moves, and nothing that §1 defines is redefined.
CELLS_PER_DONOR_BINS = ((10, 30), (30, 100), (100, 300), (300, 1000), (1000, 3000), (3000, None))

CELLS_PER_DONOR_BINS_SOURCE = (
    "half-decade log bins on cells per donor, lower edge 10 = the §1 inclusion-gate floor, ratio "
    "~sqrt(10) per step, top bin open. Chosen before looking at the frozen set's occupancy and "
    "independently of it. The quantity binned is the per-group MEDIAN cells per donor of a stratum "
    "(cells_per_donor_by_group[A|B]['median']), so each stratum contributes two values."
)

# ---------------------------------------------------------------------------
# Layer B — the pre-declared truncation (Amendment 3 Change 1).
# ---------------------------------------------------------------------------

#: If and when the sigma_donor anchor lands, the operating envelope decides which datasets survive.
#: Declared NOW, before the anchor exists, so that the surviving subset cannot be chosen afterwards
#: from among the subsets that happen to be convenient. ``dataset_ids`` is verified against the
#: manifest: a dataset is in the tier iff it holds at least one candidate stratum with at least
#: ``min_donors_per_group`` donors in BOTH groups. ``min_donors_per_group`` is verified against
#: :data:`pbcheck.gate_config.OPERATING_ENVELOPE` rather than restated as a literal.
#:
#: **All four** envelope tiers are declared, not only the two that were obviously going to fail.
#: An earlier version declared 0.5 and 0.7 alone, which meant the 0.2 and 0.35 rows of the
#: document's own table were never checked against the manifest — and the 0.35 row, the tier at
#: ``gate_config.POWER_EVAL_SIGMA`` where the instrument nominally operates, was consequently
#: published with a verdict that contradicted :data:`SPEC_DATASET_FLOOR` in this same file.
LAYER_B = (
    {
        "sigma_donor": 0.2,
        "dataset_ids": (
            "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3",
            "ac0c6561-7a48-4185-af6f-af799f699172",
            "d8da613f-e681-4c69-b463-e94f5e66847f",
            "2a498ace-872a-4935-984b-1afa70fd9886",
            "ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3",
            "f1606894-59df-4794-a37f-baa7c6fb6de1",
            "d18736c3-6292-4379-919a-d6d973204c87",
            "a12ccb9b-4fbe-457d-8590-ac78053259ef",
            "19e46756-9100-4e01-8b0e-23b557558a4c",
            "c893ddc3-f25b-45e2-8c9e-155918b4261c",
            "8e47ed12-c658-4252-b126-381df8d52a3d",
        ),
    },
    {
        "sigma_donor": 0.35,
        "dataset_ids": (
            "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3",
            "ac0c6561-7a48-4185-af6f-af799f699172",
            "2a498ace-872a-4935-984b-1afa70fd9886",
            "ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3",
            "f1606894-59df-4794-a37f-baa7c6fb6de1",
            "d18736c3-6292-4379-919a-d6d973204c87",
            "a12ccb9b-4fbe-457d-8590-ac78053259ef",
        ),
    },
    {
        "sigma_donor": 0.5,
        "dataset_ids": (
            "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3",
            "2a498ace-872a-4935-984b-1afa70fd9886",
            "f1606894-59df-4794-a37f-baa7c6fb6de1",
            "d18736c3-6292-4379-919a-d6d973204c87",
            "a12ccb9b-4fbe-457d-8590-ac78053259ef",
        ),
    },
    {
        "sigma_donor": 0.7,
        "dataset_ids": (
            "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3",
            "f1606894-59df-4794-a37f-baa7c6fb6de1",
            "a12ccb9b-4fbe-457d-8590-ac78053259ef",
        ),
    },
)

#: §1: "First pass = 8-12 datasets chosen to SPAN the outcome space". The lower end of that range is
#: the floor Layer B is measured against, and the frozen list falls below it at three of the four
#: tiers. The verdict per tier is ``below_spec_dataset_floor``, computed — never typed.
SPEC_DATASET_FLOOR = 8

# ---------------------------------------------------------------------------
# What the manifest as a whole could have supported (§6's counterfactual).
# ---------------------------------------------------------------------------

#: How many of the manifest's 68 candidate-bearing datasets hold at least one stratum at each
#: envelope tier's donors-per-group demand, keyed by that demand. §6 rests on these: they are what
#: makes the truncation a consequence of the SELECTION rather than of the public data. Declared and
#: re-derived, like everything else here.
MANIFEST_TIER_CENSUS = {4: 62, 8: 33, 13: 21, 23: 12}

#: §1 (iii) requires "some exactly 3v3". Only this many of the manifest's datasets hold an
#: exactly-3v3 stratum at all, which is why the requirement costs a dataset slot that could
#: otherwise have gone to a donor-rich one.
MANIFEST_DATASETS_WITH_EXACT_3V3 = 3

#: The counterfactual §6 states: the most datasets a TWELVE-dataset list could have kept at each
#: tier while still holding an exactly-3v3 stratum, keyed by donors-per-group. Derived from
#: :data:`MANIFEST_TIER_CENSUS` and the 3v3-bearing set, and verified against the manifest.
COUNTERFACTUAL_MAX_SURVIVORS = {4: 12, 8: 12, 13: 11, 23: 11}

#: The list that attains :data:`COUNTERFACTUAL_MAX_SURVIVORS` at every tier at once — an existence
#: witness, so the counterfactual is a construction rather than an upper bound nobody realised.
#: Rexach (#2 of the frozen twelve) is the 3v3 anchor; the other eleven are eleven of the twelve
#: datasets that clear 23v23. It is NOT a proposal and nothing may be selected from it.
COUNTERFACTUAL_WITNESS_ANCHOR = "ac0c6561-7a48-4185-af6f-af799f699172"

# ---------------------------------------------------------------------------
# Output schema.
# ---------------------------------------------------------------------------

#: Appended after the manifest's own columns, never interleaved with them, so that a reader can see
#: at a glance which values are the source artifact's and which this script computed.
DERIVED_FIELDS = (
    "role",
    "dataset_rank",
    "dataset_short",
    "expected_effect",
    "min_donors_per_group",
    "cells_per_donor_bin_A",
    "cells_per_donor_bin_B",
)

STRATUM_FIELDS = (*cs.MANIFEST_FIELDS, *DERIVED_FIELDS)


class SourceArtifactMismatch(RuntimeError):
    """The committed source artifact is not the one this freeze was derived from."""


class FrozenDeclarationMismatch(RuntimeError):
    """A figure declared in this file disagrees with the source manifest."""


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bin_label(low: int, high: int | None) -> str:
    return f"[{low},{'inf' if high is None else high})"


def bin_of(value: float) -> str:
    """The :data:`CELLS_PER_DONOR_BINS` label containing ``value``.

    A value below the lowest edge is an error rather than an "other" bucket: the §1 inclusion gate
    drops every donor with fewer than ``gate_config.MIN_CELLS`` = 10 cells in the cell type, so a
    surviving group's *median* cells per donor cannot be below 10. If one ever is, the gate and the
    bins disagree about the same number and that must stop the freeze, not fall into a bin.
    """
    for low, high in CELLS_PER_DONOR_BINS:
        if value >= low and (high is None or value < high):
            return bin_label(low, high)
    raise FrozenDeclarationMismatch(
        f"cells-per-donor median {value} falls below the lowest bin edge "
        f"{CELLS_PER_DONOR_BINS[0][0]}, which the §1 inclusion gate's >= {gc.MIN_CELLS} "
        "cells-per-donor rule should have made impossible"
    )


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FrozenDeclarationMismatch(message)


def check_pinned_file(path: Path, digest: str, size: int, what: str) -> None:
    """Refuse a hash-pinned companion file that is absent, resized or altered."""
    if not path.exists():
        raise SourceArtifactMismatch(
            f"{what} not found at {path}. Its sha256 is recorded in {DOCUMENT} and its absence is "
            "not a reason to skip the check."
        )
    actual_size = path.stat().st_size
    if actual_size != size:
        raise SourceArtifactMismatch(
            f"{path.name}: {actual_size} bytes, expected {size}"
        )
    actual = sha256_of(path)
    if actual != digest:
        raise SourceArtifactMismatch(f"{path.name}: sha256 {actual}, expected {digest}")


def load_source(path: Path = SOURCE_JSON, *, csv_path: Path | None = SOURCE_CSV) -> dict:
    """Read the source manifest, refusing anything whose bytes are not the pinned ones.

    The hash is checked before the JSON is parsed, and the header stamps are checked after: a hash
    guards against the file changing, and the stamps guard against the *constant* changing to match
    some other file. Both are needed; either alone can be walked around by one edit.

    A ``csv_path`` that is named and **absent** aborts. This document's own position is that an
    evidence file whose bytes cannot be checked is not evidence, so skipping the check because the
    file has gone missing is the one response that cannot be right: deleting the CSV would
    otherwise turn a failing hash check into a passing run. ``csv_path=None`` — used only by tests
    that hand ``load_source`` a synthetic JSON with no CSV twin — declares the absence instead.
    """
    if not path.exists():
        raise SourceArtifactMismatch(f"source manifest not found at {path}")
    size = path.stat().st_size
    if size != SOURCE_BYTES:
        raise SourceArtifactMismatch(
            f"{path.name}: {size} bytes, expected {SOURCE_BYTES} — this is not the artifact of CI "
            f"run {CI_RUN_ID}. The freeze does not run on a substitute."
        )
    digest = sha256_of(path)
    if digest != SOURCE_SHA256:
        raise SourceArtifactMismatch(
            f"{path.name}: sha256 {digest}, expected {SOURCE_SHA256}. The pre-registered numbers "
            "were derived from the pinned bytes and mean nothing against different ones."
        )
    if csv_path is not None:
        check_pinned_file(csv_path, SOURCE_CSV_SHA256, SOURCE_CSV_BYTES, "CSV twin")
    check_pinned_file(PROPOSAL_MD, PROPOSAL_SHA256, PROPOSAL_BYTES, "proposal document")

    manifest = json.loads(path.read_text(encoding="utf-8"))
    header, rows = manifest["header"], manifest["rows"]
    if header.get("generated_utc") != SOURCE_GENERATED_UTC:
        raise SourceArtifactMismatch(
            f"source generated_utc {header.get('generated_utc')!r}, expected "
            f"{SOURCE_GENERATED_UTC!r}"
        )
    if header.get("census_version") != cs.CENSUS_VERSION:
        raise SourceArtifactMismatch(
            f"source census_version {header.get('census_version')!r}, expected "
            f"{cs.CENSUS_VERSION!r} (spec §1 pin)"
        )
    if len(rows) != SOURCE_N_ROWS:
        raise SourceArtifactMismatch(f"source has {len(rows)} rows, expected {SOURCE_N_ROWS}")
    return manifest


def select_strata(rows) -> list[dict]:
    """The selection rule, in code. Nothing else in this module decides membership."""
    frozen_ids = {d["dataset_id"] for d in FROZEN_DATASETS}
    selected = [
        row for row in rows
        if row["gate_status"] == "candidate" and row["dataset_id"] in frozen_ids
    ]
    order = {d["dataset_id"]: d["rank"] for d in FROZEN_DATASETS}
    selected.sort(key=lambda r: (order[r["dataset_id"]], r["disease"], r["cell_type"]))
    return selected


def select_control_strata(rows) -> list[dict]:
    """The same rule, applied to the two siblings — they are frozen, not merely named.

    Ordered by the rank of the frozen dataset each is a sibling of, so the control block reads
    alongside the analysis set rather than in uuid order.
    """
    order = {d["dataset_id"]: d["rank"] for d in FROZEN_DATASETS}
    sibling_order = {
        s["dataset_id"]: order[s["sibling_of"]] for s in SIBLING_DATASETS
        if s["sibling_of"] in order
    }
    selected = [
        row for row in rows
        if row["gate_status"] == "candidate" and row["dataset_id"] in sibling_order
    ]
    selected.sort(key=lambda r: (sibling_order[r["dataset_id"]], r["disease"], r["cell_type"]))
    return selected


def _tier_datasets(rows, threshold: int) -> set[str]:
    """The datasets holding at least one row with ``threshold`` donors in BOTH groups."""
    return {
        row["dataset_id"] for row in rows
        if min(row["n_donors_A"], row["n_donors_B"]) >= threshold
    }


def _datasets_with_exact_3v3(rows) -> set[str]:
    return {row["dataset_id"] for row in rows
            if row["n_donors_A"] == 3 and row["n_donors_B"] == 3}


def counterfactual_max_survivors(candidates, threshold: int) -> int:
    """Most of a twelve-dataset list that could clear ``threshold``, given §1 (iii)'s 3v3 anchor.

    §1 (iii) requires "some exactly 3v3", so at least one slot of any compliant list must go to a
    dataset holding an exactly-3v3 stratum. If one of those also clears the tier, all twelve slots
    can; otherwise one slot is spent below the tier and eleven is the ceiling.
    """
    tier = _tier_datasets(candidates, threshold)
    anchors = _datasets_with_exact_3v3(candidates)
    if anchors & tier:
        return min(12, len(tier))
    return min(11, len(tier))


def _envelope_min_donors(sigma_donor: float) -> int:
    for row in gc.OPERATING_ENVELOPE:
        if float(row["sigma_donor"]) == sigma_donor:
            return int(row["min_donors_per_group"])
    raise FrozenDeclarationMismatch(
        f"sigma_donor={sigma_donor} is not a row of gate_config.OPERATING_ENVELOPE; Layer B may "
        "only be declared at sigma values the envelope itself states"
    )


def verify_declarations(manifest: dict, strata: list[dict], controls: list[dict]) -> None:
    """Recompute every declared figure from the manifest and abort on any disagreement."""
    rows = manifest["rows"]
    candidates = [r for r in rows if r["gate_status"] == "candidate"]
    _require(
        len(candidates) == SOURCE_N_CANDIDATES,
        f"source has {len(candidates)} candidate rows, declared {SOURCE_N_CANDIDATES}",
    )

    ids = [d["dataset_id"] for d in FROZEN_DATASETS]
    _require(len(set(ids)) == len(ids) == 12, "frozen_datasets must be 12 distinct dataset ids")
    _require(
        [d["rank"] for d in FROZEN_DATASETS] == list(range(1, 13)),
        "frozen_datasets ranks must be 1..12 in order",
    )
    sibling_ids = {s["dataset_id"] for s in SIBLING_DATASETS}
    _require(
        not sibling_ids & set(ids),
        "a sibling dataset is also listed as an independent dataset — the D2 denominator would "
        "count the same collection twice",
    )
    for sibling in SIBLING_DATASETS:
        _require(
            sibling["sibling_of"] in set(ids),
            f"sibling {sibling['dataset_id']} names a parent that is not in the frozen set",
        )
        own = [r for r in candidates if r["dataset_id"] == sibling["dataset_id"]]
        _require(
            bool(own),
            f"sibling {sibling['dataset_id']} holds no candidate row in the source manifest",
        )
        _require(
            len(own) == sibling["n_strata"],
            f"{sibling['short']}: manifest holds {len(own)} candidate strata, declared "
            f"{sibling['n_strata']}",
        )
        sibling_ceiling = max(min(r["n_donors_A"], r["n_donors_B"]) for r in own)
        _require(
            sibling_ceiling == sibling["ceiling_min_donors_per_group"],
            f"{sibling['short']}: ceiling min(n_donors_A, n_donors_B) is {sibling_ceiling}, "
            f"declared {sibling['ceiling_min_donors_per_group']}",
        )
        _require(
            sibling["role"] == ROLE_WITHIN_COLLECTION_CONTROL,
            f"{sibling['short']}: role is {sibling['role']!r}, and a frozen sibling may carry no "
            f"role but {ROLE_WITHIN_COLLECTION_CONTROL!r}",
        )

    declared_total = 0
    for dataset in FROZEN_DATASETS:
        own = [r for r in candidates if r["dataset_id"] == dataset["dataset_id"]]
        _require(
            len(own) == dataset["n_strata"],
            f"{dataset['short']}: manifest holds {len(own)} candidate strata, declared "
            f"{dataset['n_strata']}",
        )
        ceiling = max(min(r["n_donors_A"], r["n_donors_B"]) for r in own)
        _require(
            ceiling == dataset["ceiling_min_donors_per_group"],
            f"{dataset['short']}: ceiling min(n_donors_A, n_donors_B) is {ceiling}, declared "
            f"{dataset['ceiling_min_donors_per_group']}",
        )
        declared_total += dataset["n_strata"]

    _require(
        len(strata) == declared_total,
        f"the rule selected {len(strata)} strata; the per-dataset declarations sum to "
        f"{declared_total}",
    )
    declared_controls = sum(s["n_strata"] for s in SIBLING_DATASETS)
    _require(
        len(controls) == declared_controls,
        f"the rule selected {len(controls)} within-collection control strata; the per-sibling "
        f"declarations sum to {declared_controls}",
    )
    _require(
        not ({(r["dataset_id"], r["cell_type"], r["disease"]) for r in controls}
             & {(r["dataset_id"], r["cell_type"], r["disease"]) for r in strata}),
        "a within-collection control stratum is also in the analysis set; the two sets must stay "
        "distinguishable or the D2 denominator can acquire a sibling by accident",
    )
    _require(
        not any(r["admitted_to_sweep"] for r in rows),
        "a row of the source manifest carries admitted_to_sweep = True; the freeze is upstream of "
        "admission and cannot run on an artifact that has already admitted something",
    )

    terms = {(r["dataset_id"], r["disease"]) for r in strata}
    missing = sorted(terms - set(EXPECTED_EFFECT))
    extra = sorted(set(EXPECTED_EFFECT) - terms)
    _require(not missing, f"EXPECTED_EFFECT does not label {missing}")
    _require(not extra, f"EXPECTED_EFFECT labels {extra}, which the frozen set does not contain")
    _require(
        set(EXPECTED_EFFECT.values()) <= set(EXPECTED_EFFECT_VOCABULARY),
        f"EXPECTED_EFFECT uses labels outside {EXPECTED_EFFECT_VOCABULARY}",
    )

    envelope_sigmas = [float(row["sigma_donor"]) for row in gc.OPERATING_ENVELOPE]
    _require(
        [tier["sigma_donor"] for tier in LAYER_B] == envelope_sigmas,
        f"Layer B declares tiers {[t['sigma_donor'] for t in LAYER_B]}; "
        f"gate_config.OPERATING_ENVELOPE has {envelope_sigmas}. Every tier of the envelope must be "
        "declared, or the document publishes a row the manifest never checked.",
    )
    for tier in LAYER_B:
        threshold = _envelope_min_donors(tier["sigma_donor"])
        derived = {
            r["dataset_id"] for r in strata
            if min(r["n_donors_A"], r["n_donors_B"]) >= threshold
        }
        _require(
            derived == set(tier["dataset_ids"]),
            f"Layer B at sigma_donor={tier['sigma_donor']} (>= {threshold} donors/group): the "
            f"manifest gives {sorted(derived)}, the declaration gives "
            f"{sorted(tier['dataset_ids'])}",
        )

    for threshold, declared in MANIFEST_TIER_CENSUS.items():
        derived_n = len(_tier_datasets(candidates, threshold))
        _require(
            derived_n == declared,
            f"the manifest holds {derived_n} candidate-bearing datasets with a stratum at "
            f">= {threshold}v{threshold}, declared {declared}",
        )
    _require(
        sorted(MANIFEST_TIER_CENSUS) == sorted(
            int(row["min_donors_per_group"]) for row in gc.OPERATING_ENVELOPE
        ),
        "the manifest tier census must be keyed by the envelope's own donors-per-group demands",
    )

    anchors = _datasets_with_exact_3v3(candidates)
    _require(
        len(anchors) == MANIFEST_DATASETS_WITH_EXACT_3V3,
        f"{len(anchors)} candidate-bearing datasets hold an exactly-3v3 stratum, declared "
        f"{MANIFEST_DATASETS_WITH_EXACT_3V3}",
    )
    _require(
        COUNTERFACTUAL_WITNESS_ANCHOR in anchors,
        f"the counterfactual's 3v3 anchor {COUNTERFACTUAL_WITNESS_ANCHOR} holds no exactly-3v3 "
        "candidate stratum",
    )
    for threshold, declared in COUNTERFACTUAL_MAX_SURVIVORS.items():
        derived_n = counterfactual_max_survivors(candidates, threshold)
        _require(
            derived_n == declared,
            f"a twelve-dataset list holding an exactly-3v3 stratum could have kept {derived_n} "
            f"datasets at >= {threshold}v{threshold}, declared {declared}",
        )
        witness = set(sorted(_tier_datasets(candidates, 23))[:11]) | {
            COUNTERFACTUAL_WITNESS_ANCHOR
        }
        _require(
            len(witness) == 12,
            f"the counterfactual witness is {len(witness)} datasets, not twelve",
        )
        attained = len(witness & _tier_datasets(candidates, threshold))
        _require(
            attained == declared,
            f"the counterfactual witness keeps {attained} datasets at >= {threshold}v{threshold}, "
            f"but the maximum is {declared}; the witness no longer attains it",
        )

    for low, high in CELLS_PER_DONOR_BINS:
        _require(low < (high if high is not None else float("inf")),
                 f"bin [{low},{high}) is empty or inverted")
    edges = [low for low, _ in CELLS_PER_DONOR_BINS]
    tops = [high for _, high in CELLS_PER_DONOR_BINS]
    _require(
        tops[:-1] == edges[1:] and tops[-1] is None,
        "the cells-per-donor bins must tile [10, inf) with no gap and no overlap",
    )


def group_medians(strata: list[dict]) -> list[float]:
    return [
        float(row["cells_per_donor_by_group"][group]["median"])
        for row in strata
        for group in ("A", "B")
    ]


def bin_occupancy(strata: list[dict]) -> dict:
    medians = group_medians(strata)
    counts = {bin_label(low, high): 0 for low, high in CELLS_PER_DONOR_BINS}
    for value in medians:
        counts[bin_of(value)] += 1
    return {
        "definition": CELLS_PER_DONOR_BINS_SOURCE,
        "bins": [bin_label(low, high) for low, high in CELLS_PER_DONOR_BINS],
        "n_group_medians": len(medians),
        "min_group_median": min(medians),
        "max_group_median": max(medians),
        "occupancy": counts,
        "all_bins_occupied": all(count > 0 for count in counts.values()),
    }


def _dataset_block(dataset: dict, strata: list[dict]) -> dict:
    own = [r for r in strata if r["dataset_id"] == dataset["dataset_id"]]
    terms = sorted({r["disease"] for r in own})
    return {
        "rank": dataset["rank"],
        "dataset_id": dataset["dataset_id"],
        "short": dataset["short"],
        "doi": dataset["doi"],
        "assay": list(dataset["assay"]),
        "suspension": list(dataset["suspension"]),
        "tissue": list(dataset["tissue"]),
        "n_strata": len(own),
        "ceiling_min_donors_per_group": max(min(r["n_donors_A"], r["n_donors_B"]) for r in own),
        "n_strata_at_least_8v8": sum(
            1 for r in own if min(r["n_donors_A"], r["n_donors_B"]) >= 8
        ),
        "n_strata_exactly_3v3": sum(
            1 for r in own if r["n_donors_A"] == 3 and r["n_donors_B"] == 3
        ),
        "disease_terms": [
            {"disease": term, "expected_effect": EXPECTED_EFFECT[(dataset["dataset_id"], term)]}
            for term in terms
        ],
    }


def _sibling_block(sibling: dict, controls: list[dict]) -> dict:
    own = [r for r in controls if r["dataset_id"] == sibling["dataset_id"]]
    return {
        "dataset_id": sibling["dataset_id"],
        "short": sibling["short"],
        "sibling_of": sibling["sibling_of"],
        "role": sibling["role"],
        "role_note": sibling["role_note"],
        "n_strata": len(own),
        "ceiling_min_donors_per_group": max(min(r["n_donors_A"], r["n_donors_B"]) for r in own),
        "disease_terms": sorted({r["disease"] for r in own}),
    }


def build_header(manifest: dict, strata: list[dict], controls: list[dict]) -> dict:
    candidates = [r for r in manifest["rows"] if r["gate_status"] == "candidate"]
    layer_b = []
    for tier in LAYER_B:
        threshold = _envelope_min_donors(tier["sigma_donor"])
        layer_b.append({
            "sigma_donor": tier["sigma_donor"],
            "min_donors_per_group": threshold,
            "n_datasets": len(tier["dataset_ids"]),
            "dataset_ids": list(tier["dataset_ids"]),
            "n_strata": sum(
                1 for r in strata if min(r["n_donors_A"], r["n_donors_B"]) >= threshold
            ),
            "below_spec_dataset_floor": len(tier["dataset_ids"]) < SPEC_DATASET_FLOOR,
            "n_manifest_datasets_at_tier": len(_tier_datasets(candidates, threshold)),
            "counterfactual_max_survivors_of_twelve": counterfactual_max_survivors(
                candidates, threshold
            ),
        })

    return {
        "act": "spec §1 pre-registration of the Phase 0 stratum list, frozen before any metric is "
               "computed on these strata",
        "document": DOCUMENT,
        "frozen_date": FREEZE_DATE,
        "not_an_amendment":
            "This applies PHASE0_SPEC.md §1; it changes no threshold and supersedes no section. "
            "Once committed, every change to anything below goes through docs/AMENDMENTS.md.",
        "census_version": cs.CENSUS_VERSION,
        "source": {
            "json": SOURCE_JSON.name,
            "json_sha256": SOURCE_SHA256,
            "json_bytes": SOURCE_BYTES,
            "csv": SOURCE_CSV.name,
            "csv_sha256": SOURCE_CSV_SHA256,
            "csv_bytes": SOURCE_CSV_BYTES,
            "generated_utc": SOURCE_GENERATED_UTC,
            "ci_run_id": CI_RUN_ID,
            "ci_workflow": CI_WORKFLOW,
            "producer": "scripts/census_candidates.py over the whole pinned Census",
            "n_rows": SOURCE_N_ROWS,
            "n_candidate_rows": SOURCE_N_CANDIDATES,
        },
        "reasoning_document": {
            "file": PROPOSAL_MD.name,
            "sha256": PROPOSAL_SHA256,
            "bytes": PROPOSAL_BYTES,
            "status": PROPOSAL_STATUS,
        },
        "selection_rule": SELECTION_RULE,
        "selection_rule_notes": list(SELECTION_RULE_NOTES),
        "n_frozen_strata": len(strata),
        "n_frozen_datasets": len(FROZEN_DATASETS),
        "frozen_datasets": [_dataset_block(d, strata) for d in FROZEN_DATASETS],
        "frozen_dataset_metadata_source":
            "assay / suspension / tissue / doi are from the CELLxGENE Discover curation API index "
            "(GET /curation/v1/datasets, 2216 datasets, read 2026-08-16). They are NOT columns of "
            "the source manifest and cannot be re-derived from it; they are checkable against the "
            "public Discover record for each dataset_id.",
        "sibling_datasets": [_sibling_block(s, controls) for s in SIBLING_DATASETS],
        "n_within_collection_control_strata": len(controls),
        "sibling_datasets_note":
            "Named, FROZEN and EXCLUDED from the D2 'independent datasets' denominator. Their "
            "candidate strata are emitted in full under within_collection_control_rows, by the "
            "same rule as the analysis set, so that no runnable stratum of this Census pin is left "
            "unlisted and selectable after the fact. A result from one of them is reported as a "
            "within-collection control, never enters the D2 denominator and never counts toward a "
            "majority. Promoting one to an independent dataset is an amendment.",
        "expected_effect_source": EXPECTED_EFFECT_SOURCE,
        "expected_effect_vocabulary": list(EXPECTED_EFFECT_VOCABULARY),
        "expected_effect_on_controls":
            f"{EXPECTED_EFFECT_NOT_APPLICABLE} — §1 (i)/(ii) is a coverage claim about the twelve; "
            "labelling a control row would smuggle it into that claim.",
        "cells_per_donor_bins": bin_occupancy(strata),
        "layer_b_truncation": layer_b,
        "layer_b_note":
            "Pre-declared, before the sigma_donor anchor exists, so the surviving subset cannot be "
            "chosen later from among the convenient ones. Amendment 3: 'a negative answer is a "
            f"live outcome of this study, not a failure mode to be designed around.' §1 asks for "
            f"{SPEC_DATASET_FLOOR}-12 datasets; the frozen list falls below that floor at three of "
            "the four tiers — every tier except the most optimistic — including sigma_donor 0.35, "
            "which is gate_config.POWER_EVAL_SIGMA, the point at which the instrument nominally "
            "operates. The study in its pre-registered form is not executable there, and that is a "
            "result to report.",
        "counterfactual_note":
            "n_manifest_datasets_at_tier and counterfactual_max_survivors_of_twelve say what the "
            "manifest could have supported. The shortfall is a consequence of selecting for §1 "
            "(iii)'s coverage axes rather than for donor counts, NOT a property of the public data "
            f"at this Census pin: {MANIFEST_TIER_CENSUS[13]} of the manifest's candidate-bearing "
            f"datasets clear 13v13 and {MANIFEST_TIER_CENSUS[23]} clear 23v23. Recorded before the "
            "anchor exists so that the list cannot be reassembled on donor counts later and called "
            "the original plan.",
        "manifest_tier_census": {str(k): v for k, v in sorted(MANIFEST_TIER_CENSUS.items())},
        "manifest_datasets_with_exact_3v3": MANIFEST_DATASETS_WITH_EXACT_3V3,
        "operating_envelope": [dict(row) for row in gc.OPERATING_ENVELOPE],
        "operating_envelope_source":
            "docs/AMENDMENTS.md Amendment 3 Change 1. SYNTHETIC: sigma_donor is an unanchored "
            "simulator knob. Envelope membership is therefore UNKNOWN for every stratum here.",
        "admission": {
            "admitted_to_sweep": False,
            "n_rows_admitted": 0,
            "blockers": ["integer_check", "frozen_universe_size", "sigma_donor_estimate",
                         "envelope_membership"],
            "note":
                "Freezing the list is not admitting it. All four blockers stand on every row: "
                "integer_check and frozen_universe_size are computed at X load (§1 items 4 and 5) "
                "and sigma_donor is unanchored (Amendment 3, OPEN). The list may still SHRINK when "
                "the counts gate runs, and a shrinkage is a reported outcome — never a re-selection.",
        },
        "row_schema": {
            "manifest_fields": list(cs.MANIFEST_FIELDS),
            "derived_fields": list(DERIVED_FIELDS),
            "derived_note":
                "The manifest's own columns are carried verbatim and appear first; the seven "
                "derived columns are appended, never interleaved, so provenance is visible per "
                f"column. 'role' is {ROLE_ANALYSIS_SET!r} on the rows of the analysis set and "
                f"{ROLE_WITHIN_COLLECTION_CONTROL!r} on the sibling rows, so a row read out of "
                "context still says which set it belongs to.",
        },
        "reproduce": "python scripts/freeze_stratum_list.py --check",
    }


def _emit_row(row: dict, derived: dict) -> dict:
    merged = {**row, **derived}
    return {field: merged[field] for field in STRATUM_FIELDS}


def build(manifest: dict) -> dict:
    strata = select_strata(manifest["rows"])
    controls = select_control_strata(manifest["rows"])
    verify_declarations(manifest, strata, controls)
    ranks = {d["dataset_id"]: d["rank"] for d in FROZEN_DATASETS}
    shorts = {d["dataset_id"]: d["short"] for d in FROZEN_DATASETS}
    sibling_shorts = {s["dataset_id"]: s["short"] for s in SIBLING_DATASETS}

    rows = [
        _emit_row(row, {
            "role": ROLE_ANALYSIS_SET,
            "dataset_rank": ranks[row["dataset_id"]],
            "dataset_short": shorts[row["dataset_id"]],
            "expected_effect": EXPECTED_EFFECT[(row["dataset_id"], row["disease"])],
            "min_donors_per_group": min(row["n_donors_A"], row["n_donors_B"]),
            "cells_per_donor_bin_A": bin_of(row["cells_per_donor_by_group"]["A"]["median"]),
            "cells_per_donor_bin_B": bin_of(row["cells_per_donor_by_group"]["B"]["median"]),
        })
        for row in strata
    ]
    # ``dataset_rank`` is null on a control row rather than borrowed from its parent: a rank is a
    # position in the D2 denominator, and these are not in it.
    control_rows = [
        _emit_row(row, {
            "role": ROLE_WITHIN_COLLECTION_CONTROL,
            "dataset_rank": None,
            "dataset_short": sibling_shorts[row["dataset_id"]],
            "expected_effect": EXPECTED_EFFECT_NOT_APPLICABLE,
            "min_donors_per_group": min(row["n_donors_A"], row["n_donors_B"]),
            "cells_per_donor_bin_A": bin_of(row["cells_per_donor_by_group"]["A"]["median"]),
            "cells_per_donor_bin_B": bin_of(row["cells_per_donor_by_group"]["B"]["median"]),
        })
        for row in controls
    ]

    return {
        "header": build_header(manifest, strata, controls),
        "rows": rows,
        "within_collection_control_rows": control_rows,
    }


def render_json(frozen: dict) -> str:
    """The exact JSON text written to disk — a function so the tests can compare without writing."""
    return json.dumps(frozen, indent=1, ensure_ascii=False) + "\n"


def write(frozen: dict, out_dir: Path, *, stem: str = OUT_STEM) -> tuple[Path, Path]:
    """Write ``{stem}.json`` and ``{stem}.csv``, with line endings pinned on every platform.

    ``newline="\\n"`` on the JSON and ``newline=""`` on the CSV (which leaves csv's own ``\\r\\n``
    terminator alone) mean the bytes do not depend on the operating system. Without that, the same
    artifact would differ between the Linux, macOS and Windows CI legs and the byte-identity test
    would be a platform coin flip.

    The CSV carries the analysis set **and** the within-collection control rows, in that order,
    told apart by the ``role`` column. A CSV holding only the 251 would silently omit frozen
    content, which is the failure the control rows exist to prevent.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{stem}.json"
    csv_path = out_dir / f"{stem}.csv"

    with open(json_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(render_json(frozen))
    with open(csv_path, "w", encoding="utf-8", newline="") as fh:
        fh.write(render_csv_rows(frozen))
    return json_path, csv_path


def render_csv_rows(frozen: dict) -> str:
    """The exact CSV text written to disk — a function so ``--check`` can compare without writing."""
    buffer = io.StringIO(newline="")
    # census_select's own encoder, imported rather than copied: the two CSVs must serialise a
    # dict-valued cell the same way or they cannot be diffed against each other.
    writer = csv.DictWriter(buffer, fieldnames=list(STRATUM_FIELDS))
    writer.writeheader()
    for row in (*frozen["rows"], *frozen["within_collection_control_rows"]):
        writer.writerow({k: cs._csv_value(row.get(k)) for k in STRATUM_FIELDS})
    return buffer.getvalue()


def summarise(frozen: dict) -> str:
    header = frozen["header"]
    lines = [
        f"# Phase 0 stratum list — frozen {header['frozen_date']}",
        f"  source      : {header['source']['json']}",
        f"                sha256 {header['source']['json_sha256']}",
        f"                CI run {header['source']['ci_run_id']}, "
        f"generated {header['source']['generated_utc']}",
        f"  census      : {header['census_version']}",
        f"  datasets    : {header['n_frozen_datasets']} independent "
        f"(+{len(header['sibling_datasets'])} siblings, excluded from the D2 denominator)",
        f"  strata      : {header['n_frozen_strata']} in the analysis set, "
        f"{header['n_within_collection_control_strata']} within-collection controls",
        "",
        "  dataset                                          rank  strata  ceiling  >=8v8  3v3",
    ]
    for block in header["frozen_datasets"]:
        lines.append(
            f"  {block['short']:<47} {block['rank']:>4}  {block['n_strata']:>6}  "
            f"{block['ceiling_min_donors_per_group']:>7}  {block['n_strata_at_least_8v8']:>5}  "
            f"{block['n_strata_exactly_3v3']:>3}"
        )
    bins = header["cells_per_donor_bins"]
    lines += [
        "",
        f"  cells-per-donor bins over {bins['n_group_medians']} group medians "
        f"({bins['min_group_median']} .. {bins['max_group_median']}):",
    ]
    for label in bins["bins"]:
        lines.append(f"    {label:<12} {bins['occupancy'][label]:>4}")
    lines.append("")
    for tier in header["layer_b_truncation"]:
        floor = (
            f" BELOW §1's {SPEC_DATASET_FLOOR}-12 floor"
            if tier["below_spec_dataset_floor"] else " met"
        )
        lines.append(
            f"  Layer B at sigma_donor {tier['sigma_donor']} (>= {tier['min_donors_per_group']}"
            f"v{tier['min_donors_per_group']}): {tier['n_datasets']} datasets, "
            f"{tier['n_strata']} strata —{floor}"
            f"  [manifest-wide: {tier['n_manifest_datasets_at_tier']} datasets clear it; a "
            f"donor-count-optimised twelve would have kept "
            f"{tier['counterfactual_max_survivors_of_twelve']}]"
        )
    below = sum(1 for t in header["layer_b_truncation"] if t["below_spec_dataset_floor"])
    lines += [
        "",
        f"  {below} of {len(header['layer_b_truncation'])} envelope tiers fall below §1's "
        f"{SPEC_DATASET_FLOOR}-dataset floor.",
        "  admitted_to_sweep: False on every row. Freezing the list is not admitting it.",
    ]
    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--source", type=Path, default=SOURCE_JSON,
        help="the committed candidate manifest; its sha256 must match the pinned value or the "
             "freeze refuses to run",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=PREREG_DIR,
        help=f"where {OUT_STEM}.json / .csv are written",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="verify that both committed artifacts — JSON and CSV — are byte-identical to a fresh "
             "run, and write nothing (exit 1 if either is not)",
    )
    return parser


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)
    try:
        manifest = load_source(args.source, csv_path=SOURCE_CSV)
        frozen = build(manifest)
    except (SourceArtifactMismatch, FrozenDeclarationMismatch) as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2

    if args.check:
        # Both halves. The CSV is committed evidence too, and a check that read only the JSON
        # would have let the CSV drift silently — which is exactly the state "regenerate it and
        # diff" is supposed to make impossible.
        for suffix, expected in (
            (".json", render_json(frozen)),
            (".csv", render_csv_rows(frozen)),
        ):
            actual_path = args.output_dir / f"{OUT_STEM}{suffix}"
            if not actual_path.exists():
                print(f"MISSING: {actual_path}", file=sys.stderr)
                return 1
            # newline="" so the bytes are compared as written, not as the platform would translate
            # them (Path.read_text grew a newline argument only in 3.13; this package targets 3.12).
            with open(actual_path, encoding="utf-8", newline="") as fh:
                actual = fh.read()
            if actual != expected:
                print(f"DRIFT: {actual_path} is not what this script now produces", file=sys.stderr)
                return 1
            print(f"OK: {actual_path} matches a fresh run")
        return 0

    json_path, csv_path = write(frozen, args.output_dir)
    print(summarise(frozen))
    print(f"\n  wrote {json_path}")
    print(f"  wrote {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
