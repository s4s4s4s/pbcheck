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

The **five** same-collection siblings are frozen by the same rule and emitted under
``role = within_collection_control``, separately from the 251. They are not named by hand: a sibling
is any candidate-bearing dataset sharing a ``collection_id`` with one of the twelve, computed here
from the pinned release table. An earlier version of this freeze named two of them and missed three,
which left 79 runnable strata unlisted and therefore selectable after the fact — the exact freedom
the control block exists to remove. Three of the five carry donor-rich designs and two of those clear
every tier of the operating envelope, so the omission was not a rounding error.

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
is recomputed from the pinned sources and compared against the declaration. A disagreement is a hard
failure with the two values printed, never a silent correction: a pre-registration that quietly
re-derives its own numbers when they stop matching is not a pre-registration.

**And the document is verified too, not merely accompanied.** Three adversarial reviews of this
freeze found the same defect three times: a falsifiable claim in the prose that the artifact
committed beside it refutes. Every one of them sat in the part of the document nothing recomputed.
So :func:`verify_document` parses the tables of ``docs/PREREGISTRATION_STRATUM_LIST.md`` back out of
the file and compares them cell by cell against what has just been derived, and the freeze aborts on
a disagreement in either direction. The rule this file now enforces is: **a falsifiable number in
that document is recomputed here from a hash-pinned source, or it is labelled in the document as a
judgement, a quotation or an external fact, or it does not appear.**

**Three sources are pinned, not one.** The candidate manifest cannot answer which collection a
dataset belongs to, what assay it was run on, or what the pinned Census release contains, because it
carries none of those columns. Those questions used to be answered by *live* reads of the CELLxGENE
Discover API, which is why every collection-level claim in the document was uncheckable and why two
of them were wrong. They are now answered from two further committed, hash-pinned snapshots written
by ``scripts/fetch_preregistration_evidence.py``: the Discover curation index as read on 2026-08-16,
and the pinned release's own ``census_info/datasets`` table. The two are **not interchangeable** and
the document says which of them each claim rests on.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
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

#: The CELLxGENE Discover curation index, read once and committed. The manifest carries no assay,
#: suspension, tissue, DOI or COLLECTION column, so every claim of that kind used to rest on a live
#: API read that no reader could reproduce — and both of the blocking defects of the third review
#: were collection claims. From this commit the document's §4.2 metadata table and everything about
#: collections is recomputed from these bytes and never from a request.
#:
#: It describes **Discover on 2026-08-16**, not the pinned release. See
#: :data:`RELEASE_VS_DISCOVER` for the measured size of that gap.
DISCOVER_INDEX = PREREG_DIR / "discover_index_2026-08-16.json"
DISCOVER_SHA256 = "afc74c1c1ea8f22c9c86a7cd6a2e4eb8087b7db58c37b447f8383e76f4eaf416"
DISCOVER_BYTES = 2_343_258
DISCOVER_READ_DATE = "2026-08-16"
DISCOVER_N_DATASETS = 2216
DISCOVER_ENDPOINT = "https://api.cellxgene.cziscience.com/curation/v1/datasets"

#: The pinned release's own ``census_info/datasets`` table, read from the public TileDB array and
#: committed. This is the release enumerating itself: it is the authority on what Census
#: ``2025-01-30`` contains and on which collection each of its datasets belonged to when it was
#: built. §8 rests on it, and so does collection membership — a property of the release, not of
#: Discover, even though the two happen to agree here.
RELEASE_DATASETS = PREREG_DIR / "census_release_datasets_2025-01-30.json"
RELEASE_SHA256 = "b60e8e1920de09ef3e1a6de595d574bddca622ea64adf86498ae14bcfa26da0e"
RELEASE_BYTES = 1_367_664
RELEASE_N_DATASETS = 1573

#: What the two indexes say about each other, declared and recomputed. The point of the block is
#: that they are **different objects**: ``dataset_version_id`` matches for none of the 1573, and six
#: release datasets have no Discover record at all — which is why "absent from Discover" never
#: implied "absent from the release" and why §8 now reads the release directly.
RELEASE_VS_DISCOVER = {
    "n_release_datasets": 1573,
    "n_resolving_in_discover_by_dataset_id": 1567,
    "n_not_resolving_in_discover": 6,
    "n_dataset_version_id_matching": 0,
    "n_collection_doi_differing": 61,
    "n_collection_id_differing": 0,
}

#: The six release datasets with no Discover record, and the two collections they fall in. Named
#: rather than counted, because §8's conclusion is that none of them is Mathys and a reader has to
#: be able to check that claim against the two blind-spot collections themselves.
RELEASE_NOT_IN_DISCOVER = (
    "33911db3-f461-464b-8083-a397ab616a09",
    "78f10833-3e61-4fad-96c9-4bbd4f14bdfa",
    "9c4c8515-8f82-4c72-b0c6-f87647b00bbe",
    "bcdec5fa-a7fa-4806-92bc-0cd02f40242f",
    "da75ce6d-a395-4abd-962b-267aadb99666",
    "f7ec7bd5-04ab-453b-a8a7-c9d14812affb",
)
RELEASE_NOT_IN_DISCOVER_COLLECTIONS = (
    "0a77d4c0-d5d0-40f0-aa1a-5e1429bcbd7e",
    "0c3f148e-02ff-4c81-8946-29beaaf5fa59",
)

#: §8(d)'s needles for Mathys 2019 / ROSMAP, searched over every retained text field of BOTH
#: snapshots. Zero hits in either is the claim, and it is recomputed here rather than reported.
MATHYS_NEEDLES = (
    "mathys",
    "rosmap",
    "religious order",
    "memory and aging",
    "s41586-019-1195",
    "10.1038/s41586-019-1195-2",
)

#: The proposal document the choice of twelve was made from, committed beside the manifest and
#: pinned here for the same reason: it is referenced by §3.2, §4.2 and §8 of the document — for the
#: per-dataset rationale, the rejected datasets, the FIVE NAMED RESERVES, the third 5' candidate and
#: the Mathys search method — and a justification nobody can open is not a justification.
#:
#: It is NOT part of the binding act and nothing here derives from it: what binds is
#: :data:`FROZEN_DATASETS` and :data:`SELECTION_RULE`. It is pinned so that a substitution from its
#: reserve list would be detectable, never so that one becomes available (§9 item 7).
#:
#: **What is committed is a redacted copy, and the redaction is itself auditable.** The circulated
#: copy carried an absolute filesystem path from the author's Windows account in six places. Those
#: six occurrences of the repository root — one backslash form, five forward-slash — are replaced by
#: the literal ``<REPO>`` and **nothing else is altered**. Both hashes are recorded: the circulated
#: copy's, so anyone holding it can reproduce the substitution, and the redacted copy's, which is
#: what this script checks. The circulated copy is deliberately NOT in the repository.
PROPOSAL_MD = PREREG_DIR / "stratum_list_proposal_2026-08-16.redacted.md"
PROPOSAL_SHA256 = "5588ba845cc144b56ea27a25ca4599f3fbc33d69c5eeb8833da7772a9459f07d"
PROPOSAL_BYTES = 93_873
PROPOSAL_ORIGINAL_SHA256 = "50872414b0727c129a824b0c65ed179674ac5d6c9ecaac53327568b3eae6fb48"
PROPOSAL_ORIGINAL_BYTES = 92_589
PROPOSAL_REDACTION = (
    "Absolute filesystem paths, and only those: six occurrences of the repository root (one in the "
    "backslash form, five in the forward-slash form) replaced by the literal '<REPO>'. No prose, "
    "number, code block or line break is otherwise changed. The circulated copy's sha256 is "
    "recorded beside the redacted copy's so the substitution is checkable by anyone holding it; "
    "the circulated copy is not committed, because the point of the redaction is that the account "
    "path does not enter the repository."
)
PROPOSAL_STATUS = (
    "SUPERSEDED WORKING DOCUMENT. The reasoning behind the choice of twelve, committed for "
    "auditability. NOT part of the binding act: the binding content is the frozen list and the "
    "selection rule. It has known errors, enumerated with their correct values in §10 of the "
    "document; where the two disagree, the document governs. Committed as a redacted copy — see "
    "PROPOSAL_REDACTION — because it carried an absolute account path, and otherwise unmodified."
)

#: The banner the redacted copy must open with. Checked, so that a reader who opens the file
#: directly on GitHub cannot mistake a superseded proposal for a live request for a decision.
PROPOSAL_BANNER_MARKERS = (
    "SUPERSEDED WORKING DOCUMENT",
    "NOT A LIVE REQUEST FOR A DECISION",
    "docs/PREREGISTRATION_STRATUM_LIST.md",
    "It has known errors",
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

#: Same-collection datasets, **computed and not named by hand**. A sibling is any candidate-bearing
#: dataset of the manifest that shares a ``collection_id`` with one of the twelve, and the set is
#: re-derived from the pinned release table on every run: the declaration below is compared against
#: that derivation and a disagreement aborts the freeze in either direction — a sibling missing from
#: the declaration as loudly as one invented.
#:
#: D2 clusters evidence by dataset because same-dataset strata share donors, batch and assay; two
#: datasets of one collection share the cohort and the laboratory as well, so counting both toward
#: "majority of independent datasets" would inflate the effective n.
#:
#: They are **frozen too**, and by the same rule: every candidate row of theirs is emitted, under
#: ``role = within_collection_control``. **This is where the third review found a blocker.** An
#: earlier version of this file named two of them — SEA-AD MTG and the emphysema immune split — and
#: missed three, because collection membership was not derivable from any committed source and was
#: therefore typed from a live API read. The three that were missed hold 79 further runnable strata,
#: and two of them (Yoshida Airway at 30 donors per group, the KPMP single-cell arm at 26) clear
#: every tier of the operating envelope. Freezing them is what stops that set from being reachable
#: after the results are in.
SIBLING_DATASETS = (
    {
        "dataset_id": "c2876b1b-06d8-4d96-a56b-5304f815b99a",
        "short": "SEA-AD MTG",
        "sibling_of": "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3",
        "collection_id": "1ca90a2d-2943-483d-b678-b809bf464c30",
        "role": ROLE_WITHIN_COLLECTION_CONTROL,
        "role_note": "middle temporal gyrus instead of DLPFC: a different cortical region of the "
                     "same cohort and laboratory",
        "n_strata": 18,
        "ceiling_min_donors_per_group": 42,
    },
    {
        "dataset_id": "edc8d3fe-153c-4e3d-8be0-2108d30f8d70",
        "short": "Yoshida 2022, airway",
        "sibling_of": "2a498ace-872a-4935-984b-1afa70fd9886",
        "collection_id": "03f821b4-87be-4ff4-b65a-b5fc00061da7",
        "role": ROLE_WITHIN_COLLECTION_CONTROL,
        "role_note": "the airway compartment (bronchus, nasal cavity, trachea) of the same "
                     "SARS-CoV-2 cohort whose PBMC arm is #4",
        "n_strata": 25,
        "ceiling_min_donors_per_group": 30,
    },
    {
        "dataset_id": "8f4f8502-9170-4ac2-9707-3b6985ebfe5f",
        "short": "CAREBANK right atrium",
        "sibling_of": "f1606894-59df-4794-a37f-baa7c6fb6de1",
        "collection_id": "8c782494-01ed-491b-97b9-6f0d3b76c676",
        "role": ROLE_WITHIN_COLLECTION_CONTROL,
        "role_note": "the CAREBANK cohort of the same right-atrium collection whose PERIHEART "
                     "cohort is #6; same tissue, same assay, same suspension",
        "n_strata": 11,
        "ceiling_min_donors_per_group": 15,
    },
    {
        "dataset_id": "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d",
        "short": "KPMP adult human kidney scRNA-seq v1.5",
        "sibling_of": "a12ccb9b-4fbe-457d-8590-ac78053259ef",
        "collection_id": "0f528c8a-a25c-4840-8fa3-d156fa11086f",
        "role": ROLE_WITHIN_COLLECTION_CONTROL,
        "role_note": "the single-CELL arm of the same KPMP kidney atlas whose single-NUCLEUS arm "
                     "is #8 — same consortium, same donors' tissue, different dissociation",
        "n_strata": 43,
        "ceiling_min_donors_per_group": 26,
    },
    {
        "dataset_id": "1e5bd3b8-6a0e-4959-8d69-cafed30fe814",
        "short": "Emphysema Cell Atlas, immune cells",
        "sibling_of": "4b6af54a-4a21-46e0-bc8d-673c0561a836",
        "collection_id": "03cdc7f4-bd08-49d0-a395-4487c0e5a168",
        "role": ROLE_WITHIN_COLLECTION_CONTROL,
        "role_note": "immune compartment of the same 3v3 cohort whose non-immune compartment "
                     "is #12",
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
# What the manifest as a whole could have supported (§6's trade-off).
# ---------------------------------------------------------------------------

#: How many of the manifest's 68 candidate-bearing datasets hold at least one stratum at each
#: envelope tier's donors-per-group demand, keyed by that demand. §6 rests on these: they are what
#: makes the truncation a consequence of the SELECTION rather than of the public data.
MANIFEST_TIER_CENSUS = {4: 62, 8: 33, 13: 21, 23: 12}

#: The same census counted in **collections** rather than datasets, which is the honest unit: two
#: datasets of one collection are not two independent choices, and a dataset count silently offers
#: slots that D2 would refuse. This column is why the earlier constructed witness was wrong — it
#: put both SEA-AD datasets and both KPMP datasets in one twelve and double-counted two collections.
MANIFEST_TIER_COLLECTIONS = {4: 46, 8: 25, 13: 15, 23: 10}

#: The denominators: the manifest's candidate-bearing datasets, and the collections they fall in.
MANIFEST_N_CANDIDATE_DATASETS = 68
MANIFEST_N_CANDIDATE_COLLECTIONS = 50

#: §1 (iii) requires "some exactly 3v3". Only this many of the manifest's datasets — in this many
#: collections — hold an exactly-3v3 stratum at all, which is why the requirement costs a slot that
#: could otherwise have gone to a donor-rich dataset.
MANIFEST_DATASETS_WITH_EXACT_3V3 = 3
MANIFEST_COLLECTIONS_WITH_EXACT_3V3 = 2

#: **There is deliberately no constructed witness here any more.** The previous version of §6 built
#: an explicit twelve-dataset list to prove the counterfactual attainable; it double-counted two
#: collections, its headline figure was wrong by one, and it leaned on a dataset the shipped
#: covariate-less arm cannot analyse. The retention figures above and in :data:`LAYER_B` say
#: everything the section needs, and every one of them is a count over the manifest rather than a
#: construction — so there is nothing left in §6 for a later reader to select from, and nothing left
#: that has to be re-derived to stay true.
COUNTERFACTUAL_REMOVED_NOTE = (
    "§6 states counts, not a constructed list. An earlier version built a twelve-dataset witness "
    "to show the counterfactual attainable; it contained both SEA-AD datasets and both KPMP "
    "datasets, so it double-counted two collections and overstated its own >=23v23 figure. It is "
    "deleted rather than repaired: the plain tier counts, in datasets AND in collections, carry "
    "the argument and cannot be selected from."
)

# ---------------------------------------------------------------------------
# Figures the document states in prose, recomputed here so that none of them is unchecked.
# ---------------------------------------------------------------------------

#: Every load-bearing number §3.2, §4.2, §5, §7, §9 and §10 state about the frozen set, declared and
#: recomputed. They were correct when a reviewer checked them by hand; the point of moving them here
#: is that nothing else in this document's history stayed correct without being recomputed.
#:
#: ``median_counts_per_cell`` needs its population stated, because three plausible populations give
#: three different numbers: it is the median of the 502 values of
#: ``median_counts_per_cell_by_group[A|B]`` over the 251 frozen strata, two per stratum. The
#: per-stratum mean of A and B gives 3126.75 and the whole candidate set gives 3454.0.
ATTESTED = {
    "n_cells": 4_609_595,
    "n_dataset_cell_type_strata": 182,
    "n_cell_type_labels": 124,
    "n_disease_terms": 15,
    "residual_df_min": 4,
    "residual_df_max": 108,
    "permutation_count_min": 20,
    "permutation_count_max": 727_646_193_812_764_637_422_200,
    "n_strata_below_1000_permutations": 45,
    "cells_per_donor_min": 10,
    "cells_per_donor_max": 16_383,
    "group_median_cells_per_donor_min": 11.0,
    "group_median_cells_per_donor_max": 6671.5,
    "counts_per_cell_min": 284.0,
    "counts_per_cell_max": 56_841.5,
    "median_counts_per_cell": 3081.0,
    "n_contrasts_in_frozen_datasets": 306,
    "n_inclusion_gate_failures_in_frozen_datasets": 55,
    "n_confound_exclusions_in_frozen_datasets": 0,
    "n_tagged_sequencing_depth_partial": 231,
    "n_tagged_assay_partial": 24,
    "n_tagged_near_confound": 8,
    "n_strata_exactly_3v3": 10,
    "n_strata_at_least_8v8": 150,
    "n_datasets_at_least_8v8": 7,
    "n_strata_min_donors_3": 24,
    "n_strata_min_donors_4_to_7": 77,
    "n_strata_below_8v8": 101,
    "largest_design": [39, 44],
    "min_bins_occupied_by_a_dataset": 4,
    "n_datasets_occupying_all_bins": 3,
    "kpmp_strata_at_least_23v23": 6,
    "rexach_assay_v_n_strata": 27,
    "rexach_assay_v_n_non_null": 24,
    "rexach_assay_v_min": 0.029,
    "rexach_assay_v_max": 0.085,
    "seaad_assay_v_n_strata": 18,
    "seaad_assay_v_n_non_null": 0,
    "n_control_strata_exactly_3v3": 9,
    # The corrected side of §10's list of discrepancies in the superseded proposal. The wrong value
    # in each is a quotation from a hash-pinned file; the right value is one of these, recomputed.
    "rexach_strata_exactly_3v3": 2,
    "n_rows_failing_donor_nesting": 284,
    "n_rows_with_a_gate_status_blocker": 993,
    "n_candidate_rows_tagged_sequencing_depth": 1172,
    "n_candidate_rows_tagged_assay": 443,
    "n_candidate_rows_tagged_tissue_general": 269,
    "n_candidate_rows_tagged_suspension_type": 173,
    "candidate_group_median_cells_per_donor_median": 96.25,
    "n_datasets_beyond_the_twelve_at_13v13": 16,
    "seaad_dlpfc_n_donors_B_min": 40,
    "seaad_dlpfc_n_donors_B_max": 44,
    "seaad_mtg_n_donors_B_min": 34,
    "seaad_mtg_n_donors_B_max": 46,
    "yoshida_sequencing_depth_v_min": 0.175,
    "yoshida_sequencing_depth_v_max": 0.816,
    "smallest_frozen_dataset_rank": 12,
    "smallest_frozen_dataset_cell_count": 18_386,
    "second_smallest_frozen_dataset_rank": 11,
    "second_smallest_frozen_dataset_cell_count": 22_502,
    "n_candidate_datasets_with_drop_seq": 2,
    "n_candidate_datasets_with_seq_well": 2,
    # Prose figures elsewhere in the document that nothing else recomputed.
    "n_strata_in_datasets_at_least_4v4": 243,
    "combat_covid_permutation_count_max": 46_897_636_623_981,
    "combat_max_n_donors_A": 100,
    "combat_max_n_donors_B": 10,
    "melms_max_n_donors_A": 20,
    "melms_max_n_donors_B": 7,
    "wang_near_confound_strata": 3,
}

#: The KPMP dataset ids, singled out because both of them — #8 and its sibling — carry
#: ``is_primary_data = [False]`` in the 2026-08-16 Discover snapshot while the manifest's own Census
#: query filtered on ``is_primary_data == True``. Recorded as a measured discrepancy between the
#: two sources, not explained away: ``obs`` at the pin returned cells for both under that filter,
#: and Discover today says neither holds a primary cell. Nothing here decides which is right.
DISCOVER_NOT_PRIMARY = (
    "a12ccb9b-4fbe-457d-8590-ac78053259ef",
    "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d",
)
MANIFEST_VALUE_FILTER = "is_primary_data == True and disease != 'na'"

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
    check_proposal_redaction()

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


def check_proposal_redaction() -> None:
    """The redacted proposal must carry its banner and no absolute path at all.

    Two failures are worth catching here rather than at review time. A redaction that missed an
    occurrence would put the account path back in the repository, and the check for that is
    mechanical: no drive-letter path may survive anywhere in the file. And a superseded proposal
    with no banner reads, to anyone who opens it on GitHub, as a live request for a decision that
    was in fact taken.
    """
    text = PROPOSAL_MD.read_text(encoding="utf-8")
    leaked = re.findall(r"[A-Za-z]:[\\/](?:Users|home)[\\/]\S*", text)
    if leaked:
        raise SourceArtifactMismatch(
            f"{PROPOSAL_MD.name} still contains {len(leaked)} absolute filesystem path(s) — the "
            f"redaction is incomplete: {sorted(set(leaked))}"
        )
    missing = [marker for marker in PROPOSAL_BANNER_MARKERS if marker not in text[:4000]]
    if missing:
        raise SourceArtifactMismatch(
            f"{PROPOSAL_MD.name} does not open with the superseded-document banner; missing "
            f"{missing}. A reader opening it directly must not mistake it for a live proposal."
        )


def load_discover(path: Path = DISCOVER_INDEX) -> dict:
    """The pinned Discover snapshot. Read for collections, assay, suspension, tissue and DOI.

    Hash-checked before it is parsed, exactly like the manifest, and for the same reason: until this
    file existed those five properties came from a live API read, which is why no reader could
    check them and why two of them were wrong.
    """
    check_pinned_file(path, DISCOVER_SHA256, DISCOVER_BYTES, "Discover index snapshot")
    snapshot = json.loads(path.read_text(encoding="utf-8"))
    header = snapshot["header"]
    if header.get("read_date") != DISCOVER_READ_DATE:
        raise SourceArtifactMismatch(
            f"Discover snapshot read_date {header.get('read_date')!r}, expected "
            f"{DISCOVER_READ_DATE!r}"
        )
    if header.get("endpoint") != DISCOVER_ENDPOINT:
        raise SourceArtifactMismatch(
            f"Discover snapshot endpoint {header.get('endpoint')!r}, expected {DISCOVER_ENDPOINT!r}"
        )
    if header.get("n_datasets_in_index") != DISCOVER_N_DATASETS:
        raise SourceArtifactMismatch(
            f"Discover index held {header.get('n_datasets_in_index')} datasets, declared "
            f"{DISCOVER_N_DATASETS}"
        )
    if len(snapshot["datasets"]) != DISCOVER_N_DATASETS:
        raise SourceArtifactMismatch(
            f"Discover snapshot carries {len(snapshot['datasets'])} records, expected "
            f"{DISCOVER_N_DATASETS}"
        )
    return snapshot


def load_release(path: Path = RELEASE_DATASETS) -> dict:
    """The pinned release's own dataset table — the authority on what Census 2025-01-30 contains."""
    check_pinned_file(path, RELEASE_SHA256, RELEASE_BYTES, "Census release dataset table")
    snapshot = json.loads(path.read_text(encoding="utf-8"))
    header = snapshot["header"]
    if header.get("census_version") != cs.CENSUS_VERSION:
        raise SourceArtifactMismatch(
            f"release table census_version {header.get('census_version')!r}, expected "
            f"{cs.CENSUS_VERSION!r} (spec §1 pin)"
        )
    if header.get("n_datasets") != RELEASE_N_DATASETS or len(snapshot["datasets"]) != (
        RELEASE_N_DATASETS
    ):
        raise SourceArtifactMismatch(
            f"release table holds {len(snapshot['datasets'])} datasets (header says "
            f"{header.get('n_datasets')}), expected {RELEASE_N_DATASETS}"
        )
    return snapshot


def collections_of(release: dict) -> dict[str, str]:
    """``{dataset_id: collection_id}`` from the pinned release, which is the collection at the pin.

    Discover carries the same mapping and agrees with it for every dataset the two share — that
    agreement is asserted in :func:`verify_external_sources` rather than assumed, because the two
    snapshots demonstrably disagree about other fields.
    """
    return {row["dataset_id"]: row["collection_id"] for row in release["datasets"]}


def derive_siblings(candidates, collections: dict[str, str]) -> dict[str, str]:
    """``{sibling_dataset_id: frozen_dataset_id}`` — the sibling set, computed, never typed.

    A sibling is a candidate-bearing dataset of the manifest that shares a ``collection_id`` with
    one of the twelve and is not itself one of the twelve. This function is the whole definition;
    :data:`SIBLING_DATASETS` is a declaration checked against it.
    """
    frozen_ids = [d["dataset_id"] for d in FROZEN_DATASETS]
    frozen_by_collection = {}
    for dataset_id in frozen_ids:
        if dataset_id not in collections:
            raise FrozenDeclarationMismatch(
                f"frozen dataset {dataset_id} is absent from the pinned release table; the freeze "
                "cannot establish which collection it belongs to"
            )
        frozen_by_collection.setdefault(collections[dataset_id], []).append(dataset_id)

    siblings = {}
    for dataset_id in sorted({row["dataset_id"] for row in candidates}):
        if dataset_id in frozen_ids:
            continue
        parents = frozen_by_collection.get(collections.get(dataset_id))
        if parents:
            siblings[dataset_id] = parents[0]
    return siblings


def _needle_hits(text: str) -> list[str]:
    lowered = text.lower()
    return [needle for needle in MATHYS_NEEDLES if needle in lowered]


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
    """The same rule, applied to the five siblings — they are frozen, not merely named.

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


def assert_sets_disjoint(strata: list[dict], controls: list[dict]) -> None:
    """The analysis set and the control set may not share a stratum.

    A named function rather than an inline check so that the guard itself is directly testable: one
    set is the D2 denominator and the other is explicitly outside it, and a row in both would let
    the denominator acquire a sibling by accident.
    """
    key = lambda row: (row["dataset_id"], row["cell_type"], row["disease"])  # noqa: E731
    shared = sorted({key(r) for r in controls} & {key(r) for r in strata})
    _require(
        not shared,
        f"{len(shared)} within-collection control strata are also in the analysis set "
        f"({shared[:3]}); the two sets must stay distinguishable or the D2 denominator can "
        "acquire a sibling by accident",
    )


def _tier_collections(candidates, threshold: int, collections: dict[str, str]) -> set[str]:
    """The **collections** holding a candidate stratum at ``threshold`` donors in both groups.

    The unit that matters. A dataset count offers slots D2 would refuse: two datasets of one
    collection share the cohort and the laboratory, so a list that took both would not be two
    independent choices. Counting collections is what makes §6's trade-off statement honest, and
    the absence of this function is what let the previous version of §6 publish a twelve-dataset
    witness containing two SEA-AD datasets and two KPMP datasets.
    """
    return {collections[dataset_id] for dataset_id in _tier_datasets(candidates, threshold)}


def _envelope_min_donors(sigma_donor: float) -> int:
    for row in gc.OPERATING_ENVELOPE:
        if float(row["sigma_donor"]) == sigma_donor:
            return int(row["min_donors_per_group"])
    raise FrozenDeclarationMismatch(
        f"sigma_donor={sigma_donor} is not a row of gate_config.OPERATING_ENVELOPE; Layer B may "
        "only be declared at sigma values the envelope itself states"
    )


def verify_external_sources(manifest: dict, discover: dict, release: dict) -> dict[str, str]:
    """Check the two pinned snapshots against each other and against every claim drawn from them.

    Returns the ``{dataset_id: collection_id}`` map the rest of the freeze uses. Everything the
    document says about collections, assay, suspension, tissue, DOI, the pinned release's contents
    and the Mathys search is settled here, from committed bytes.
    """
    discover_by_id = {row["dataset_id"]: row for row in discover["datasets"]}
    release_by_id = {row["dataset_id"]: row for row in release["datasets"]}
    _require(
        len(discover_by_id) == DISCOVER_N_DATASETS,
        f"the Discover snapshot holds {len(discover_by_id)} distinct dataset ids, expected "
        f"{DISCOVER_N_DATASETS}",
    )
    _require(
        len(release_by_id) == RELEASE_N_DATASETS,
        f"the release table holds {len(release_by_id)} distinct dataset ids, expected "
        f"{RELEASE_N_DATASETS}",
    )

    # --- the two indexes are different objects, and the document rests on how different ---
    resolving = [row for row in release["datasets"] if row["dataset_id"] in discover_by_id]
    discover_version_ids = {row["dataset_version_id"] for row in discover["datasets"]}
    measured = {
        "n_release_datasets": len(release["datasets"]),
        "n_resolving_in_discover_by_dataset_id": len(resolving),
        "n_not_resolving_in_discover": len(release["datasets"]) - len(resolving),
        "n_dataset_version_id_matching": sum(
            1 for row in release["datasets"] if row["dataset_version_id"] in discover_version_ids
        ),
        "n_collection_doi_differing": sum(
            1 for row in resolving
            if (row["collection_doi"] or None)
            != (discover_by_id[row["dataset_id"]].get("collection_doi") or None)
        ),
        "n_collection_id_differing": sum(
            1 for row in resolving
            if row["collection_id"] != discover_by_id[row["dataset_id"]]["collection_id"]
        ),
    }
    for key, declared in RELEASE_VS_DISCOVER.items():
        _require(
            measured[key] == declared,
            f"release vs Discover: {key} is {measured[key]}, declared {declared}",
        )
    not_resolving = sorted(
        row["dataset_id"] for row in release["datasets"] if row["dataset_id"] not in discover_by_id
    )
    _require(
        tuple(not_resolving) == RELEASE_NOT_IN_DISCOVER,
        f"the release datasets with no Discover record are {not_resolving}, declared "
        f"{list(RELEASE_NOT_IN_DISCOVER)}",
    )
    blind_spot_collections = sorted(
        {release_by_id[dataset_id]["collection_id"] for dataset_id in not_resolving}
    )
    _require(
        tuple(blind_spot_collections) == RELEASE_NOT_IN_DISCOVER_COLLECTIONS,
        f"the Discover blind spot falls in collections {blind_spot_collections}, declared "
        f"{list(RELEASE_NOT_IN_DISCOVER_COLLECTIONS)}",
    )

    # --- §8: Mathys is in neither index, recomputed over both rather than reported ---
    for label, snapshot in (("Discover snapshot", discover), ("release table", release)):
        hits = {
            row["dataset_id"]: found
            for row in snapshot["datasets"]
            if (found := _needle_hits(json.dumps(row, ensure_ascii=False)))
        }
        _require(
            not hits,
            f"the {label} matches a Mathys/ROSMAP needle on {sorted(hits)}; §8's conclusion that "
            "the oracle is absent no longer holds",
        )

    # --- the manifest's datasets exist in both, which is what makes the comparison meaningful ---
    manifest_ids = {row["dataset_id"] for row in manifest["rows"]}
    missing_release = sorted(manifest_ids - set(release_by_id))
    missing_discover = sorted(manifest_ids - set(discover_by_id))
    _require(
        not missing_release,
        f"manifest datasets absent from the pinned release table: {missing_release}",
    )
    _require(
        not missing_discover,
        f"manifest datasets absent from the Discover snapshot: {missing_discover}",
    )

    collections = collections_of(release)
    disagreeing = sorted(
        dataset_id for dataset_id in manifest_ids
        if collections[dataset_id] != discover_by_id[dataset_id]["collection_id"]
    )
    _require(
        not disagreeing,
        "the release table and the Discover snapshot disagree about the collection of "
        f"{disagreeing}; collection membership cannot be asserted while they do",
    )

    # --- §4.2's metadata table: declared per dataset, checked against Discover ---
    for dataset in FROZEN_DATASETS:
        record = discover_by_id[dataset["dataset_id"]]
        for field, discover_field in (
            ("assay", "assay"), ("suspension", "suspension_type"), ("tissue", "tissue"),
        ):
            declared = tuple(sorted(dataset[field]))
            actual = tuple(sorted(record[discover_field]))
            _require(
                declared == actual,
                f"{dataset['short']}: Discover gives {discover_field} {list(actual)}, declared "
                f"{list(declared)}",
            )
        actual_doi = record.get("collection_doi") or None
        _require(
            (dataset["doi"] or None) == actual_doi,
            f"{dataset['short']}: Discover gives collection_doi {actual_doi!r}, declared "
            f"{dataset['doi']!r}",
        )

    # --- the recorded is_primary_data discrepancy, verified rather than remembered ---
    frozen_and_siblings = [d["dataset_id"] for d in FROZEN_DATASETS] + [
        s["dataset_id"] for s in SIBLING_DATASETS
    ]
    not_primary = tuple(
        dataset_id for dataset_id in frozen_and_siblings
        if discover_by_id[dataset_id].get("is_primary_data") == [False]
    )
    _require(
        sorted(not_primary) == sorted(DISCOVER_NOT_PRIMARY),
        f"Discover marks {sorted(not_primary)} as holding no primary data, declared "
        f"{sorted(DISCOVER_NOT_PRIMARY)}",
    )
    _require(
        manifest["header"].get("value_filter") == MANIFEST_VALUE_FILTER,
        f"the manifest's value_filter is {manifest['header'].get('value_filter')!r}, declared "
        f"{MANIFEST_VALUE_FILTER!r}; the is_primary_data discrepancy is stated against that filter",
    )
    return collections


def verify_declarations(
    manifest: dict, strata: list[dict], controls: list[dict], collections: dict[str, str]
) -> None:
    """Recompute every declared figure from the pinned sources and abort on any disagreement."""
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
    _require(
        len(sibling_ids) == len(SIBLING_DATASETS),
        "SIBLING_DATASETS repeats a dataset_id",
    )

    # The twelve occupy twelve distinct collections. This is what makes the D2 denominator clean:
    # if two of them shared a collection, "independent datasets" would already be counting one
    # cohort twice before any sibling was considered.
    frozen_collections = [collections[dataset_id] for dataset_id in ids]
    _require(
        len(set(frozen_collections)) == 12,
        "the twelve frozen datasets occupy "
        f"{len(set(frozen_collections))} distinct collections, not twelve; the D2 denominator "
        "would count a collection twice",
    )

    # The sibling set is DERIVED and the declaration is compared against it, in both directions.
    # A sibling missing from the declaration is the defect that shipped: three were missed, and the
    # 79 strata they hold were runnable, unlisted and therefore selectable after the fact.
    derived_siblings = derive_siblings(candidates, collections)
    declared_siblings = {s["dataset_id"]: s["sibling_of"] for s in SIBLING_DATASETS}
    _require(
        derived_siblings == declared_siblings,
        "the sibling set derived from the pinned release table is "
        f"{sorted(derived_siblings)} (parents {sorted(derived_siblings.values())}), the "
        f"declaration is {sorted(declared_siblings)}. Every candidate-bearing dataset sharing a "
        "collection with one of the twelve must be frozen as a within-collection control.",
    )
    for sibling in SIBLING_DATASETS:
        _require(
            collections[sibling["dataset_id"]] == sibling["collection_id"],
            f"{sibling['short']}: the release table puts it in collection "
            f"{collections[sibling['dataset_id']]}, declared {sibling['collection_id']}",
        )
        _require(
            collections[sibling["sibling_of"]] == sibling["collection_id"],
            f"{sibling['short']}: its declared parent is not in collection "
            f"{sibling['collection_id']}",
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
    assert_sets_disjoint(strata, controls)
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
        sorted(MANIFEST_TIER_CENSUS) == sorted(MANIFEST_TIER_COLLECTIONS) == sorted(
            int(row["min_donors_per_group"]) for row in gc.OPERATING_ENVELOPE
        ),
        "the manifest tier census must be keyed by the envelope's own donors-per-group demands",
    )
    for threshold, declared in MANIFEST_TIER_COLLECTIONS.items():
        derived_n = len(_tier_collections(candidates, threshold, collections))
        _require(
            derived_n == declared,
            f"the manifest's datasets with a stratum at >= {threshold}v{threshold} fall in "
            f"{derived_n} distinct collections, declared {declared}",
        )
    _require(
        len({r["dataset_id"] for r in candidates}) == MANIFEST_N_CANDIDATE_DATASETS,
        f"the manifest holds {len({r['dataset_id'] for r in candidates})} candidate-bearing "
        f"datasets, declared {MANIFEST_N_CANDIDATE_DATASETS}",
    )
    candidate_collections = {collections[r["dataset_id"]] for r in candidates}
    _require(
        len(candidate_collections) == MANIFEST_N_CANDIDATE_COLLECTIONS,
        f"the candidate-bearing datasets fall in {len(candidate_collections)} collections, "
        f"declared {MANIFEST_N_CANDIDATE_COLLECTIONS}",
    )

    anchors = _datasets_with_exact_3v3(candidates)
    _require(
        len(anchors) == MANIFEST_DATASETS_WITH_EXACT_3V3,
        f"{len(anchors)} candidate-bearing datasets hold an exactly-3v3 stratum, declared "
        f"{MANIFEST_DATASETS_WITH_EXACT_3V3}",
    )
    anchor_collections = {collections[dataset_id] for dataset_id in anchors}
    _require(
        len(anchor_collections) == MANIFEST_COLLECTIONS_WITH_EXACT_3V3,
        f"the exactly-3v3 datasets fall in {len(anchor_collections)} collections, declared "
        f"{MANIFEST_COLLECTIONS_WITH_EXACT_3V3}",
    )
    for threshold in (13, 23):
        _require(
            not (anchors & _tier_datasets(candidates, threshold)),
            f"a dataset holding an exactly-3v3 stratum also clears {threshold}v{threshold}; §6's "
            "statement that the 3v3 anchor costs a slot at the hard tiers no longer holds",
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


def _partial_tag(row: dict, covariate: str) -> bool:
    """A ``partial_confound`` tag specifically — not a near-confound, not a set-aside note."""
    return any(
        flag.startswith(f"{covariate}:partial_confound") for flag in row["confound_flags"]
    )


def _flagged(row: dict, covariate: str) -> bool:
    """Any confound flag naming ``covariate``, at any level.

    Distinguished from :func:`_partial_tag` because the superseded proposal conflated the two:
    it counted flag *occurrences* at any level and printed them under a heading that said rows.
    """
    return any(flag.startswith(f"{covariate}:") for flag in row["confound_flags"])


def measure_attested(
    manifest: dict, strata: list[dict], controls: list[dict], discover: dict
) -> dict:
    """Recompute every figure of :data:`ATTESTED` from the pinned sources.

    §10 used to *attest* these — a human saying they had checked. Eleven of them appeared nowhere
    else, so nothing would have noticed them drifting. They are measurements now.
    """
    rows = manifest["rows"]
    candidates = [r for r in rows if r["gate_status"] == "candidate"]
    frozen_ids = {d["dataset_id"] for d in FROZEN_DATASETS}
    in_frozen_datasets = [r for r in rows if r["dataset_id"] in frozen_ids]
    discover_by_id = {row["dataset_id"]: row for row in discover["datasets"]}

    def group_values(source, key):
        return [row[key][group] for row in source for group in ("A", "B")]

    medians = [v["median"] for v in group_values(strata, "cells_per_donor_by_group")]
    counts_per_cell = group_values(strata, "median_counts_per_cell_by_group")
    permutations = [r["permutation_count"] for r in strata]
    min_donors = [min(r["n_donors_A"], r["n_donors_B"]) for r in strata]

    bins_per_dataset: dict[str, set] = {}
    for row in strata:
        for group in ("A", "B"):
            bins_per_dataset.setdefault(row["dataset_id"], set()).add(
                bin_of(row["cells_per_donor_by_group"][group]["median"])
            )

    rexach = [r for r in strata if r["dataset_id"] == "ac0c6561-7a48-4185-af6f-af799f699172"]
    rexach_v = [r["confound_cramers_v"]["assay"] for r in rexach]
    rexach_non_null = [v for v in rexach_v if v is not None]
    seaad = [r for r in strata if r["dataset_id"] == "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3"]
    mtg = [r for r in controls if r["dataset_id"] == "c2876b1b-06d8-4d96-a56b-5304f815b99a"]
    yoshida = [r for r in strata if r["dataset_id"] == "2a498ace-872a-4935-984b-1afa70fd9886"]
    yoshida_depth = [
        r["confound_cramers_v"]["sequencing_depth_bin"] for r in yoshida
        if r["confound_cramers_v"]["sequencing_depth_bin"] is not None
    ]
    kpmp = [r for r in strata if r["dataset_id"] == "a12ccb9b-4fbe-457d-8590-ac78053259ef"]
    combat = [r for r in strata if r["dataset_id"] == "ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3"]
    melms = [r for r in strata if r["dataset_id"] == "d8da613f-e681-4c69-b463-e94f5e66847f"]
    wang = [r for r in strata if r["dataset_id"] == "4b6af54a-4a21-46e0-bc8d-673c0561a836"]

    candidate_medians = sorted(
        v["median"] for v in group_values(candidates, "cells_per_donor_by_group")
    )
    middle = len(candidate_medians) // 2
    candidate_median = (
        candidate_medians[middle] if len(candidate_medians) % 2
        else (candidate_medians[middle - 1] + candidate_medians[middle]) / 2
    )

    by_cells = sorted(
        ((discover_by_id[d["dataset_id"]]["cell_count"], d["rank"]) for d in FROZEN_DATASETS)
    )
    candidate_dataset_ids = {r["dataset_id"] for r in candidates}

    def datasets_with_assay(label: str) -> int:
        return sum(
            1 for dataset_id in candidate_dataset_ids
            if label in discover_by_id[dataset_id]["assay"]
        )

    return {
        "n_cells": sum(r["n_cells"] for r in strata),
        "n_dataset_cell_type_strata": len({(r["dataset_id"], r["cell_type"]) for r in strata}),
        "n_cell_type_labels": len({r["cell_type"] for r in strata}),
        "n_disease_terms": len({r["disease"] for r in strata}),
        "residual_df_min": min(r["residual_df"] for r in strata),
        "residual_df_max": max(r["residual_df"] for r in strata),
        "permutation_count_min": min(permutations),
        "permutation_count_max": max(permutations),
        "n_strata_below_1000_permutations": sum(1 for v in permutations if v < 1000),
        "cells_per_donor_min": min(
            v["min"] for v in group_values(strata, "cells_per_donor_by_group")
        ),
        "cells_per_donor_max": max(
            v["max"] for v in group_values(strata, "cells_per_donor_by_group")
        ),
        "group_median_cells_per_donor_min": min(medians),
        "group_median_cells_per_donor_max": max(medians),
        "counts_per_cell_min": min(counts_per_cell),
        "counts_per_cell_max": max(counts_per_cell),
        "median_counts_per_cell": sorted(counts_per_cell)[len(counts_per_cell) // 2],
        "n_contrasts_in_frozen_datasets": len(in_frozen_datasets),
        "n_inclusion_gate_failures_in_frozen_datasets": sum(
            1 for r in in_frozen_datasets if r["gate_status"] == "excluded_inclusion_gate"
        ),
        "n_confound_exclusions_in_frozen_datasets": sum(
            1 for r in in_frozen_datasets if r["gate_status"] == "excluded_confound"
        ),
        "n_tagged_sequencing_depth_partial": sum(
            1 for r in strata if _partial_tag(r, "sequencing_depth_bin")
        ),
        "n_tagged_assay_partial": sum(1 for r in strata if _partial_tag(r, "assay")),
        "n_tagged_near_confound": sum(
            1 for r in strata
            if any(v is not None and v >= 0.8 for v in r["confound_cramers_v"].values())
        ),
        "n_strata_exactly_3v3": sum(
            1 for r in strata if r["n_donors_A"] == 3 and r["n_donors_B"] == 3
        ),
        "n_strata_at_least_8v8": sum(1 for v in min_donors if v >= 8),
        "n_datasets_at_least_8v8": len(
            {r["dataset_id"] for r in strata if min(r["n_donors_A"], r["n_donors_B"]) >= 8}
        ),
        "n_strata_min_donors_3": sum(1 for v in min_donors if v == 3),
        "n_strata_min_donors_4_to_7": sum(1 for v in min_donors if 4 <= v <= 7),
        "n_strata_below_8v8": sum(1 for v in min_donors if v < 8),
        "largest_design": list(
            max((min(r["n_donors_A"], r["n_donors_B"]), r["n_donors_A"], r["n_donors_B"])
                for r in strata)[1:]
        ),
        "min_bins_occupied_by_a_dataset": min(len(v) for v in bins_per_dataset.values()),
        "n_datasets_occupying_all_bins": sum(
            1 for v in bins_per_dataset.values() if len(v) == len(CELLS_PER_DONOR_BINS)
        ),
        "kpmp_strata_at_least_23v23": sum(
            1 for r in kpmp if min(r["n_donors_A"], r["n_donors_B"]) >= 23
        ),
        "rexach_assay_v_n_strata": len(rexach_v),
        "rexach_assay_v_n_non_null": len(rexach_non_null),
        "rexach_assay_v_min": min(rexach_non_null),
        "rexach_assay_v_max": max(rexach_non_null),
        "seaad_assay_v_n_strata": len(seaad),
        "seaad_assay_v_n_non_null": sum(
            1 for r in seaad if r["confound_cramers_v"]["assay"] is not None
        ),
        "n_control_strata_exactly_3v3": sum(
            1 for r in controls if r["n_donors_A"] == 3 and r["n_donors_B"] == 3
        ),
        "rexach_strata_exactly_3v3": sum(
            1 for r in rexach if r["n_donors_A"] == 3 and r["n_donors_B"] == 3
        ),
        "n_rows_failing_donor_nesting": sum(
            1 for r in rows if any("not nested" in f or "nested within" in f
                                   for f in r["gate_failures"])
        ),
        "n_rows_with_a_gate_status_blocker": sum(
            1 for r in rows if r["gate_status"] != "candidate"
        ),
        "n_candidate_rows_tagged_sequencing_depth": sum(
            1 for r in candidates if _flagged(r, "sequencing_depth_bin")
        ),
        "n_candidate_rows_tagged_assay": sum(1 for r in candidates if _flagged(r, "assay")),
        "n_candidate_rows_tagged_tissue_general": sum(
            1 for r in candidates if _flagged(r, "tissue_general")
        ),
        "n_candidate_rows_tagged_suspension_type": sum(
            1 for r in candidates if _flagged(r, "suspension_type")
        ),
        "candidate_group_median_cells_per_donor_median": candidate_median,
        "n_datasets_beyond_the_twelve_at_13v13": len(
            _tier_datasets(candidates, 13) - frozen_ids
        ),
        "seaad_dlpfc_n_donors_B_min": min(r["n_donors_B"] for r in seaad),
        "seaad_dlpfc_n_donors_B_max": max(r["n_donors_B"] for r in seaad),
        "seaad_mtg_n_donors_B_min": min(r["n_donors_B"] for r in mtg),
        "seaad_mtg_n_donors_B_max": max(r["n_donors_B"] for r in mtg),
        "yoshida_sequencing_depth_v_min": min(yoshida_depth),
        "yoshida_sequencing_depth_v_max": max(yoshida_depth),
        "smallest_frozen_dataset_rank": by_cells[0][1],
        "smallest_frozen_dataset_cell_count": by_cells[0][0],
        "second_smallest_frozen_dataset_rank": by_cells[1][1],
        "second_smallest_frozen_dataset_cell_count": by_cells[1][0],
        "n_candidate_datasets_with_drop_seq": datasets_with_assay("Drop-seq"),
        "n_candidate_datasets_with_seq_well": datasets_with_assay("Seq-Well"),
        "n_strata_in_datasets_at_least_4v4": sum(
            1 for r in strata if r["dataset_id"] in {
                row["dataset_id"] for row in strata
                if min(row["n_donors_A"], row["n_donors_B"]) >= 4
            }
        ),
        "combat_covid_permutation_count_max": max(
            r["permutation_count"] for r in combat if r["disease"] == "COVID-19"
        ),
        "combat_max_n_donors_A": max(r["n_donors_A"] for r in combat),
        "combat_max_n_donors_B": max(r["n_donors_B"] for r in combat),
        "melms_max_n_donors_A": max(r["n_donors_A"] for r in melms),
        "melms_max_n_donors_B": max(r["n_donors_B"] for r in melms),
        "wang_near_confound_strata": sum(
            1 for r in wang
            if any(v is not None and v >= 0.8 for v in r["confound_cramers_v"].values())
        ),
    }


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


def verify_attested(measured: dict) -> None:
    """Compare :data:`ATTESTED` against :func:`measure_attested` and abort on any disagreement."""
    unknown = sorted(set(ATTESTED) - set(measured))
    unmeasured = sorted(set(measured) - set(ATTESTED))
    _require(not unknown, f"ATTESTED declares figures nothing measures: {unknown}")
    _require(not unmeasured, f"measured figures nothing attests: {unmeasured}")
    wrong = {
        key: (measured[key], declared)
        for key, declared in ATTESTED.items()
        if measured[key] != declared
    }
    _require(
        not wrong,
        "attested figures disagree with the manifest "
        f"(measured, declared): {json.dumps(wrong, sort_keys=True)}",
    )


def build_header(
    manifest: dict,
    strata: list[dict],
    controls: list[dict],
    collections: dict[str, str],
    measured: dict,
) -> dict:
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
            "n_manifest_collections_at_tier": len(
                _tier_collections(candidates, threshold, collections)
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
        "external_sources": {
            "discover_index": {
                "file": DISCOVER_INDEX.name,
                "sha256": DISCOVER_SHA256,
                "bytes": DISCOVER_BYTES,
                "endpoint": DISCOVER_ENDPOINT,
                "read_date": DISCOVER_READ_DATE,
                "n_datasets": DISCOVER_N_DATASETS,
                "scope":
                    "Discover as of the read date, NOT the pinned release. Assay, suspension, "
                    "tissue, publication DOI and COLLECTION membership come from here; none of "
                    "them is a column of the candidate manifest.",
            },
            "census_release_datasets": {
                "file": RELEASE_DATASETS.name,
                "sha256": RELEASE_SHA256,
                "bytes": RELEASE_BYTES,
                "n_datasets": RELEASE_N_DATASETS,
                "scope":
                    "The pinned release's own census_info/datasets table. The authority on what "
                    "Census 2025-01-30 contains and on collection membership at the pin.",
            },
            "release_vs_discover": dict(RELEASE_VS_DISCOVER),
            "release_not_in_discover": list(RELEASE_NOT_IN_DISCOVER),
            "release_not_in_discover_collections": list(RELEASE_NOT_IN_DISCOVER_COLLECTIONS),
            "mathys_needles": list(MATHYS_NEEDLES),
            "mathys_hits": {"discover_index": 0, "census_release_datasets": 0},
            "limitation":
                "The Discover snapshot was read on 2026-08-16 and the Census is pinned to "
                f"{cs.CENSUS_VERSION}. They are not the same object: dataset_version_id matches "
                f"for {RELEASE_VS_DISCOVER['n_dataset_version_id_matching']} of "
                f"{RELEASE_N_DATASETS} release datasets, collection_doi differs for "
                f"{RELEASE_VS_DISCOVER['n_collection_doi_differing']}, and "
                f"{RELEASE_VS_DISCOVER['n_not_resolving_in_discover']} release datasets have no "
                "Discover record at all. Pinning the snapshot makes the Discover claims "
                "checkable; it does not make them claims about the release.",
            "discover_is_primary_data_discrepancy": {
                "dataset_ids": list(DISCOVER_NOT_PRIMARY),
                "manifest_value_filter": MANIFEST_VALUE_FILTER,
                "note":
                    "Both are KPMP: the frozen #8 (single-nucleus) and its within-collection "
                    "control (single-cell). Discover records is_primary_data = [False] for each, "
                    "while the manifest's Census query filtered on is_primary_data == True and "
                    "returned cells for both. Recorded as a measured disagreement between a 2026 "
                    "Discover read and the 2025-01-30 release; nothing here decides which is "
                    "right, and no stratum is dropped or kept on the strength of it.",
            },
        },
        "reasoning_document": {
            "file": PROPOSAL_MD.name,
            "sha256": PROPOSAL_SHA256,
            "bytes": PROPOSAL_BYTES,
            "circulated_sha256": PROPOSAL_ORIGINAL_SHA256,
            "circulated_bytes": PROPOSAL_ORIGINAL_BYTES,
            "redaction": PROPOSAL_REDACTION,
            "status": PROPOSAL_STATUS,
        },
        "selection_rule": SELECTION_RULE,
        "selection_rule_notes": list(SELECTION_RULE_NOTES),
        "n_frozen_strata": len(strata),
        "n_frozen_datasets": len(FROZEN_DATASETS),
        "frozen_datasets": [_dataset_block(d, strata) for d in FROZEN_DATASETS],
        "frozen_dataset_metadata_source":
            f"assay / suspension / tissue / doi are recomputed from {DISCOVER_INDEX.name}, the "
            "hash-pinned Discover snapshot, and the freeze aborts if any of them disagrees with "
            "the declaration. They are NOT columns of the source manifest and cannot be "
            "re-derived from it.",
        "frozen_dataset_collections": {
            dataset["dataset_id"]: collections[dataset["dataset_id"]]
            for dataset in FROZEN_DATASETS
        },
        "n_frozen_dataset_collections": len(
            {collections[d["dataset_id"]] for d in FROZEN_DATASETS}
        ),
        "frozen_dataset_collections_note":
            "Twelve datasets in twelve distinct collections, verified against the pinned release "
            "table. That is what makes the D2 denominator clean: no two entries of it share a "
            "cohort and a laboratory.",
        "sibling_datasets": [_sibling_block(s, controls) for s in SIBLING_DATASETS],
        "n_within_collection_control_strata": len(controls),
        "sibling_datasets_note":
            "COMPUTED, FROZEN and EXCLUDED from the D2 'independent datasets' denominator. A "
            "sibling is any candidate-bearing dataset of the manifest sharing a collection_id with "
            "one of the twelve; the set is derived from the pinned release table on every run and "
            "the declaration is compared against it, so one cannot be missed by hand. All of their "
            "candidate strata are emitted under within_collection_control_rows by the same rule as "
            "the analysis set, so no candidate stratum of any collection represented in the "
            "analysis set is left unlisted and selectable after the fact. That is the scope of the "
            "claim: it is NOT a claim about the whole Census pin, where 1197 candidate strata over "
            "68 datasets exist and only these are frozen. A result from a control never enters the "
            "D2 denominator and never counts toward a majority; promoting one to an independent "
            "dataset is an amendment.",
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
            "n_manifest_datasets_at_tier and n_manifest_collections_at_tier say what the manifest "
            "could have supported. The shortfall is a consequence of selecting for §1 (iii)'s "
            "coverage axes rather than for donor counts, NOT a property of the public data at this "
            f"Census pin: {MANIFEST_TIER_CENSUS[13]} candidate-bearing datasets in "
            f"{MANIFEST_TIER_COLLECTIONS[13]} collections clear 13v13 and "
            f"{MANIFEST_TIER_CENSUS[23]} in {MANIFEST_TIER_COLLECTIONS[23]} clear 23v23, against "
            "the frozen list's 5 and 3. Recorded before the anchor exists so that the list cannot "
            "be reassembled on donor counts later and called the original plan. "
            + COUNTERFACTUAL_REMOVED_NOTE,
        "manifest_tier_census": {str(k): v for k, v in sorted(MANIFEST_TIER_CENSUS.items())},
        "manifest_tier_collections": {
            str(k): v for k, v in sorted(MANIFEST_TIER_COLLECTIONS.items())
        },
        "manifest_n_candidate_datasets": MANIFEST_N_CANDIDATE_DATASETS,
        "manifest_n_candidate_collections": MANIFEST_N_CANDIDATE_COLLECTIONS,
        "manifest_datasets_with_exact_3v3": MANIFEST_DATASETS_WITH_EXACT_3V3,
        "manifest_collections_with_exact_3v3": MANIFEST_COLLECTIONS_WITH_EXACT_3V3,
        "attested_figures": dict(sorted(measured.items())),
        "attested_figures_note":
            "Every figure §3.2, §4.2, §5, §7, §9 and §10 of the document state in prose about the "
            "frozen set, recomputed here on every run. §10 used to attest eleven of these and "
            "nothing else in the repository held them, so nothing would have caught them drifting. "
            "median_counts_per_cell is the median of the 502 per-group medians, two per stratum — "
            "the per-stratum mean of A and B gives 3126.75 and the whole candidate set 3454.0, so "
            "the population has to be stated to make the number mean anything.",
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


def build(manifest: dict, discover: dict | None = None, release: dict | None = None) -> dict:
    discover = load_discover() if discover is None else discover
    release = load_release() if release is None else release
    collections = verify_external_sources(manifest, discover, release)
    strata = select_strata(manifest["rows"])
    controls = select_control_strata(manifest["rows"])
    verify_declarations(manifest, strata, controls, collections)
    measured = measure_attested(manifest, strata, controls, discover)
    verify_attested(measured)
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
        "header": build_header(manifest, strata, controls, collections, measured),
        "rows": rows,
        "within_collection_control_rows": control_rows,
    }


# ---------------------------------------------------------------------------
# The document. Parsed back out and compared, because three reviews found the same defect three
# times: a falsifiable claim in the prose that the artifact committed beside it refutes.
# ---------------------------------------------------------------------------

DOCUMENT_PATH = REPO / DOCUMENT
SPEC_PATH = REPO / "docs" / "PHASE0_SPEC.md"


def _cells(line: str) -> list[str]:
    """The cells of one markdown table row, with emphasis markers and padding stripped."""
    return [cell.replace("*", "").replace("`", "").strip()
            for cell in line.strip().strip("|").split("|")]


def markdown_tables(text: str) -> list[list[list[str]]]:
    """Every pipe table in ``text``, as ``[header, *rows]`` of stripped cells."""
    tables, current = [], []
    lines = text.splitlines()
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            current.append(_cells(stripped))
        else:
            if len(current) >= 2:
                tables.append(current)
            current = []
    if len(current) >= 2:
        tables.append(current)
    # Drop the |---|---| separator row of each table.
    return [
        [row for row in table if not all(set(cell) <= set("-: ") for cell in row)]
        for table in tables
    ]


def _find_table(tables, *markers: str) -> list[list[str]]:
    """The first table containing all of ``markers`` anywhere in it.

    Matching on the whole table rather than the header row, because two of the tables this has to
    find — §1's summary and §2's provenance blocks — are key/value tables with an empty header.
    The markers are chosen to be unique to one table each.
    """
    for table in tables:
        blob = " | ".join(cell for row in table for cell in row)
        if all(marker in blob for marker in markers):
            return table
    raise FrozenDeclarationMismatch(f"{DOCUMENT}: no table containing all of {markers}")


def _int(cell: str) -> int:
    match = re.search(r"-?\d[\d   ]*", cell.replace(",", ""))
    if not match:
        raise FrozenDeclarationMismatch(f"{DOCUMENT}: no integer in table cell {cell!r}")
    return int(re.sub(r"[^\d-]", "", match.group(0)))


def _float(cell: str) -> float:
    match = re.search(r"-?\d[\d   ]*(?:\.\d+)?", cell.replace(",", ""))
    if not match:
        raise FrozenDeclarationMismatch(f"{DOCUMENT}: no number in table cell {cell!r}")
    return float(re.sub(r"[^\d.-]", "", match.group(0)))


def verify_document(frozen: dict, text: str | None = None) -> None:
    """Parse the document's tables and abort if any cell disagrees with what was just derived.

    This is the guard the three previous rounds did not have. Every table below carries numbers a
    reader will quote; every one of them is now compared against the artifact rather than trusted,
    in both directions — a document that understates the freeze fails as loudly as one that
    overstates it.
    """
    text = DOCUMENT_PATH.read_text(encoding="utf-8") if text is None else text
    header = frozen["header"]
    tables = markdown_tables(text)

    # --- §1: what the act binds -------------------------------------------------------------
    summary = _find_table(tables, "Independent datasets")
    facts = {row[0]: row[1] for row in summary if len(row) == 2}
    _require(
        _int(facts["Independent datasets"]) == header["n_frozen_datasets"],
        f"§1 says {facts['Independent datasets']!r} independent datasets, the artifact has "
        f"{header['n_frozen_datasets']}",
    )
    strata_row = next(k for k in facts if k.startswith("Frozen strata"))
    _require(
        _int(facts[strata_row]) == header["n_frozen_strata"],
        f"§1 says {facts[strata_row]!r} frozen strata, the artifact has "
        f"{header['n_frozen_strata']}",
    )
    sibling_row = next(k for k in facts if "sibling" in k.lower())
    numbers = [int(n) for n in re.findall(r"\d+", facts[sibling_row])]
    _require(
        numbers == [len(header["sibling_datasets"]),
                    header["n_within_collection_control_strata"]],
        f"§1's sibling row reads {facts[sibling_row]!r}; the artifact has "
        f"{len(header['sibling_datasets'])} sibling datasets and "
        f"{header['n_within_collection_control_strata']} control strata",
    )

    # --- §3.2: the twelve --------------------------------------------------------------------
    twelve = _find_table(tables, "dataset_id", "Ceiling", "8v8")
    body = [row for row in twelve[1:] if row[0].isdigit()]
    _require(len(body) == 12, f"§3.2's table has {len(body)} numbered rows, expected 12")
    blocks = {block["dataset_id"]: block for block in header["frozen_datasets"]}
    for row in body:
        rank, dataset_id = int(row[0]), row[1]
        block = blocks.get(dataset_id)
        _require(block is not None, f"§3.2 row {rank} names {dataset_id}, which is not frozen")
        where = f"§3.2 row {rank} ({block['short']})"
        _require(block["rank"] == rank, f"{where}: the artifact ranks it {block['rank']}")
        for column, key in ((3, "n_strata"), (4, "ceiling_min_donors_per_group"),
                            (5, "n_strata_at_least_8v8"), (6, "n_strata_exactly_3v3")):
            _require(
                _int(row[column]) == block[key],
                f"{where}: column {column} reads {row[column]!r}, the artifact's {key} is "
                f"{block[key]}",
            )
        terms = [t.strip() for t in row[7].split(";")]
        declared_terms = [t["disease"] for t in block["disease_terms"]]
        _require(
            sorted(terms) == sorted(declared_terms),
            f"{where}: disease terms {terms} against the artifact's {declared_terms}",
        )
    total = next(row for row in twelve[1:] if "Total" in " ".join(row))
    _require(
        _int(total[3]) == header["n_frozen_strata"],
        f"§3.2's total row reads {total[3]!r} strata, the artifact has "
        f"{header['n_frozen_strata']}",
    )

    # --- §3.3: the siblings ------------------------------------------------------------------
    siblings = _find_table(tables, "Sibling of", "Frozen strata")
    sibling_rows = [row for row in siblings[1:] if row[0]]
    emitted = {block["dataset_id"]: block for block in header["sibling_datasets"]}
    _require(
        len(sibling_rows) == len(emitted),
        f"§3.3's table has {len(sibling_rows)} rows, the artifact freezes {len(emitted)} siblings",
    )
    for row in sibling_rows:
        dataset_id = row[0]
        block = emitted.get(dataset_id)
        _require(block is not None, f"§3.3 names {dataset_id}, which is not a frozen sibling")
        _require(
            _int(row[3]) == block["n_strata"],
            f"§3.3 ({block['short']}): {row[3]!r} frozen strata, the artifact has "
            f"{block['n_strata']}",
        )
        _require(
            _int(row[4]) == block["ceiling_min_donors_per_group"],
            f"§3.3 ({block['short']}): ceiling {row[4]!r}, the artifact has "
            f"{block['ceiling_min_donors_per_group']}",
        )
        _require(
            block["sibling_of"] in row[2],
            f"§3.3 ({block['short']}): the 'sibling of' cell {row[2]!r} does not name "
            f"{block['sibling_of']}",
        )

    # --- §4.1: the effect-size labels -------------------------------------------------------
    effects = _find_table(tables, "Expected", "Arms")
    ranks = {d["dataset_id"]: d["rank"] for d in FROZEN_DATASETS}
    by_label: dict[str, list[tuple[int, str]]] = {}
    for (dataset_id, disease), label in EXPECTED_EFFECT.items():
        by_label.setdefault(label, []).append((ranks[dataset_id], disease))
    parsed_labels = []
    for row in effects[1:]:
        label = row[0]
        if label not in EXPECTED_EFFECT_VOCABULARY:
            continue
        parsed_labels.append(label)
        expected_arms = sorted(by_label.get(label, []))
        expected_ranks = sorted({rank for rank, _ in expected_arms})
        _require(
            _int(row[1]) == len(expected_ranks),
            f"§4.1 row {label!r}: {row[1]!r} datasets, the artifact labels {len(expected_ranks)}",
        )
        _require(
            [int(n) for n in re.findall(r"#(\d+)", row[1])] == expected_ranks,
            f"§4.1 row {label!r}: dataset list {row[1]!r} against ranks {expected_ranks}",
        )
        arms = sorted(
            (int(match.group(1)), match.group(2).strip())
            for chunk in row[2].split(";")
            if (match := re.match(r"\s*#(\d+)\s+(.+)", chunk))
        )
        _require(
            arms == expected_arms,
            f"§4.1 row {label!r}: arms {arms} against the artifact's {expected_arms}",
        )
    _require(
        sorted(parsed_labels) == sorted(EXPECTED_EFFECT_VOCABULARY),
        f"§4.1's table carries labels {parsed_labels}, the vocabulary is "
        f"{list(EXPECTED_EFFECT_VOCABULARY)}",
    )

    # --- §5: bin occupancy -------------------------------------------------------------------
    occupancy_table = _find_table(tables, "Bin", "Group medians")
    occupancy = dict(header["cells_per_donor_bins"]["occupancy"])
    seen = {}
    for row in occupancy_table[1:]:
        for label, count in zip(row[0::3], row[1::3], strict=False):
            if label:
                seen[label.replace(" ", "").replace("∞", "inf")] = _int(count)
    normalised = {k.replace(" ", ""): v for k, v in occupancy.items()}
    _require(
        seen == normalised,
        f"§5's occupancy table reads {seen}, the artifact has {normalised}",
    )

    # --- §6: Layer B, and the tier census ----------------------------------------------------
    layer_b_table = _find_table(tables, "Surviving datasets", "Surviving strata")
    layer_rows = [row for row in layer_b_table[1:] if row[0].startswith("≈")]
    tiers = header["layer_b_truncation"]
    _require(
        len(layer_rows) == len(tiers),
        f"§6's Layer B table has {len(layer_rows)} tier rows, the artifact has {len(tiers)}",
    )
    for row, tier in zip(layer_rows, tiers, strict=True):
        where = f"§6 Layer B row {row[0]!r}"
        _require(_float(row[0]) == tier["sigma_donor"], f"{where}: sigma mismatch")
        _require(
            _int(row[1]) == tier["min_donors_per_group"],
            f"{where}: demands {row[1]!r}, the envelope demands {tier['min_donors_per_group']}",
        )
        _require(
            _int(row[2]) == tier["n_datasets"],
            f"{where}: {row[2]!r} surviving datasets, the artifact has {tier['n_datasets']}",
        )
        _require(
            _int(row[3]) == tier["n_strata"],
            f"{where}: {row[3]!r} surviving strata, the artifact has {tier['n_strata']}",
        )
        _require(
            (row[4].upper() == "BELOW") is tier["below_spec_dataset_floor"],
            f"{where}: verdict {row[4]!r} against below_spec_dataset_floor="
            f"{tier['below_spec_dataset_floor']}",
        )

    census_table = _find_table(tables, "Envelope demand", "collections")
    census_rows = [row for row in census_table[1:] if row[0].startswith("≥")]
    _require(
        len(census_rows) == len(tiers),
        f"§6's tier-census table has {len(census_rows)} rows, expected {len(tiers)}",
    )
    for row, tier in zip(census_rows, tiers, strict=True):
        where = f"§6 tier census row {row[0]!r}"
        _require(_int(row[0]) == tier["min_donors_per_group"], f"{where}: demand mismatch")
        _require(
            _int(row[1]) == tier["n_manifest_datasets_at_tier"],
            f"{where}: {row[1]!r} datasets, the artifact has "
            f"{tier['n_manifest_datasets_at_tier']}",
        )
        _require(
            _int(row[2]) == tier["n_manifest_collections_at_tier"],
            f"{where}: {row[2]!r} collections, the artifact has "
            f"{tier['n_manifest_collections_at_tier']}",
        )
        _require(
            _int(row[3]) == tier["n_datasets"],
            f"{where}: the frozen list retains {row[3]!r}, the artifact says {tier['n_datasets']}",
        )

    # --- the two claims the document makes ABOUT its own tables ------------------------------
    normalised = " ".join(text.split())
    cells_claim = (
        f"every cell of both tables — {len(layer_rows) * len(layer_b_table[0])} in the first, "
        f"{len(census_rows) * len(census_table[0])} in the second — is recomputed"
    )
    _require(
        cells_claim in normalised,
        f"§6 must state the size of the tables it claims to recompute: expected {cells_claim!r}",
    )
    sibling_sum = " + ".join(str(block["n_strata"]) for block in header["sibling_datasets"])
    sum_claim = f"{sibling_sum} = **{header['n_within_collection_control_strata']} strata**"
    _require(
        sum_claim in normalised,
        f"§3.3 must sum the per-sibling counts in the artifact's own order: expected {sum_claim!r}",
    )

    verify_document_cross_references(text)


def verify_document_cross_references(text: str) -> None:
    """Every § in the document resolves, and says which document it resolves in.

    The convention is stated in the document and enforced here: ``spec §N`` is PHASE0_SPEC.md and a
    bare or dotted ``§N`` is this document. It is enforced because it was violated — ``§7.1`` meant
    the spec's §7 item 1 and named nothing at all in either file, and eight bare references switched
    document relative to their neighbours.
    """
    # Headings are written both ways in this repository: "## 1. Title" and "### 3.1 Title".
    heading = r"(?m)^#{2,4}\s+(\d+(?:\.\d+)?)[.\s]"
    own_sections = {match.group(1) for match in re.finditer(heading, text)}
    spec_sections = {
        match.group(1)
        for match in re.finditer(heading, SPEC_PATH.read_text(encoding="utf-8"))
    }
    _require(bool(own_sections) and bool(spec_sections), "no numbered sections found to check")

    for match in re.finditer(r"(?:((?i:spec))\s+)?§\s*(\d+(?:\.\d+)?)", text):
        is_spec = bool(match.group(1))
        number = match.group(2)
        target, universe = (
            ("spec", spec_sections) if is_spec else (DOCUMENT, own_sections)
        )
        _require(
            number in universe,
            f"{DOCUMENT} cites {match.group(0)!r}, which is not a section of {target}. Bare and "
            "dotted § are this document; a reference to PHASE0_SPEC.md must be written 'spec §N'.",
        )


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
            f"  [manifest-wide: {tier['n_manifest_datasets_at_tier']} datasets in "
            f"{tier['n_manifest_collections_at_tier']} collections clear it]"
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
        # The document is checked here and not inside build(), so that a test exercising a
        # declaration guard is not also asserting the prose. Both halves must agree to ship.
        verify_document(frozen)
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
