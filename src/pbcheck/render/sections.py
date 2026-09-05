"""The report section model: a frozen intermediate representation of a pbcheck audit report.

``build_sections`` turns a validated ``pbcheck-audit/1`` payload into the ten report sections of
the v0.1.0 plan, in order, for every ``status`` (``complete``, ``naive_only``, ``design_only``): a
section whose arm did not run at this status prints a short "not run" note instead of its usual
blocks. ``pbcheck.render.markdown`` (and later ``pbcheck.render.html``) walk this same list, so
the two output formats never diverge in content, only in how a block is drawn.

Every fixed sentence template lives in :mod:`pbcheck.render.text` (``CAVEATS``, ``SENTENCES``,
``GLOSSARY``); this module only picks payload values and fills them in. Caveats are the one
exception: ``pbcheck.audit`` (built in parallel, elsewhere) already formats ``payload["caveats"]``
into finished text, so this module prints ``caveats[].text`` verbatim and never re-formats a
caveat template itself.
"""

from __future__ import annotations

from dataclasses import dataclass

from pbcheck import gate_config
from pbcheck.render import text


@dataclass(frozen=True)
class Paragraph:
    """A block of plain prose. Scanned by the forbidden-pattern test."""

    text: str


@dataclass(frozen=True)
class Table:
    """A block with a header row and data rows. Gene symbols and other payload values live here,
    never in a :class:`Paragraph` or :class:`Callout`, so the forbidden-pattern scan (which reads
    only prose blocks) never has to special-case a gene name that happens to contain a banned
    substring."""

    headers: tuple[str, ...]
    rows: tuple[tuple[object, ...], ...]
    caption: str | None = None


@dataclass(frozen=True)
class KeyValues:
    """A block of field/value pairs, rendered as a two-column table."""

    items: tuple[tuple[str, object], ...]
    caption: str | None = None


@dataclass(frozen=True)
class Callout:
    """A block of prose set off from the surrounding text: a caveat or the action sentence."""

    text: str
    kind: str = "note"


Block = Paragraph | Table | KeyValues | Callout


@dataclass(frozen=True)
class Section:
    id: str
    title: str
    blocks: tuple[Block, ...]


#: ``status_reason`` values in plain words, used wherever a section explains why an arm did not
#: run (plan section 1.3's status table).
STATUS_REASON_WORDS: dict[str, str] = {
    "design_only_requested": "a metadata-only run was requested",
    "donor_spans_conditions": (
        "a donor appears in both conditions, which the design gate does not allow"
    ),
    "too_few_donors": "fewer than the minimum donors were present in at least one group",
    "non_integer_counts": "no counts matrix passed the raw-count check",
    "universe_too_small": "the frozen gene universe was too small to proceed",
    "too_few_profiles_after_thin_filter": (
        "too few pseudobulk profiles remained after the thin-donor filter"
    ),
}


def _fmt(value: object) -> str:
    """Render a payload scalar for a table cell or a sentence: ``None`` as ``"n/a"``, a bool as
    ``"yes"``/``"no"``, a float to four significant figures, everything else via ``str``."""
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def _caveat_text(payload: dict, caveat_id: str) -> str | None:
    """The verbatim ``text`` of ``caveat_id`` in ``payload["caveats"]``, or ``None`` if absent."""
    for item in payload["caveats"]:
        if item["id"] == caveat_id:
            return item["text"]
    return None


def _status_reason_words(payload: dict) -> str:
    reason = payload["status_reason"]
    return STATUS_REASON_WORDS.get(reason, "no reason recorded")


def _build_sentences(payload: dict) -> list[str]:
    """The plain-language read-out lines: ``text.SENTENCES`` filled from ``payload["readout"]``
    and its neighbouring blocks. Built here, not read from ``payload["readout"]["sentences"]``,
    because the ``SENTENCES`` templates are this work package's contract, not the audit engine's."""
    lines: list[str] = []
    perm = payload["permutation_null"]
    if perm is None:
        return lines

    design = payload["design"]
    inputb = payload["input"]
    universe = payload["universe"]
    readout = payload["readout"]
    real = payload["real_label"]

    naive_perm = perm["naive"]
    floor = naive_perm["floor_solo"]
    coarse_note = ""
    if _caveat_text(payload, "C9") is not None:
        coarse_note = ", coarse because few distinct donor splits exist"
    lines.append(
        text.SENTENCES["floor_solo"].format_map(
            {
                "median_count": _fmt(floor["median_count"]),
                "median_frac_pct": _fmt(floor["median_frac"] * 100),
                "G": universe["size"],
                "mc_se": _fmt(floor["mc_se"]),
                "n_perm_achieved": readout["n_perm_naive_achieved"],
                "coarse_note": coarse_note,
                "real_solo": _fmt(readout["naive_real_solo"]),
            }
        )
    )
    lines.append(
        text.SENTENCES["lambda_naive"].format_map(
            {
                "lambda": _fmt(naive_perm["lambda"]),
                "iqr": _fmt(naive_perm["lambda_iqr"]),
                "class_word": _fmt(readout["lambda_naive_class"]),
                "band_lo": gate_config.LAMBDA_BAND[0],
                "band_hi": gate_config.LAMBDA_BAND[1],
            }
        )
    )

    pb_perm = perm["pseudobulk"]
    if pb_perm is not None:
        pb_real = real["pseudobulk"] if real is not None else None
        real_paired = _fmt(pb_real["n_significant_paired"]) if pb_real is not None else "n/a"
        lines.append(
            text.SENTENCES["pseudobulk"].format_map(
                {
                    "lambda": _fmt(pb_perm["lambda"]),
                    "class_word": _fmt(readout["lambda_pseudobulk_class"]),
                    "fp_rate": _fmt(pb_perm["fp_rate"]),
                    "se": _fmt(pb_perm["fp_rate_mc_se"]),
                    "median_count": _fmt(pb_perm["floor"]["median_count"]),
                    "real_paired": real_paired,
                }
            )
        )
    else:
        lines.append(
            text.SENTENCES["pseudobulk_not_run"].format_map(
                {"reason_text": _status_reason_words(payload)}
            )
        )

    lines.append(
        text.SENTENCES["donors"].format_map(
            {
                "n_test": _fmt(design["donors_per_group"].get(inputb["test_level"])),
                "test_level": inputb["test_level"],
                "n_ref": _fmt(design["donors_per_group"].get(inputb["ref_level"])),
                "ref_level": inputb["ref_level"],
                "n_distinct_splits": perm["n_distinct_splits"],
            }
        )
    )

    profiles = universe["profiles_per_group_after_thin_filter"]
    if profiles is not None:
        protocol = payload["settings"]["protocol_constants"]
        lines.append(
            text.SENTENCES["profiles"].format_map(
                {
                    "min_cells": protocol["min_cells"],
                    "min_counts": protocol["min_counts"],
                    "p_test": _fmt(profiles.get(inputb["test_level"])),
                    "p_ref": _fmt(profiles.get(inputb["ref_level"])),
                }
            )
        )
    return lines


def _basename(path: str | None) -> str:
    if not path:
        return "(no file)"
    return path.replace("\\", "/").rsplit("/", 1)[-1]


def _section_header(payload: dict) -> Section:
    version = payload["pbcheck_version"]
    status = payload["status"]
    reason_words = STATUS_REASON_WORDS.get(payload["status_reason"] or "", "")
    status_line = f"Status: {status}"
    if reason_words:
        status_line += f" ({reason_words})"
    blocks: list[Block] = [
        Paragraph(f"pbcheck audit - {_basename(payload['input']['path'])}"),
        Paragraph(f"pbcheck {version}, generated {payload['generated_utc']}."),
        Paragraph(status_line),
    ]
    c1 = _caveat_text(payload, "C1")
    if c1 is not None:
        blocks.append(Callout(c1))
    return Section(id="1", title="Header", blocks=tuple(blocks))


def _section_readout(payload: dict) -> Section:
    perm = payload["permutation_null"]
    readout = payload["readout"]
    blocks: list[Block] = []

    if perm is None:
        blocks.append(Paragraph(f"No detection arms were run at this status: {_status_reason_words(payload)}."))
    else:
        blocks.extend(Paragraph(line) for line in _build_sentences(payload))
        ratio_line = (
            "Real-label calls over the permutation floor: naive "
            f"{_fmt(readout['naive_real_over_floor_solo'])} times the floor"
        )
        if readout["pseudobulk_real_over_floor"] is not None:
            ratio_line += f"; pseudobulk {_fmt(readout['pseudobulk_real_over_floor'])} times the floor"
        ratio_line += "."
        blocks.append(Paragraph(ratio_line))

    a1 = _caveat_text(payload, "A1")
    if a1 is not None:
        blocks.append(Callout(a1))
    for caveat_id in ("C2", "C3", "C4", "C5", "C6", "C8", "C9"):
        caveat = _caveat_text(payload, caveat_id)
        if caveat is not None:
            blocks.append(Callout(caveat))

    return Section(id="2", title="Read-out", blocks=tuple(blocks))


def _section_glossary() -> Section:
    blocks = tuple(Paragraph(f"{term}: {definition}") for term, definition in text.GLOSSARY)
    return Section(id="3", title="What these words mean", blocks=blocks)


def _section_design(payload: dict) -> Section:
    design = payload["design"]
    inputb = payload["input"]
    universe = payload["universe"]
    blocks: list[Block] = []

    profiles = universe["profiles_per_group_after_thin_filter"] or {}
    donor_rows = tuple(
        (level, n_donors, _fmt(profiles.get(level)))
        for level, n_donors in design["donors_per_group"].items()
    )
    blocks.append(
        Table(
            headers=("group", "donors (pre-filter)", "pseudobulk profiles (post-filter)"),
            rows=donor_rows,
            caption="Donors per group",
        )
    )

    ordered = sorted(design["cells_per_donor"].items(), key=lambda kv: kv[1], reverse=True)
    top20 = ordered[:20]
    rest_count = len(ordered) - len(top20)
    blocks.append(
        Table(
            headers=("donor", "cells"),
            rows=tuple(top20),
            caption=f"Cells per donor (top {len(top20)} of {len(ordered)}; {rest_count} more not shown)",
        )
    )

    batch_cols = sorted(set(design["batch_confounded"]) | set(design["batch_separates_condition"]))
    if batch_cols:
        batch_rows = tuple(
            (
                col,
                _fmt(design["batch_confounded"].get(col)),
                _fmt(design["batch_separates_condition"].get(col)),
            )
            for col in batch_cols
        )
        blocks.append(
            Table(
                headers=("batch column", "confounding", "separates condition"),
                rows=batch_rows,
                caption="Batch confounding",
            )
        )
    else:
        blocks.append(Paragraph("No batch columns were provided."))

    blocks.append(
        KeyValues(
            items=(
                ("donor nests in condition", design["donor_nests_in_condition"]),
                ("imbalance ratio", _fmt(design["imbalance_ratio"])),
                ("usable for pseudobulk", design["usable_for_pseudobulk"]),
                ("flags", ", ".join(design["flags"]) if design["flags"] else "none"),
            ),
            caption="Design flags",
        )
    )

    blocks.append(
        KeyValues(
            items=(
                ("dropped, other condition", inputb["n_cells_dropped_other_condition"]),
                ("dropped, other cell type", inputb["n_cells_dropped_other_celltype"]),
                ("dropped, missing donor", inputb["n_cells_dropped_missing_donor"]),
                ("dropped, missing condition", inputb["n_cells_dropped_missing_condition"]),
                ("dropped, missing cell type", inputb["n_cells_dropped_missing_celltype"]),
            ),
            caption="Dropped cells",
        )
    )

    return Section(id="4", title="Design audit", blocks=tuple(blocks))


def _section_counts(payload: dict) -> Section:
    counts_check = payload["counts_check"]
    if counts_check is None:
        return Section(
            id="5",
            title="Counts check",
            blocks=(Paragraph(f"Counts check: not run ({_status_reason_words(payload)})."),),
        )
    examples = counts_check["examples"]
    examples_text = ", ".join(str(example) for example in examples) if examples else "none"
    blocks: list[Block] = [
        KeyValues(
            items=(
                ("source", _fmt(payload["input"]["counts_source"])),
                ("passed", counts_check["passed"]),
                ("reason", _fmt(counts_check["reason"])),
                ("dtype", counts_check["dtype"]),
                ("values scanned", counts_check["n_values_checked"]),
            ),
            caption="Counts check",
        ),
        Paragraph(f"Examples of the values checked: {examples_text}."),
    ]
    return Section(id="5", title="Counts check", blocks=tuple(blocks))


def _section_naive(payload: dict) -> Section:
    perm = payload["permutation_null"]
    real = payload["real_label"]
    readout = payload["readout"]
    if perm is None or real is None:
        return Section(
            id="6",
            title="Naive per-cell arm",
            blocks=(Paragraph(f"Naive per-cell arm: not run: {_status_reason_words(payload)}."),),
        )

    naive_perm = perm["naive"]
    naive_real = real["naive"]
    blocks: list[Block] = []

    count_items: list[tuple[str, object]] = [("real-label calls, solo BH", naive_real["n_significant_solo"])]
    if readout["paired_floor_shown"]:
        count_items.append(("real-label calls, paired BH", _fmt(naive_real["n_significant_paired"])))
    blocks.append(KeyValues(items=tuple(count_items), caption="Real-label counts"))

    top_rows = tuple(
        (
            item["gene"],
            _fmt(item["pval"]),
            _fmt(item["padj"]),
            _fmt(item["log2fc"]),
            _fmt(item["pct_group"]),
            _fmt(item["pct_reference"]),
        )
        for item in naive_real["top"]
    )
    blocks.append(
        Table(
            headers=("gene", "pval", "padj", "log2fc", "pct group", "pct reference"),
            rows=top_rows,
            caption="Top genes, naive arm",
        )
    )

    blocks.append(
        KeyValues(
            items=(
                ("lambda", _fmt(naive_perm["lambda"])),
                ("lambda IQR", _fmt(naive_perm["lambda_iqr"])),
                ("class", _fmt(readout["lambda_naive_class"])),
            ),
            caption="Naive inflation factor",
        )
    )

    floor_solo = naive_perm["floor_solo"]
    blocks.append(
        KeyValues(
            items=(
                ("median count", _fmt(floor_solo["median_count"])),
                ("median fraction", _fmt(floor_solo["median_frac"])),
                ("IQR count", _fmt(floor_solo["iqr_count"])),
                ("bh mode", floor_solo["bh_mode"]),
                ("Monte-Carlo SE", _fmt(floor_solo["mc_se"])),
            ),
            caption="Solo permutation floor",
        )
    )

    if readout["paired_floor_shown"] and naive_perm["floor_paired"] is not None:
        floor_paired = naive_perm["floor_paired"]
        blocks.append(
            KeyValues(
                items=(
                    ("median count", _fmt(floor_paired["median_count"])),
                    ("median fraction", _fmt(floor_paired["median_frac"])),
                    ("IQR count", _fmt(floor_paired["iqr_count"])),
                    ("bh mode", floor_paired["bh_mode"]),
                    ("Monte-Carlo SE", _fmt(floor_paired["mc_se"])),
                ),
                caption="Paired permutation floor",
            )
        )
    else:
        na_count = real["paired_bh"]["n_na_pseudobulk"]
        blocks.append(
            Paragraph(
                f"The paired floor is not shown: the pseudobulk arm left {na_count} genes without "
                "a value, so the paired series is not comparable; the solo floor above stands alone."
            )
        )

    blocks.append(
        KeyValues(
            items=(
                ("requested", readout["n_perm_naive_requested"]),
                ("achieved", _fmt(readout["n_perm_naive_achieved"])),
            ),
            caption="Permutations, naive arm",
        )
    )

    blocks.append(
        Callout(
            "B5 machinery check, not a calibration criterion: empirical-permutation-p lambda "
            f"{_fmt(naive_perm['b5_lambda_empirical'])}."
        )
    )

    return Section(id="6", title="Naive per-cell arm", blocks=tuple(blocks))


def _section_pseudobulk(payload: dict) -> Section:
    perm = payload["permutation_null"]
    real = payload["real_label"]
    readout = payload["readout"]
    pb_perm = perm["pseudobulk"] if perm is not None else None
    pb_real = real["pseudobulk"] if real is not None else None

    if pb_perm is None or pb_real is None:
        return Section(
            id="7",
            title="Donor-pseudobulk arm",
            blocks=(Paragraph(f"Donor-pseudobulk arm: not run: {_status_reason_words(payload)}."),),
        )

    blocks: list[Block] = [
        KeyValues(
            items=(("real-label calls, paired BH", pb_real["n_significant_paired"]),),
            caption="Real-label counts",
        )
    ]

    top_rows = tuple(
        (item["gene"], _fmt(item["pval"]), _fmt(item["padj"]), _fmt(item["log2fc"]))
        for item in pb_real["top"]
    )
    blocks.append(
        Table(
            headers=("gene", "pval", "padj", "log2fc"),
            rows=top_rows,
            caption="Top genes, donor-pseudobulk arm",
        )
    )

    blocks.append(
        KeyValues(
            items=(
                ("lambda", _fmt(pb_perm["lambda"])),
                ("lambda IQR", _fmt(pb_perm["lambda_iqr"])),
                ("class", _fmt(readout["lambda_pseudobulk_class"])),
                ("false-positive rate", _fmt(pb_perm["fp_rate"])),
                ("false-positive rate Monte-Carlo SE", _fmt(pb_perm["fp_rate_mc_se"])),
            ),
            caption="Pseudobulk inflation and false-positive rate",
        )
    )

    floor = pb_perm["floor"]
    blocks.append(
        KeyValues(
            items=(
                ("median count", _fmt(floor["median_count"])),
                ("median fraction", _fmt(floor["median_frac"])),
                ("bh mode", floor["bh_mode"]),
                ("Monte-Carlo SE", _fmt(floor["mc_se"])),
            ),
            caption="Pseudobulk permutation floor",
        )
    )

    blocks.append(
        KeyValues(
            items=(
                ("requested", readout["n_perm_pb_requested"]),
                ("achieved", _fmt(readout["n_perm_pb_achieved"])),
            ),
            caption="Permutations, pseudobulk arm",
        )
    )

    blocks.append(
        Callout(
            "B5 machinery check, not a calibration criterion: empirical-permutation-p lambda "
            f"{_fmt(pb_perm['b5_lambda_empirical'])}."
        )
    )

    moderation = pb_real["moderation"]
    tech_rows = (
        ("d0", _fmt(moderation.get("d0")), "moderated eBayes prior degrees of freedom"),
        (
            "shrinkage factor",
            _fmt(moderation.get("shrinkage_factor")),
            "how strongly a gene's own variance is pulled toward the prior",
        ),
        (
            "complete pooling",
            _fmt(moderation.get("complete_pooling")),
            "whether every gene's variance was replaced by the prior outright",
        ),
        (
            "residual df",
            _fmt(moderation.get("residual_df")),
            "residual degrees of freedom feeding the moderated test",
        ),
    )
    blocks.append(
        Table(
            headers=("quantity", "value", "meaning"),
            rows=tech_rows,
            caption="Moderated eBayes technical detail",
        )
    )

    return Section(id="7", title="Donor-pseudobulk arm", blocks=tuple(blocks))


def _section_universe(payload: dict) -> Section:
    universe = payload["universe"]
    real = payload["real_label"]
    perm = payload["permutation_null"]
    blocks: list[Block] = [
        KeyValues(
            items=(
                ("size", universe["size"]),
                ("minimum size", universe["min_size"]),
                ("builder", _fmt(universe["builder"])),
            ),
            caption="Frozen gene universe",
        ),
        Paragraph(f"Builder rule: {universe['builder_rule']}"),
    ]
    thin = universe["thin_donor_filter"]
    if thin:
        blocks.append(
            KeyValues(items=tuple((str(k), _fmt(v)) for k, v in thin.items()), caption="Thin-donor filter")
        )
    else:
        blocks.append(Paragraph("Thin-donor filter: not run."))

    if real is not None:
        paired_bh = real["paired_bh"]
        blocks.append(
            KeyValues(
                items=(
                    ("common tested", paired_bh["n_tested_common"]),
                    ("dropped for fairness", paired_bh["n_dropped_for_fairness"]),
                    ("no value in naive arm", paired_bh["n_na_naive"]),
                    ("no value in pseudobulk arm", paired_bh["n_na_pseudobulk"]),
                    ("pseudobulk arm complete for every gene", paired_bh["pseudobulk_na_free"]),
                ),
                caption="Shared universe, real-label bookkeeping",
            )
        )

    if perm is not None:
        blocks.append(
            KeyValues(
                items=(
                    ("real split inside permutation range", _fmt(perm["real_split_inside_perm_range"])),
                    ("real split percentile in permutations", _fmt(perm["real_split_percentile_in_perms"])),
                ),
                caption="Real split against the permutation range",
            )
        )
        if perm["monte_carlo"]:
            blocks.append(
                KeyValues(
                    items=tuple((str(k), _fmt(v)) for k, v in perm["monte_carlo"].items()),
                    caption="Monte-Carlo detail",
                )
            )

    return Section(id="8", title="Shared universe bookkeeping", blocks=tuple(blocks))


def _section_settings(payload: dict) -> Section:
    tool = payload["settings"]["tool"]
    protocol = payload["settings"]["protocol_constants"]
    provenance = payload["provenance"]
    blocks: list[Block] = [
        Table(
            headers=("setting", "value"),
            rows=tuple((key, _fmt(value)) for key, value in tool.items()),
            caption="Tool settings",
        ),
        Table(
            headers=("protocol constant", "value"),
            rows=tuple((key, _fmt(value)) for key, value in protocol.items()),
            caption="Protocol constants",
        ),
        Table(
            headers=("sigma_donor", "min donors per group", "grid support"),
            rows=tuple(
                (_fmt(row["sigma_donor"]), row["min_donors_per_group"], row["grid_support"])
                for row in provenance["operating_envelope"]
            ),
            caption="Operating envelope",
        ),
        Table(
            headers=("package", "version"),
            rows=tuple((key, value) for key, value in provenance["packages"].items()),
            caption="Package versions",
        ),
        KeyValues(
            items=(("platform", provenance["platform"]), ("python", provenance["python"])),
            caption="Platform",
        ),
        Table(
            headers=("stage", "seconds"),
            rows=tuple((key, _fmt(value)) for key, value in payload["runtime_by_stage_seconds"].items()),
            caption="Runtime by stage",
        ),
    ]
    c7 = _caveat_text(payload, "C7")
    if c7 is not None:
        blocks.append(Callout(c7))
    return Section(id="9", title="Settings and provenance", blocks=tuple(blocks))


def _section_footer(payload: dict) -> Section:
    version = payload["pbcheck_version"]
    line = (
        f"Generated by pbcheck {version}. Protocol: docs/PHASE0_SPEC.md and docs/AMENDMENTS.md "
        "in the pbcheck repository."
    )
    return Section(id="10", title="Footer", blocks=(Paragraph(line),))


def build_sections(payload: dict) -> list[Section]:
    """The ten report sections of the plan, in order, filled from ``payload``.

    Works at every ``status``: a section whose data did not run at this status (the DE arms in
    ``design_only``, the pseudobulk arm in ``naive_only``) prints a short "not run" note in place
    of its usual blocks, rather than being omitted.
    """
    return [
        _section_header(payload),
        _section_readout(payload),
        _section_glossary(),
        _section_design(payload),
        _section_counts(payload),
        _section_naive(payload),
        _section_pseudobulk(payload),
        _section_universe(payload),
        _section_settings(payload),
        _section_footer(payload),
    ]


def prose_blocks(sections: list[Section]) -> list[str]:
    """The text of every :class:`Paragraph` and :class:`Callout` block, in section order.

    What the forbidden-pattern test scans: table cells (gene symbols, raw values) are excluded on
    purpose, so a gene name like ``GOLGA8A`` never trips a prose rule meant for sentences.
    """
    texts: list[str] = []
    for section in sections:
        for block in section.blocks:
            if isinstance(block, (Paragraph, Callout)):
                texts.append(block.text)
    return texts


def summary_lines(payload: dict) -> list[str]:
    """The CLI's end-of-run summary: the header and read-out paragraphs, reused verbatim from
    :func:`build_sections` rather than composed again from the payload."""
    lines: list[str] = []
    for section in build_sections(payload):
        if section.id in ("1", "2"):
            lines.extend(block.text for block in section.blocks if isinstance(block, Paragraph))
    return lines
