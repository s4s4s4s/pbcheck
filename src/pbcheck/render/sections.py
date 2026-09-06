"""The report section model: a frozen intermediate representation of a pbcheck audit report.

``build_sections`` turns a validated ``pbcheck-audit/1`` payload into the ten report sections of
the v0.1.0 plan, in order, for every ``status`` (``complete``, ``naive_only``, ``design_only``): a
section whose arm did not run at this status prints a short "not run" note instead of its usual
blocks. ``pbcheck.render.markdown`` and ``pbcheck.render.html`` walk this same list, so the two
output formats never diverge in content, only in how a block is drawn.

This module composes no sentence of its own. Every fixed line it prints comes from
:mod:`pbcheck.render.text` (``CAVEATS``, ``SENTENCES``, ``GLOSSARY``, ``REPORT_LINES``,
``STATUS_REASON_WORDS``); the only strings written here are the captions that label a table. Two
kinds of finished prose are printed verbatim rather than filled in here: ``payload["caveats"][].text``
and ``payload["readout"]["sentences"]``, both assembled by :mod:`pbcheck.audit` from those same
templates, so the sentences a reader gets are the sentences the JSON payload carries.

Numbers are drawn by :func:`_fmt`: four significant figures for a float, ``"n/a"`` for a missing
value, ``"yes"``/``"no"`` for a flag. That convention is the report's own (the protocol states
none) and is shared with ``pbcheck.audit.format_scalar``, which formats the payload's sentences.
"""

from __future__ import annotations

from dataclasses import dataclass

from pbcheck.render import text
from pbcheck.render.text import NO_REASON_RECORDED, NOT_REACHED, STATUS_REASON_WORDS

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


def _fmt(value: object) -> str:
    """Render a payload scalar for a table cell: ``None`` as ``"n/a"``, a bool as ``"yes"``/``"no"``,
    a float to four significant figures, everything else via ``str``.

    The four-figure convention is the report's, not the protocol's: it is applied here and in
    ``pbcheck.audit.format_scalar`` (which formats the read-out sentences the payload carries), and
    ``tests/test_audit.py`` pins the two against each other so one number never appears in two
    shapes in one report.
    """
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def _fmt_achieved(value: object) -> str:
    """An achieved permutation count: ``None`` means the arm never reached the null, which is a
    fact about the run rather than a missing cell, so it is said in words."""
    return NOT_REACHED if value is None else _fmt(value)


def _caveat_text(payload: dict, caveat_id: str) -> str | None:
    """The verbatim ``text`` of ``caveat_id`` in ``payload["caveats"]``, or ``None`` if absent."""
    for item in payload["caveats"]:
        if item["id"] == caveat_id:
            return item["text"]
    return None


def _status_reason_words(payload: dict) -> str:
    reason = payload["status_reason"]
    return STATUS_REASON_WORDS.get(reason, NO_REASON_RECORDED)


def _basename(path: str | None) -> str:
    if not path:
        return "(no file)"
    return path.replace("\\", "/").rsplit("/", 1)[-1]


def _section_header(payload: dict) -> Section:
    reason = payload["status_reason"]
    reason_words = STATUS_REASON_WORDS.get(reason) if reason else None
    status_line = (
        text.report_line("header_status", status=payload["status"])
        if reason_words is None
        else text.report_line(
            "header_status_with_reason", status=payload["status"], reason_text=reason_words
        )
    )
    blocks: list[Block] = [
        Paragraph(text.report_line("header_title", file=_basename(payload["input"]["path"]))),
        Paragraph(text.report_line(
            "header_version",
            version=payload["pbcheck_version"],
            generated_utc=payload["generated_utc"],
        )),
        Paragraph(status_line),
    ]
    n1 = _caveat_text(payload, "N1")
    if n1 is not None:
        blocks.append(Callout(n1))
    return Section(id="1", title="Header", blocks=tuple(blocks))


def _ratio_line(payload: dict) -> str:
    """The read-out ratio line: real-label calls over the permutation floor, per arm.

    Each ratio names the BH convention of the two counts it divides, and the pseudobulk ratio is
    printed only when the payload shows the paired floor, so the line never sets a solo-corrected
    count against a paired-corrected one. Below the donor threshold the ratios are the quantity
    change 2 of the fifth amendment (docs/AMENDMENTS.md) calls leak-contaminated, and the line
    says so instead of leaving the number to be read as a signal-to-floor factor.
    """
    readout = payload["readout"]
    line = text.report_line(
        "ratio_solo", naive_ratio=_fmt(readout["naive_real_over_floor_solo"])
    )
    if readout["paired_floor_shown"] and readout["pseudobulk_real_over_floor"] is not None:
        line += text.report_line(
            "ratio_paired", pseudobulk_ratio=_fmt(readout["pseudobulk_real_over_floor"])
        )
    if readout["few_donors"]:
        line += text.report_line("ratio_leak_contaminated")
    return line


def _section_readout(payload: dict) -> Section:
    """Section 2: the payload's own read-out sentences, printed verbatim, then the notes.

    ``payload["readout"]["sentences"]`` is the single source of these lines (``pbcheck.audit``
    fills it from :data:`pbcheck.render.text.SENTENCES`): the renderer does not re-derive them, so
    the report and the JSON a reader gets cannot disagree.
    """
    blocks: list[Block] = []

    if payload["permutation_null"] is None:
        blocks.append(Paragraph(
            text.report_line("no_arms_run", reason_text=_status_reason_words(payload))
        ))
    else:
        blocks.extend(Paragraph(line) for line in payload["readout"]["sentences"])
        blocks.append(Paragraph(_ratio_line(payload)))

    for caveat_id in ("R1", "N2", "N3", "N4", "N5", "N6", "N8", "N9"):
        caveat = _caveat_text(payload, caveat_id)
        if caveat is not None:
            blocks.append(Callout(caveat))

    return Section(id="2", title="Read-out", blocks=tuple(blocks))


def _section_glossary() -> Section:
    blocks = tuple(
        Paragraph(text.report_line("glossary_entry", term=term, definition=definition))
        for term, definition in text.GLOSSARY
    )
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
        blocks.append(Paragraph(text.report_line("no_batch_columns")))

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
            blocks=(Paragraph(text.report_line(
                "counts_check_not_run", reason_text=_status_reason_words(payload)
            )),),
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
        Paragraph(text.report_line("counts_examples", examples=examples_text)),
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
            blocks=(Paragraph(text.report_line(
                "arm_not_run", arm="Naive per-cell arm",
                reason_text=_status_reason_words(payload),
            )),),
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

    floor_paired = naive_perm["floor_paired"]
    if not readout["paired_floor_shown"]:
        blocks.append(Paragraph(text.report_line(
            "paired_floor_not_comparable",
            n_na_pseudobulk=real["paired_bh"]["n_na_pseudobulk"],
        )))
    elif floor_paired is None:
        blocks.append(Paragraph(text.report_line("paired_floor_not_measured")))
    else:
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
    blocks.append(
        KeyValues(
            items=(
                ("requested", readout["n_perm_naive_requested"]),
                ("achieved", _fmt_achieved(readout["n_perm_naive_achieved"])),
            ),
            caption="Permutations, naive arm",
        )
    )

    blocks.append(
        Callout(text.report_line(
            "machinery_check",
            b5_lambda_empirical=_fmt(naive_perm["b5_lambda_empirical"]),
        ))
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
            blocks=(Paragraph(text.report_line(
                "arm_not_run", arm="Donor-pseudobulk arm",
                reason_text=_status_reason_words(payload),
            )),),
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
                ("achieved", _fmt_achieved(readout["n_perm_pb_achieved"])),
            ),
            caption="Permutations, pseudobulk arm",
        )
    )

    blocks.append(
        Callout(text.report_line(
            "machinery_check",
            b5_lambda_empirical=_fmt(pb_perm["b5_lambda_empirical"]),
        ))
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
        Paragraph(text.report_line("builder_rule", rule=universe["builder_rule"])),
    ]
    thin = universe["thin_donor_filter"]
    if thin:
        blocks.append(
            KeyValues(items=tuple((str(k), _fmt(v)) for k, v in thin.items()), caption="Thin-donor filter")
        )
    else:
        blocks.append(Paragraph(text.report_line("thin_filter_not_run")))

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
    n7 = _caveat_text(payload, "N7")
    if n7 is not None:
        blocks.append(Callout(n7))
    return Section(id="9", title="Settings and provenance", blocks=tuple(blocks))


def _section_footer(payload: dict) -> Section:
    line = text.report_line("footer", version=payload["pbcheck_version"])
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
    """Every string of a section that is pbcheck's own prose rather than a payload value.

    What the forbidden-pattern gate scans: the text of every :class:`Paragraph` and
    :class:`Callout`, and the caption of every :class:`Table` and :class:`KeyValues` (a caption is
    written by this module and reaches both outputs, so it is prose too). Table cells are excluded
    on purpose: they hold payload values, and a gene symbol like ``GOLGA8A`` must not trip a rule
    meant for sentences.
    """
    texts: list[str] = []
    for section in sections:
        for block in section.blocks:
            if isinstance(block, (Paragraph, Callout)):
                texts.append(block.text)
            elif isinstance(block, (Table, KeyValues)) and block.caption:
                texts.append(block.caption)
    return texts


def summary_lines(payload: dict) -> list[str]:
    """The CLI's end-of-run summary: the header and read-out sections, reused verbatim from
    :func:`build_sections` rather than composed again from the payload.

    Callouts are included, so the summary a user sees in the terminal carries the always-on scope
    note (N1) and the read-out paragraph (R1) rather than the numbers alone.
    """
    lines: list[str] = []
    for section in build_sections(payload):
        if section.id in ("1", "2"):
            lines.extend(
                block.text
                for block in section.blocks
                if isinstance(block, (Paragraph, Callout))
            )
    return lines
