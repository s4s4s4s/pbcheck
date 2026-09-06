"""The fixed prose of the pbcheck report: caveats, read-out sentences, glossary, forbidden patterns.

These texts are a protocol surface (the wording-rules section of the v0.1.0 implementation plan):
they are checked against the forbidden patterns listed below and are changed only together with
their tests. This module owns every fixed sentence the report can print, in six blocks:

* :data:`CAVEATS` (the notes N1..N9 and the read-out paragraph R1) with :func:`caveat_text` and
  :func:`readout_caveat_text`, which assemble R1's conditional clauses;
* :data:`SENTENCES`, the read-out line templates the audit fills and the renderer prints verbatim;
* :data:`GLOSSARY`, the plain-language section;
* :data:`REPORT_LINES` with :func:`report_line`, the report's remaining fixed lines (header,
  "not run" notes, the read-out ratio line, the paired-floor explanations, the footer), and
  :data:`STATUS_REASON_WORDS`, the plain words for each ``status_reason``;
* :data:`FORBIDDEN_PATTERNS` with :func:`compiled_forbidden_patterns`, :func:`prose_for_pattern_check`
  and :func:`forbidden_pattern_hits`, the single place those patterns are compiled and applied.

Nothing here is improvised at render time: :mod:`pbcheck.render.sections` composes sections out of
these constants, and :mod:`pbcheck.audit_schema` re-applies the pattern check to every caveat text
in a payload, so a report cannot reach a reader past this gate.

The module imports ``pbcheck.gate_config`` and ``pbcheck.product_constants`` and nothing else of
the package: every protocol number in the prose comes from the first, the donor threshold from the
second, and no number is typed by hand here.
"""

from __future__ import annotations

import re
from types import MappingProxyType
from typing import Any

from pbcheck import gate_config
from pbcheck.product_constants import FEW_DONORS_THRESHOLD

#: The committed synthetic-gate run the provenance note (N7) cites as the record of the one
#: measurement the engine has: an artifact path, not a verdict.
GATE_ARTIFACT = "pilot/gate/synthetic_gate_2026-08-15.json"

#: Longest run of characters a single-quoted span may hold and still be masked by
#: :func:`prose_for_pattern_check`. Values longer than :data:`MAX_QUOTED_VALUE_CHARS` are shortened
#: by :func:`quoted` so a user-controlled name can never grow past the mask.
QUOTED_SPAN_LIMIT = 120
MAX_QUOTED_VALUE_CHARS = 100

#: Placeholders whose value comes from the user's file (column names, layer names, the raw-count
#: check's reason string) rather than from the payload's own numbers. :func:`caveat_text` puts
#: every one of them through :func:`quoted`, so the wording gate treats them as identifiers rather
#: than as prose and a column named after a forbidden word cannot forge a sentence.
USER_VALUE_PLACEHOLDERS = frozenset(
    {"col", "cols", "reason", "test_level", "ref_level", "file"}
)


def quoted(value: object) -> str:
    """``value`` as a single-quoted identifier that :func:`prose_for_pattern_check` can mask.

    Three things are normalised so the quoting cannot be escaped by the value itself: single
    quotes become double quotes (an embedded quote would close the span early), newlines and runs
    of whitespace become one space (the mask does not cross a newline), and a value longer than
    :data:`MAX_QUOTED_VALUE_CHARS` is shortened with a trailing ellipsis (the mask does not span
    more than :data:`QUOTED_SPAN_LIMIT` characters). The result is a display form; the payload's
    own fields keep the value as it was read.
    """
    text = re.sub(r"\s+", " ", str(value)).replace("'", '"').strip()
    if len(text) > MAX_QUOTED_VALUE_CHARS:
        text = text[: MAX_QUOTED_VALUE_CHARS - 3] + "..."
    return f"'{text}'"


#: A quoted identifier: an opening quote no letter precedes, a closing quote no letter follows,
#: and no newline or second quote in between. The letter guards are what keeps an English
#: possessive (the arm's band) from opening a span and masking the sentence that follows it, which
#: would turn the mask itself into the loophole it exists to close.
_QUOTED_SPAN = re.compile(rf"(?<![A-Za-z])'[^'\n]{{1,{QUOTED_SPAN_LIMIT}}}'(?![A-Za-z])")


def prose_for_pattern_check(prose: str) -> str:
    """``prose`` with every quoted identifier emptied, for the forbidden-pattern check.

    The wording rules bind pbcheck's own sentences. A user's column name is quoted where it is
    interpolated (see :data:`USER_VALUE_PLACEHOLDERS`), and a quoted identifier is not the tool
    making a claim, so it is masked here before the patterns are applied. Possessive apostrophes
    are not delimiters (see :data:`_QUOTED_SPAN`); nothing else is changed.
    """
    return _QUOTED_SPAN.sub("''", prose)


#: Prose patterns the rendered report must never contain (the plan's wording rules). There is no
#: whitelist: a sentence that would need one is reworded instead. Compiled in exactly one place,
#: :func:`compiled_forbidden_patterns`, case-insensitively, and applied to prose that has been put
#: through :func:`prose_for_pattern_check` first.
FORBIDDEN_PATTERNS: tuple[str, ...] = (
    r"validit|validat",
    r"\bvalid\b",
    r"\bfindings?\b",
    r"false discover",
    r"\bpublished\b",
    r"\bproves?\b",
    r"\bGO\b",
    r"no[- ]?go",
    r"risk[ _-]score",
    r"INSTRUMENT VALID",
    r"\bTier\b",
    r"\breportab",
    r"\bconcordance\b",
    r"\bJaccard\b",
    r"\boverlap\b",
)

_COMPILED_FORBIDDEN_PATTERNS: tuple[re.Pattern[str], ...] = tuple(
    re.compile(pattern, re.IGNORECASE) for pattern in FORBIDDEN_PATTERNS
)


def compiled_forbidden_patterns() -> tuple[re.Pattern[str], ...]:
    """The forbidden patterns, compiled once, case-insensitively. The only compilation site."""
    return _COMPILED_FORBIDDEN_PATTERNS


def forbidden_pattern_hits(prose: str) -> tuple[str, ...]:
    """The patterns of :data:`FORBIDDEN_PATTERNS` that ``prose`` matches, masked spans excluded.

    Empty when the prose passes. Callers (the schema validator, the renderers, the prose tests)
    use this rather than compiling the list themselves, so all of them apply the same mask and the
    same flags.
    """
    masked = prose_for_pattern_check(prose)
    return tuple(
        pattern.pattern
        for pattern in _COMPILED_FORBIDDEN_PATTERNS
        if pattern.search(masked) is not None
    )


#: The prose phrase for each ``readout.lambda_*_class`` value of the payload schema. The class
#: names are neutral by decision (they describe a position against a band, they do not award the
#: protocol's reserved words), and the band itself is the donor-pseudobulk arm's; the naive arm's
#: sentence in :data:`SENTENCES` says so where it uses one of these phrases.
LAMBDA_CLASS_WORDS = MappingProxyType({
    "in_band": "inside the band",
    "above_band": "above the band",
    "below_band": "below the band",
})

#: One operating-envelope row, filled from ``gate_config.OPERATING_ENVELOPE`` and the pre-registered
#: oracle by :func:`envelope_rows`. The row's own ``grid_support`` is quoted because it is a
#: verbatim quotation of the frozen constant, not a sentence of this report.
ENVELOPE_ROW_TEMPLATE = (
    "sigma_donor {sigma_donor}: at least {min_donors_per_group} donors per group (power at least "
    "{power_target} at log2FC {oracle_log2fc} in {oracle_k} genes; grid support {grid_support})"
)


def envelope_rows() -> str:
    """The operating envelope, one row per line, built from ``pbcheck.gate_config``.

    No envelope number is typed here: the donor counts, the sigmas and each row's grid support come
    from ``gate_config.OPERATING_ENVELOPE``, the power target and the oracle from the pre-registered
    constants next to it. The rows say what the envelope means (a minimum donor count per group at
    a donor-effect size, for a stated power at the pre-registered oracle) and what the committed
    grid says about each point, so a row that is derived or extrapolated says so in its own words.
    """
    return "\n".join(
        ENVELOPE_ROW_TEMPLATE.format(
            sigma_donor=row["sigma_donor"],
            min_donors_per_group=row["min_donors_per_group"],
            power_target=gate_config.POWER_TARGET,
            oracle_log2fc=gate_config.ORACLE_LOG2FC,
            oracle_k=gate_config.ORACLE_K,
            grid_support=quoted(row["grid_support"]),
        )
        for row in gate_config.OPERATING_ENVELOPE
    )


#: Values every caveat template may interpolate without the caller supplying them: the protocol
#: constants of the one regime the engine was measured at, the donor threshold, and the artifact
#: that records the measurement. Keeping them here is what makes the "no number typed by hand"
#: rule checkable over the whole module.
_TEMPLATE_CONSTANTS: dict[str, Any] = {
    "calibration_sigma": gate_config.CALIBRATION_EVAL_SIGMA,
    "oracle_donors": gate_config.ORACLE_SIM["n_donors_per_group"],
    "oracle_genes": gate_config.ORACLE_SIM["n_genes"],
    "gate_artifact": GATE_ARTIFACT,
    "threshold": FEW_DONORS_THRESHOLD,
}

#: Caveat and read-out templates, formatted by :func:`caveat_text` and :func:`readout_caveat_text`.
#: N1, N2, N3 and N7 are rendered in every report; N4, N5, N6, N8, N9 fire under their own
#: condition; R1 is the read-out paragraph of a complete run. No template holds a single quote:
#: user-controlled values are quoted by :func:`quoted` when they are interpolated.
CAVEATS = MappingProxyType({
    "N1": (
        "This is a pbcheck v{version} audit of one stratum of one file. It measures how a naive "
        "per-cell test and a donor-pseudobulk test behave on this data under a donor-permutation "
        "null. It is run outside pbcheck's pre-registered Phase 0 protocol, is not a Phase 0 "
        "measurement, and makes no claim about any publication or about whether any reported "
        "result is true."
    ),
    "N2": (
        "Amendment 3 declares an operating envelope for the pseudobulk arm on synthetic oracles. "
        "For each donor-variance point it states the donor count per group at which the power "
        "target is reached, and each row says what the committed grid supports there: a point "
        "measured on the grid, or a count derived or extrapolated from it.\n{envelope_rows}\n"
        "The arm's calibration was evaluated at one hard regime (sigma_donor {calibration_sigma}, "
        "{oracle_donors} against {oracle_donors} donors) and nowhere else. pbcheck does not "
        "estimate sigma_donor for real data, so whether this stratum lies inside that envelope is "
        "not determined here."
    ),
    "N3": (
        "Everything in this report is a diagnostic of this file: the naive arm's inflation factor "
        "and permutation floor, the pseudobulk arm's inflation factor, floor and false-positive "
        "rate as its negative control, and the bookkeeping of the shared gene universe. None of "
        "it is a Phase 0 result, and pbcheck does not rate this file against the Phase 0 decision "
        "rule."
    ),
    "N4": (
        "At least one group has fewer than {threshold} donors: the permutation null has few "
        "distinct donor splits and the floors are coarse. Amendment 5 Change 2 of pbcheck's "
        "protocol admits floor-based quantities outside the envelope only when every group has at "
        "least {threshold} donors, and this run is below that in at least one group, so its "
        "floor-based numbers describe this run alone. Do not compare them with a run on another "
        "file."
    ),
    "N5": (
        "The design audit found a batch column that separates the two conditions ({cols}): "
        "condition and batch cannot be told apart in this stratum, and neither arm's real-label "
        "result can be read as a condition effect."
    ),
    "N6": (
        "The count matrix failed the raw-count check ({reason}). The pseudobulk arm was dropped, "
        "never rounded, so this run has no negative control. The naive arm was run on the matrix "
        "as found, with the naive pipeline's own normalisation and log transform applied on top "
        "of it; its numbers describe that pipeline on this matrix and are not comparable to a run "
        "on raw counts."
    ),
    "N7": (
        "The engine was measured on one synthetic oracle point (sigma_donor {calibration_sigma}, "
        "{oracle_donors} against {oracle_donors} donors, {oracle_genes} genes), recorded in "
        "{gate_artifact}; that measurement is not repeated on this file. Of the settings below, "
        "only the ones listed as protocol constants are pre-registered values taken from "
        "pbcheck.gate_config ({protocol_constant_names}); the permutation counts, the universe "
        "filter parameters, the fallback universe rule and the display settings are the tool's "
        "own and are not protocol values."
    ),
    "N8": (
        "No cell type was selected, so all cells were pooled into one stratum; the file has a "
        "column that looks like a cell-type annotation ({col}). Pooling mixes composition shifts "
        "between conditions into the contrast; rerun with --celltype {col} --celltype-value "
        "<level> for a per-cell-type audit."
    ),
    "N9": (
        "Only {n} distinct donor splits exist for this design, so the permutation null has {n} "
        "draws, not the {requested} requested; the floor's Monte-Carlo standard error ({se}) is "
        "correspondingly large."
    ),
    "R1": (
        "On this file, with condition labels shuffled between donors and therefore no real signal "
        "to find, the per-cell test still calls a median of {floor_solo} of {universe_size} genes "
        "at FDR {alpha} ({floor_pct}%), corrected over the whole gene universe on its own (solo "
        "BH); on the real labels, corrected the same way, it calls {real_solo}. "
        "{separation_clause} The donor is the replication unit this design supports."
        "{pseudobulk_clause}"
    ),
})

#: R1's conditional clauses. ``separable`` and ``leak_contaminated`` are alternatives, chosen by
#: whether every group has at least :data:`~pbcheck.product_constants.FEW_DONORS_THRESHOLD` donors:
#: below that the protocol calls the real-over-floor ratio leak-contaminated, so the categorical
#: sentence is replaced by a statement about this run's numbers. ``pseudobulk`` is rendered only
#: when the paired floor is shown, and names the BH convention of the two numbers it carries.
R1_CLAUSES = MappingProxyType({
    "separable": (
        "A gene list produced by a per-cell test on this data cannot be separated from that floor."
    ),
    "leak_contaminated": (
        "With fewer than {threshold} donors in a group, the ratio of the real-label count to that "
        "floor is contaminated by the per-cell leak at this donor count and is not interpretable "
        "on this run."
    ),
    "pseudobulk": (
        " The donor-pseudobulk test called {pb_real} genes on the real labels against a "
        "permutation median of {pb_floor}, both corrected across the two arms together over the "
        "genes they have in common (paired BH)."
    ),
})


def _format(template: str, values: dict[str, Any]) -> str:
    merged = dict(_TEMPLATE_CONSTANTS)
    for key, value in values.items():
        merged[key] = quoted(value) if key in USER_VALUE_PLACEHOLDERS else value
    return template.format_map(merged)


def caveat_text(caveat_id: str, **values: object) -> str:
    """Format ``CAVEATS[caveat_id]`` with ``values`` and the module's own constants.

    Placeholders that name a protocol constant, the donor threshold or the gate artifact are
    filled from :data:`_TEMPLATE_CONSTANTS`; ``{envelope_rows}`` is filled by :func:`envelope_rows`
    when the caller does not pass it; a placeholder listed in :data:`USER_VALUE_PLACEHOLDERS` is
    put through :func:`quoted`. Raises ``KeyError`` naming the missing placeholder when a template
    references a name none of those sources supplies.
    """
    if caveat_id == "R1":
        raise ValueError("R1 is assembled by readout_caveat_text, not by caveat_text")
    if "{envelope_rows}" in CAVEATS[caveat_id] and "envelope_rows" not in values:
        values = {**values, "envelope_rows": envelope_rows()}
    return _format(CAVEATS[caveat_id], values)


def envelope_sentence() -> str:
    """The N2 caveat with its envelope rows filled in."""
    return caveat_text("N2")


def readout_caveat_text(
    *,
    floor_solo: object,
    universe_size: object,
    alpha: object,
    floor_pct: object,
    real_solo: object,
    few_donors: bool,
    paired_floor_shown: bool,
    pb_real: object = None,
    pb_floor: object = None,
) -> str:
    """The read-out paragraph R1, with the clauses this run is entitled to.

    ``few_donors`` is true when at least one group has fewer than
    :data:`~pbcheck.product_constants.FEW_DONORS_THRESHOLD` donors: the categorical clause is then
    replaced by the leak-contaminated one. ``paired_floor_shown`` is the payload's own flag (true
    only when the paired correction left the pseudobulk arm free of missing values), and the
    pseudobulk clause is rendered only then, with ``pb_real`` and ``pb_floor`` required.
    """
    if paired_floor_shown and (pb_real is None or pb_floor is None):
        raise ValueError("paired_floor_shown requires both pb_real and pb_floor")
    separation = R1_CLAUSES["leak_contaminated" if few_donors else "separable"]
    pseudobulk = R1_CLAUSES["pseudobulk"] if paired_floor_shown else ""
    return _format(
        CAVEATS["R1"],
        {
            "floor_solo": floor_solo,
            "universe_size": universe_size,
            "alpha": alpha,
            "floor_pct": floor_pct,
            "real_solo": real_solo,
            "separation_clause": _format(separation, {}),
            "pseudobulk_clause": _format(pseudobulk, {"pb_real": pb_real, "pb_floor": pb_floor}),
        },
    )


#: What a table cell says instead of a number when an arm never reached the permutation null: a
#: fact about the run, not a missing value, so it is said in words.
NOT_REACHED = "not reached"

#: What the report says when a payload carries a ``status_reason`` this module has no words for.
#: A schema-valid payload never reaches it (``status_reason`` is an enum), so it is a fallback,
#: not a sentence the report is expected to print.
NO_REASON_RECORDED = "no reason recorded"

#: ``status_reason`` values in plain words, interpolated wherever the report says why an arm did
#: not run (the status table of the plan's section 1.3). ``donor_spans_conditions`` carries the
#: whole reason rather than its first clause: a donor measured under both conditions makes the
#: design paired, pbcheck v0.1.0 implements no paired or mixed model, and the donor-permutation
#: null would treat one donor's two halves as independent.
#: :mod:`pbcheck.render.sections` renders these phrases and :mod:`pbcheck.audit` interpolates one
#: of them into the read-out line of a run whose pseudobulk arm never started.
STATUS_REASON_WORDS = MappingProxyType({
    "design_only_requested": "a metadata-only run was requested",
    "donor_spans_conditions": (
        "at least one donor was measured under both conditions, which makes this a paired design; "
        "a paired design needs a paired or mixed model, which pbcheck does not implement in this "
        "release, and the donor-permutation null would treat the two halves of one donor as "
        "independent"
    ),
    "too_few_donors": "fewer than the minimum donors were present in at least one group",
    "non_integer_counts": "no counts matrix passed the raw-count check",
    "universe_too_small": "the frozen gene universe was too small to proceed",
    "too_few_profiles_after_thin_filter": (
        "too few pseudobulk profiles remained after the thin-donor filter"
    ),
})

#: The report's fixed lines that are not caveats, read-out sentences or glossary entries: the
#: header, the "not run" notes, the read-out ratio line, the two paired-floor explanations, the
#: machinery-check callout and the footer. They live here, with the rest of the report's fixed
#: prose, so one module holds every sentence pbcheck can print and one test sweeps all of them
#: through the wording gate; :mod:`pbcheck.render.sections` fills them by :func:`report_line` and
#: composes no sentence of its own.
REPORT_LINES = MappingProxyType({
    "header_title": "pbcheck audit of {file}",
    "header_version": "pbcheck {version}, generated {generated_utc}.",
    "header_status": "Status: {status}",
    "header_status_with_reason": "Status: {status} ({reason_text})",
    "no_arms_run": "No detection arms were run at this status: {reason_text}.",
    "arm_not_run": "{arm}: not run, because {reason_text}.",
    "counts_check_not_run": "Counts check: not run, because {reason_text}.",
    "ratio_solo": (
        "Real-label calls over the permutation floor: the per-cell arm calls "
        "{naive_ratio} times its own floor, both counts corrected over the whole universe on "
        "their own (solo BH)."
    ),
    "ratio_paired": (
        " The donor-pseudobulk arm calls {pseudobulk_ratio} times its floor, both counts "
        "corrected across the two arms together (paired BH)."
    ),
    "ratio_leak_contaminated": (
        " At fewer than {threshold} donors in a group these ratios are contaminated by the "
        "per-cell leak and describe this run alone."
    ),
    "counts_examples": "Examples of the values checked: {examples}.",
    "no_batch_columns": "No batch columns were provided.",
    "thin_filter_not_run": "Thin-donor filter: not run.",
    "builder_rule": "Builder rule: {rule}",
    "glossary_entry": "{term}: {definition}",
    "paired_floor_not_comparable": (
        "The paired floor is not shown: the pseudobulk arm left {n_na_pseudobulk} genes without a "
        "value, so a paired correction over the two arms would not cover the same genes; the solo "
        "floor above stands alone."
    ),
    "paired_floor_not_measured": (
        "The paired floor is not shown: the paired correction covers every gene of this run, but "
        "the permutation null recorded no paired floor for the per-cell arm; the solo floor above "
        "stands alone."
    ),
    "machinery_check": (
        "Machinery check of the permutation engine, not a criterion of any kind: the inflation "
        "factor of the empirical permutation p-values is {b5_lambda_empirical}."
    ),
    "footer": (
        "Generated by pbcheck {version}. Protocol: docs/PHASE0_SPEC.md and docs/AMENDMENTS.md in "
        "the pbcheck repository."
    ),
})


def report_line(line_id: str, **values: object) -> str:
    """Format ``REPORT_LINES[line_id]`` the way :func:`caveat_text` formats a caveat.

    The same quoting applies: a placeholder listed in :data:`USER_VALUE_PLACEHOLDERS` (the audited
    file's name among them) is rendered as a quoted identifier, so a file or column named after a
    forbidden word is masked by :func:`prose_for_pattern_check` instead of forging a sentence.
    """
    return _format(REPORT_LINES[line_id], values)


#: Plain-language glossary, section 3 of the report ("What these words mean"), one sentence per
#: term, in the order the plan lists them. ``pbcheck.render.sections`` renders one paragraph per
#: entry; nothing here is interpolated, so the wording is fixed at import time.
GLOSSARY: tuple[tuple[str, str], ...] = (
    (
        "lambda",
        "the genomic inflation factor: the ratio of observed test statistics to the null "
        "expectation, with 1.0 meaning no inflation.",
    ),
    (
        "permutation floor",
        "the number of genes a test calls under donor permutation, when condition labels carry "
        "no real signal; it is the test's own false-positive baseline on this data.",
    ),
    (
        "donor-permutation null",
        "the null distribution built by reassigning condition labels between donors, keeping "
        "every donor's cells together, and rerunning the test on each reassignment.",
    ),
    (
        "gene universe",
        "the fixed set of genes both arms are tested and corrected over, frozen before either "
        "arm sees the real labels.",
    ),
    (
        "replication unit",
        "the unit whose independent draws the statistics assume; for donor data that unit is "
        "the donor, not the cell.",
    ),
    (
        "thin-donor filter",
        "the rule that drops a donor's pseudobulk profile when it is built from too few cells "
        "or too few counts, rather than keeping a noisy profile.",
    ),
    (
        "operating envelope",
        "the region of donor count and donor-to-donor variability that Amendment 3 declares for "
        "the pseudobulk arm on synthetic data; each of its points states the donor count per "
        "group at which the power target is reached, on the committed grid or by the derivation "
        "the grid support of that point names.",
    ),
    (
        "sigma_donor",
        "a knob of pbcheck's synthetic simulator for how much donors differ from each other; it "
        "cannot be measured on your data, which is why the envelope question is left open.",
    ),
    (
        "Monte-Carlo SE",
        "the standard error of a quantity estimated from a finite number of permutations; it "
        "shrinks as more permutations are drawn.",
    ),
    (
        "BH convention (solo vs paired)",
        "solo BH corrects an arm's p-values over the whole gene universe on its own; paired BH "
        "corrects both arms together over the genes common to both, so their real-label counts "
        "are directly comparable.",
    ),
)

#: The read-out sentence templates, formatted with ``str.format_map`` by :mod:`pbcheck.audit` from
#: the payload's ``readout`` block and written into ``readout.sentences``; the renderer prints
#: those strings verbatim. Fixed here; not improvised at render time.
SENTENCES: dict[str, str] = {
    "floor_solo": (
        "Under donor permutation, with no true signal to find, the naive per-cell test calls a "
        "median of {median_count} genes ({median_frac_pct}% of the {universe_size}-gene universe; "
        "Monte-Carlo SE {mc_se}) over {n_perm_achieved} permutations{coarse_note}; on the real "
        "labels it calls {real_solo}. Both counts are corrected over the whole universe on their "
        "own (solo BH)."
    ),
    "lambda_naive": (
        "The naive arm's inflation factor lambda is {lambda} (IQR {iqr}): {class_word} "
        "{band_lo} to {band_hi}, which is the donor-pseudobulk arm's band, shown here to describe "
        "the naive number and not as the naive arm's own criterion."
    ),
    "pseudobulk": (
        "The donor-pseudobulk arm's lambda is {lambda} ({class_word} {band_lo} to {band_hi}); its "
        "permutation false-positive rate is {fp_rate} (MC SE {se}) and its floor a median of "
        "{median_count} genes; on the real labels it calls {real_paired} genes, corrected across "
        "both arms together (paired BH)."
    ),
    "pseudobulk_not_run": "The donor-pseudobulk arm was not run: {reason_text}.",
    "donors": (
        "{n_test} donors in {test_level}, {n_ref} in {ref_level}; {n_distinct_splits} distinct "
        "donor splits exist."
    ),
    "profiles": (
        "After the thin-donor filter (fewer than {min_cells} cells or {min_counts} counts), "
        "{p_test} and {p_ref} pseudobulk profiles remain."
    ),
}


def sentence_text(sentence_id: str, **values: object) -> str:
    """Format ``SENTENCES[sentence_id]`` the way :func:`caveat_text` formats a caveat.

    The read-out lines interpolate the user's condition levels, so they go through the same
    quoting of :data:`USER_VALUE_PLACEHOLDERS`: a level named after a forbidden word is rendered
    as a quoted identifier and masked by :func:`prose_for_pattern_check`, exactly as in a caveat.
    """
    return _format(SENTENCES[sentence_id], values)
