"""The fixed prose of the pbcheck report: caveats, the action sentence, forbidden patterns.

These texts are a protocol surface (the wording-rules section of the v0.1.0 implementation plan):
they are tested for the forbidden patterns listed below and are changed only together with those
tests. ``pbcheck.audit`` (built in
parallel by another work package, in this same release) imports ``CAVEATS`` and
``LAMBDA_CLASS_WORDS`` to fill the payload's ``caveats`` and ``readout`` blocks;
``pbcheck.render`` (also built in parallel) renders them. Neither of those two modules is
imported here, so this module can be finished and tested independently of both.
"""

from __future__ import annotations

from pbcheck import gate_config

#: Caveat and action-sentence templates, formatted with ``str.format_map`` by :func:`caveat_text`.
#: C1, C2, C3, C7 are rendered in every report; the rest are conditional (the plan's wording rules).
CAVEATS: dict[str, str] = {
    "C1": (
        "This is a pbcheck v{version} audit of one stratum of one file. It measures how a naive "
        "per-cell test and a donor-pseudobulk test behave on this data under a donor-permutation "
        "null. It is run outside pbcheck's pre-registered Phase 0 protocol, is not a Phase 0 "
        "measurement, and makes no claim about any publication or about whether any reported "
        "result is true."
    ),
    "C2": (
        "The pseudobulk arm's calibration and power were established on synthetic oracles only "
        "inside the operating envelope declared in Amendment 3 ({envelope_rows}). pbcheck does "
        "not estimate sigma_donor for real data, so whether this stratum lies inside that "
        "envelope is not determined here."
    ),
    "C3": (
        "Everything in this report is a diagnostic of this file: the naive arm's inflation factor "
        "and permutation floor, the pseudobulk arm's inflation factor, floor and false-positive "
        "rate as its negative control, and the bookkeeping of the shared gene universe. None of "
        "it is a Phase 0 result, and pbcheck does not rate this file against the Phase 0 decision "
        "rule."
    ),
    "C4": (
        "Fewer than 8 donors in at least one group: the permutation null has few distinct donor "
        "splits, floors are coarse, and pbcheck's own protocol treats floor-based comparisons "
        "below 8 donors per group as contaminated by the per-cell leak. Do not compare these "
        "numbers with a run on another file."
    ),
    "C5": (
        "The design audit found a batch column that separates the two conditions ({cols}): "
        "condition and batch cannot be told apart in this stratum, and neither arm's real-label "
        "result can be read as a condition effect."
    ),
    "C6": (
        "The count matrix failed the raw-count check ({reason}). The pseudobulk arm was dropped, "
        "never rounded, so this run has no negative control. The naive arm was run on the matrix "
        "as found, with the naive pipeline's own normalisation and log transform applied on top "
        "of it; its numbers describe that pipeline on this matrix and are not comparable to a run "
        "on raw counts."
    ),
    "C7": (
        "This instrument was calibrated on synthetic oracles "
        "(pilot/gate/synthetic_gate_2026-08-15.json). Of the settings below, only those listed "
        "under 'protocol constants' are pre-registered values taken from pbcheck.gate_config "
        "({protocol_constant_names}); the permutation counts, the universe filter parameters, the "
        "fallback universe rule and the display settings are the tool's own and are not protocol "
        "values."
    ),
    "C8": (
        "No cell type was selected, so all cells were pooled into one stratum; the file has a "
        "column that looks like a cell-type annotation ({col}). Pooling mixes composition shifts "
        "between conditions into the contrast; rerun with --celltype {col} --celltype-value "
        "<level> for a per-cell-type audit."
    ),
    "C9": (
        "Only {n} distinct donor splits exist for this design, so the permutation null has {n} "
        "draws, not the {requested} requested; the floor's Monte-Carlo standard error ({se}) is "
        "correspondingly large."
    ),
    "A1": (
        "On this file, with condition labels shuffled between donors and therefore no real signal "
        "to find, the per-cell test still calls a median of {floor_solo} of {G} genes at FDR "
        "{alpha} ({floor_pct}%); on the real labels it calls {real_solo}. A gene list produced by "
        "a per-cell test on this data cannot be separated from that floor. The donor is the "
        "replication unit this design supports; the donor-pseudobulk test called {pb_real} genes "
        "on the real labels against a permutation median of {pb_floor}."
    ),
}

#: Prose patterns the rendered report must never contain (the plan's wording rules). Applied by
#: ``test_render_prose_has_no_forbidden_patterns`` and ``test_cli_report_prose_has_no_forbidden_patterns``
#: (both in the rendering work package) with ``re.search``, case-sensitive where capitalised.
FORBIDDEN_PATTERNS: tuple[str, ...] = (
    r"\bvalid\b(?! within)",
    r"\bfindings?\b",
    r"false discover",
    r"\bthe published\b",
    r"\bproves?\b",
    r"\bGO\b",
    r"NO-GO",
    r"risk[ _]score",
    r"INSTRUMENT VALID",
    r"\bTier\b",
    r"\breportab",
    r"\bconcordance\b",
    r"\bJaccard\b",
    r"\boverlap\b",
)

#: The three words a lambda class is ever rendered as (plan section 1.5, ``readout.lambda_*_class``).
LAMBDA_CLASS_WORDS = ("calibrated", "inflated", "under")


def envelope_rows() -> str:
    """The C2 ``{envelope_rows}`` fragment, built from ``gate_config.OPERATING_ENVELOPE``.

    No envelope number is typed here: every value comes from the constant, in the order it is
    declared there.
    """
    donors = " / ".join(str(row["min_donors_per_group"]) for row in gate_config.OPERATING_ENVELOPE)
    sigmas = " / ".join(str(row["sigma_donor"]) for row in gate_config.OPERATING_ENVELOPE)
    return f"minimum donors per group {donors} at sigma_donor {sigmas}"


def envelope_sentence() -> str:
    """The C2 caveat filled in with :func:`envelope_rows`."""
    return caveat_text("C2", envelope_rows=envelope_rows())


def caveat_text(caveat_id: str, **values: object) -> str:
    """Format ``CAVEATS[caveat_id]`` with ``values``.

    Raises ``KeyError`` naming the missing placeholder when a template references a name not
    present in ``values`` (the default behaviour of ``str.format_map`` on a plain ``dict``).
    """
    return CAVEATS[caveat_id].format_map(values)


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
        "the region of donor count and donor-to-donor variability where the pseudobulk arm's "
        "calibration and power were established on synthetic data.",
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

#: The read-out sentence templates (plan section, "Read-out templates (SENTENCES), fixed here"),
#: formatted with ``str.format_map`` by :mod:`pbcheck.render.sections` from the payload's
#: ``readout`` block and neighbouring fields. Fixed verbatim; not improvised at render time.
SENTENCES: dict[str, str] = {
    "floor_solo": (
        "Under donor permutation, with no true signal to find, the naive per-cell test calls a "
        "median of {median_count} genes ({median_frac_pct}% of the {G}-gene universe; "
        "Monte-Carlo SE {mc_se}) over {n_perm_achieved} permutations{coarse_note}; on the real "
        "labels it calls {real_solo}."
    ),
    "lambda_naive": (
        "The naive arm's inflation factor lambda is {lambda} (IQR {iqr}): {class_word} against "
        "the band {band_lo}-{band_hi}."
    ),
    "pseudobulk": (
        "The donor-pseudobulk arm's lambda is {lambda} ({class_word}); its permutation "
        "false-positive rate is {fp_rate} (MC SE {se}) and its floor a median of {median_count} "
        "genes; on the real labels it calls {real_paired} genes."
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
