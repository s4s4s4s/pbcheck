"""Checks on the shipped documentation for the v0.1.0 release (WP5, fix pass CORE-4).

Every assertion is made against whitespace-normalised text (``re.sub(r"\\s+", " ", text)``), so
these tests are insensitive to how the source markdown happens to wrap a line.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
REPO_BLOB_PREFIX = "https://github.com/s4s4s4s/pbcheck/blob/main/"


def _normalised(path: Path) -> str:
    """Whitespace-collapsed text, with markdown blockquote markers stripped first.

    A quoted block like ``> line one\\n> line two`` must read the same as the plain prose it
    quotes; without stripping the ``>`` prefixes, collapsing whitespace alone would leave a
    stray ``>`` glued into the middle of the sentence at every wrapped line.
    """
    raw = path.read_text(encoding="utf-8")
    unquoted = re.sub(r"(?m)^\s*>\s?", "", raw)
    return re.sub(r"\s+", " ", unquoted)


def _no_em_dash(path: Path) -> bool:
    em_dash = chr(0x2014)
    return em_dash not in path.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# README.md
# ---------------------------------------------------------------------------


def test_readme_contains_install_and_quickstart_commands():
    text = _normalised(ROOT / "README.md")
    assert "pip install pbcheck" in text
    assert "pbcheck example" in text


def test_readme_does_not_call_itself_unreleased():
    text = _normalised(ROOT / "README.md")
    assert "not a released tool" not in text


def test_readme_states_no_risk_score_and_nowhere_else():
    text = _normalised(ROOT / "README.md")
    assert "No `risk_score` exists" in text
    occurrences = text.count("risk_score")
    # The one occurrence is inside "No `risk_score` exists and none is promised."; no other
    # mention of risk_score anywhere else in the file.
    assert occurrences == 1, f"expected exactly one 'risk_score' mention, found {occurrences}"


def test_readme_quotes_envelope_sentence():
    from pbcheck.render.text import envelope_sentence

    text = _normalised(ROOT / "README.md")
    quoted = re.sub(r"\s+", " ", envelope_sentence())
    assert quoted in text


def test_usage_quotes_envelope_sentence():
    from pbcheck.render.text import envelope_sentence

    text = _normalised(ROOT / "docs" / "USAGE.md")
    quoted = re.sub(r"\s+", " ", envelope_sentence())
    assert quoted in text


def test_readme_states_the_headline_numbers_correctly():
    """r_wp5 major 2: floor is a count, the ratio is the real-over-floor field, lambda is separate."""
    text = _normalised(ROOT / "README.md")
    assert "how many genes a test still calls significant under that shuffling" in text
    assert "real_over_floor" in text
    assert "genomic inflation" in text
    assert "a separate summary of the shape of" in text


def test_readme_quickstart_command_has_no_line_continuation():
    """r_wp5 minor 8: a trailing backslash is not a line continuation in cmd or PowerShell."""
    text = ROOT.joinpath("README.md").read_text(encoding="utf-8")
    quickstart = text.split("## Quickstart", 1)[1].split("## What the numbers mean", 1)[0]
    assert "\\\n" not in quickstart


def _readout_numeric_tokens(readout: object, tokens: set[str]) -> None:
    """Recursively collect every numeric leaf of a readout-shaped structure, 4 significant figures.

    r_wp5 blocker 1: the guard walked ``readout.values()`` one level deep only, missing every
    number nested one level down (``naive_floor_solo`` is itself a dict) and missing
    ``permutation_null`` entirely. This walks both, recursively, and rounds to 4 significant
    figures so a README author's rounded quote (54.57) still matches the raw payload value
    (54.5678).
    """
    if isinstance(readout, dict):
        for value in readout.values():
            _readout_numeric_tokens(value, tokens)
    elif isinstance(readout, (list, tuple)):
        for value in readout:
            _readout_numeric_tokens(value, tokens)
    elif isinstance(readout, bool):
        return
    elif isinstance(readout, (int, float)):
        if readout == 0:
            tokens.add("0")
            return
        from math import floor, log10

        digits = 4
        magnitude = floor(log10(abs(readout)))
        rounded = round(readout, -magnitude + (digits - 1))
        tokens.add(str(rounded))
        if isinstance(rounded, float) and rounded == int(rounded):
            tokens.add(str(int(rounded)))


def test_readme_contains_no_numeric_token_from_demo_readouts():
    demo_dir = ROOT / "demo"
    if not demo_dir.exists():
        pytest.skip("demo/ is produced by a later work package; nothing to check yet.")

    numeric_tokens: set[str] = set()
    for payload_path in demo_dir.glob("**/pbcheck_audit.json"):
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
        _readout_numeric_tokens(payload.get("readout"), numeric_tokens)
        _readout_numeric_tokens(payload.get("permutation_null"), numeric_tokens)

    text = _normalised(ROOT / "README.md")
    found_tokens = set(re.findall(r"[-+]?\d[\d.]*", text))
    overlap = numeric_tokens & found_tokens
    assert not overlap, (
        f"README contains numeric token(s) {sorted(overlap)!r} that appear in a demo readout"
    )


# ---------------------------------------------------------------------------
# README.md links: PyPI renders README.md verbatim, so a relative link 404s there.
# ---------------------------------------------------------------------------

_MARKDOWN_LINK = re.compile(r"\[[^\]]*\]\(([^)]+)\)")


def test_readme_has_no_relative_markdown_link():
    text = ROOT.joinpath("README.md").read_text(encoding="utf-8")
    for target in _MARKDOWN_LINK.findall(text):
        assert target.startswith(("http://", "https://")), (
            f"README.md has a relative link {target!r}; PyPI renders README.md verbatim and "
            "would 404 on it (see docs/USAGE.md and CONTRIBUTING.md for the release rationale)"
        )


def test_readme_absolute_repo_links_resolve_to_existing_files():
    # demo/ is produced by a later work package (see the demo-readout guard above); README links
    # to it in advance, so its absence here is expected and not a broken-link regression.
    not_yet_landed = {"demo/README.md"}

    text = ROOT.joinpath("README.md").read_text(encoding="utf-8")
    checked_any = False
    for target in _MARKDOWN_LINK.findall(text):
        if not target.startswith(REPO_BLOB_PREFIX):
            continue
        checked_any = True
        relative = target[len(REPO_BLOB_PREFIX):]
        if relative in not_yet_landed and not (ROOT / relative).exists():
            continue
        assert (ROOT / relative).exists(), f"README.md links {target!r}, which does not exist"
    assert checked_any, "expected at least one absolute repository link in README.md"


# ---------------------------------------------------------------------------
# CHANGELOG.md
# ---------------------------------------------------------------------------


def test_changelog_has_v010_and_prerelease_history_headings():
    text = _normalised(ROOT / "CHANGELOG.md")
    assert "## [0.1.0]" in text
    assert "## Pre-release history" in text


# ---------------------------------------------------------------------------
# CONTRIBUTING.md
# ---------------------------------------------------------------------------


def test_contributing_drops_pypi_release_workflow_omission():
    text = _normalised(ROOT / "CONTRIBUTING.md")
    assert "PyPI release workflow" not in text


def test_contributing_has_dated_release_decision():
    text = _normalised(ROOT / "CONTRIBUTING.md")
    assert "Decision of 2026-09-05" in text


def test_contributing_states_the_release_procedure():
    """c_release_gaps 'release order / CITATION.cff' + r_s4 MINOR-6: the pre-tag gate is documented."""
    text = _normalised(ROOT / "CONTRIBUTING.md")
    assert "## Release procedure" in text
    assert "check_version_consistency.py --tag" in text
    assert "--require-date" in text
    assert "protocol_safety_check.py --with-gate" in text
    assert "Python 3.12" in text


# ---------------------------------------------------------------------------
# SECURITY.md
# ---------------------------------------------------------------------------


def test_security_states_supported_version():
    text = _normalised(ROOT / "SECURITY.md")
    assert "0.1.x" in text
    assert "No supported releases exist" not in text


# ---------------------------------------------------------------------------
# pilot/README.md
# ---------------------------------------------------------------------------


def test_pilot_readme_distinguishes_protocol_from_release():
    text = _normalised(ROOT / "pilot" / "README.md")
    assert "produced under this protocol" in text
    assert "outside the protocol" in text


def test_pilot_readme_at_most_one_deleted_line_against_main():
    """r_wp5 minor: the frozen-adjacent narrative must not be rewritten wholesale."""
    import subprocess

    result = subprocess.run(
        ["git", "diff", "main...HEAD", "--numstat", "--", "pilot/README.md"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0 or not result.stdout.strip():
        pytest.skip("no 'main' ref reachable in this checkout to diff against")
    added, deleted, _ = result.stdout.split(maxsplit=2)
    assert int(deleted) <= 1, f"pilot/README.md has {deleted} deleted lines against main, expected <= 1"


# ---------------------------------------------------------------------------
# .claude/skills/pbcheck-map/SKILL.md
# ---------------------------------------------------------------------------


def test_skill_map_mentions_the_v010_release_surface():
    text = _normalised(ROOT / ".claude" / "skills" / "pbcheck-map" / "SKILL.md")
    for token in ("audit.py", "cli.py", "release.yml", "render/"):
        assert token in text, f"SKILL.md is missing a mention of {token!r}"


def test_skill_map_lists_the_new_test_files():
    """r_wp5 major 5: the map lists every new module and test file, with a real, matching count."""
    text = _normalised(ROOT / ".claude" / "skills" / "pbcheck-map" / "SKILL.md")
    for token in (
        "test_audit.py", "test_audit_schema.py", "test_render.py", "test_render_output.py",
        "test_render_text.py", "test_cli.py", "test_example.py", "test_packaging.py",
        "test_docs.py", "test_checklist_scripts.py", "test_protocol_safety_check.py",
        "test_measure_audit_runtime.py", "markdown.py", "html.py", "sections.py",
        "product_constants.py",
    ):
        assert token in text, f"SKILL.md is missing a mention of {token!r}"


def test_skill_map_test_count_matches_collection():
    """The map's own Done criterion: its 'N passed' line matches a fresh collection count."""
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "--collect-only"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    assert lines, f"pytest --collect-only produced no output; stderr: {result.stderr}"
    match = re.match(r"(\d+) tests? collected", lines[-1])
    assert match, f"unexpected --collect-only output: {lines[-1]!r}"
    collected = int(match.group(1))
    text = _normalised(ROOT / ".claude" / "skills" / "pbcheck-map" / "SKILL.md")
    assert f"**{collected} passed**" in text, (
        f"SKILL.md's passing-test count does not match a fresh collection of {collected}"
    )


# ---------------------------------------------------------------------------
# No em dash anywhere in the release-facing prose
# ---------------------------------------------------------------------------


def test_no_em_dash_in_release_docs():
    for relative in ("README.md", "docs/USAGE.md", "SECURITY.md"):
        path = ROOT / relative
        assert _no_em_dash(path), f"{relative} contains an em dash (U+2014)"

    # CONTRIBUTING.md keeps its pre-existing prose untouched; only the new Releases section
    # (added for v0.1.0) is held to the no-em-dash rule.
    contributing = (ROOT / "CONTRIBUTING.md").read_text(encoding="utf-8")
    releases_section = contributing.split("## Releases", 1)[1]
    em_dash = chr(0x2014)
    assert em_dash not in releases_section, (
        "CONTRIBUTING.md's Releases section contains an em dash (U+2014)"
    )


# ---------------------------------------------------------------------------
# Doc-level forbidden-wording scan (r_wp5 major 4)
# ---------------------------------------------------------------------------
#
# FORBIDDEN_PATTERNS (src/pbcheck/render/text.py) binds pbcheck's *rendered report* prose: a
# v0.1.0 audit must never read as a Phase 0 verdict. These project documents, in contrast, are
# *about* the Phase 0 protocol and the tool that measures it, and legitimately discuss validity,
# the GO/NO-GO decision rule, publications and the literal `risk_score`/`INSTRUMENT VALID` strings
# pbcheck's own report is forbidden from asserting. Each whitelist entry below is one such
# reviewed, accepted sentence, quoted verbatim with the reason it is not a claim pbcheck's product
# surface is forbidden from making; nothing else in these files may match. A whitelist entry that
# stops matching (a stale quote) fails loudly rather than silently widening the pass.

_DOC_FORBIDDEN_WHITELIST: dict[str, tuple[tuple[str, str], ...]] = {
    "README.md": (
        (
            "None of this certifies a published result.",
            "negation: pbcheck's own report does not certify a publication",
        ),
        (
            "A large fraction of published single-cell RNA-seq",
            "describes the literature pbcheck audits, not a claim pbcheck makes",
        ),
        (
            "audit a published study from its metadata alone",
            "describes what the design auditor can do, not a claim of validity",
        ),
        (
            "Instrument validity | **established only within the operating envelope**",
            "quotes the Phase 0 gate's own scoped claim, already qualified by the envelope",
        ),
        (
            "the instrument now states the region it is valid in instead of being judged at "
            "an arbitrary simulator setting",
            "describes the Amendment 3 envelope narrowing, a Phase 0 status statement",
        ),
        (
            'reports `INSTRUMENT VALID WITHIN THE STATED OPERATING ENVELOPE` - never an '
            'unqualified "valid".',
            "quotes the gate script's own literal, already-qualified verdict string",
        ),
        (
            "No `risk_score` exists and none is promised.",
            "negation: states the feature does not exist",
        ),
        (
            "No GO/NO-GO verdict, no envelope-membership verdict, no `sigma_donor` estimate "
            "for real data.",
            "negation: states what the v0.1.0 report does not contain",
        ),
        (
            "| GO / NO-GO | not taken |",
            "Phase 0 status table: the decision has explicitly not been taken",
        ),
        (
            '**GO** - inflation is large and consistent',
            "names the Phase 0 pre-registered decision rule's own two outcomes",
        ),
        (
            'publish the "map of false discoveries".',
            "names the Phase 0 study's stated goal if the GO branch is reached",
        ),
        (
            "**NO-GO** - inflation is small or erratic",
            "names the Phase 0 pre-registered decision rule's own two outcomes",
        ),
        (
            "Confronting false discoveries in single-cell differential expression",
            "a citation title (Squair et al. 2021)",
        ),
    ),
    "docs/USAGE.md": (
        (
            "it was validated on. It states the range and stops there.",
            "negation: pbcheck cannot tell whether a file was inside the validated range",
        ),
        (
            "`pbcheck_audit.json` validates against `pbcheck.audit_schema.validate`",
            "names the Python function `validate`, a code reference, not a claim",
        ),
        (
            "it does not certify any published result.",
            "negation: pbcheck's own report does not certify a publication",
        ),
        (
            "No `risk_score` exists and none is promised.",
            "negation: states the feature does not exist",
        ),
    ),
    "CONTRIBUTING.md": (
        (
            "YAML validity, a large-file guard",
            "describes pre-commit's own hygiene hooks, unrelated to Phase 0 validity",
        ),
    ),
    "pilot/README.md": (
        (
            "The instrument now passes its validity gate",
            "states the Phase 0 gate's own scoped result",
        ),
        (
            "validation on known truth",
            "describes the synthetic-oracle calibration method, not a claim about real data",
        ),
        (
            "does not gate the pseudobulk validity claim",
            "cross-references the Phase 0 gate's own scope",
        ),
        (
            "GO/NO-GO rule, clustered by dataset (D2), pseudobulk validity gate first.",
            "names the Phase 0 pre-registered decision rule module (`decision.py`, not built)",
        ),
        (
            "What remains for the actual GO/NO-GO gate is to run the",
            "names the Phase 0 pre-registered decision rule, not a claim pbcheck makes",
        ),
        (
            "no result here may be read as a finding;",
            "negation: no result under this protocol may be read as a finding",
        ),
        (
            "GO if, on a majority of **independent datasets** at matched cells-per-donor: "
            "pseudobulk is calibrated (perm-null rejects ~0) **and** powered (synthetic "
            "positive), naive λ ≥ 2, and the naive permutation floor is far above the "
            "BH complete-null expectation (~0). NO-GO if pseudobulk fails its validity gate, or "
            "naive inflation vanishes once cells-per-donor is matched (i.e. it was a "
            "depth/cell-count artifact).",
            "the Phase 0 pre-registered decision rule itself, quoted in full, short form",
        ),
    ),
}


def _apply_doc_whitelist(text: str, relative: str) -> str:
    masked = text
    for sentence, _reason in _DOC_FORBIDDEN_WHITELIST.get(relative, ()):
        found = masked.count(sentence)
        assert found >= 1, f"{relative}: whitelisted sentence no longer found: {sentence!r}"
        masked = masked.replace(sentence, "", 1)
    return masked


@pytest.mark.parametrize(
    "relative", ["README.md", "docs/USAGE.md", "CONTRIBUTING.md", "pilot/README.md"]
)
def test_doc_prose_matches_no_forbidden_pattern_outside_the_whitelist(relative: str):
    from pbcheck.render.text import forbidden_pattern_hits

    text = _normalised(ROOT / relative)
    masked = _apply_doc_whitelist(text, relative)
    hits = forbidden_pattern_hits(masked)
    assert not hits, f"{relative} matches forbidden pattern(s) {hits!r} outside the whitelist"


def test_doc_prose_forbidden_pattern_scan_covers_demo_readme_when_present():
    demo_readme = ROOT / "demo" / "README.md"
    if not demo_readme.exists():
        pytest.skip("demo/README.md is produced by a later work package; nothing to check yet.")

    from pbcheck.render.text import forbidden_pattern_hits

    text = _normalised(demo_readme)
    masked = _apply_doc_whitelist(text, "demo/README.md")
    hits = forbidden_pattern_hits(masked)
    assert not hits, f"demo/README.md matches forbidden pattern(s) {hits!r} outside the whitelist"


def test_demonstrations_sentence_drops_bare_as_a_finding():
    text = _normalised(ROOT / "README.md")
    assert "as a finding about the underlying publications" not in text


# ---------------------------------------------------------------------------
# CLI defaults documented independently of the runtime table (r_wp5 major 3)
# ---------------------------------------------------------------------------


def test_usage_states_the_product_permutation_and_alpha_defaults():
    from pbcheck import audit, gate_config

    text = _normalised(ROOT / "docs" / "USAGE.md")
    assert f"default {audit.PRODUCT_N_PERM} (pbcheck.audit.PRODUCT_N_PERM" in text
    assert f"default {audit.PRODUCT_N_PERM_PB} (pbcheck.audit.PRODUCT_N_PERM_PB" in text
    assert f"default {gate_config.ALPHA} (gate_config.ALPHA" in text


def test_usage_states_the_example_command_defaults():
    from pbcheck.example import REFERENCE_SHAPE, SMALL_SHAPE  # noqa: F401 - existence check

    text = _normalised(ROOT / "docs" / "USAGE.md")
    assert "--seed default 0" in text
    assert "--shape default small" in text


def test_usage_lists_the_three_extra_product_constants():
    """r_wp1 minor: the three product constants outside the contract list are documented."""
    from pbcheck import audit

    text = _normalised(ROOT / "docs" / "USAGE.md")
    for name in (
        "PRODUCT_MIN_DONORS_PER_GROUP",
        "PRODUCT_UNIVERSE_MIN_TOTAL_COUNT",
        "PRODUCT_UNIVERSE_MIN_PROP",
    ):
        assert name in text, f"docs/USAGE.md does not name {name}"
    assert "min_donors_per_group" in text
    assert "universe_min_total_count" in text
    assert "universe_min_prop" in text
    assert audit.PRODUCT_MIN_DONORS_PER_GROUP is not None  # constants actually exist


def test_usage_lists_every_status_reason_value():
    from pbcheck.audit_schema import STATUS_REASON_VALUES

    text = _normalised(ROOT / "docs" / "USAGE.md")
    for reason in STATUS_REASON_VALUES:
        assert f"`{reason}`" in text, f"docs/USAGE.md does not list status_reason {reason!r}"


def test_usage_uses_the_neutral_lambda_class_names():
    """X1: the old words never appear as the class descriptors in USAGE."""
    text = _normalised(ROOT / "docs" / "USAGE.md")
    assert "`in_band`" in text
    assert "`above_band`" in text
    assert "`below_band`" in text
