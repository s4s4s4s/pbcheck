"""Checks on the shipped documentation for the v0.1.0 release (WP5).

Every assertion is made against whitespace-normalised text (``re.sub(r"\\s+", " ", text)``), so
these tests are insensitive to how the source markdown happens to wrap a line.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


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


def test_readme_contains_no_numeric_token_from_demo_readouts():
    demo_dir = ROOT / "demo"
    if not demo_dir.exists():
        return  # demo/ is produced by a later work package; nothing to check yet.

    numeric_tokens: set[str] = set()
    for payload_path in demo_dir.glob("**/pbcheck_audit.json"):
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
        readout = payload.get("readout")
        if not readout:
            continue
        for value in readout.values():
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                numeric_tokens.add(str(value))

    text = _normalised(ROOT / "README.md")
    for token in numeric_tokens:
        assert token not in re.findall(r"[-+]?\d[\d.]*", text), (
            f"README contains a numeric token {token!r} that appears in a demo readout"
        )


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


# ---------------------------------------------------------------------------
# .claude/skills/pbcheck-map/SKILL.md
# ---------------------------------------------------------------------------


def test_skill_map_mentions_the_v010_release_surface():
    text = _normalised(ROOT / ".claude" / "skills" / "pbcheck-map" / "SKILL.md")
    for token in ("audit.py", "cli.py", "release.yml", "render/"):
        assert token in text, f"SKILL.md is missing a mention of {token!r}"


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
