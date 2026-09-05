"""Packaging, versioning, release workflow and Zenodo metadata (WP4)."""

from __future__ import annotations

import json
import re
import sys
from importlib import metadata
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import check_version_consistency  # noqa: E402
import changelog_section  # noqa: E402

import pbcheck  # noqa: E402


def test_version_matches_citation_cff() -> None:
    cff_text = (REPO / "CITATION.cff").read_text(encoding="utf-8")
    cff_version = check_version_consistency.extract_cff_version(cff_text)
    assert cff_version is not None, "CITATION.cff has no 'version:' line"
    assert pbcheck.__version__ == cff_version


def test_console_script_entry_point_registered() -> None:
    try:
        metadata.distribution("pbcheck")
    except metadata.PackageNotFoundError:
        pytest.skip("pbcheck is not installed in this environment")

    names = {ep.name for ep in metadata.entry_points(group="console_scripts")}
    if "pbcheck" not in names:
        pytest.skip(
            "no 'pbcheck' console_scripts entry point registered yet: [project.scripts] is "
            "owned by the CLI work package, excluded from WP4's pyproject.toml scope, and has "
            "not landed in this branch; this assertion activates once it does"
        )
    assert "pbcheck" in names


def test_check_version_consistency_passes_with_matching_changelog(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        f"# Changelog\n\n## [{pbcheck.__version__}] - 2026-09-05\n\nFirst release.\n",
        encoding="utf-8",
        newline="\n",
    )
    exit_code = check_version_consistency.main(["--changelog", str(changelog)])
    assert exit_code == 0


def test_check_version_consistency_fails_without_changelog_heading(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text("# Changelog\n\nNothing released yet.\n", encoding="utf-8", newline="\n")
    exit_code = check_version_consistency.main(["--changelog", str(changelog)])
    assert exit_code == 1


def test_check_version_consistency_tag_flag(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        f"## [{pbcheck.__version__}] - 2026-09-05\n",
        encoding="utf-8",
        newline="\n",
    )
    ok = check_version_consistency.main(
        ["--changelog", str(changelog), "--tag", f"v{pbcheck.__version__}"]
    )
    assert ok == 0

    mismatched = check_version_consistency.main(
        ["--changelog", str(changelog), "--tag", "v9.9.9"]
    )
    assert mismatched == 1


def test_check_version_consistency_require_date(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        f"## [{pbcheck.__version__}] - 2026-09-05\n",
        encoding="utf-8",
        newline="\n",
    )
    # CITATION.cff at this point in the repo has no date-released yet (it is set by the
    # integrator in the release commit, per s5_release.md 5.3), so this must fail.
    exit_code = check_version_consistency.main(
        ["--changelog", str(changelog), "--require-date"]
    )
    assert exit_code == 1


def test_zenodo_json_parses_and_has_expected_fields() -> None:
    data = json.loads((REPO / ".zenodo.json").read_text(encoding="utf-8"))
    assert data["upload_type"] == "software"
    assert "doi" not in data
    assert data["creators"] == [{"name": "Poliakov, Alexander"}]
    assert data["access_right"] == "open"
    assert data["license"] == "BSD-3-Clause"
    assert data["related_identifiers"][0]["resource_type"] == "software"


def test_release_workflow_parses_and_has_required_shape() -> None:
    workflow = yaml.safe_load((REPO / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8"))
    jobs = workflow["jobs"]

    for job_name in ("publish-testpypi", "publish-pypi"):
        job = jobs[job_name]
        assert job["permissions"]["id-token"] == "write"
        assert "environment" in job

    check_steps = jobs["check"]["steps"]
    tag_check = next(
        step
        for step in check_steps
        if "--require-date" in step.get("run", "")
    )
    assert tag_check["if"] == "startsWith(github.ref, 'refs/tags/v')"


def test_changelog_section_extracts_the_released_version(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        "# Changelog\n\n"
        "## [Unreleased]\n\nNothing yet.\n\n"
        "## [0.1.0] - 2026-09-05\n\nFirst release.\n\n"
        "## [0.0.1] - 2026-01-01\n\nInitial scaffold.\n",
        encoding="utf-8",
        newline="\n",
    )
    section = changelog_section.extract_section(
        changelog.read_text(encoding="utf-8"), "0.1.0"
    )
    assert section == "First release."


def test_changelog_section_missing_version_exits_1(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text("# Changelog\n\n## [0.1.0]\n\nFirst release.\n", encoding="utf-8", newline="\n")
    exit_code = changelog_section.main(["9.9.9", "--changelog", str(changelog)])
    assert exit_code == 1


@pytest.mark.parametrize("field", ["title", "description", "keywords"])
def test_zenodo_json_fields_are_copied_from_citation_cff(field: str) -> None:
    cff_text = (REPO / "CITATION.cff").read_text(encoding="utf-8")
    zenodo = json.loads((REPO / ".zenodo.json").read_text(encoding="utf-8"))
    if field == "title":
        match = re.search(r'^title:\s*"(.+)"\s*$', cff_text, flags=re.MULTILINE)
        assert match is not None
        assert zenodo["title"] == match.group(1)
    elif field == "keywords":
        keywords_block = re.search(r"^keywords:\n((?:\s+-\s+.+\n)+)", cff_text, flags=re.MULTILINE)
        assert keywords_block is not None
        cff_keywords = [line.strip("- ").strip() for line in keywords_block.group(1).splitlines()]
        assert zenodo["keywords"] == cff_keywords
    else:
        assert zenodo["description"].startswith("pbcheck estimates how much")
