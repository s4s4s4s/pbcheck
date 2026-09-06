"""Packaging, versioning, release workflow and Zenodo metadata (WP4)."""

from __future__ import annotations

import json
import re
import sys
import tomllib
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
        distribution = metadata.distribution("pbcheck")
    except metadata.PackageNotFoundError:
        pytest.skip("pbcheck is not installed in this environment")

    names = {ep.name for ep in distribution.entry_points if ep.group == "console_scripts"}
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


def _write_citation(tmp_path: Path, *, date_released: str | None) -> Path:
    citation = tmp_path / "CITATION.cff"
    body = (
        "cff-version: 1.2.0\n"
        'title: "pbcheck"\n'
        "type: software\n"
        f"version: {pbcheck.__version__}\n"
    )
    if date_released is not None:
        body += f"date-released: {date_released}\n"
    citation.write_text(body, encoding="utf-8", newline="\n")
    return citation


def _write_zenodo(tmp_path: Path, *, version: str) -> Path:
    zenodo = tmp_path / ".zenodo.json"
    zenodo.write_text(json.dumps({"version": version}), encoding="utf-8", newline="\n")
    return zenodo


def test_check_version_consistency_require_date_missing_is_exit_1(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        f"## [{pbcheck.__version__}] - YYYY-MM-DD\n",
        encoding="utf-8",
        newline="\n",
    )
    citation = _write_citation(tmp_path, date_released=None)
    zenodo = _write_zenodo(tmp_path, version=pbcheck.__version__)
    exit_code = check_version_consistency.main(
        [
            "--changelog",
            str(changelog),
            "--citation",
            str(citation),
            "--zenodo",
            str(zenodo),
            "--require-date",
        ]
    )
    assert exit_code == 1


def test_check_version_consistency_require_date_valid_is_exit_0(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        f"## [{pbcheck.__version__}] - 2026-09-10\n",
        encoding="utf-8",
        newline="\n",
    )
    citation = _write_citation(tmp_path, date_released="2026-09-10")
    zenodo = _write_zenodo(tmp_path, version=pbcheck.__version__)
    exit_code = check_version_consistency.main(
        [
            "--changelog",
            str(changelog),
            "--citation",
            str(citation),
            "--zenodo",
            str(zenodo),
            "--require-date",
        ]
    )
    assert exit_code == 0


def test_check_version_consistency_require_date_mismatched_dates_is_exit_1(tmp_path: Path) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        f"## [{pbcheck.__version__}] - 2026-09-10\n",
        encoding="utf-8",
        newline="\n",
    )
    citation = _write_citation(tmp_path, date_released="2026-09-11")
    zenodo = _write_zenodo(tmp_path, version=pbcheck.__version__)
    exit_code = check_version_consistency.main(
        [
            "--changelog",
            str(changelog),
            "--citation",
            str(citation),
            "--zenodo",
            str(zenodo),
            "--require-date",
        ]
    )
    assert exit_code == 1


@pytest.mark.parametrize("malformed_date", ["2026-9-5", "05.09.2026"])
def test_check_version_consistency_require_date_malformed_is_exit_1(
    tmp_path: Path, malformed_date: str
) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        f"## [{pbcheck.__version__}] - 2026-09-10\n",
        encoding="utf-8",
        newline="\n",
    )
    citation = _write_citation(tmp_path, date_released=malformed_date)
    zenodo = _write_zenodo(tmp_path, version=pbcheck.__version__)
    exit_code = check_version_consistency.main(
        [
            "--changelog",
            str(changelog),
            "--citation",
            str(citation),
            "--zenodo",
            str(zenodo),
            "--require-date",
        ]
    )
    assert exit_code == 1


def test_check_version_consistency_mismatched_version_between_files_is_exit_1(
    tmp_path: Path,
) -> None:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        f"## [{pbcheck.__version__}] - 2026-09-05\n",
        encoding="utf-8",
        newline="\n",
    )
    citation = _write_citation(tmp_path, date_released=None)
    zenodo = _write_zenodo(tmp_path, version="9.9.9")
    exit_code = check_version_consistency.main(
        [
            "--changelog",
            str(changelog),
            "--citation",
            str(citation),
            "--zenodo",
            str(zenodo),
        ]
    )
    assert exit_code == 1


def test_check_version_consistency_missing_citation_file_is_exit_2(tmp_path: Path) -> None:
    exit_code = check_version_consistency.main(
        ["--citation", str(tmp_path / "does-not-exist.cff")]
    )
    assert exit_code == 2


def test_check_version_consistency_package_not_found_is_exit_2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def _raise(_name: str) -> str:
        raise metadata.PackageNotFoundError("pbcheck")

    monkeypatch.setattr(check_version_consistency.metadata, "version", _raise)
    citation = _write_citation(tmp_path, date_released=None)
    zenodo = _write_zenodo(tmp_path, version=pbcheck.__version__)
    exit_code = check_version_consistency.main(
        ["--citation", str(citation), "--zenodo", str(zenodo)]
    )
    assert exit_code == 2


def test_check_version_consistency_main_passes_against_the_repository_files() -> None:
    """The Done criterion the plan requires: this must be true on the merged branch."""
    assert check_version_consistency.main([]) == 0


def test_zenodo_json_parses_and_has_expected_fields() -> None:
    data = json.loads((REPO / ".zenodo.json").read_text(encoding="utf-8"))
    assert data["upload_type"] == "software"
    assert "doi" not in data
    assert data["creators"] == [{"name": "Poliakov, Alexander"}]
    assert data["access_right"] == "open"
    # Zenodo's own vocabulary id is lowercase; a differently-cased string is silently dropped
    # or rejected by Zenodo's legacy .zenodo.json reader.
    assert data["license"] == "bsd-3-clause"
    assert data["version"] == pbcheck.__version__
    assert data["related_identifiers"][0]["resource_type"] == "software"


def test_release_workflow_parses_and_has_required_shape() -> None:
    workflow = yaml.safe_load((REPO / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8"))
    jobs = workflow["jobs"]
    workflow_text = (REPO / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")

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
    assert tag_check["if"] == "github.event_name == 'push' && startsWith(github.ref, 'refs/tags/v')"

    # The tag-format guard: only v<major>.<minor>.<patch> reaches the release pipeline.
    tag_guard = next(
        step
        for step in check_steps
        if "v[0-9]+.[0-9]+.[0-9]+" in step.get("run", "")
    )
    assert tag_guard["if"] == "github.event_name == 'push' && startsWith(github.ref, 'refs/tags/v')"

    # publish-pypi: only a real tag push, never workflow_dispatch (a dispatch against a tag ref
    # would otherwise satisfy a bare startsWith(github.ref, ...) guard).
    assert jobs["publish-pypi"]["if"] == "github.event_name == 'push' && startsWith(github.ref, 'refs/tags/v')"

    # publish-testpypi: only workflow_dispatch, and only when the operator chose 'testpypi'.
    assert jobs["publish-testpypi"]["if"] == "github.event_name == 'workflow_dispatch' && inputs.target == 'testpypi'"

    # No "secrets." reference anywhere: trusted publishing only (id-token: write), the only
    # token in the whole workflow is the ambient github.token.
    assert "secrets." not in workflow_text

    # The needs chain: check -> build -> {publish-testpypi, publish-pypi} -> github-release.
    assert jobs["build"]["needs"] == "check"
    assert jobs["publish-testpypi"]["needs"] == "build"
    assert jobs["publish-pypi"]["needs"] == "build"
    assert jobs["github-release"]["needs"] == "publish-pypi"

    # The GitHub Release is created from the built artifacts, not re-derived.
    release_step = next(
        step
        for step in jobs["github-release"]["steps"]
        if "gh release create" in step.get("run", "")
    )
    assert "dist/*" in release_step["run"]

    # The publish action is pinned by commit SHA (40 hex chars), not a moving branch/tag ref.
    for job_name in ("publish-testpypi", "publish-pypi"):
        publish_step = next(
            step
            for step in jobs[job_name]["steps"]
            if step.get("uses", "").startswith("pypa/gh-action-pypi-publish@")
        )
        pinned_ref = publish_step["uses"].split("@", 1)[1]
        assert re.fullmatch(r"[0-9a-f]{40}", pinned_ref), pinned_ref


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


def _fold_cff_abstract(cff_text: str) -> str:
    abstract_block = re.search(r"^abstract:\s*>-\n((?:[ \t]+.+\n)+)", cff_text, flags=re.MULTILINE)
    assert abstract_block is not None
    return " ".join(line.strip() for line in abstract_block.group(1).splitlines())


@pytest.mark.parametrize("field", ["title", "description", "keywords", "license", "creators", "repository"])
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
    elif field == "description":
        assert zenodo["description"] == _fold_cff_abstract(cff_text)
    elif field == "license":
        match = re.search(r"^license:\s*(\S+)\s*$", cff_text, flags=re.MULTILINE)
        assert match is not None
        # CITATION.cff carries the SPDX form (BSD-3-Clause); .zenodo.json carries Zenodo's
        # lowercase vocabulary id for the same licence.
        assert zenodo["license"] == match.group(1).lower()
    elif field == "creators":
        match = re.search(
            r"^\s*-\s*family-names:\s*(\S+)\n\s*given-names:\s*(\S+)",
            cff_text,
            flags=re.MULTILINE,
        )
        assert match is not None
        family_names, given_names = match.group(1), match.group(2)
        assert zenodo["creators"][0]["name"] == f"{family_names}, {given_names}"
    else:  # repository
        match = re.search(r'^repository-code:\s*"?(\S+?)"?\s*$', cff_text, flags=re.MULTILINE)
        assert match is not None
        assert zenodo["related_identifiers"][0]["identifier"] == match.group(1)


def test_pyproject_keywords_match_citation_and_zenodo() -> None:
    pyproject = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
    pyproject_keywords = pyproject["project"]["keywords"]

    cff_text = (REPO / "CITATION.cff").read_text(encoding="utf-8")
    keywords_block = re.search(r"^keywords:\n((?:\s+-\s+.+\n)+)", cff_text, flags=re.MULTILINE)
    assert keywords_block is not None
    cff_keywords = [line.strip("- ").strip() for line in keywords_block.group(1).splitlines()]

    zenodo = json.loads((REPO / ".zenodo.json").read_text(encoding="utf-8"))

    assert pyproject_keywords == cff_keywords == zenodo["keywords"]


def test_pyproject_dependency_lower_bounds_match_the_tested_stack() -> None:
    pyproject = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
    dependencies = {
        re.match(r"^[A-Za-z0-9_.-]+", dep).group(0): dep
        for dep in pyproject["project"]["dependencies"]
    }
    lock_versions = {}
    for line in (REPO / "requirements.lock").read_text(encoding="utf-8").splitlines():
        match = re.match(r"^([A-Za-z0-9_.-]+)==(\S+)$", line.strip())
        if match:
            lock_versions[match.group(1)] = match.group(2)

    expected_lower_bounds = {
        "numpy": "2.5",
        "pandas": "3.0",
        "scipy": "1.18",
        "statsmodels": "0.14",
    }
    for name, bound in expected_lower_bounds.items():
        assert dependencies[name] == f"{name}>={bound}"
        assert lock_versions[name].startswith(bound.rsplit(".", 1)[0])


def test_pyproject_uses_pep639_license_expression() -> None:
    pyproject_text = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    pyproject = tomllib.loads(pyproject_text)
    assert pyproject["project"]["license"] == "BSD-3-Clause"
    assert pyproject["project"]["license-files"] == ["LICENSE"]
    # PEP 639: a project using the license-expression form does not also carry the classifier.
    assert "License :: OSI Approved :: BSD License" not in pyproject["project"]["classifiers"]


def test_pyproject_sdist_excludes_claude_directory() -> None:
    pyproject = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
    sdist_exclude = pyproject["tool"]["hatch"]["build"]["targets"]["sdist"]["exclude"]
    assert ".claude/**" in sdist_exclude


def test_documentation_project_url_target_exists_on_the_merged_branch() -> None:
    pyproject_text = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    pyproject = tomllib.loads(pyproject_text)
    for name, url in pyproject["project"]["urls"].items():
        prefix = "https://github.com/s4s4s4s/pbcheck/blob/main/"
        if not url.startswith(prefix):
            continue
        relative_path = url[len(prefix) :]
        assert (REPO / relative_path).exists(), f"{name} Project-URL points at a missing file: {url}"
