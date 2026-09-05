"""Tests for scripts/protocol_safety_check.py: the executable version of the
implementation plan's section 4 protocol-safety checklist."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT_PATH = _REPO_ROOT / "scripts" / "protocol_safety_check.py"

_spec = importlib.util.spec_from_file_location("protocol_safety_check", _SCRIPT_PATH)
protocol_safety_check = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(protocol_safety_check)


def _git(repo: Path, *args: str) -> subprocess.CompletedProcess:
    cp = subprocess.run(
        ["git", "-C", str(repo), *args],
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    assert cp.returncode == 0, f"git {args} failed: {cp.stderr}"
    return cp


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def _make_repo(tmp_path: Path) -> Path:
    """A minimal two-commit git repository with the frozen tree shape the checklist expects."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "test@example.com")
    _git(repo, "config", "user.name", "Test")

    _write(repo / "docs" / "PHASE0_SPEC.md", "spec\n")
    _write(repo / "docs" / "AMENDMENTS.md", "amendments\n")
    _write(repo / "docs" / "PILOT_FINDINGS.md", "findings\n")
    _write(repo / "docs" / "PREREGISTRATION_STRATUM_LIST.md", "strata\n")
    _write(repo / "docs" / "ENV_NOTES.md", "env\n")
    _write(repo / "pilot" / "preregistration" / "list.csv", "a,b\n")
    _write(repo / "pilot" / "gate" / "gate.json", "{}\n")
    _write(repo / "pilot" / "testsel" / "sel.json", "{}\n")
    _write(repo / "pilot" / "README.md", "line1\nline2\nline3\n")
    _write(repo / "src" / "pbcheck" / "gate_config.py", "PRE_REGISTERED = True\n")
    _write(repo / "src" / "pbcheck" / "design.py", "# design\n")
    _write(repo / "src" / "pbcheck" / "gene_universe.py", "# universe\n")
    _write(repo / "src" / "pbcheck" / "permutation.py", "# permutation\n")
    _write(repo / "src" / "pbcheck" / "mtc.py", "# mtc\n")
    _write(repo / "src" / "pbcheck" / "io_counts.py", "# io\n")
    _write(repo / "src" / "pbcheck" / "census_select.py", "# census\n")
    _write(repo / "src" / "pbcheck" / "methods" / "__init__.py", "# methods\n")
    _write(repo / "synthetic" / "oracles.py", "# oracles\n")
    _write(repo / "scripts" / "synthetic_gate.py", "# gate\n")
    _write(repo / "tests" / "conftest.py", "# conftest\n")
    _write(repo / "tests" / "test_something.py", "def test_ok():\n    assert True\n")
    _write(repo / ".pre-commit-config.yaml", "repos: []\n")
    _write(repo / ".gitattributes", "* text=auto\n")

    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "base")
    _git(repo, "branch", "main")

    return repo


def test_pass_when_frozen_paths_untouched(tmp_path):
    repo = _make_repo(tmp_path)
    # a change well outside every frozen row: a new, non-frozen file
    _write(repo / "README.md", "hello\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: unrelated addition")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    fails = [r.line() for r in results if r.status == "FAIL"]
    assert fails == []
    assert ok is True
    assert any(r.status == "PASS" and r.name == "docs/** frozen" for r in results)


def test_fail_when_a_frozen_doc_is_edited(tmp_path):
    repo = _make_repo(tmp_path)
    (repo / "docs" / "AMENDMENTS.md").write_text("edited\n", encoding="utf-8", newline="\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: touches a frozen doc")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    assert ok is False
    by_name = {r.name: r for r in results}
    assert by_name["docs/** frozen"].status == "FAIL"
    assert by_name["docs/AMENDMENTS.md untouched"].status == "FAIL"


def test_skip_when_a_row_input_is_absent(tmp_path):
    repo = _make_repo(tmp_path)
    _git(repo, "commit", "-q", "--allow-empty", "-m", "head: no-op")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    # tests/test_docs.py, scripts/compare_gate_scalars.py etc. do not exist in this repo
    assert by_name["pilot/README.md"].status == "SKIP"
    assert "tests/test_docs.py" in by_name["pilot/README.md"].detail
    assert by_name["metrics.py docstring-only exception"].status == "SKIP"
    assert by_name["Gate numbers do not move"].status == "SKIP"
    assert by_name["Gate numbers do not move"].detail == "requires --with-gate"
    # a SKIPped row must never be reported as PASS
    assert all(r.status != "PASS" for r in results if r.name == "pilot/README.md")
    assert ok is True


def test_engine_module_change_fails(tmp_path):
    repo = _make_repo(tmp_path)
    (repo / "src" / "pbcheck" / "design.py").write_text("# changed\n", encoding="utf-8", newline="\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: touches an engine module")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    assert by_name["Engine modules"].status == "FAIL"
    assert ok is False


def test_conftest_deletion_fails_existing_tests_row(tmp_path):
    repo = _make_repo(tmp_path)
    (repo / "tests" / "conftest.py").write_text("", encoding="utf-8", newline="\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: deletes a line from conftest.py")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    assert by_name["Existing tests unchanged"].status == "FAIL"
    assert ok is False


def test_gate_row_requires_with_gate_flag(tmp_path):
    repo = _make_repo(tmp_path)
    _git(repo, "commit", "-q", "--allow-empty", "-m", "head: no-op")

    results, _ = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=True
    )
    by_name = {r.name: r for r in results}
    # compare_gate_scalars.py does not exist in this synthetic repo
    assert by_name["Gate numbers do not move"].status == "SKIP"
    assert "compare_gate_scalars.py" in by_name["Gate numbers do not move"].detail


def test_cli_exits_nonzero_on_failure(tmp_path):
    repo = _make_repo(tmp_path)
    (repo / "docs" / "AMENDMENTS.md").write_text("edited\n", encoding="utf-8", newline="\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: touches a frozen doc")

    cp = subprocess.run(
        [sys.executable, str(_SCRIPT_PATH), "--base", "main", "--head", "HEAD", "--repo", str(repo)],
        capture_output=True,
        text=True,
    )
    assert cp.returncode == 1, cp.stdout
    assert "docs/** frozen: FAIL" in cp.stdout
    assert "summary:" in cp.stdout


def test_cli_exits_zero_on_success(tmp_path):
    repo = _make_repo(tmp_path)
    _git(repo, "commit", "-q", "--allow-empty", "-m", "head: no-op")

    cp = subprocess.run(
        [sys.executable, str(_SCRIPT_PATH), "--base", "main", "--head", "HEAD", "--repo", str(repo)],
        capture_output=True,
        text=True,
    )
    assert cp.returncode == 0, cp.stdout


def test_real_repository_run_exits_zero():
    """This must exit 0 at the commit this release ships from: rows whose inputs
    (product code, demo docs, new test files) do not exist yet at HEAD report SKIP,
    never PASS or FAIL."""
    cp = subprocess.run(
        [sys.executable, str(_SCRIPT_PATH), "--base", "main", "--head", "HEAD"],
        cwd=_REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert cp.returncode == 0, cp.stdout
    lines = [ln for ln in cp.stdout.splitlines() if ln.strip()]
    # exactly one line per table row plus the trailing summary line, nothing else
    assert len(lines) == len(protocol_safety_check._ROWS) + 1
    assert lines[-1].startswith("summary:")
