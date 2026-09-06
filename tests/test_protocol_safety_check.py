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


# 16 checklist table rows (r_s4 MINOR-7): pinned to the literal count and the
# exact set of names, so deleting a row from _ROWS cannot keep this green the
# way comparing len(lines) against len(_ROWS) itself could.
_EXPECTED_ROW_NAMES = {
    "docs/** frozen",
    "pilot/preregistration, pilot/gate, pilot/testsel frozen",
    "pilot/README.md",
    "gate_config.PRE_REGISTERED and the file",
    "Stratum list freeze",
    "Engine modules",
    "metrics.py docstring-only exception",
    "Existing tests unchanged",
    "Product constants are not protocol constants",
    "No protocol language in product output",
    "No Census in the product path",
    "Demo datasets outside the freeze",
    "Gate numbers do not move",
    "README carries no demo number",
    "Frozen files excluded from tooling",
    "docs/AMENDMENTS.md untouched",
}


def test_row_reported_names_match_the_plan_table_exactly(tmp_path):
    repo = _make_repo(tmp_path)
    _git(repo, "commit", "-q", "--allow-empty", "-m", "head: no-op")
    results, _ = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    assert len(results) == 16
    assert {r.name for r in results} == _EXPECTED_ROW_NAMES


def test_grep_row_skips_rather_than_passes_when_a_target_is_partly_missing(tmp_path):
    """r_s4 BLOCKER-1: with one of several grep targets present and clean, and
    the rest missing, the row must SKIP, never report a false-clean PASS."""
    repo = _make_repo(tmp_path)
    _write(repo / "src" / "pbcheck" / "render" / "__init__.py", "# render package\n")
    # src/pbcheck/audit.py, cli.py, example.py stay absent
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: add one product target, not all")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    row = by_name["Product constants are not protocol constants"]
    assert row.status == "SKIP"
    assert "missing:" in row.detail
    assert "audit.py" in row.detail
    assert ok is True


def test_grep_row_fails_on_a_real_hit_even_with_other_targets_missing(tmp_path):
    """FAIL still wins over the BLOCKER-1 SKIP fold: a genuine hit is never masked."""
    repo = _make_repo(tmp_path)
    _write(
        repo / "src" / "pbcheck" / "render" / "__init__.py",
        "N_PERM = gate_config.N_PERM  # a hand-typed protocol constant reference\n",
    )
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: leak a protocol constant reference")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    row = by_name["Product constants are not protocol constants"]
    assert row.status == "FAIL"
    assert "render/__init__.py" in row.detail
    assert ok is False


def test_grep_row_scans_non_python_files_under_a_directory_target(tmp_path):
    """r_s4 MAJOR-2: the plan's grep -rnE reads every file, not only *.py."""
    repo = _make_repo(tmp_path)
    for rel in ("audit.py", "cli.py", "example.py"):
        _write(repo / "src" / "pbcheck" / rel, "# stub\n")
    _write(repo / "src" / "pbcheck" / "render" / "__init__.py", "# render package\n")
    _write(repo / "src" / "pbcheck" / "render" / "report.txt.j2", "verdict: NO-GO\n")
    _write(repo / "tests" / "test_render.py", "def test_forbidden():\n    assert True\n")
    _write(repo / "tests" / "test_cli.py", "def test_forbidden():\n    assert True\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: forbidden word in a non-python template")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    row = by_name["No protocol language in product output"]
    assert row.status == "FAIL"
    assert "report.txt.j2" in row.detail
    assert ok is False


def test_grep_row_skips_binary_files_and_reports_them(tmp_path):
    """r_s4 MAJOR-2: a binary file (NUL byte in the first 8 KiB) is skipped
    rather than grepped as text, and the skip is visible in the row's detail."""
    repo = _make_repo(tmp_path)
    for rel in ("audit.py", "cli.py", "example.py"):
        _write(repo / "src" / "pbcheck" / rel, "# stub\n")
    _write(repo / "src" / "pbcheck" / "render" / "__init__.py", "# render package\n")
    _write(repo / "tests" / "test_render.py", "def test_forbidden():\n    assert True\n")
    _write(repo / "tests" / "test_cli.py", "def test_forbidden():\n    assert True\n")
    (repo / "src" / "pbcheck" / "render" / "asset.bin").parent.mkdir(parents=True, exist_ok=True)
    (repo / "src" / "pbcheck" / "render" / "asset.bin").write_bytes(b"NO-GO\x00binary payload")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: forbidden word only inside a binary file")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    row = by_name["No protocol language in product output"]
    assert row.status == "PASS"
    assert "binary skipped" in row.detail
    assert "asset.bin" in row.detail
    assert ok is True


def test_existing_tests_row_exempts_only_files_added_in_the_range(tmp_path):
    """r_s4 MAJOR-3: the exemption is derived from git diff --diff-filter=A,
    not a hard-coded name list; a genuinely new test file is exempt, an edit
    to a pre-existing one still fails the row."""
    repo = _make_repo(tmp_path)
    _write(repo / "tests" / "test_new_thing.py", "def test_new():\n    assert True\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: add a brand-new test file")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    assert by_name["Existing tests unchanged"].status == "PASS"
    assert ok is True


def test_existing_tests_row_fails_on_edits_to_pre_existing_test_files(tmp_path):
    repo = _make_repo(tmp_path)
    (repo / "tests" / "test_something.py").write_text(
        "def test_ok():\n    assert False  # edited\n", encoding="utf-8", newline="\n"
    )
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: edit a pre-existing test file")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    assert by_name["Existing tests unchanged"].status == "FAIL"
    assert ok is False


def test_binary_rewrite_of_pilot_readme_fails_instead_of_zero_deletions(tmp_path):
    """r_s4 MINOR-4: git diff --numstat prints '-' for a binary diff; that
    must be a FAIL, not silently counted as zero deleted lines."""
    repo = _make_repo(tmp_path)
    (repo / "pilot" / "README.md").write_bytes(b"\x00\x01\x02\x03binary now")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: replace pilot/README.md with binary bytes")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    row = by_name["pilot/README.md"]
    assert row.status == "FAIL"
    assert "binary" in row.detail
    assert ok is False


def test_binary_rewrite_of_conftest_fails_existing_tests_row(tmp_path):
    repo = _make_repo(tmp_path)
    (repo / "tests" / "conftest.py").write_bytes(b"\x00\x01\x02\x03binary now")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: replace conftest.py with binary bytes")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    row = by_name["Existing tests unchanged"]
    assert row.status == "FAIL"
    assert "binary" in row.detail
    assert ok is False


def test_pytest_row_reports_an_unmatched_selector_explicitly(tmp_path):
    """r_s4 MINOR-9: a -k selector matching nothing is an explicit FAIL that
    names the selector, not an opaque pytest-exit-5 message."""
    repo = _make_repo(tmp_path)
    _write(
        repo / "tests" / "test_docs.py",
        "def test_unrelated():\n    assert True\n",
    )
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "head: add test_docs.py without a pilot_readme test")

    results, ok = protocol_safety_check.run_checklist(
        repo, "main", "HEAD", tmp_path / "scratch", with_gate=False
    )
    by_name = {r.name: r for r in results}
    row = by_name["pilot/README.md"]
    assert row.status == "FAIL"
    assert "no test matched" in row.detail
    assert "pilot_readme" in row.detail
    assert ok is False


def test_gate_row_fails_on_platform_or_interpreter_mismatch(tmp_path, monkeypatch):
    """r_s4 MINOR-6: the gate artifact was recorded on win32, Python 3.12;
    a reproduction attempted elsewhere is a FAIL, not a silent PASS or SKIP."""
    repo = _make_repo(tmp_path)
    _git(repo, "commit", "-q", "--allow-empty", "-m", "head: no-op")

    monkeypatch.setattr(protocol_safety_check.sys, "platform", "linux")
    row = protocol_safety_check.row_gate_numbers(
        repo=repo, base="main", head="HEAD", scratch=tmp_path / "scratch", with_gate=True
    )
    assert row.status == "FAIL"
    assert "win32" in row.detail
    assert "Python 3.12" in row.detail


def test_invalid_ref_reports_a_fail_row_without_a_traceback(tmp_path):
    """r_s4 MINOR-10: an invalid --base/--head is a FAIL row, not an uncaught
    traceback; exit status stays non-zero either way."""
    repo = _make_repo(tmp_path)
    _git(repo, "commit", "-q", "--allow-empty", "-m", "head: no-op")

    cp = subprocess.run(
        [
            sys.executable, str(_SCRIPT_PATH),
            "--base", "this-ref-does-not-exist",
            "--head", "HEAD",
            "--repo", str(repo),
        ],
        capture_output=True,
        text=True,
    )
    assert cp.returncode == 1
    assert "Traceback" not in cp.stdout
    assert "Traceback" not in cp.stderr
    assert "FAIL" in cp.stdout
    assert "summary:" in cp.stdout


def test_no_em_dash_in_this_areas_own_source():
    """Acceptance ("Every area: no em dash in added lines") checked in Python
    rather than with the shell one-liner the acceptance run uses
    (``git diff ... | grep -P '^\\+.*\\x{2014}'``): on this machine's Git Bash,
    that ``grep -P`` needs an explicit UTF-8 locale (``LANG`` or ``LC_ALL``) to
    evaluate a ``\\x{2014}`` escape at all; with neither set, GNU grep 3.0
    reports "character value in \\x{} or \\o{} is too large" and the pipeline's
    exit status is the trailing ``head``'s, which is 0 regardless, masking a
    real hit rather than proving its absence. Reading the files with Python's
    own text decoding sidesteps the shell locale entirely, so this check is
    the same acceptance rule made to hold regardless of the invoking shell's
    environment. It is scoped to the two files this stage owns; the whole-diff
    check that "Every area" runs is unaffected and still applies at the level
    the repair spec sets it at.
    """
    em_dash = chr(0x2014)
    for path in (_SCRIPT_PATH, _REPO_ROOT / "tests" / "test_protocol_safety_check.py"):
        text = path.read_text(encoding="utf-8")
        assert em_dash not in text, f"em dash (U+2014) found in {path}"
