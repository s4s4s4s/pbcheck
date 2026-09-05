"""Tests for the two release-checklist scripts in ``scripts/``.

Both answer a yes/no question the checklist asks before a release, and both are only useful if
their "no" is trustworthy: a docstring-only check that also passed a renamed variable, or a gate
comparison that also passed a moved number, would sign off exactly the change it exists to catch.
So every test here comes in pairs: the case that must pass, and the neighbouring case that must
fail.

The scripts are stdlib-only command-line tools run in the release environment, so they are
exercised as such: through ``subprocess`` on the real interpreter, with the real exit status.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCSTRING_SCRIPT = REPO_ROOT / "scripts" / "check_docstring_only_diff.py"
GATE_SCRIPT = REPO_ROOT / "scripts" / "compare_gate_scalars.py"
PILOT_GATE_JSON = REPO_ROOT / "pilot" / "gate" / "synthetic_gate_2026-08-15.json"

MODULE_V1 = '''"""Module docstring, first version."""


CONSTANT = 3


def f(x, y=2):
    """Add two numbers."""
    total = x + y
    return total


class C:
    """A class."""

    def method(self):
        """Only a docstring in this body."""
'''

#: Same code, every docstring rewritten (module, function, class, docstring-only body).
MODULE_DOCSTRINGS_CHANGED = '''"""Module docstring, second version, saying something else entirely."""


CONSTANT = 3


def f(x, y=2):
    """Add two numbers, and note that pbcheck.audit reads this too."""
    total = x + y
    return total


class C:
    """A class, redescribed."""

    def method(self):
        """A different docstring in this body."""
'''

#: Same docstrings, one renamed local: the change the script exists to catch.
MODULE_CODE_CHANGED = MODULE_V1.replace("    total = x + y\n    return total",
                                        "    result = x + y\n    return result")


def _run(script: Path, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    """Run one checklist script the way the checklist runs it."""
    return subprocess.run(
        [sys.executable, str(script), *args],
        cwd=str(cwd) if cwd is not None else None,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )


def _git(repo: Path, *args: str) -> None:
    """Run git in ``repo``, with identity on the command line so no global config is required."""
    completed = subprocess.run(
        ["git", "-c", "user.name=pbcheck test", "-c", "user.email=test@example.invalid", *args],
        cwd=str(repo), capture_output=True, text=True, encoding="utf-8",
    )
    assert completed.returncode == 0, completed.stderr


@pytest.fixture
def git_repo(tmp_path: Path) -> Path:
    """A throwaway repository with ``mod.py`` committed once, at revision ``base``."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "--quiet")
    (repo / "mod.py").write_text(MODULE_V1, encoding="utf-8", newline="\n")
    _git(repo, "add", "mod.py")
    _git(repo, "commit", "--quiet", "-m", "first")
    _git(repo, "tag", "base")
    return repo


def _commit_version(repo: Path, source: str) -> None:
    (repo / "mod.py").write_text(source, encoding="utf-8", newline="\n")
    _git(repo, "add", "mod.py")
    _git(repo, "commit", "--quiet", "-m", "second")


def test_docstring_only_change_passes(git_repo: Path):
    """Every docstring in the file rewritten, no statement moved: the checklist's allowed edit."""
    _commit_version(git_repo, MODULE_DOCSTRINGS_CHANGED)

    result = _run(DOCSTRING_SCRIPT, "base", "HEAD", "mod.py", cwd=git_repo)

    assert result.returncode == 0, result.stdout + result.stderr
    assert "docstring-only" in result.stdout


def test_unchanged_file_passes(git_repo: Path):
    """No change at all is a docstring-only change: the trees are identical."""
    (git_repo / "other.py").write_text("x = 1\n", encoding="utf-8", newline="\n")
    _git(git_repo, "add", "other.py")
    _git(git_repo, "commit", "--quiet", "-m", "unrelated")

    result = _run(DOCSTRING_SCRIPT, "base", "HEAD", "mod.py", cwd=git_repo)

    assert result.returncode == 0, result.stdout + result.stderr


def test_renamed_local_variable_fails(git_repo: Path):
    """A rename inside a function body leaves the diff looking cosmetic; the AST does not."""
    _commit_version(git_repo, MODULE_CODE_CHANGED)

    result = _run(DOCSTRING_SCRIPT, "base", "HEAD", "mod.py", cwd=git_repo)

    assert result.returncode == 1
    assert "not only docstrings" in result.stdout


def test_added_statement_fails(git_repo: Path):
    """A statement appended after the docstring is not a docstring change."""
    _commit_version(git_repo, MODULE_V1.replace('    """Add two numbers."""\n',
                                                '    """Add two numbers."""\n    x = abs(x)\n'))

    result = _run(DOCSTRING_SCRIPT, "base", "HEAD", "mod.py", cwd=git_repo)

    assert result.returncode == 1


def test_docstring_removed_from_body_that_becomes_empty(git_repo: Path):
    """Deleting the only statement of a body is a code change, not a docstring change.

    ``C.method``'s body is nothing but its docstring. Stripping it leaves an empty body, which the
    script replaces with ``pass`` on both sides; a version that dropped the method entirely must
    still come out different.
    """
    _commit_version(git_repo, MODULE_V1.replace(
        '    def method(self):\n        """Only a docstring in this body."""\n', ""))

    result = _run(DOCSTRING_SCRIPT, "base", "HEAD", "mod.py", cwd=git_repo)

    assert result.returncode == 1


def test_docstring_only_usage_and_git_errors(git_repo: Path):
    """Two arguments, or a revision that does not exist, are exit 2, not a verdict."""
    usage = _run(DOCSTRING_SCRIPT, "base", "HEAD", cwd=git_repo)
    assert usage.returncode == 2
    assert "usage" in usage.stderr

    unknown_revision = _run(DOCSTRING_SCRIPT, "no-such-rev", "HEAD", "mod.py", cwd=git_repo)
    assert unknown_revision.returncode == 2
    assert unknown_revision.stderr.strip()

    unknown_path = _run(DOCSTRING_SCRIPT, "base", "HEAD", "not_a_file.py", cwd=git_repo)
    assert unknown_path.returncode == 2


def test_docstring_only_reports_unparsable_revision(git_repo: Path):
    """A revision whose file does not parse is an error, never a silent "code changed"."""
    _commit_version(git_repo, "def f(:\n    pass\n")

    result = _run(DOCSTRING_SCRIPT, "base", "HEAD", "mod.py", cwd=git_repo)

    assert result.returncode == 2
    assert "does not parse" in result.stderr


# ---------------------------------------------------------------------------
# compare_gate_scalars.py
# ---------------------------------------------------------------------------

#: A pilot-shaped artifact: the eight compared scalars at the precisions the real pilot file
#: records them at (a full double repr for the lambdas, two decimals for the power sensitivity),
#: plus keys the comparison must ignore.
PILOT_JSON = {
    "generated_utc": "2026-08-15T18:43:10+00:00",
    "runtime_seconds": 511.8,
    "lambda_naive": 54.571780406237494,
    "lambda_pseudobulk": 1.0065633454721827,
    "pb_perm_null_fp_rate": 0.035,
    "naive_perm_floor": {"median_count": 1161.5, "iqr_count": 22.25},
    "pb_perm_floor": {"median_count": 0.0, "iqr_count": 0.0},
    "power_sensitivity": 0.86,
    "power_empirical_fdr": 0.044444444444444446,
    "verdict": "the pilot's verdict string, copied as recorded",
}


def _write_json(path: Path, payload: dict) -> Path:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8", newline="\n")
    return path


@pytest.fixture
def pilot_file(tmp_path: Path) -> Path:
    return _write_json(tmp_path / "pilot.json", PILOT_JSON)


def test_gate_scalars_identical_files_pass(pilot_file: Path, tmp_path: Path):
    """The same numbers on both sides, and every one of them named in the table."""
    new = _write_json(tmp_path / "new.json", PILOT_JSON)

    result = _run(GATE_SCRIPT, str(new), str(pilot_file))

    assert result.returncode == 0, result.stdout + result.stderr
    assert "all 8 compared scalars are unchanged" in result.stdout
    for name in ("lambda_naive", "lambda_pseudobulk", "pb_perm_null_fp_rate",
                 "naive_perm_floor.median_count", "pb_perm_floor.median_count",
                 "power_sensitivity", "power_empirical_fdr", "verdict"):
        assert name in result.stdout


def test_gate_scalars_ignore_key_set_differences(pilot_file: Path, tmp_path: Path):
    """A rerun adds keys (``bh_mode``, ``naive_perm_floor_solo``, ``naive_engine``) and drops none
    of the compared ones; that is not a difference."""
    payload = json.loads(json.dumps(PILOT_JSON))
    payload["naive_engine"] = "fast"
    payload["naive_perm_floor_solo"] = {"median_count": 1161.5}
    payload["naive_perm_floor"]["bh_mode"] = "paired"
    payload["generated_utc"] = "2026-09-05T00:00:00+00:00"
    del payload["runtime_seconds"]
    new = _write_json(tmp_path / "new.json", payload)

    result = _run(GATE_SCRIPT, str(new), str(pilot_file))

    assert result.returncode == 0, result.stdout + result.stderr


def test_gate_scalars_difference_below_recorded_precision_passes(pilot_file: Path, tmp_path: Path):
    """``power_sensitivity`` is recorded as 0.86: two decimals is what the pilot claims.

    A rerun that reports 0.8600000000000001 has not moved the recorded number, and a comparison
    that called it a difference would fail every release for a float repr.
    """
    payload = json.loads(json.dumps(PILOT_JSON))
    payload["power_sensitivity"] = 0.8600000000000001
    new = _write_json(tmp_path / "new.json", payload)

    result = _run(GATE_SCRIPT, str(new), str(pilot_file))

    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    "key, value",
    [
        ("lambda_naive", 54.57178040623750),
        ("lambda_pseudobulk", 1.0165633454721827),
        ("pb_perm_null_fp_rate", 0.04),
        ("power_sensitivity", 0.85),
        ("power_empirical_fdr", 0.05),
        ("verdict", "a different verdict string"),
    ],
)
def test_gate_scalars_moved_value_fails(pilot_file: Path, tmp_path: Path, key: str, value):
    """Each compared scalar, moved on its own, must be reported and must fail the run."""
    payload = json.loads(json.dumps(PILOT_JSON))
    payload[key] = value
    new = _write_json(tmp_path / "new.json", payload)

    result = _run(GATE_SCRIPT, str(new), str(pilot_file))

    assert result.returncode == 1, result.stdout
    assert "DIFFERS" in result.stdout
    assert "1 of 8 compared scalars moved" in result.stdout


def test_gate_scalars_nested_floor_moved_fails(pilot_file: Path, tmp_path: Path):
    """The two floors are nested one level down and are compared by path, not by name."""
    payload = json.loads(json.dumps(PILOT_JSON))
    payload["naive_perm_floor"]["median_count"] = 1160.5
    payload["pb_perm_floor"]["median_count"] = 1.0
    new = _write_json(tmp_path / "new.json", payload)

    result = _run(GATE_SCRIPT, str(new), str(pilot_file))

    assert result.returncode == 1
    assert "2 of 8 compared scalars moved" in result.stdout


def test_gate_scalars_missing_key_in_new_file_fails(pilot_file: Path, tmp_path: Path):
    """A compared scalar the new run did not produce is a difference, not something to skip."""
    payload = json.loads(json.dumps(PILOT_JSON))
    del payload["pb_perm_floor"]
    new = _write_json(tmp_path / "new.json", payload)

    result = _run(GATE_SCRIPT, str(new), str(pilot_file))

    assert result.returncode == 1
    assert "(missing)" in result.stdout


def test_gate_scalars_usage_and_unreadable_files(pilot_file: Path, tmp_path: Path):
    """One argument, a file that is not there, and a file that is not JSON are all exit 2."""
    usage = _run(GATE_SCRIPT, str(pilot_file))
    assert usage.returncode == 2
    assert "usage" in usage.stderr

    absent = _run(GATE_SCRIPT, str(tmp_path / "nope.json"), str(pilot_file))
    assert absent.returncode == 2

    broken = tmp_path / "broken.json"
    broken.write_text("{not json", encoding="utf-8", newline="\n")
    assert _run(GATE_SCRIPT, str(broken), str(pilot_file)).returncode == 2


def test_gate_scalars_frozen_pilot_artifact_compares_to_itself():
    """The real artifact against itself: the compared paths exist in the file as shipped.

    This is what makes the fixture above more than a description of the script's own conventions:
    if the frozen gate JSON ever stopped carrying one of these eight paths, the checklist row would
    silently compare nothing, and this test fails instead.
    """
    result = _run(GATE_SCRIPT, str(PILOT_GATE_JSON), str(PILOT_GATE_JSON))

    assert result.returncode == 0, result.stdout + result.stderr
    assert "(missing)" not in result.stdout
