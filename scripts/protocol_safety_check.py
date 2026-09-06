"""Executable version of the protocol-safety checklist (implementation plan section 4).

Runs every row of that table as its own check and prints one PASS / FAIL / SKIP
line per row, then a summary. Exit status is 0 only when no row FAILs; a row
that SKIPs (because one of its inputs does not exist yet at --head) never
counts as PASS and never fails the run on its own.

Uses only the standard library and subprocess (git is invoked with a list
argv, never through a shell).
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# Priority used when a single checklist row bundles more than one machine
# check (as several rows in the table do): FAIL beats SKIP beats PASS, so a
# row that references a not-yet-existing input can never report PASS.
_STATUS_RANK = {"FAIL": 2, "SKIP": 1, "PASS": 0}


class RowResult:
    __slots__ = ("name", "status", "detail")

    def __init__(self, name: str, status: str, detail: str = "") -> None:
        assert status in ("PASS", "FAIL", "SKIP")
        self.name = name
        self.status = status
        self.detail = detail

    def line(self) -> str:
        if self.status == "PASS":
            return f"{self.name}: PASS"
        if self.status == "SKIP":
            detail = f" {self.detail}" if self.detail else ""
            return f"{self.name}: SKIP{detail}"
        detail = f" ({self.detail})" if self.detail else ""
        return f"{self.name}: FAIL{detail}"


def _combine(name: str, parts: list[RowResult]) -> RowResult:
    """Reduce several sub-checks of one table row into a single RowResult."""
    winner = max(parts, key=lambda r: _STATUS_RANK[r.status])
    details = [p.detail for p in parts if p.status == winner.status and p.detail]
    return RowResult(name, winner.status, "; ".join(details))


def _run_git(repo: Path, args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", "-C", str(repo), *args],
        capture_output=True,
        text=True,
        encoding="utf-8",
    )


def _diff_quiet(repo: Path, base: str, head: str, paths: list[str]) -> bool:
    """True when the two revisions are identical on the given paths."""
    cp = _run_git(repo, ["diff", "--quiet", f"{base}..{head}", "--", *paths])
    if cp.returncode not in (0, 1):
        raise RuntimeError(f"git diff --quiet failed: {cp.stderr.strip()}")
    return cp.returncode == 0


def _numstat(repo: Path, base: str, head: str, path: str) -> list[tuple[int, int, str]]:
    cp = _run_git(repo, ["diff", "--numstat", f"{base}..{head}", "--", path])
    if cp.returncode != 0:
        raise RuntimeError(f"git diff --numstat failed: {cp.stderr.strip()}")
    rows = []
    for line in cp.stdout.splitlines():
        if not line.strip():
            continue
        added, deleted, name = line.split("\t", 2)
        deleted_n = -1 if deleted == "-" else int(deleted)
        added_n = -1 if added == "-" else int(added)
        rows.append((added_n, deleted_n, name))
    return rows


def _log_empty(repo: Path, base: str, head: str, path: str) -> bool:
    cp = _run_git(repo, ["log", f"{base}..{head}", "--", path])
    if cp.returncode != 0:
        raise RuntimeError(f"git log failed: {cp.stderr.strip()}")
    return cp.stdout.strip() == ""


def _diff_row(repo: Path, base: str, head: str, paths: list[str], name: str) -> RowResult:
    if _diff_quiet(repo, base, head, paths):
        return RowResult(name, "PASS")
    return RowResult(name, "FAIL", f"diff on {', '.join(paths)}")


def _pytest_row(
    repo: Path,
    name: str,
    target: str,
    extra_args: list[str] | None = None,
) -> RowResult:
    path = repo / target
    if not path.exists():
        return RowResult(name, "SKIP", f"missing: {target}")
    cmd = [sys.executable, "-m", "pytest", "-q", target, *(extra_args or [])]
    cp = subprocess.run(cmd, cwd=repo, capture_output=True, text=True)
    if cp.returncode != 0:
        tail = "\n".join(cp.stdout.strip().splitlines()[-5:])
        return RowResult(name, "FAIL", tail or f"pytest exit {cp.returncode}")
    return RowResult(name, "PASS")


def _grep_row(
    repo: Path,
    name: str,
    pattern: str,
    targets: list[str],
    exclude_files: frozenset[str] = frozenset(),
) -> RowResult:
    """Fail on any regex hit in the named files/directories, skipping absent ones.

    ``exclude_files`` names files (repo-relative posix paths) that are skipped even
    when they live under a matched directory, for files whose whole documented role
    is to hold the very words being screened for (e.g. ``render/text.py``'s
    ``FORBIDDEN_PATTERNS`` constant, per implementation-plan section 2.1).
    """
    regex = re.compile(pattern)
    existing: list[Path] = []
    for target in targets:
        p = repo / target
        if p.exists():
            existing.append(p)
    if not existing:
        return RowResult(name, "SKIP", f"missing: {', '.join(targets)}")
    hits: list[str] = []
    for p in existing:
        files = [p] if p.is_file() else sorted(p.rglob("*.py"))
        files = [f for f in files if f.relative_to(repo).as_posix() not in exclude_files]
        for f in files:
            text = f.read_text(encoding="utf-8", errors="replace")
            for lineno, line in enumerate(text.splitlines(), start=1):
                if regex.search(line):
                    hits.append(f"{f.relative_to(repo).as_posix()}:{lineno}")
    if hits:
        return RowResult(name, "FAIL", "; ".join(hits[:5]))
    return RowResult(name, "PASS")


# --- individual checklist rows -------------------------------------------------

def row_docs_frozen(repo: Path, base: str, head: str, **_: object) -> RowResult:
    paths = [
        "docs/PHASE0_SPEC.md",
        "docs/AMENDMENTS.md",
        "docs/PILOT_FINDINGS.md",
        "docs/PREREGISTRATION_STRATUM_LIST.md",
        "docs/ENV_NOTES.md",
    ]
    return _diff_row(repo, base, head, paths, "docs/** frozen")


def row_pilot_subtrees_frozen(repo: Path, base: str, head: str, **_: object) -> RowResult:
    parts = [_diff_row(
        repo, base, head,
        ["pilot/preregistration", "pilot/gate", "pilot/testsel"],
        "pilot/preregistration, pilot/gate, pilot/testsel frozen",
    )]
    if sys.platform.startswith("win"):
        cp = _run_git(repo, ["ls-files", "--eol", "pilot/gate", "pilot/testsel"])
        if cp.returncode != 0:
            parts.append(RowResult(
                "pilot/preregistration, pilot/gate, pilot/testsel frozen",
                "FAIL",
                f"git ls-files --eol failed: {cp.stderr.strip()}",
            ))
        else:
            bad = [
                line for line in cp.stdout.splitlines()
                if line.strip() and not line.split()[0] == "i/lf"
            ]
            if bad:
                parts.append(RowResult(
                    "pilot/preregistration, pilot/gate, pilot/testsel frozen",
                    "FAIL",
                    f"not i/lf: {'; '.join(bad[:5])}",
                ))
    return _combine("pilot/preregistration, pilot/gate, pilot/testsel frozen", parts)


def row_pilot_readme(repo: Path, base: str, head: str, **_: object) -> RowResult:
    name = "pilot/README.md"
    parts = []
    rows = _numstat(repo, base, head, "pilot/README.md")
    deleted = sum(d for _, d, _ in rows if d >= 0)
    if deleted > 1:
        parts.append(RowResult(name, "FAIL", f"{deleted} deleted lines in pilot/README.md"))
    else:
        parts.append(RowResult(name, "PASS"))
    parts.append(_pytest_row(repo, name, "tests/test_docs.py", ["-k", "pilot_readme"]))
    return _combine(name, parts)


def row_gate_config_frozen(repo: Path, base: str, head: str, **_: object) -> RowResult:
    name = "gate_config.PRE_REGISTERED and the file"
    parts = [
        _diff_row(repo, base, head, ["src/pbcheck/gate_config.py"], name),
        _pytest_row(repo, name, "tests/test_gate_config.py"),
    ]
    return _combine(name, parts)


def row_stratum_list_freeze(repo: Path, base: str, head: str, **_: object) -> RowResult:
    return _pytest_row(repo, "Stratum list freeze", "tests/test_stratum_list_freeze.py")


def row_engine_modules(repo: Path, base: str, head: str, **_: object) -> RowResult:
    paths = [
        "src/pbcheck/design.py",
        "src/pbcheck/gene_universe.py",
        "src/pbcheck/permutation.py",
        "src/pbcheck/mtc.py",
        "src/pbcheck/io_counts.py",
        "src/pbcheck/census_select.py",
        "src/pbcheck/methods",
        "synthetic",
        "scripts/synthetic_gate.py",
    ]
    return _diff_row(repo, base, head, paths, "Engine modules")


def row_metrics_docstring(repo: Path, base: str, head: str, **_: object) -> RowResult:
    name = "metrics.py docstring-only exception"
    parts = []
    script = repo / "scripts" / "check_docstring_only_diff.py"
    if not script.exists():
        parts.append(RowResult(name, "SKIP", "missing: scripts/check_docstring_only_diff.py"))
    else:
        cp = subprocess.run(
            [sys.executable, str(script), base, head, "src/pbcheck/metrics.py"],
            cwd=repo, capture_output=True, text=True,
        )
        if cp.returncode != 0:
            tail = "\n".join(cp.stdout.strip().splitlines()[-5:] or [cp.stderr.strip()])
            parts.append(RowResult(name, "FAIL", tail))
        else:
            parts.append(RowResult(name, "PASS"))
    parts.append(_pytest_row(repo, name, "tests/test_metrics.py"))
    return _combine(name, parts)


_EXISTING_TESTS_EXCLUDE = [
    ":!tests/test_audit*.py",
    ":!tests/test_render.py",
    ":!tests/test_cli.py",
    ":!tests/test_example.py",
    ":!tests/test_packaging.py",
    ":!tests/test_docs.py",
    ":!tests/test_demo_scripts.py",
    ":!tests/fixtures",
    ":!tests/conftest.py",
    # new files shipped in this release, not part of the frozen test suite
    ":!tests/test_render_text.py",
    ":!tests/test_checklist_scripts.py",
    ":!tests/test_protocol_safety_check.py",
    ":!tests/test_measure_audit_runtime.py",
]


def row_existing_tests_unchanged(repo: Path, base: str, head: str, **_: object) -> RowResult:
    name = "Existing tests unchanged"
    parts = []
    if _diff_quiet(repo, base, head, ["tests/", *_EXISTING_TESTS_EXCLUDE]):
        parts.append(RowResult(name, "PASS"))
    else:
        parts.append(RowResult(name, "FAIL", "diff outside the excluded tests/ paths"))
    conf_rows = _numstat(repo, base, head, "tests/conftest.py")
    conf_deleted = sum(d for _, d, _ in conf_rows if d >= 0)
    if conf_deleted > 0:
        parts.append(RowResult(name, "FAIL", f"{conf_deleted} deleted lines in tests/conftest.py"))
    else:
        parts.append(RowResult(name, "PASS"))
    return _combine(name, parts)


_PRODUCT_TARGETS = [
    "src/pbcheck/audit.py",
    "src/pbcheck/cli.py",
    "src/pbcheck/render",
    "src/pbcheck/example.py",
]


def row_product_constants_not_protocol(repo: Path, base: str, head: str, **_: object) -> RowResult:
    return _grep_row(
        repo,
        "Product constants are not protocol constants",
        r"gate_config\.N_PERM",
        _PRODUCT_TARGETS,
    )


_FORBIDDEN_PROTOCOL_LANGUAGE = (
    r"Amendment 5|\bTier\b|reportab|jaccard|concordance|ratio_naive|"
    r"INSTRUMENT VALID|NO-GO|risk_score"
)


def row_no_protocol_language(repo: Path, base: str, head: str, **_: object) -> RowResult:
    name = "No protocol language in product output"
    parts = [
        _grep_row(
            repo, name, _FORBIDDEN_PROTOCOL_LANGUAGE, _PRODUCT_TARGETS,
            # render/text.py's FORBIDDEN_PATTERNS constant is the canonical list of
            # these very words (plan section 2.1); it is not rendered product output.
            exclude_files=frozenset({"src/pbcheck/render/text.py"}),
        ),
        _pytest_row(repo, name, "tests/test_render.py", ["-k", "forbidden"]),
        _pytest_row(repo, name, "tests/test_cli.py", ["-k", "forbidden"]),
    ]
    return _combine(name, parts)


def row_no_census(repo: Path, base: str, head: str, **_: object) -> RowResult:
    demo_scripts = sorted((repo / "scripts").glob("demo_*.py"))
    targets = _PRODUCT_TARGETS + [p.relative_to(repo).as_posix() for p in demo_scripts]
    return _grep_row(
        repo,
        "No Census in the product path",
        r"cellxgene_census|census_select",
        targets,
    )


def row_demo_datasets_outside_freeze(repo: Path, base: str, head: str, **_: object) -> RowResult:
    return _pytest_row(
        repo, "Demo datasets outside the freeze", "tests/test_demo_scripts.py", ["-k", "membership"]
    )


def row_gate_numbers(repo: Path, base: str, head: str, scratch: Path, with_gate: bool, **_: object) -> RowResult:
    name = "Gate numbers do not move"
    if not with_gate:
        return RowResult(name, "SKIP", "requires --with-gate")
    compare = repo / "scripts" / "compare_gate_scalars.py"
    if not compare.exists():
        return RowResult(name, "SKIP", "missing: scripts/compare_gate_scalars.py")
    baseline = repo / "pilot" / "gate" / "synthetic_gate_2026-08-15.json"
    if not baseline.exists():
        return RowResult(name, "SKIP", f"missing: {baseline.relative_to(repo).as_posix()}")
    scratch.mkdir(parents=True, exist_ok=True)
    out_path = scratch / "gate_release.json"
    cp = subprocess.run(
        [sys.executable, "scripts/synthetic_gate.py", "--out", str(out_path)],
        cwd=repo, capture_output=True, text=True,
    )
    if cp.returncode != 0:
        tail = "\n".join(cp.stdout.strip().splitlines()[-5:] or [cp.stderr.strip()])
        return RowResult(name, "FAIL", f"synthetic_gate.py: {tail}")
    cp2 = subprocess.run(
        [sys.executable, str(compare), str(out_path), str(baseline)],
        cwd=repo, capture_output=True, text=True,
    )
    if cp2.returncode != 0:
        tail = "\n".join(cp2.stdout.strip().splitlines()[-5:] or [cp2.stderr.strip()])
        return RowResult(name, "FAIL", f"compare_gate_scalars.py: {tail}")
    return RowResult(name, "PASS")


def row_readme_no_demo_number(repo: Path, base: str, head: str, **_: object) -> RowResult:
    return _pytest_row(
        repo, "README carries no demo number", "tests/test_docs.py", ["-k", "demo_readouts"]
    )


def row_frozen_tooling(repo: Path, base: str, head: str, **_: object) -> RowResult:
    return _diff_row(
        repo, base, head,
        [".pre-commit-config.yaml", ".gitattributes"],
        "Frozen files excluded from tooling",
    )


def row_amendments_untouched(repo: Path, base: str, head: str, **_: object) -> RowResult:
    name = "docs/AMENDMENTS.md untouched"
    parts = [_diff_row(repo, base, head, ["docs/AMENDMENTS.md"], name)]
    if not _log_empty(repo, base, head, "docs/AMENDMENTS.md"):
        parts.append(RowResult(name, "FAIL", "git log shows commits touching docs/AMENDMENTS.md"))
    else:
        parts.append(RowResult(name, "PASS"))
    return _combine(name, parts)


_ROWS = [
    row_docs_frozen,
    row_pilot_subtrees_frozen,
    row_pilot_readme,
    row_gate_config_frozen,
    row_stratum_list_freeze,
    row_engine_modules,
    row_metrics_docstring,
    row_existing_tests_unchanged,
    row_product_constants_not_protocol,
    row_no_protocol_language,
    row_no_census,
    row_demo_datasets_outside_freeze,
    row_gate_numbers,
    row_readme_no_demo_number,
    row_frozen_tooling,
    row_amendments_untouched,
]


def run_checklist(
    repo: Path, base: str, head: str, scratch: Path, with_gate: bool
) -> tuple[list[RowResult], bool]:
    results = [row(repo=repo, base=base, head=head, scratch=scratch, with_gate=with_gate) for row in _ROWS]
    ok = not any(r.status == "FAIL" for r in results)
    return results, ok


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default="main", help="base git ref (default: main)")
    parser.add_argument("--head", default="HEAD", help="head git ref (default: HEAD)")
    parser.add_argument(
        "--repo",
        default=None,
        help="repository root (default: the repository containing this script)",
    )
    parser.add_argument(
        "--with-gate",
        action="store_true",
        help="also run the gate-reproduction row (off by default)",
    )
    parser.add_argument(
        "--scratch",
        default=None,
        help="gate output directory, created if absent (default: a temp dir)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    repo = Path(args.repo).resolve() if args.repo else Path(__file__).resolve().parent.parent
    scratch = Path(args.scratch).resolve() if args.scratch else Path(tempfile.mkdtemp(prefix="pbcheck_gate_"))

    results, ok = run_checklist(repo, args.base, args.head, scratch, args.with_gate)
    for r in results:
        print(r.line())

    passed = sum(1 for r in results if r.status == "PASS")
    failed = sum(1 for r in results if r.status == "FAIL")
    skipped = sum(1 for r in results if r.status == "SKIP")
    print(f"summary: {passed} passed, {failed} failed, {skipped} skipped, {len(results)} total")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
