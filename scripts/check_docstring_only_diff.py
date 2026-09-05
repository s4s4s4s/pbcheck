#!/usr/bin/env python
"""Check that a file changed only in its docstrings between two git revisions.

    python scripts/check_docstring_only_diff.py BASE HEAD PATH

The release checklist allows exactly one edit to ``src/pbcheck/metrics.py``: a sentence appended to
the docstring of ``signal_above_floor``. "Only a docstring changed" is not something a reviewer can
assert by eye on a diff, so it is decided here by comparing the two abstract syntax trees with
every docstring removed: if the stripped trees are identical, no statement, expression, signature,
decorator or default value moved.

Exit status:

    0   the two revisions differ only in docstrings (or not at all)
    1   the stripped trees differ: something other than a docstring changed
    2   wrong usage, a git error, or a revision whose version of the file does not parse

Standard library only, so it runs in the release environment without the project installed.
"""

from __future__ import annotations

import ast
import subprocess
import sys

USAGE = "usage: check_docstring_only_diff.py BASE HEAD PATH"

#: Nodes whose first statement may be a docstring.
_DOCSTRING_OWNERS = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)


def git_show(revision: str, path: str) -> str:
    """The file's text at ``revision``. Exits 2 with git's own message when git refuses."""
    target = f"{revision}:{path}"
    completed = subprocess.run(
        ["git", "show", target],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or f"git show {target} failed"
        sys.stderr.write(f"{message}\n")
        raise SystemExit(2)
    return completed.stdout


def strip_docstrings(tree: ast.AST) -> ast.AST:
    """Remove every docstring from ``tree`` in place and return it.

    A body left empty by the removal gets an explicit ``pass`` so the tree stays a well-formed
    program; both sides of the comparison are treated the same way, so this cannot hide a change.
    """
    for node in ast.walk(tree):
        if not isinstance(node, _DOCSTRING_OWNERS):
            continue
        body = node.body
        if (body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)):
            node.body = body[1:] or [ast.Pass()]
    return tree


def parse_stripped(source: str, label: str) -> str:
    """``ast.dump`` of ``source`` with every docstring removed. Exits 2 when it does not parse."""
    try:
        tree = ast.parse(source)
    except SyntaxError as exc:
        sys.stderr.write(f"{label} does not parse: {exc}\n")
        raise SystemExit(2) from exc
    return ast.dump(strip_docstrings(tree))


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 3:
        sys.stderr.write(USAGE + "\n")
        return 2
    base, head, path = args

    base_dump = parse_stripped(git_show(base, path), f"{base}:{path}")
    head_dump = parse_stripped(git_show(head, path), f"{head}:{path}")

    if base_dump == head_dump:
        print(f"{path}: docstring-only difference between {base} and {head}")
        return 0
    print(f"{path}: code changed between {base} and {head}, not only docstrings")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
