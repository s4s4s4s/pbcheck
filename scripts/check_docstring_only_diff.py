#!/usr/bin/env python
"""Check that a file changed only in its docstrings between two git revisions.

    python scripts/check_docstring_only_diff.py BASE HEAD PATH

The release checklist allows exactly one edit to ``src/pbcheck/metrics.py``: a sentence appended to
the docstring of ``signal_above_floor``. "Only a docstring changed" is not something a reviewer can
assert by eye on a diff, so it is decided here by comparing the two files as token streams with
every docstring token removed: if the streams are identical, no statement, expression, signature,
decorator, default value or comment moved.

Tokens, not syntax trees: a comment is not a node of the abstract syntax tree, so an AST
comparison answers "no code changed" to a revision that rewrote every comment in the file. The
stream compared here keeps comments and drops only the string literal that is a docstring (and the
newline token that terminated it), so the printed sentence says exactly what was checked.

Exit status:

    0   the two revisions differ only in docstrings (or not at all)
    1   the streams differ: something other than a docstring changed
    2   wrong usage, a git error, or a revision whose version of the file does not tokenize

Standard library only, so it runs in the release environment without the project installed.
"""

from __future__ import annotations

import ast
import io
import subprocess
import sys
import token as token_module
import tokenize

USAGE = "usage: check_docstring_only_diff.py BASE HEAD PATH"

#: Nodes whose first statement may be a docstring.
_DOCSTRING_OWNERS = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)

#: Token types carrying no information about the program: the encoding marker, the non-logical
#: newlines that end blank and comment-only lines, and the end marker. Dropping them keeps the
#: comparison from failing on a blank line, while COMMENT itself is kept.
_IGNORED_TOKENS = frozenset({
    token_module.ENCODING,
    token_module.NL,
    token_module.ENDMARKER,
})


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


def docstring_spans(source: str, label: str) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    """The ``(start, end)`` positions of every docstring in ``source``, as ``(line, column)``.

    Located on the syntax tree, because "a string that is the first statement of a module, class or
    function" is a statement position, not something a token stream can tell apart from any other
    expression statement. A span rather than a single position, so that a docstring written as
    several adjacent string literals is removed whole. Exits 2 when the source does not parse.
    """
    try:
        tree = ast.parse(source)
    except SyntaxError as exc:
        sys.stderr.write(f"{label} does not parse: {exc}\n")
        raise SystemExit(2) from exc

    spans: list[tuple[tuple[int, int], tuple[int, int]]] = []
    for node in ast.walk(tree):
        if not isinstance(node, _DOCSTRING_OWNERS):
            continue
        body = node.body
        if (body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)):
            first = body[0].value
            spans.append((
                (first.lineno, first.col_offset),
                (first.end_lineno, first.end_col_offset),
            ))
    return spans


def token_stream(source: str, label: str) -> list[tuple[int, str]]:
    """``source`` as a list of ``(token type, text)`` with the docstrings removed.

    A removed docstring takes the NEWLINE that terminated its statement with it, so a file with a
    docstring and the same file without one produce the same stream; a removed statement, a
    renamed name, a changed default and a rewritten comment all survive. Exits 2 when the source
    does not tokenize.
    """
    spans = docstring_spans(source, label)
    stream: list[tuple[int, str]] = []
    drop_next_newline = False

    def inside_a_docstring(start: tuple[int, int]) -> bool:
        return any(first <= start <= last for first, last in spans)
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
    except (tokenize.TokenError, IndentationError, SyntaxError) as exc:
        sys.stderr.write(f"{label} does not tokenize: {exc}\n")
        raise SystemExit(2) from exc

    for item in tokens:
        if item.type in _IGNORED_TOKENS:
            continue
        if item.type == token_module.STRING and inside_a_docstring(item.start):
            drop_next_newline = True
            continue
        if drop_next_newline and item.type == token_module.NEWLINE:
            drop_next_newline = False
            continue
        drop_next_newline = False
        stream.append((item.type, item.string))
    return stream


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 3:
        sys.stderr.write(USAGE + "\n")
        return 2
    base, head, path = args

    base_stream = token_stream(git_show(base, path), f"{base}:{path}")
    head_stream = token_stream(git_show(head, path), f"{head}:{path}")

    if base_stream == head_stream:
        print(f"{path}: docstring-only difference between {base} and {head}")
        return 0
    print(f"{path}: code or comments changed between {base} and {head}, not only docstrings")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
