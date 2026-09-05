"""pbcheck's report renderer.

``render_markdown`` and ``render_html`` build the Markdown and HTML reports from a validated
``pbcheck-audit/1`` payload (``pbcheck.audit_schema.validate``), via the shared section model in
:mod:`pbcheck.render.sections`. ``write_outputs`` writes the three files a pbcheck run produces
(the raw JSON payload plus both reports) to a directory.
"""

from __future__ import annotations

import json
from pathlib import Path

from pbcheck.render.html import render_html
from pbcheck.render.markdown import render_markdown

__all__ = ["render_markdown", "render_html", "write_outputs"]

_OUTPUT_NAMES: dict[str, str] = {
    "json": "pbcheck_audit.json",
    "md": "pbcheck_report.md",
    "html": "pbcheck_report.html",
}


def _render(fmt: str, payload: dict) -> str:
    if fmt == "json":
        return json.dumps(payload, indent=2, sort_keys=False, ensure_ascii=False)
    if fmt == "md":
        return render_markdown(payload)
    if fmt == "html":
        return render_html(payload)
    raise ValueError(f"unknown output format: {fmt!r}")


def write_outputs(
    payload: dict,
    out_dir: str | Path,
    *,
    formats: tuple[str, ...] = ("json", "md", "html"),
    overwrite: bool = False,
) -> dict[str, Path]:
    """Write ``formats`` of ``payload`` under ``out_dir``, creating it if needed.

    ``formats`` selects among ``"json"`` (``pbcheck_audit.json``, ``indent=2``,
    ``sort_keys=False``, ``ensure_ascii=False``), ``"md"`` (``pbcheck_report.md``) and ``"html"``
    (``pbcheck_report.html``). Every file is written with ``encoding="utf-8"`` and
    ``newline="\\n"``. Unless ``overwrite`` is true, an existing target file raises
    ``FileExistsError`` naming the file rather than being silently replaced.

    Returns a mapping from each requested format to the ``Path`` written.
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    targets: dict[str, Path] = {}
    for fmt in formats:
        name = _OUTPUT_NAMES.get(fmt)
        if name is None:
            raise ValueError(f"unknown output format: {fmt!r}")
        target = out_path / name
        if target.exists() and not overwrite:
            raise FileExistsError(f"refusing to overwrite existing file: {target}")
        targets[fmt] = target

    written: dict[str, Path] = {}
    for fmt, target in targets.items():
        content = _render(fmt, payload)
        with open(target, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        written[fmt] = target
    return written
