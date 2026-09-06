"""The HTML renderer: walks :func:`pbcheck.render.sections.build_sections` and draws the same
section model as :mod:`pbcheck.render.markdown`, as a single self-contained HTML document. Every
value drawn from the payload passes through :func:`html.escape`; the document carries one inline
``<style>`` block, no ``<script>``, and no external font or asset, so it opens and prints correctly
offline. Every table sits inside a container that scrolls horizontally rather than wrapping or
truncating, since column counts (donors, genes) are not bounded by this module.
"""

from __future__ import annotations

from html import escape

from pbcheck.render.sections import Block, Callout, KeyValues, Paragraph, Table, build_sections

_STYLE = """
body { font-family: Georgia, "Times New Roman", serif; max-width: 60rem; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; background: #fff; line-height: 1.5; }
h1 { font-size: 1.6rem; margin-bottom: 0.25rem; }
h2 { font-size: 1.2rem; margin-top: 2rem; border-bottom: 1px solid #ccc; padding-bottom: 0.2rem; }
p { margin: 0.75rem 0; }
.table-scroll { overflow-x: auto; margin: 0.75rem 0; }
table { border-collapse: collapse; width: 100%; min-width: max-content; }
caption { caption-side: top; text-align: left; font-weight: bold; margin-bottom: 0.25rem; }
th, td { border: 1px solid #ccc; padding: 0.35rem 0.6rem; text-align: left; white-space: nowrap; }
th { background: #f2f2f2; }
blockquote { border-left: 4px solid #888; margin: 0.75rem 0; padding: 0.25rem 1rem; color: #333; background: #f7f7f7; }
"""


def _escape(value: object) -> str:
    """A payload scalar as escaped HTML text: ``None`` renders as an empty string, everything else
    goes through ``str`` before ``html.escape``."""
    cell = "" if value is None else str(value)
    return escape(cell)


def _escape_prose(text: str) -> str:
    """A block of prose as escaped HTML, keeping the line breaks its author put in it.

    Some of the report's prose is a list rather than a run of sentences (the operating-envelope
    rows of note N2 are one row per line), and HTML would otherwise collapse those newlines into
    spaces and run the rows together.
    """
    return "<br>\n".join(_escape(line) for line in text.split("\n"))


def _render_table(headers: tuple[str, ...], rows: tuple[tuple[object, ...], ...], caption: str | None) -> str:
    parts = ['<div class="table-scroll"><table>']
    if caption:
        parts.append(f"<caption>{_escape(caption)}</caption>")
    parts.append("<thead><tr>" + "".join(f"<th>{_escape(header)}</th>" for header in headers) + "</tr></thead>")
    parts.append("<tbody>")
    for row in rows:
        parts.append("<tr>" + "".join(f"<td>{_escape(cell)}</td>" for cell in row) + "</tr>")
    parts.append("</tbody></table></div>")
    return "".join(parts)


def _render_block(block: Block) -> str:
    if isinstance(block, Paragraph):
        return f"<p>{_escape_prose(block.text)}</p>"
    if isinstance(block, Callout):
        return f"<blockquote>{_escape_prose(block.text)}</blockquote>"
    if isinstance(block, Table):
        return _render_table(block.headers, block.rows, block.caption)
    if isinstance(block, KeyValues):
        rows = tuple((key, value) for key, value in block.items)
        return _render_table(("Field", "Value"), rows, block.caption)
    raise TypeError(f"unknown block type: {type(block).__name__}")


def render_html(payload: dict) -> str:
    """The full HTML report for a ``pbcheck-audit/1`` payload, as a single self-contained string:
    one ``<html>`` document, one inline ``<style>``, no ``<script>``, no external asset."""
    parts = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        "<title>pbcheck audit report</title>",
        f"<style>{_STYLE}</style>",
        "</head>",
        "<body>",
        "<h1>pbcheck audit report</h1>",
    ]
    for section in build_sections(payload):
        parts.append(f"<h2>{_escape(section.id)}. {_escape(section.title)}</h2>")
        parts.extend(_render_block(block) for block in section.blocks)
    parts.append("</body></html>")
    return "\n".join(parts) + "\n"
