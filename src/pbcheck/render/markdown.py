"""The Markdown renderer: walks :func:`pbcheck.render.sections.build_sections` and draws every
block as GitHub-flavoured Markdown. Tables are always pipe tables (``| a | b |`` with a
``| --- | --- |`` separator row), including the two-column tables :class:`~pbcheck.render.sections.KeyValues`
blocks are drawn as.
"""

from __future__ import annotations

from pbcheck.render.sections import Block, Callout, KeyValues, Paragraph, Table, build_sections


def _escape_cell(value: object) -> str:
    """A table cell as Markdown text: pipes and newlines would otherwise break the row."""
    cell = "" if value is None else str(value)
    return cell.replace("|", "\\|").replace("\n", " ")


def _render_table(headers: tuple[str, ...], rows: tuple[tuple[object, ...], ...]) -> str:
    lines = [
        "| " + " | ".join(_escape_cell(header) for header in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_escape_cell(cell) for cell in row) + " |")
    return "\n".join(lines)


def _render_block(block: Block) -> str:
    if isinstance(block, Paragraph):
        return block.text
    if isinstance(block, Callout):
        return "> " + block.text
    if isinstance(block, Table):
        parts = [f"**{block.caption}**"] if block.caption else []
        parts.append(_render_table(block.headers, block.rows))
        return "\n\n".join(parts)
    if isinstance(block, KeyValues):
        parts = [f"**{block.caption}**"] if block.caption else []
        rows = tuple((key, value) for key, value in block.items)
        parts.append(_render_table(("Field", "Value"), rows))
        return "\n\n".join(parts)
    raise TypeError(f"unknown block type: {type(block).__name__}")


def render_markdown(payload: dict) -> str:
    """The full Markdown report for a ``pbcheck-audit/1`` payload, as a single string."""
    parts = ["# pbcheck audit report"]
    for section in build_sections(payload):
        parts.append(f"## {section.id}. {section.title}")
        parts.extend(_render_block(block) for block in section.blocks)
    return "\n\n".join(parts) + "\n"
