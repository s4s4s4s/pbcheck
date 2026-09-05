"""pbcheck's report renderer.

``render_markdown`` builds the Markdown report from a validated ``pbcheck-audit/1`` payload
(``pbcheck.audit_schema.validate``), via the section model in :mod:`pbcheck.render.sections`.
``render_html`` and ``write_outputs`` are added by the next pass of this work package.
"""

from __future__ import annotations

from pbcheck.render.markdown import render_markdown

__all__ = ["render_markdown"]
