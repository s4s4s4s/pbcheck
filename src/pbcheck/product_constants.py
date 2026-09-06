"""Product constants shared by more than one module of the v0.1.0 release.

Every value here is a PRODUCT value: this tool's own, changeable by an ordinary engineering
change, and never a pre-registered threshold of ``pbcheck.gate_config``. The module exists so a
constant used both by the prose (:mod:`pbcheck.render.text`) and by the audit
(:mod:`pbcheck.audit`) has exactly one definition and neither module has to import the other
(``pbcheck.audit`` imports ``pbcheck.render.text``, so the reverse import would be a cycle).

It imports nothing from the package, so it can be imported from anywhere in it.
"""

from __future__ import annotations

#: Donors per group below which this tool stops treating a floor comparison as free of the
#: per-cell leak: floors are called coarse, cross-file comparison is warned against, and the
#: read-out paragraph drops its categorical clause.
#:
#: PRODUCT VALUE for display; it gates nothing in the engine. Its origin is correction A1 of
#: pbcheck's frozen protocol as narrowed by change 2 of the fifth amendment
#: (``docs/AMENDMENTS.md``), which admits floor-based quantities outside the operating envelope
#: only when every group has at least this many donors. The number is reused here as the point
#: where the tool stops warning; the protocol's own use of it is quoted, not re-decided, in the
#: caveat text that interpolates it.
FEW_DONORS_THRESHOLD = 8
