"""Guard against the version drifting across its two remaining independent copies.

``pyproject.toml`` no longer carries its own copy: R3 made it ``dynamic = ["version"]``, read by
hatchling out of ``src/pbcheck/__init__.py`` (:data:`pbcheck.__version__`), which is the single
source of truth. ``CITATION.cff`` intentionally keeps a separate ``version:`` field — CFF is
consumed by tools that read the file standalone (GitHub's "Cite this repository", Zenodo) and
should not require importing the package — so it can still drift silently. This script is the
one check that would catch that; run it locally or from CI (`.github/workflows/tests.yml`).
"""

from __future__ import annotations

import re
import sys
from importlib import metadata
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def main() -> int:
    pkg_version = metadata.version("pbcheck")

    cff_text = (REPO / "CITATION.cff").read_text(encoding="utf-8")
    m = re.search(r"^version:\s*(\S+)\s*$", cff_text, flags=re.MULTILINE)
    if not m:
        print("check_version_consistency: no 'version:' line found in CITATION.cff", file=sys.stderr)
        return 2
    cff_version = m.group(1)

    if pkg_version != cff_version:
        print(
            f"version drift: pbcheck.__version__ (via src/pbcheck/__init__.py) = {pkg_version!r} "
            f"but CITATION.cff version = {cff_version!r}",
            file=sys.stderr,
        )
        return 1

    print(f"version consistent: {pkg_version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
