"""Guard against the version drifting across its independent copies.

``pyproject.toml`` no longer carries its own copy: R3 made it ``dynamic = ["version"]``, read by
hatchling out of ``src/pbcheck/__init__.py`` (:data:`pbcheck.__version__`), which is the single
source of truth. ``CITATION.cff`` intentionally keeps a separate ``version:`` field - CFF is
consumed by tools that read the file standalone (GitHub's "Cite this repository", Zenodo) and
should not require importing the package - so it can still drift silently. ``CHANGELOG.md`` keeps
a third independent copy, as a ``## [<version>]`` heading. This script is the one check that would
catch any of the three drifting apart; run it locally or from CI (`.github/workflows/tests.yml`,
`.github/workflows/release.yml`).
"""

from __future__ import annotations

import argparse
import re
import sys
from importlib import metadata
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

_CFF_VERSION_RE = re.compile(r"^version:\s*(\S+)\s*$", flags=re.MULTILINE)


def extract_cff_version(cff_text: str) -> str | None:
    """Return the ``version:`` value from a CITATION.cff document, or ``None`` if absent.

    The single parser used both by this script's own drift check and by
    ``tests/test_packaging.py``, so the two cannot silently diverge.
    """
    match = _CFF_VERSION_RE.search(cff_text)
    return match.group(1) if match else None


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--changelog",
        type=Path,
        default=REPO / "CHANGELOG.md",
        help="path to the changelog file to check for a '## [<version>]' heading "
        "(default: CHANGELOG.md at the repository root)",
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="release tag (e.g. 'v0.1.0'); asserts tag[1:] == pbcheck.__version__",
    )
    parser.add_argument(
        "--require-date",
        action="store_true",
        help="also assert that CITATION.cff has an ISO 'date-released' value",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)

    pkg_version = metadata.version("pbcheck")

    cff_text = (REPO / "CITATION.cff").read_text(encoding="utf-8")
    cff_version = extract_cff_version(cff_text)
    if cff_version is None:
        print("check_version_consistency: no 'version:' line found in CITATION.cff", file=sys.stderr)
        return 2

    if pkg_version != cff_version:
        print(
            f"version drift: pbcheck.__version__ (via src/pbcheck/__init__.py) = {pkg_version!r} "
            f"but CITATION.cff version = {cff_version!r}",
            file=sys.stderr,
        )
        return 1

    if not args.changelog.exists():
        print(
            f"check_version_consistency: changelog not found at {args.changelog}",
            file=sys.stderr,
        )
        return 2

    changelog_text = args.changelog.read_text(encoding="utf-8")
    heading_pattern = re.compile(
        r"^##\s*\[" + re.escape(pkg_version) + r"\]", flags=re.MULTILINE
    )
    if not heading_pattern.search(changelog_text):
        print(
            f"version drift: no '## [{pkg_version}]' heading found in {args.changelog}",
            file=sys.stderr,
        )
        return 1

    if args.tag is not None:
        if not args.tag.startswith("v"):
            print(
                f"check_version_consistency: --tag {args.tag!r} does not start with 'v'",
                file=sys.stderr,
            )
            return 2
        tag_version = args.tag[1:]
        if tag_version != pkg_version:
            print(
                f"version drift: tag {args.tag!r} implies version {tag_version!r} "
                f"but pbcheck.__version__ = {pkg_version!r}",
                file=sys.stderr,
            )
            return 1

    if args.require_date:
        date_match = re.search(
            r"^date-released:\s*(\d{4}-\d{2}-\d{2})\s*$", cff_text, flags=re.MULTILINE
        )
        if not date_match:
            print(
                "check_version_consistency: CITATION.cff has no ISO 'date-released' value",
                file=sys.stderr,
            )
            return 1

    print(f"version consistent: {pkg_version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
