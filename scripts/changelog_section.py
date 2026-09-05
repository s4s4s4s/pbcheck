"""Print the CHANGELOG.md body for one released version.

Extracts the text between a ``## [<version>]`` heading and the next line starting with ``## ``
(or end of file). Used by `.github/workflows/release.yml` to fill the GitHub Release notes from
the CHANGELOG entry so the two never drift apart.

Usage: python scripts/changelog_section.py <version> [--changelog PATH]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def extract_section(changelog_text: str, version: str) -> str | None:
    heading_pattern = re.compile(r"^##\s*\[" + re.escape(version) + r"\].*$", flags=re.MULTILINE)
    match = heading_pattern.search(changelog_text)
    if match is None:
        return None

    start = match.end()
    next_heading = re.search(r"^##\s+", changelog_text[start:], flags=re.MULTILINE)
    end = start + next_heading.start() if next_heading else len(changelog_text)
    return changelog_text[start:end].strip("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", help="version to extract, without a leading 'v' (e.g. '0.1.0')")
    parser.add_argument(
        "--changelog",
        type=Path,
        default=REPO / "CHANGELOG.md",
        help="path to the changelog file (default: CHANGELOG.md at the repository root)",
    )
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    if not args.changelog.exists():
        print(f"changelog_section: changelog not found at {args.changelog}", file=sys.stderr)
        return 1

    changelog_text = args.changelog.read_text(encoding="utf-8")
    section = extract_section(changelog_text, args.version)
    if section is None:
        print(
            f"changelog_section: no '## [{args.version}]' heading found in {args.changelog}",
            file=sys.stderr,
        )
        return 1

    print(section)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
