"""Guard against the version drifting across its independent copies.

``src/pbcheck/__init__.py`` (:data:`pbcheck.__version__`) is the single source of truth: this
script parses it directly from source rather than importing the package or asking installed
package metadata, so it works against a checkout that has never been installed. Four other copies
can still drift silently and are each compared against the source value:

* ``CITATION.cff`` (``version:``) - consumed by tools that read the file standalone (GitHub's
  "Cite this repository", Zenodo).
* ``.zenodo.json`` (``version``) - Zenodo's own metadata sidecar; without this key Zenodo falls
  back to the release tag name instead of the package version.
* installed package metadata (``importlib.metadata.version("pbcheck")``), when the package
  happens to be installed - this is one more source to catch, not a requirement to be installed.
* ``CHANGELOG.md`` (a ``## [<version>]`` heading).

Run it locally or from CI (`.github/workflows/tests.yml`, `.github/workflows/release.yml`).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from importlib import metadata
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

_SOURCE_VERSION_RE = re.compile(r'^__version__\s*=\s*"([^"]+)"\s*$', flags=re.MULTILINE)
_CFF_VERSION_RE = re.compile(r"^version:\s*(\S+)\s*$", flags=re.MULTILINE)
_CFF_DATE_RELEASED_RE = re.compile(r"^date-released:\s*(\S+)\s*$", flags=re.MULTILINE)
_ISO_DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def extract_source_version(init_text: str) -> str | None:
    """Return the ``__version__`` value from ``src/pbcheck/__init__.py``, or ``None`` if absent."""
    match = _SOURCE_VERSION_RE.search(init_text)
    return match.group(1) if match else None


def extract_cff_version(cff_text: str) -> str | None:
    """Return the ``version:`` value from a CITATION.cff document, or ``None`` if absent.

    The single parser used both by this script's own drift check and by
    ``tests/test_packaging.py``, so the two cannot silently diverge.
    """
    match = _CFF_VERSION_RE.search(cff_text)
    return match.group(1) if match else None


def extract_cff_date_released(cff_text: str) -> str | None:
    """Return the raw ``date-released:`` value from a CITATION.cff document, unvalidated."""
    match = _CFF_DATE_RELEASED_RE.search(cff_text)
    return match.group(1) if match else None


def extract_changelog_heading_date(changelog_text: str, version: str) -> str | None:
    """Return the raw date token after ``## [<version>] - ``, unvalidated, or ``None``."""
    pattern = re.compile(r"^##\s*\[" + re.escape(version) + r"\]\s*-\s*(\S+)", flags=re.MULTILINE)
    match = pattern.search(changelog_text)
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
        "--citation",
        type=Path,
        default=REPO / "CITATION.cff",
        help="path to CITATION.cff (default: CITATION.cff at the repository root)",
    )
    parser.add_argument(
        "--zenodo",
        type=Path,
        default=REPO / ".zenodo.json",
        help="path to .zenodo.json (default: .zenodo.json at the repository root)",
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="release tag (e.g. 'v0.1.0'); asserts tag[1:] == pbcheck.__version__",
    )
    parser.add_argument(
        "--require-date",
        action="store_true",
        help="also assert that CITATION.cff has a real ISO 'date-released' value equal to the "
        "CHANGELOG heading's date",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)

    init_path = REPO / "src" / "pbcheck" / "__init__.py"
    init_text = init_path.read_text(encoding="utf-8")
    pkg_version = extract_source_version(init_text)
    if pkg_version is None:
        print(
            f"check_version_consistency: no '__version__' assignment found in {init_path}",
            file=sys.stderr,
        )
        return 2

    if not args.citation.exists():
        print(
            f"check_version_consistency: citation file not found at {args.citation}",
            file=sys.stderr,
        )
        return 2
    cff_text = args.citation.read_text(encoding="utf-8")
    cff_version = extract_cff_version(cff_text)
    if cff_version is None:
        print(
            f"check_version_consistency: no 'version:' line found in {args.citation}",
            file=sys.stderr,
        )
        return 2

    if pkg_version != cff_version:
        print(
            f"version drift: pbcheck.__version__ (via src/pbcheck/__init__.py) = {pkg_version!r} "
            f"but {args.citation} version = {cff_version!r}",
            file=sys.stderr,
        )
        return 1

    if not args.zenodo.exists():
        print(
            f"check_version_consistency: zenodo metadata not found at {args.zenodo}",
            file=sys.stderr,
        )
        return 2
    zenodo_data = json.loads(args.zenodo.read_text(encoding="utf-8"))
    zenodo_version = zenodo_data.get("version")
    if zenodo_version is None:
        print(
            f"check_version_consistency: no 'version' key found in {args.zenodo}",
            file=sys.stderr,
        )
        return 2
    if zenodo_version != pkg_version:
        print(
            f"version drift: pbcheck.__version__ (via src/pbcheck/__init__.py) = {pkg_version!r} "
            f"but {args.zenodo} version = {zenodo_version!r}",
            file=sys.stderr,
        )
        return 1

    try:
        installed_version = metadata.version("pbcheck")
    except metadata.PackageNotFoundError:
        print(
            "check_version_consistency: pbcheck is not installed in this environment "
            "(pip install -e . first)",
            file=sys.stderr,
        )
        return 2
    if installed_version != pkg_version:
        print(
            f"version drift: pbcheck.__version__ (via src/pbcheck/__init__.py) = {pkg_version!r} "
            f"but the installed distribution metadata reports {installed_version!r} "
            "(reinstall with 'pip install -e .')",
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
        date_released = extract_cff_date_released(cff_text)
        if date_released is None or not _ISO_DATE_RE.match(date_released):
            print(
                f"check_version_consistency: {args.citation} has no valid ISO 'date-released' "
                f"value (got {date_released!r})",
                file=sys.stderr,
            )
            return 1
        heading_date = extract_changelog_heading_date(changelog_text, pkg_version)
        if heading_date is None or not _ISO_DATE_RE.match(heading_date):
            print(
                f"check_version_consistency: {args.changelog} heading for [{pkg_version}] has no "
                f"valid ISO date (got {heading_date!r})",
                file=sys.stderr,
            )
            return 1
        if heading_date != date_released:
            print(
                f"version drift: {args.changelog} heading date {heading_date!r} != "
                f"{args.citation} date-released {date_released!r}",
                file=sys.stderr,
            )
            return 1

    print(f"version consistent: {pkg_version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
