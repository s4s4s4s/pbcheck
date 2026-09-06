"""``python -m pbcheck``: the same entry point as the ``pbcheck`` console script."""

from __future__ import annotations

import sys

from pbcheck.cli import main

if __name__ == "__main__":
    sys.exit(main())
