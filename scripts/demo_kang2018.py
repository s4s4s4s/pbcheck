"""Demonstration 1: the design gate on Kang et al. 2018 (WP6, part 1).

Runs ``pbcheck.audit.run_audit`` (the same code path the ``pbcheck`` command line uses, not a
work-around of it) on the largest cell type of the Kang et al. 2018 IFN-beta stimulation dataset,
and writes the three report files a ``pbcheck audit`` run produces to
``demo/kang2018_design_gate/``.

This is a demonstration outside the pre-registered Phase 0 protocol: Kang et al. 2018 is not among
the 17 frozen datasets or 357 strata of ``pilot/preregistration/stratum_list_2026-08-16.csv``
(``tests/test_demo_scripts.py -k membership`` checks this by machine), and nothing this script
prints or writes is a Phase 0 measurement. Kang et al. 2018 pairs each of 8 donors under both
conditions (control and IFN-beta stimulation): the design gate refuses a paired design before
either DE arm runs, and that refusal, not a DE result, is what this demonstration shows.

Data source: pertpy's ``kang_2018()`` loader downloads the same file from the figshare API
(article file id 34464122); the id was confirmed on 2026-09-06 by fetching
``https://raw.githubusercontent.com/scverse/pertpy/main/pertpy/data/_datasets.py``. That fetch
landed on a redirected path (the file has since moved to
``src/pertpy/data/_datasets.py`` in pertpy's repository) and the current source there downloads
from ``https://exampledata.scverse.org/pertpy/kang_2018.h5ad`` instead of figshare, meaning
pertpy has migrated its own default source off figshare since this file id was pinned. The
figshare API endpoint (``https://api.figshare.com/v2/file/download/<id>``) still serves the
identical file (verified by sha256 against the value below) and is kept here as the pinned
source, per the task's instruction to keep the id when it cannot be freshly confirmed against the
current default. Do not use ``https://figshare.com/ndownloader/files/34464122``: on 2026-09-06 it
answered HTTP 202 with an empty body to both HEAD and a ranged GET.

Licence: not confirmed by this script. The candidate figshare article ids probed via
``https://api.figshare.com/v2/articles/<id>`` search on 2026-09-06 did not resolve to the item
that serves this file id; see the figshare item page for the file's own licence terms. Raw data
is not redistributed by this repository (``data/`` is git-ignored); only the derived report files
under ``demo/kang2018_design_gate/`` are committed.

Citation: Kang, H., Subramaniam, M., Targ, S. et al. Multiplexed droplet single-cell
RNA-sequencing using natural genetic variation. Nat Biotechnol 36, 89-94 (2018).
doi:10.1038/nbt.4042
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from demo_common import DEMO_DATASETS, download, load_and_check, pick_celltype  # noqa: E402

from pbcheck.audit import AuditSettings, run_audit  # noqa: E402
from pbcheck.render import write_outputs  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent

#: The figshare API "file download" endpoint for pertpy's kang_2018 artifact (file id 34464122).
#: Follows a redirect to a presigned S3 URL internally (``urllib.request`` does this by default);
#: see the module docstring for why this URL, not ``figshare.com/ndownloader/...``, is used.
KANG_2018_URL = "https://api.figshare.com/v2/file/download/34464122"

#: Pinned after the first verified run of this script (2026-09-06); a mismatch aborts the download.
KANG_2018_SHA256 = "e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830"

#: The file's size in bytes, recorded alongside the pinned digest as a second, cheap check a
#: reader can compare without recomputing a hash.
KANG_2018_SIZE_BYTES = 38356412

KANG_2018_PATH = REPO_ROOT / "data" / "kang_2018.h5ad"
OUT_DIR = REPO_ROOT / "demo" / "kang2018_design_gate"

#: Column names as pertpy's ``kang_2018()`` artifact defines them (verified against the file on
#: 2026-09-06): ``label`` holds ``"ctrl"``/``"stim"``, ``replicate`` holds one value per donor,
#: ``cell_type`` holds the 8 cell-type levels. No adaptation was needed: these match the plan's
#: assumed names exactly.
DONOR_COL = "replicate"
CONDITION_COL = "label"
CELLTYPE_COL = "cell_type"
TEST_LEVEL = "stim"
REF_LEVEL = "ctrl"


def main() -> None:
    dataset = DEMO_DATASETS["kang_2018"]
    print(f"Downloading/verifying {dataset.display_name} ({KANG_2018_URL}) ...")
    path = download(KANG_2018_URL, KANG_2018_PATH, KANG_2018_SHA256)
    size = path.stat().st_size
    if size != KANG_2018_SIZE_BYTES:
        raise AssertionError(
            f"{path}: size {size} bytes does not match the pinned {KANG_2018_SIZE_BYTES} bytes"
        )
    print(f"{path}: {size} bytes, sha256 {KANG_2018_SHA256}")

    adata = load_and_check(path, (DONOR_COL, CONDITION_COL, CELLTYPE_COL))
    celltype_value = pick_celltype(adata, CELLTYPE_COL)
    print(f"Largest cell type: {celltype_value!r}")

    settings = AuditSettings(
        donor_col=DONOR_COL,
        condition_col=CONDITION_COL,
        test_level=TEST_LEVEL,
        ref_level=REF_LEVEL,
        celltype_col=CELLTYPE_COL,
        celltype_value=celltype_value,
        design_only=False,
    )

    started = time.perf_counter()
    payload = run_audit(adata, settings)
    elapsed = time.perf_counter() - started
    print(f"run_audit finished in {elapsed:.1f}s: status={payload['status']!r} "
          f"status_reason={payload['status_reason']!r}")

    written = write_outputs(payload, OUT_DIR, overwrite=True)
    for fmt, target in written.items():
        print(f"wrote {fmt}: {target}")


if __name__ == "__main__":
    main()
