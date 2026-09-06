"""Demonstration 2: a two-arm audit on an unpaired dataset (WP6, part 2).

Runs ``pbcheck.audit.run_audit`` (the same code path the ``pbcheck`` command line uses) on the
largest cell type of a real, unpaired public dataset (independent donors per condition, no design
gate expected to stop the run), and writes the three report files a ``pbcheck audit`` run produces
to ``demo/<dataset>_two_arm/``.

Candidate order (fixed by the release plan, section 6): Stephenson et al. 2021 first; if its matrix
carries no raw counts anywhere, or its design fails a gate, or a size/time bound below is hit, this
script exits 3 naming the reason and the caller falls back to the synthetic ``example_reference``
dataset, which always succeeds.

This is a demonstration outside the pre-registered Phase 0 protocol: neither candidate dataset is
among the 17 frozen datasets or 357 strata of ``pilot/preregistration/stratum_list_2026-08-16.csv``
(``tests/test_demo_scripts.py -k membership`` checks this by machine), and nothing this script
prints or writes is a Phase 0 measurement.

Data source (``stephenson2021``): pertpy's ``stephenson_2021_subsampled()`` loader now downloads
from ``https://exampledata.scverse.org/pertpy/stephenson_2021_subsampled.h5ad`` (verified against
``pertpy/data/_datasets.py`` on 2026-09-06), not figshare; pertpy migrated its own default source
off figshare, the same migration ``demo_kang2018.py`` documents for the Kang 2018 artifact. The
figshare API endpoint for the same file (article file id 38171703) is recovered from pertpy's git
history, the last commit before that migration
(``https://github.com/scverse/pertpy/blob/c14dd9b1b505246e4b894581583168c13a89800f/pertpy/data/_datasets.py``,
commit before "Replace figshare & zenodo URLs with scverse URLs (#849)"), and is kept here as the
pinned source, the same reasoning ``demo_kang2018.py`` documents. Its ``Content-Length`` (probed via
a ranged request on 2026-09-06, before any download) is 700549709 bytes, about 0.65 GiB, well under
this script's 3 GB size bound.

Column names (verified against the downloaded file on 2026-09-06, since pertpy's loader itself
names no ``.obs`` columns): ``patient_id`` (119 donors), ``Status`` (``Covid``/``Healthy``/``LPS``;
this demo contrasts ``Covid`` versus ``Healthy``, dropping the 10 ``LPS`` stimulation-control rows
the way ``pbcheck.audit.prepare_stratum`` drops any condition value outside the two named levels),
``cell_type`` (Cell Ontology labels).

Result on 2026-09-06: ``adata.X`` is ``float32`` but **not** integer-valued (log-normalized, values
up to ~7.5 on ~62.5k cells, e.g. cell 0 gene 0 onward), the file has no named ``.layers`` and no
``.raw``, so no raw integer counts exist anywhere in it. ``check_integer_counts`` confirms this by
machine at runtime and this script exits 3 (``no_raw_counts_anywhere``) before ``run_audit`` is
ever called; the committed demonstration under ``demo/`` is therefore the ``example_reference``
fallback (see ``demo/README.md`` for the exact reason recorded from this run).

Licence: not confirmed by this script; the figshare item's own page states its terms. Raw data is
not redistributed by this repository (``data/`` is git-ignored); only the derived report files
under ``demo/<dataset>_two_arm/`` are committed, and only for whichever dataset actually produced
them (``example_reference`` here).

Citation: Stephenson, E., Reynolds, G., Botting, R.A. et al. Single-cell multi-omics analysis of
the immune response in COVID-19. Nat Med 27, 904-916 (2021). doi:10.1038/s41591-021-01329-2
"""

from __future__ import annotations

import argparse
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from demo_common import (  # noqa: E402
    DEMO_DATASETS,
    download,
    load_and_check,
    pick_celltype,
    probe_content_length,
)

from pbcheck.audit import AuditSettings, run_audit  # noqa: E402
from pbcheck.example import REFERENCE_SHAPE, example_adata  # noqa: E402
from pbcheck.io_counts import check_integer_counts  # noqa: E402
from pbcheck.render import write_outputs  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent

#: This script's own exit code for "stopped early, use the fallback dataset instead"; distinct from
#: the ``pbcheck`` CLI's exit codes (``cli.py``), which this script never sets, since it drives the
#: audit through the Python API (``run_audit`` + ``write_outputs``), not the CLI process.
EXIT_FALLBACK = 3

#: A download this size or larger is not attempted; the caller falls back to ``example_reference``
#: instead (release plan section 6). Chosen as a demo-script budget for a developer's own machine,
#: not a product constant.
MAX_DOWNLOAD_BYTES = 3 * 1024 * 1024 * 1024

#: An audit run this long or longer is abandoned in favor of the fallback (release plan section 6).
MAX_AUDIT_SECONDS = 90 * 60

#: The figshare API "file download" endpoint recovered from pertpy's git history for the
#: ``stephenson_2021_subsampled`` artifact (file id 38171703); see the module docstring.
STEPHENSON_2021_URL = "https://api.figshare.com/v2/file/download/38171703"

#: Pinned after the first verified fetch of this script (2026-09-06); a mismatch aborts the
#: download rather than trusting an unverified file.
STEPHENSON_2021_SHA256 = "2345f3c512552892bb8327f4779dc69972faccac3b00dcf90c0e9d9ae5ab738e"

STEPHENSON_2021_PATH = REPO_ROOT / "data" / "stephenson_2021.h5ad"

DONOR_COL = "patient_id"
CONDITION_COL = "Status"
CELLTYPE_COL = "cell_type"
TEST_LEVEL = "Covid"
REF_LEVEL = "Healthy"


def _exit_fallback(reason: str) -> None:
    print(f"demo_two_arm.py: stopping, reason={reason!r}; use --dataset example_reference instead",
          file=sys.stderr)
    raise SystemExit(EXIT_FALLBACK)


def _run_with_wall_time_bound(fn, *, budget_seconds: float):
    """Call ``fn()`` in a daemon thread, giving up after ``budget_seconds``.

    A daemon thread that is still running when the budget expires is abandoned (not killed:
    CPython gives no safe way to force-stop another thread); the process exits without waiting for
    it, which is what "abandoned" means for a demonstration script. On success within the budget,
    returns ``fn()``'s result.
    """
    result: dict = {}
    error: dict = {}

    def _target() -> None:
        try:
            result["value"] = fn()
        except BaseException as exc:  # re-raised on the caller's side below
            error["exc"] = exc

    thread = threading.Thread(target=_target, daemon=True)
    started = time.perf_counter()
    thread.start()
    thread.join(timeout=budget_seconds)
    if thread.is_alive():
        _exit_fallback(
            f"audit exceeded the {budget_seconds:.0f}s wall-time bound "
            f"(still running after {time.perf_counter() - started:.0f}s)"
        )
    if "exc" in error:
        raise error["exc"]
    return result["value"]


def _check_raw_counts_everywhere(adata) -> str | None:
    """Return the name of the first layer of ``adata`` (``"X"`` or a ``.layers`` key) that passes
    :func:`pbcheck.io_counts.check_integer_counts`, or ``None`` if none does.
    """
    x_check = check_integer_counts(adata.X)
    if x_check.passed:
        return None  # None here means "use X", AuditSettings.counts_layer's own default
    print(f"X is not raw integer counts: reason={x_check.reason!r}")
    # ``adata.layers.keys()`` yields a spurious ``None`` entry on an empty ``.layers`` group with
    # this installed anndata version (0.13.2), and ``adata.layers[None]`` then returns ``adata.X``
    # itself rather than raising; skip it so an empty ``.layers`` group is reported as "no layers",
    # not as a second, identical check of X under a fake name.
    for layer_name in adata.layers.keys():
        if layer_name is None:
            continue
        layer_check = check_integer_counts(adata.layers[layer_name])
        if layer_check.passed:
            return layer_name
        print(f"layer {layer_name!r} is not raw integer counts: reason={layer_check.reason!r}")
    return "__none_found__"


def _run_audit_and_write(adata, settings: AuditSettings, out_dir: Path) -> dict:
    started = time.perf_counter()
    payload = run_audit(adata, settings)
    elapsed = time.perf_counter() - started
    print(f"run_audit finished in {elapsed:.1f}s: status={payload['status']!r} "
          f"status_reason={payload['status_reason']!r}")
    if payload["status"] == "design_only":
        _exit_fallback(f"design gate stopped the run: {payload['status_reason']}")
    written = write_outputs(payload, out_dir, overwrite=True)
    for fmt, target in written.items():
        print(f"wrote {fmt}: {target}")
    return payload


def run_stephenson2021(out_dir: Path) -> None:
    dataset = DEMO_DATASETS["stephenson_2021"]
    print(f"Probing size of {dataset.display_name} ({STEPHENSON_2021_URL}) ...")
    size = probe_content_length(STEPHENSON_2021_URL)
    if size is None:
        _exit_fallback("could not determine the download size before fetching it")
    print(f"Content-Length: {size} bytes")
    if size >= MAX_DOWNLOAD_BYTES:
        _exit_fallback(f"download size {size} bytes is at or above the {MAX_DOWNLOAD_BYTES} byte bound")

    print(f"Downloading/verifying {dataset.display_name} ...")
    path = download(STEPHENSON_2021_URL, STEPHENSON_2021_PATH, STEPHENSON_2021_SHA256)
    print(f"{path}: {path.stat().st_size} bytes, sha256 {STEPHENSON_2021_SHA256}")

    adata = load_and_check(path, (DONOR_COL, CONDITION_COL, CELLTYPE_COL))

    counts_layer = _check_raw_counts_everywhere(adata)
    if counts_layer == "__none_found__":
        _exit_fallback("no_raw_counts_anywhere: X and every layer failed check_integer_counts")
    if counts_layer is not None:
        print(f"raw integer counts found in layer {counts_layer!r}")

    celltype_value = pick_celltype(adata, CELLTYPE_COL)
    print(f"Largest cell type: {celltype_value!r}")

    settings = AuditSettings(
        donor_col=DONOR_COL,
        condition_col=CONDITION_COL,
        test_level=TEST_LEVEL,
        ref_level=REF_LEVEL,
        celltype_col=CELLTYPE_COL,
        celltype_value=celltype_value,
        counts_layer=counts_layer,
        design_only=False,
    )
    _run_with_wall_time_bound(
        lambda: _run_audit_and_write(adata, settings, out_dir),
        budget_seconds=MAX_AUDIT_SECONDS,
    )


def run_example_reference(out_dir: Path) -> None:
    adata = example_adata(
        seed=0,
        n_genes=REFERENCE_SHAPE["n_genes"],
        n_donors_per_group=REFERENCE_SHAPE["n_donors_per_group"],
        n_cells_per_donor=REFERENCE_SHAPE["n_cells_per_donor"],
    )
    celltype_col = "cell_type"
    celltype_value = pick_celltype(adata, celltype_col)
    print(f"Largest cell type: {celltype_value!r} (synthetic, one level by construction)")

    settings = AuditSettings(
        donor_col="donor",
        condition_col="condition",
        test_level="disease",
        ref_level="ctrl",
        celltype_col=celltype_col,
        celltype_value=celltype_value,
        design_only=False,
    )
    _run_audit_and_write(adata, settings, out_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", required=True, choices=("stephenson2021", "example_reference"))
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    if args.dataset == "stephenson2021":
        run_stephenson2021(args.out)
    else:
        run_example_reference(args.out)


if __name__ == "__main__":
    main()
