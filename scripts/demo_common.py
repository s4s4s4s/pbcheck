"""Shared helpers for the WP6 demonstration scripts (``scripts/demo_*.py``).

Both demonstrations run outside the pre-registered Phase 0 protocol, on datasets that are not
among the 17 frozen datasets or 357 strata of ``pilot/preregistration/stratum_list_2026-08-16.csv``
(machine-checked by ``tests/test_demo_scripts.py -k membership``). This module owns the pieces both
scripts share: a small, dependency-free downloader with checksum verification, an ``.obs`` column
assertion that fails loudly with the available columns, a "largest cell type" picker, and the
dataset registry the membership check reads.

Nothing here imports the pbcheck engine outside ``pbcheck.io_counts``: the demo scripts, not this
module, call ``pbcheck.audit.run_audit`` and ``pbcheck.render.write_outputs``, the same code path
the ``pbcheck`` command line uses.
"""

from __future__ import annotations

import hashlib
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

#: How many times ``download`` retries a transient figshare response (HTTP 202 or an empty body)
#: before giving up. The figshare API redirects a file-download request to a presigned S3 URL that
#: is not always ready on the first hit; retrying a few times with a short wait clears this without
#: masking a genuine outage (5 tries, 5 s apart, is a demo-script budget, not a protocol constant).
DOWNLOAD_MAX_RETRIES = 5

#: Seconds to wait between retries of a transient figshare response.
DOWNLOAD_RETRY_WAIT_SECONDS = 5.0

#: Bytes read per chunk while streaming a download to disk.
DOWNLOAD_CHUNK_SIZE = 1 << 20


@dataclass(frozen=True)
class DemoDataset:
    """One entry of the demo dataset registry.

    ``search_tokens`` and ``known_identifiers`` are the strings the protocol-safety checklist's
    frozen-list membership test (``tests/test_demo_scripts.py -k membership``) looks for in
    ``pilot/preregistration/stratum_list_2026-08-16.csv`` and in ``docs/`` / ``pilot/`` text, case
    insensitively: the demonstration is only meaningful if neither dataset is a dataset the
    pre-registered protocol could ever measure.
    """

    display_name: str
    citation: str
    search_tokens: tuple[str, ...]
    known_identifiers: tuple[str, ...] = field(default_factory=tuple)


#: The two demonstration datasets (section 6 of the release plan). Every token here is what the
#: frozen-list membership check searches for; add a dataset only with its tokens filled in.
DEMO_DATASETS: dict[str, DemoDataset] = {
    "kang_2018": DemoDataset(
        display_name="Kang et al. 2018",
        citation=(
            "Kang, H., Subramaniam, M., Targ, S. et al. Multiplexed droplet single-cell "
            "RNA-sequencing using natural genetic variation. Nat Biotechnol 36, 89-94 (2018). "
            "doi:10.1038/nbt.4042"
        ),
        search_tokens=("Kang", "GSE96583"),
        known_identifiers=(),
    ),
    "stephenson_2021": DemoDataset(
        display_name="Stephenson et al. 2021",
        citation=(
            "Stephenson, E., Reynolds, G., Botting, R.A. et al. Single-cell multi-omics analysis "
            "of the immune response in COVID-19. Nat Med 27, 904-916 (2021). "
            "doi:10.1038/s41591-021-01329-2"
        ),
        search_tokens=("Stephenson",),
        known_identifiers=(),
    ),
}


def _is_transient_figshare_response(status: int, body_length: int) -> bool:
    """Whether a figshare API response should be retried rather than treated as a failure.

    HTTP 202 with an empty body is the observed transient shape (the API has accepted the request
    for a presigned S3 redirect that is not ready yet); an empty body on any other 2xx status is
    treated the same way, since a 0-byte "download" is never a valid h5ad file.
    """
    return status == 202 or body_length == 0


def download(url: str, dest: str | Path, sha256: str) -> Path:
    """Download ``url`` to ``dest``, verifying its sha256, skipping a download that already matches.

    Streams to a temporary file in ``dest``'s directory, then renames it onto ``dest`` only after
    the digest matches, so a crashed or interrupted download never leaves a half-written file at
    the final path. A digest mismatch raises ``ValueError`` and deletes the temporary file; the
    file at ``dest`` (if any) is left untouched.

    ``url`` is expected to be a figshare API "file download" endpoint, which replies with a
    redirect to a presigned S3 URL that expires in roughly 10 seconds. ``urllib.request`` follows
    the redirect inside the same request by default, so the presigned URL itself is never stored
    or reused; a transient HTTP 202 or empty-body response (the API not yet ready to redirect) is
    retried up to :data:`DOWNLOAD_MAX_RETRIES` times, :data:`DOWNLOAD_RETRY_WAIT_SECONDS` apart,
    before ``download`` aborts naming the URL and the last status seen.
    """
    dest_path = Path(dest)
    expected = sha256.lower()

    if dest_path.exists() and _sha256_of(dest_path) == expected:
        return dest_path

    dest_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = dest_path.with_name(dest_path.name + ".part")

    last_status: int | None = None
    last_body_length: int | None = None
    for attempt in range(1, DOWNLOAD_MAX_RETRIES + 1):
        digest = hashlib.sha256()
        body_length = 0
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                status = getattr(response, "status", None) or 200
                with open(tmp_path, "wb") as handle:
                    while True:
                        chunk = response.read(DOWNLOAD_CHUNK_SIZE)
                        if not chunk:
                            break
                        digest.update(chunk)
                        body_length += len(chunk)
                        handle.write(chunk)
        except urllib.error.HTTPError as exc:
            status = exc.code
            body_length = 0

        last_status = status
        last_body_length = body_length

        if _is_transient_figshare_response(status, body_length):
            tmp_path.unlink(missing_ok=True)
            if attempt < DOWNLOAD_MAX_RETRIES:
                time.sleep(DOWNLOAD_RETRY_WAIT_SECONDS)
                continue
            raise RuntimeError(
                f"download of {url!r} kept returning a transient response "
                f"(last status {last_status}, body length {last_body_length}) after "
                f"{DOWNLOAD_MAX_RETRIES} attempts; aborting"
            )

        if status >= 400:
            tmp_path.unlink(missing_ok=True)
            raise RuntimeError(f"download of {url!r} failed with HTTP {status}")

        got = digest.hexdigest()
        if got != expected:
            tmp_path.unlink(missing_ok=True)
            raise ValueError(
                f"downloaded {url!r} but its sha256 {got} does not match the expected {expected}; "
                "aborting rather than trusting an unverified file"
            )

        tmp_path.replace(dest_path)
        return dest_path

    # Unreachable: the loop above always either returns, raises, or continues.
    raise RuntimeError(f"download of {url!r} did not complete (last status {last_status})")


def _sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(DOWNLOAD_CHUNK_SIZE)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def load_and_check(path: str | Path, columns: tuple[str, ...]):
    """Load an ``.h5ad`` file and assert it has every column of ``columns`` in ``.obs``.

    On a missing column, raises ``AssertionError`` naming the missing column, the full list of
    ``adata.obs.columns``, and, for every column that *is* present among ``columns``, its distinct
    level values, so a mismatch between a secondary source's assumed column names and the file's
    actual names is diagnosable from the error message alone rather than a notebook session.
    """
    import anndata as ad  # local import: neither pbcheck nor its tests otherwise need anndata here

    adata = ad.read_h5ad(str(path))
    missing = [c for c in columns if c not in adata.obs.columns]
    if missing:
        present_levels = {
            c: sorted(map(str, adata.obs[c].dropna().unique().tolist()))
            for c in columns
            if c in adata.obs.columns
        }
        raise AssertionError(
            f"{path}: obs is missing column(s) {missing!r}; obs.columns = "
            f"{list(adata.obs.columns)!r}; levels of the columns that are present: {present_levels!r}"
        )
    return adata


def pick_celltype(adata, col: str) -> str:
    """Return the most populous level of ``adata.obs[col]``.

    Raises ``AssertionError`` if ``col`` is not a column of ``adata.obs`` or has no non-null value,
    naming the available columns in the first case and the column in the second.
    """
    if col not in adata.obs.columns:
        raise AssertionError(
            f"obs has no column {col!r}; available columns: {list(adata.obs.columns)!r}"
        )
    counts = adata.obs[col].value_counts()
    if counts.empty:
        raise AssertionError(f"obs[{col!r}] has no non-null value to pick a largest level from")
    return str(counts.idxmax())
