"""A2 cell-count stratification: the tau-ladder feasibility bracket (correction A2, plan sec 6 T1).

**What this is.** A2 (Amendment 6, not yet written) proposes restricting the donor-permutation null
to assignments whose per-group cell-count split lands within ``t* = floor(tau * T)`` of the real
split, where ``T`` is the stratum's total cell count. Before that tolerance ``tau`` can be
pre-registered, someone has to know what it would cost: how many strata would become too thin to
build a restricted set at all (VOID under the sec 1.2 thin-set rule), across a ladder of candidate
tau values, on both the frozen 150-stratum analysis set and the wider 554-stratum manifest tier.
That is what this script computes and commits. **It is evidence, not a decision** — it does not pick
``A2_TAU`` and is not itself an amendment.

**What it is not.** It cannot compute the real assignment distribution, because the committed data
does not contain it — see "The central limitation" below. Every number this script emits is a
*bracket* between two deliberately extreme constructions, not a measurement of the true donor-level
cell-count distribution. Re-run against the same two pinned source files, it must reproduce the
committed artifact byte-for-byte (``--check``); run against anything else, it refuses.

**Sources, pinned by sha256 before parsing** (the same discipline as ``scripts/freeze_stratum_list.py``):

* ``pilot/preregistration/stratum_list_2026-08-16.json`` — the frozen sec 1 stratum list. Its
  ``rows`` (``role = analysis_set``, 251 strata) filtered to ``min(n_donors_A, n_donors_B) >= 8``
  give the **frozen** tier.
* ``pilot/preregistration/census_candidates_run31910799023_2026-08-15.json`` — the whole-Census
  candidate manifest the frozen list itself derives from. Its ``gate_status == "candidate"`` rows
  filtered the same way give the **manifest** tier. The frozen tier is a strict subset of the
  manifest tier (same source data, twelve datasets out of the manifest's 68 candidate-bearing ones);
  the two tiers are reported side by side deliberately, not because they are disjoint.

**The central limitation — read this before reading any number below.** A committed row carries, for
each of a stratum's two groups, only five summary statistics of that group's per-donor cell counts:
``n_donors``, ``n_cells`` (the group total), ``min``, ``median``, ``max``. It does **not** carry the
per-donor vector itself, and the true vector cannot be recovered from what is committed — many
distinct integer vectors share the same five statistics, and they do not all induce the same
donor-permutation assignment-total distribution. This script builds, from those five numbers alone,
**two** deliberately extreme per-donor vectors per group — one of minimum plausible dispersion, one
of maximum — concatenates each into a whole-stratum vector, and reports the interval the two extremes
bracket. Where the two disagree, that disagreement is the honest answer: the true figure is somewhere
between them and this script does not know where.

**The bracket construction, in full** (also carried in the emitted artifact's own header, so a reader
never has to come back to this file to audit a number in it):

* Two donors are always pinned exactly at the group's reported ``min`` and ``max`` — every real
  per-donor vector must contain a donor at each, so this part is not an approximation.
* *Low dispersion* (minimum plausible variance): the remaining ``n_donors - 2`` donors are set as
  close to one shared value as integer arithmetic allows — ``n_cells - min - max`` split into
  ``n_donors - 2`` parts differing by at most one cell. **This is where "clamping the vector to the
  recorded total perturbs the median" happens, exactly as the plan requires disclosed**: the shared
  value lands near ``(n_cells - min - max) / (n_donors - 2)``, which is *not* the reported median in
  general (only coincidentally, when the group's true distribution is itself close to uniform). No
  attempt is made to force the reported median onto this vector at all; forcing it and then
  re-clamping to hit the exact total is the perturbation being described, so this construction takes
  the already-clamped value directly instead of computing and then discarding an intermediate one.
* *High dispersion* (maximum plausible variance): one interior donor (two, for an even donor count)
  is pinned at the reported ``median`` — so the value is present in the vector — and the rest of the
  ``n_donors - 2`` interior donors are filled two-point (as many as possible at ``max``, the rest at
  ``min``) to hit the exact remaining total, the maximum-variance construction on a bounded support
  for a fixed sum. At most one "leftover" donor takes an intermediate integer value to close the
  final rounding gap; it is proven (and checked by ``tests/test_a2_feasibility.py``) to always land
  back inside ``[min, max]``, so this construction never needs to clip. Because the fill targets the
  exact total rather than balancing how many donors land above versus below the pinned median, the
  *realised* order-statistic median of the finished, sorted vector can still differ from the reported
  one — this is measured per stratum, per group, and carried in the artifact as
  ``median_drift_A`` / ``median_drift_B`` under each dispersion, rather than assumed away.
* Donor counts of 4 or fewer per group need no approximation at all: for ``n_donors in {1, 2, 3, 4}``
  the five summary statistics determine the group's per-donor multiset *exactly* (for ``n_donors
  <= 3`` the sorted vector literally **is** ``(min, median, max)``; for ``n_donors == 4`` the median's
  own definition means the two middle values sum to ``2 * median`` as an identity). The frozen and
  manifest tiers analysed here never see this case (every group holds >= 8 donors by the tier
  filter), but the construction functions handle it exactly and it is covered by
  ``tests/test_a2_feasibility.py``.

**Classification against the sec 1.2 thin-set rule differs by whether a stratum's design count is
enumerable**, and the two paths are reported honestly as different kinds of number, not blended into
one:

* ``permutation_count = C(n_donors_A + n_donors_B, n_donors_A) <= A2_ENUM_CAP`` (2e6): every design's
  assignment total is computed **exactly** via a subset-sum dynamic program (one cell per (count
  chosen, running sum) pair; no sampling). The window count ``m`` is exact, and is compared directly
  against ``A2_FULL_MIN`` (1000) / ``A2_COARSE_MIN`` (200) — full / coarse / VOID, matching sec 1.2's
  enumerated branch exactly.
* ``permutation_count > A2_ENUM_CAP``: full enumeration is not attempted (nor is it in the real
  machinery — sec 1.2's own sampled branch draws a fixed budget rather than enumerating). This script
  estimates the survival mass ``q`` from ``MC_DRAWS`` = 20 000 balanced draws (seeded — see below) and
  classifies against the *acceptance-rate* form of the same rule sec 1.2 states as equivalent to its
  sample-budget thresholds: full at ``q >= 1e-3``, coarse at ``2e-5 <= q < 1e-3``, VOID below. This is
  a materially smaller budget than the real machinery's own (``A2_SAMPLE_BUDGET`` = 1e6, escalating to
  1e7) — deliberately: this script exists to bracket VOID counts cheaply and reproducibly across a
  five-point tau ladder and two tiers (1108 stratum-dispersion pairs), not to reproduce the real
  runtime's budget. The Monte-Carlo standard error of every sampled ``q`` is carried alongside it
  (``q_mc_se``), and is not small near the coarse boundary — at ``q = 2e-5`` the expected hit count in
  20 000 draws is 0.4, so a projected "VOID" versus "coarse" call this close to the boundary should be
  read as noisy, not decisive. This is disclosed rather than hidden behind a single number.

**What is exact regardless of branch.** The assignment total's mean and variance under balanced
sampling without replacement have closed forms that depend only on the constructed vector and the
group sizes (finite-population sampling theory: ``mean = n_test/D * T``,
``var = n_test*(D-n_test)/(D-1) * population_variance(vector)``) — no enumeration or sampling is
needed for these, and ``tests/test_a2_feasibility.py`` cross-checks the formula against the exact
enumerated distribution on real committed strata. ``z`` (the real split's position in units of that
exact standard deviation) is therefore exact under both dispersion extremes; only the *percentile*
and the tau-ladder ``q`` values are enumerated-exact or MC-estimated depending on branch.

**Determinism.** No timestamp, environment, or package version appears anywhere in the emitted
artifact. The Monte-Carlo draws are seeded per stratum, per dispersion, from
``numpy.random.SeedSequence([MC_SEED, sha256-derived key, dispersion index])`` — independent of
iteration order, so reordering strata in a future revision of this script cannot change any other
stratum's draws. Re-running ``python scripts/a2_feasibility.py`` must reproduce the committed JSON
and CSV bytes exactly; ``--check`` verifies this without writing anything.

**Provenance of the four tier facts** the artifact states before computing anything else — frozen tier
150 strata, minimum ``permutation_count`` 24 310; manifest tier 554 strata, minimum 12 870; 32 of 150
and 88 of 554 fully enumerable at the 2e6 cap — is arithmetic over the two pinned sources with no
judgement involved, and is asserted against hard-coded expectations at load time (``main`` refuses to
write an artifact on a mismatch, the same discipline ``freeze_stratum_list.py`` applies to its own
declared figures).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
PREREG_DIR = REPO / "pilot" / "preregistration"

# ---------------------------------------------------------------------------
# Pinned sources. Hash-checked before parsing; refuses to run on a mismatch.
# ---------------------------------------------------------------------------

STRATUM_LIST_JSON = PREREG_DIR / "stratum_list_2026-08-16.json"
STRATUM_LIST_CSV = PREREG_DIR / "stratum_list_2026-08-16.csv"
STRATUM_LIST_SHA256 = "b39334e861020a2cfda503f576e2a7bd0b2197df9002af404bc44b1acbc1535e"
STRATUM_LIST_BYTES = 786_952
STRATUM_LIST_CSV_SHA256 = "72089ce8b9c96edfee250dc4efd69ede3f046f3f89812c65d88bf626cce7626e"
STRATUM_LIST_CSV_BYTES = 362_394
STRATUM_LIST_FROZEN_DATE = "2026-08-16"
STRATUM_LIST_N_FROZEN_STRATA = 251
STRATUM_LIST_N_CONTROL_STRATA = 106

CI_RUN_ID = "31910799023"
MANIFEST_JSON = PREREG_DIR / f"census_candidates_run{CI_RUN_ID}_2026-08-15.json"
MANIFEST_CSV = PREREG_DIR / f"census_candidates_run{CI_RUN_ID}_2026-08-15.csv"
MANIFEST_SHA256 = "33f8a800229dccc5f58f311e7d0c493655068d43563b31ff53fdaebb3b44e4b4"
MANIFEST_BYTES = 6_630_446
MANIFEST_CSV_SHA256 = "09eb110dd308155f64e10b2b05beff36854f7125b1434699935707d3551f12d6"
MANIFEST_CSV_BYTES = 4_513_660
MANIFEST_GENERATED_UTC = "2026-08-15T22:18:37+00:00"
MANIFEST_N_ROWS = 2190
MANIFEST_N_CANDIDATE_ROWS = 1197

# ---------------------------------------------------------------------------
# Constants. NOT sourced from gate_config — A2 has no gate_config entry yet (that is plan T4,
# gated on the Amendment 6 text of T2, which is itself gated on this file). Defined locally so this
# script has zero dependency on unwritten code; a future T4 that centralises these values must not
# change what this committed artifact already says without a fresh, dated run.
# ---------------------------------------------------------------------------

#: The candidate tolerances being bracketed (plan sec 0 / sec 1.1).
TAU_LADDER = (0.01, 0.025, 0.05, 0.075, 0.1)

#: Enumerate exactly at or below this design count; sample above it (plan sec 1.2).
A2_ENUM_CAP = 2_000_000

#: Enumerated-branch thin-set thresholds on the exact window count (plan sec 1.2).
A2_FULL_MIN = 1000
A2_COARSE_MIN = 200

#: Sampled-branch thresholds on acceptance rate, algebraically equivalent to sec 1.2's sample-budget
#: thresholds (1000 / 1e6, 200 / 1e6) — see the module docstring for why this script does not try to
#: match the real machinery's own 1e6/1e7 budget.
Q_FULL_THRESHOLD = 1000 / 1_000_000
Q_COARSE_THRESHOLD = 200 / 1_000_000

#: This script's own Monte-Carlo budget for the sampled branch — an estimation shortcut, smaller than
#: the real machinery's A2_SAMPLE_BUDGET (1e6) / A2_ESCALATED_BUDGET (1e7). See module docstring.
MC_DRAWS = 20_000
MC_SEED = 7

STATUS_FULL = "full"
STATUS_COARSE = "coarse"
STATUS_VOID = "void"

FROZEN_TIER = "frozen"
MANIFEST_TIER = "manifest"

FREEZE_DATE = "2026-08-16"
OUT_STEM = f"a2_feasibility_{FREEZE_DATE}"


class SourceArtifactMismatch(RuntimeError):
    """A pinned source file is missing, resized, or does not hash to the pinned value."""


class DeclaredFigureMismatch(RuntimeError):
    """A figure this script asserts about its own sources does not hold."""


# ---------------------------------------------------------------------------
# Loading, hash-pinned.
# ---------------------------------------------------------------------------


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_pinned_file(path: Path, digest: str, size: int, what: str) -> None:
    if not path.exists():
        raise SourceArtifactMismatch(f"{what} not found at {path}")
    actual_size = path.stat().st_size
    if actual_size != size:
        raise SourceArtifactMismatch(
            f"{what} ({path.name}): {actual_size} bytes, expected {size} — this is not the pinned "
            "artifact and the bracket below would describe different bytes."
        )
    actual = sha256_of(path)
    if actual != digest:
        raise SourceArtifactMismatch(f"{what} ({path.name}): sha256 {actual}, expected {digest}")


def load_stratum_list(
    path: Path = STRATUM_LIST_JSON, *, csv_path: Path | None = STRATUM_LIST_CSV
) -> dict:
    check_pinned_file(path, STRATUM_LIST_SHA256, STRATUM_LIST_BYTES, "frozen stratum list")
    if csv_path is not None:
        check_pinned_file(csv_path, STRATUM_LIST_CSV_SHA256, STRATUM_LIST_CSV_BYTES,
                           "frozen stratum list CSV twin")
    data = json.loads(path.read_text(encoding="utf-8"))
    header = data["header"]
    if header.get("frozen_date") != STRATUM_LIST_FROZEN_DATE:
        raise SourceArtifactMismatch(
            f"stratum list frozen_date {header.get('frozen_date')!r}, expected "
            f"{STRATUM_LIST_FROZEN_DATE!r}"
        )
    if header.get("n_frozen_strata") != STRATUM_LIST_N_FROZEN_STRATA:
        raise SourceArtifactMismatch(
            f"stratum list n_frozen_strata {header.get('n_frozen_strata')}, expected "
            f"{STRATUM_LIST_N_FROZEN_STRATA}"
        )
    if len(data["rows"]) != STRATUM_LIST_N_FROZEN_STRATA:
        raise SourceArtifactMismatch(
            f"stratum list carries {len(data['rows'])} analysis-set rows, expected "
            f"{STRATUM_LIST_N_FROZEN_STRATA}"
        )
    if len(data["within_collection_control_rows"]) != STRATUM_LIST_N_CONTROL_STRATA:
        raise SourceArtifactMismatch(
            f"stratum list carries {len(data['within_collection_control_rows'])} control rows, "
            f"expected {STRATUM_LIST_N_CONTROL_STRATA}"
        )
    return data


def load_manifest(path: Path = MANIFEST_JSON, *, csv_path: Path | None = MANIFEST_CSV) -> dict:
    check_pinned_file(path, MANIFEST_SHA256, MANIFEST_BYTES, "census candidate manifest")
    if csv_path is not None:
        check_pinned_file(csv_path, MANIFEST_CSV_SHA256, MANIFEST_CSV_BYTES,
                           "census candidate manifest CSV twin")
    data = json.loads(path.read_text(encoding="utf-8"))
    header = data["header"]
    if header.get("generated_utc") != MANIFEST_GENERATED_UTC:
        raise SourceArtifactMismatch(
            f"manifest generated_utc {header.get('generated_utc')!r}, expected "
            f"{MANIFEST_GENERATED_UTC!r}"
        )
    if len(data["rows"]) != MANIFEST_N_ROWS:
        raise SourceArtifactMismatch(f"manifest has {len(data['rows'])} rows, expected {MANIFEST_N_ROWS}")
    n_candidate = sum(1 for r in data["rows"] if r["gate_status"] == "candidate")
    if n_candidate != MANIFEST_N_CANDIDATE_ROWS:
        raise SourceArtifactMismatch(
            f"manifest has {n_candidate} candidate rows, expected {MANIFEST_N_CANDIDATE_ROWS}"
        )
    return data


# ---------------------------------------------------------------------------
# Tier selection.
# ---------------------------------------------------------------------------


def _min_donors(row: dict) -> int:
    return min(row["n_donors_A"], row["n_donors_B"])


def select_frozen_tier(stratum_list: dict, *, min_donors: int = 8) -> list[dict]:
    """The frozen 150: analysis-set rows (``role`` is always ``analysis_set`` here) at >= 8v8."""
    return [r for r in stratum_list["rows"] if _min_donors(r) >= min_donors]


def select_manifest_tier(manifest: dict, *, min_donors: int = 8) -> list[dict]:
    """The manifest 554: every candidate row of the whole Census pin at >= 8v8."""
    return [
        r for r in manifest["rows"]
        if r["gate_status"] == "candidate" and _min_donors(r) >= min_donors
    ]


def stratum_key(row: dict) -> str:
    return f"{row['dataset_id']}|{row['cell_type']}|{row['disease']}|{row.get('reference', '')}"


# ---------------------------------------------------------------------------
# Bracket vector construction (see module docstring for the full derivation).
# ---------------------------------------------------------------------------


def spread_evenly(total: int, n: int) -> list[int]:
    """``n`` integers, differing by at most one, summing exactly to ``total``."""
    base, rem = divmod(total, n)
    return [base + 1] * rem + [base] * (n - rem)


def build_low_group_vector(n_donors: int, n_cells: int, lo: int, median: float, hi: int) -> list[int]:
    """Minimum-plausible-variance per-donor vector consistent with the four hard facts.

    ``median`` is accepted for signature symmetry with :func:`build_high_group_vector` and is not
    used: this construction targets the achievable total directly rather than starting from the
    median and re-clamping, which is exactly the perturbation the module docstring discloses for
    the *conceptual* clamp-then-correct route. The realised median of the vector this returns is
    computed by the caller and reported as a drift regardless.
    """
    del median
    if n_donors == 1:
        if n_cells != lo or n_cells != hi:
            raise DeclaredFigureMismatch(
                f"single-donor group: n_cells={n_cells} but min={lo}, max={hi} disagree"
            )
        return [n_cells]
    interior_n = n_donors - 2
    if interior_n == 0:
        if lo + hi != n_cells:
            raise DeclaredFigureMismatch(
                f"two-donor group: min+max={lo + hi} != n_cells={n_cells}"
            )
        return sorted([lo, hi])
    target = n_cells - lo - hi
    interior = spread_evenly(target, interior_n)
    full = [lo, hi, *interior]
    if sum(full) != n_cells:  # pragma: no cover - arithmetic invariant of spread_evenly
        raise DeclaredFigureMismatch(f"low-dispersion vector sums to {sum(full)}, expected {n_cells}")
    return sorted(full)


def _mid_values(mid_slots: int, median: float) -> list[int]:
    """The one or two integer donor values a reported median implies, exactly.

    For an odd donor count the median *is* one real integer value. For an even count the two
    middle real values sum to ``2 * median`` by the definition of an even-length median — an
    identity, not an approximation — so ``floor(median), ceil(median)`` reproduces them whenever
    they are adjacent integers (true whenever ``median`` has a fractional part of exactly 0.5) and
    reproduces a repeated integer median exactly when ``median`` is a whole number.
    """
    if mid_slots == 1:
        if float(median) != math.floor(median):
            raise DeclaredFigureMismatch(
                f"odd donor count but median {median} is not an integer"
            )
        return [int(median)]
    lo_m, hi_m = math.floor(median), math.ceil(median)
    if lo_m + hi_m != 2 * median:  # pragma: no cover - arithmetic identity of floor/ceil
        raise DeclaredFigureMismatch(f"floor+ceil of median {median} do not sum to 2*median")
    return [lo_m, hi_m]


def build_high_group_vector(n_donors: int, n_cells: int, lo: int, median: float, hi: int) -> list[int]:
    """Maximum-plausible-variance per-donor vector consistent with all five reported statistics.

    Pins ``min``/``max`` at the endpoints and the median at one (odd) or two (even) interior slots,
    then piles the remaining interior donors two-point at ``max``/``min`` to hit the exact total,
    with at most one intermediate "leftover" donor absorbing the integer remainder. The leftover is
    proven to always land in ``[lo, hi]`` (see module docstring; checked by the test suite) so this
    function never needs to clip.
    """
    if n_donors == 1:
        if n_cells != lo or n_cells != hi:
            raise DeclaredFigureMismatch(
                f"single-donor group: n_cells={n_cells} but min={lo}, max={hi} disagree"
            )
        return [n_cells]
    if n_donors == 2:
        # no room for a median-designated slot distinct from the two anchors: for D=2 the real
        # vector *is* (lo, hi), and mid_slots/free_n below would go negative if this were skipped.
        if lo + hi != n_cells:
            raise DeclaredFigureMismatch(
                f"two-donor group: min+max={lo + hi} != n_cells={n_cells}"
            )
        return sorted([lo, hi])
    mid_slots = 1 if n_donors % 2 == 1 else 2
    mids = _mid_values(mid_slots, median)
    interior_n = n_donors - 2
    free_n = interior_n - mid_slots
    remaining = n_cells - lo - hi - sum(mids)
    if free_n == 0:
        if remaining != 0:
            raise DeclaredFigureMismatch(
                f"donor count {n_donors}: min+max+median-slots do not sum to n_cells "
                f"(off by {remaining}) — the reported statistics are inconsistent"
            )
        free_vals: list[int] = []
    elif hi == lo:
        if remaining != lo * free_n:
            raise DeclaredFigureMismatch(
                f"degenerate group (min==max=={lo}) but remaining total {remaining} != {lo * free_n}"
            )
        free_vals = [lo] * free_n
    else:
        k_real = (remaining - lo * free_n) / (hi - lo)
        if k_real >= free_n - 1e-9:
            free_vals = [hi] * free_n
        elif k_real <= 1e-9:
            free_vals = [lo] * free_n
        else:
            k = min(max(int(math.floor(k_real + 1e-9)), 0), free_n - 1)
            leftover = remaining - (k * hi + (free_n - k - 1) * lo)
            if not (lo <= leftover <= hi):  # pragma: no cover - proven unreachable, kept as a guard
                raise DeclaredFigureMismatch(
                    f"high-dispersion leftover donor value {leftover} escaped [{lo}, {hi}] for "
                    f"n_donors={n_donors}, n_cells={n_cells}, min={lo}, median={median}, max={hi}"
                )
            free_vals = [hi] * k + [lo] * (free_n - k - 1) + [leftover]
    full = [lo, hi, *mids, *free_vals]
    if sum(full) != n_cells:  # pragma: no cover - arithmetic invariant of the branches above
        raise DeclaredFigureMismatch(f"high-dispersion vector sums to {sum(full)}, expected {n_cells}")
    return sorted(full)


def realised_median(values: list[int]) -> float:
    s = sorted(values)
    n = len(s)
    mid = n // 2
    return float(s[mid]) if n % 2 else (s[mid - 1] + s[mid]) / 2.0


def build_group_vector(gs: dict, dispersion: str) -> tuple[list[int], float]:
    """Returns ``(vector, median_drift)`` for one group under one dispersion extreme."""
    n_donors, n_cells = gs["n_donors"], gs["n_cells"]
    lo, median, hi = gs["min"], gs["median"], gs["max"]
    builder = build_low_group_vector if dispersion == "low" else build_high_group_vector
    vector = builder(n_donors, n_cells, lo, median, hi)
    drift = abs(realised_median(vector) - median)
    return vector, drift


def build_combined_vector(row: dict, dispersion: str) -> tuple[list[int], float, float]:
    """The whole-stratum vector (group A donors, then group B), plus each group's median drift."""
    vec_a, drift_a = build_group_vector(row["cells_per_donor_by_group"]["A"], dispersion)
    vec_b, drift_b = build_group_vector(row["cells_per_donor_by_group"]["B"], dispersion)
    return [*vec_a, *vec_b], drift_a, drift_b


# ---------------------------------------------------------------------------
# Exact enumeration (subset-sum dynamic program) — the enumerated branch.
# ---------------------------------------------------------------------------


def enumerate_subset_sum_counts(values: list[int], k: int) -> np.ndarray:
    """``counts[s]`` = number of size-``k`` subsets of ``values`` (donors distinct by position,
    even when two donors tie in value) summing to ``s``, for ``s`` in ``[0, sum(values)]``.

    A standard 0/1-knapsack-with-cardinality dynamic program, vectorised over the sum axis with
    numpy so the Python-level loop is only ``len(values) * k`` iterations. Counts fit easily in
    int64: they are bounded by ``math.comb(len(values), k)``, which is <= ``A2_ENUM_CAP`` (2e6) for
    every stratum this function is called on.
    """
    total = sum(values)
    dp = np.zeros((k + 1, total + 1), dtype=np.int64)
    dp[0, 0] = 1
    for val in values:
        for j in range(k, 0, -1):
            if val == 0:
                dp[j] += dp[j - 1]
            else:
                dp[j, val:] += dp[j - 1, : total + 1 - val]
    return dp[k]


def finite_population_mean_var(values: list[int], k: int) -> tuple[float, float]:
    """Exact mean/variance of a size-``k`` sample-without-replacement total (no enumeration needed).

    Closed form from finite-population sampling theory: for population ``values`` of size ``D``
    with population variance ``sigma2`` (divisor ``D``, not ``D-1``), the sum of a random size-``k``
    subset has ``mean = k/D * T`` and ``var = k*(D-k)/(D-1) * sigma2``. Derived and checked against
    the exact enumerated distribution in ``tests/test_a2_feasibility.py``.
    """
    arr = np.asarray(values, dtype=np.float64)
    d = arr.size
    total = arr.sum()
    mean = k / d * total
    if d == 1:
        return float(mean), 0.0
    sigma2 = float(np.mean((arr - arr.mean()) ** 2))
    var = k * (d - k) / (d - 1) * sigma2
    return float(mean), float(var)


# ---------------------------------------------------------------------------
# Monte-Carlo estimation — the sampled branch (permutation_count > A2_ENUM_CAP).
# ---------------------------------------------------------------------------


def stratum_seed_material(row: dict, dispersion: str) -> list[int]:
    key = stratum_key(row).encode("utf-8")
    digest = hashlib.sha256(key).digest()
    derived = int.from_bytes(digest[:8], "big")
    return [MC_SEED, derived, 0 if dispersion == "low" else 1]


def mc_assignment_totals(values: list[int], k: int, seed_material: list[int], n_draws: int) -> np.ndarray:
    """``n_draws`` balanced-assignment totals, drawn without replacement, seeded deterministically."""
    arr = np.asarray(values, dtype=np.int64)
    d = arr.size
    rng = np.random.default_rng(seed_material)
    order = np.argsort(rng.random((n_draws, d)), axis=1)
    picks = order[:, :k]
    return arr[picks].sum(axis=1)


# ---------------------------------------------------------------------------
# Per-stratum, per-dispersion analysis.
# ---------------------------------------------------------------------------


def classify(status_kind: str, *, m: int | None = None, q: float | None = None) -> str:
    if status_kind == "enumerated":
        if m >= A2_FULL_MIN:
            return STATUS_FULL
        if m >= A2_COARSE_MIN:
            return STATUS_COARSE
        return STATUS_VOID
    if q >= Q_FULL_THRESHOLD:
        return STATUS_FULL
    if q >= Q_COARSE_THRESHOLD:
        return STATUS_COARSE
    return STATUS_VOID


def analyze_dispersion(row: dict, dispersion: str, tau_ladder: tuple[float, ...]) -> dict:
    vector, drift_a, drift_b = build_combined_vector(row, dispersion)
    d = len(vector)
    n_test = row["n_donors_A"]
    total_cells = row["n_cells"]
    real_n = row["cells_per_donor_by_group"]["A"]["n_cells"]
    permutation_count = math.comb(d, n_test)

    mean, var = finite_population_mean_var(vector, n_test)
    std = math.sqrt(var) if var > 0 else 0.0
    z = (real_n - mean) / std if std > 0 else None

    enumerated = permutation_count <= A2_ENUM_CAP
    per_tau = []
    if enumerated:
        counts = enumerate_subset_sum_counts(vector, n_test)
        if int(counts.sum()) != permutation_count:  # pragma: no cover - DP invariant
            raise DeclaredFigureMismatch(
                f"enumerated distribution sums to {int(counts.sum())}, expected "
                f"{permutation_count} designs"
            )
        percentile = float(counts[: real_n + 1].sum()) / permutation_count
        for tau in tau_ladder:
            t_star = math.floor(tau * total_cells)
            lo_w = max(0, real_n - t_star)
            hi_w = min(total_cells, real_n + t_star)
            m = int(counts[lo_w : hi_w + 1].sum())
            per_tau.append({
                "tau": tau, "t_star": t_star, "q": m / permutation_count,
                "n_survivors": m, "survivors_denominator": permutation_count,
                "q_mc_se": None, "status": classify("enumerated", m=m),
            })
    else:
        draws = mc_assignment_totals(vector, n_test, stratum_seed_material(row, dispersion), MC_DRAWS)
        percentile = float(np.mean(draws <= real_n))
        for tau in tau_ladder:
            t_star = math.floor(tau * total_cells)
            hits = int(np.sum(np.abs(draws - real_n) <= t_star))
            q_hat = hits / MC_DRAWS
            se = math.sqrt(q_hat * (1 - q_hat) / MC_DRAWS)
            per_tau.append({
                "tau": tau, "t_star": t_star, "q": q_hat,
                "n_survivors": hits, "survivors_denominator": MC_DRAWS,
                "q_mc_se": se, "status": classify("sampled", q=q_hat),
            })

    for entry in per_tau:
        entry["void"] = entry["status"] == STATUS_VOID

    return {
        "resolution": "enumerated" if enumerated else "sampled",
        "permutation_count": permutation_count,
        "mean": mean, "std": std, "z": z, "percentile": percentile,
        "median_drift_A": drift_a, "median_drift_B": drift_b,
        "per_tau": per_tau,
    }


def analyze_stratum(row: dict, tau_ladder: tuple[float, ...] = TAU_LADDER) -> dict:
    total_cells = row["n_cells"]
    real_n = row["cells_per_donor_by_group"]["A"]["n_cells"]
    real_b = row["cells_per_donor_by_group"]["B"]["n_cells"]
    if real_n + real_b != total_cells:
        raise DeclaredFigureMismatch(
            f"{row['dataset_id']}/{row['cell_type']}/{row['disease']}: group cell totals "
            f"{real_n} + {real_b} != stratum n_cells {total_cells}"
        )

    low = analyze_dispersion(row, "low", tau_ladder)
    high = analyze_dispersion(row, "high", tau_ladder)
    if low["permutation_count"] != high["permutation_count"]:  # pragma: no cover - invariant
        raise DeclaredFigureMismatch("low/high dispersion disagree on permutation_count")
    declared_permutation_count = row["permutation_count"]
    if low["permutation_count"] != declared_permutation_count:
        raise DeclaredFigureMismatch(
            f"{row['dataset_id']}/{row['cell_type']}/{row['disease']}: computed permutation_count "
            f"{low['permutation_count']}, row declares {declared_permutation_count}"
        )

    return {
        "dataset_id": row["dataset_id"],
        "dataset_short": row.get("dataset_short"),
        "cell_type": row["cell_type"],
        "disease": row["disease"],
        "reference": row["reference"],
        "role": row.get("role", "candidate"),
        "n_donors_A": row["n_donors_A"],
        "n_donors_B": row["n_donors_B"],
        "n_donors_total": row["n_donors_A"] + row["n_donors_B"],
        "n_cells": total_cells,
        "n_cells_test_real": real_n,
        "f_real": real_n / total_cells,
        "permutation_count": declared_permutation_count,
        "resolution": low["resolution"],
        "mean_assignment_total": low["mean"],
        "low": {
            "std": low["std"], "z": low["z"], "percentile": low["percentile"],
            "median_drift_A": low["median_drift_A"], "median_drift_B": low["median_drift_B"],
            "per_tau": low["per_tau"],
        },
        "high": {
            "std": high["std"], "z": high["z"], "percentile": high["percentile"],
            "median_drift_A": high["median_drift_A"], "median_drift_B": high["median_drift_B"],
            "per_tau": high["per_tau"],
        },
    }


# ---------------------------------------------------------------------------
# Tier-level aggregation.
# ---------------------------------------------------------------------------


def _median(values: list[float]) -> float:
    s = sorted(values)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2.0


def summarise_tier(strata: list[dict], tau_ladder: tuple[float, ...]) -> list[dict]:
    n = len(strata)
    rows = []
    for i, tau in enumerate(tau_ladder):
        row = {"tau": tau}
        for disp in ("low", "high"):
            statuses = [s[disp]["per_tau"][i]["status"] for s in strata]
            qs = [s[disp]["per_tau"][i]["q"] for s in strata]
            row[f"void_{disp}"] = statuses.count(STATUS_VOID)
            row[f"coarse_{disp}"] = statuses.count(STATUS_COARSE)
            row[f"full_{disp}"] = statuses.count(STATUS_FULL)
            row[f"void_{disp}_pct"] = round(100.0 * statuses.count(STATUS_VOID) / n, 4) if n else 0.0
            row[f"q_median_{disp}"] = round(_median(qs), 6) if qs else None
        rows.append(row)
    return rows


def extreme_f_summary(strata: list[dict], *, lo: float = 0.1, hi: float = 0.9) -> dict:
    outside = [s for s in strata if not (lo <= s["f_real"] <= hi)]
    return {
        "band": [lo, hi],
        "n_outside_band": len(outside),
        "dataset_ids": sorted({s["dataset_id"] for s in outside}),
    }


def void_index(strata: list[dict], tau_ladder: tuple[float, ...]) -> dict:
    index: dict[str, dict[str, list[dict]]] = {}
    for i, tau in enumerate(tau_ladder):
        tau_key = f"{tau:g}"
        index[tau_key] = {"low": [], "high": []}
        for s in strata:
            for disp in ("low", "high"):
                if s[disp]["per_tau"][i]["void"]:
                    index[tau_key][disp].append({
                        "dataset_id": s["dataset_id"], "dataset_short": s["dataset_short"],
                        "cell_type": s["cell_type"], "disease": s["disease"],
                    })
    return index


# ---------------------------------------------------------------------------
# Artifact assembly.
# ---------------------------------------------------------------------------

BRACKET_CONSTRUCTION = {
    "inputs": (
        "each committed row carries, per group (A, B), only n_donors, n_cells (group total), min, "
        "median and max of per-donor cell counts -- never the per-donor vector. Two integer "
        "vectors per group are built from those five numbers, one per dispersion extreme, and "
        "concatenated (group A donors, then group B) into one whole-stratum vector of length "
        "n_donors_A + n_donors_B. The true vector is unrecoverable from committed data; these are "
        "disclosed, deliberately extreme constructions, not estimates of it."
    ),
    "low_dispersion": (
        "one donor is pinned at the reported min, one at the reported max. The remaining "
        "(n_donors - 2) donors are all set as close to a single shared value as integer "
        "arithmetic allows: (n_cells - min - max) split into (n_donors - 2) parts differing by at "
        "most one cell. This is the minimum-variance vector consistent with n_donors, n_cells, min "
        "and max. It does not target the reported median at all -- see 'clamping_perturbs_median'."
    ),
    "high_dispersion": (
        "the same min/max pins, plus one (odd n_donors) or two (even n_donors) interior donors "
        "pinned exactly at the reported median, plus a two-point fill of the remaining interior "
        "donors (as many as possible at max, the rest at min, at most one intermediate 'leftover' "
        "donor absorbing the integer remainder) chosen to hit the exact total. This is the "
        "maximum-variance vector on a bounded support consistent with all five reported "
        "statistics. The leftover donor is proven to always land inside [min, max] -- it is never "
        "clipped -- but because the fill targets the exact total rather than balancing how many "
        "donors sit above versus below the pinned median, the finished vector's own sorted "
        "(realised) median can still differ from the reported one. This drift is measured per "
        "group and reported as median_drift_A / median_drift_B under 'high' in every stratum "
        "record; it is not assumed away."
    ),
    "clamping_perturbs_median": (
        "the plan's own note is that forcing a vector's interior donors to hit the exact recorded "
        "total generally moves them off the reported median. The low-dispersion construction "
        "exhibits this directly: its interior donors land near "
        "(n_cells - min - max) / (n_donors - 2), which need not equal, and generally does not "
        "equal, the reported median. This is measured per group and reported as "
        "median_drift_A / median_drift_B under 'low' in every stratum record."
    ),
    "exact_reconstruction_below_5_donors": (
        "for n_donors in {1, 2, 3, 4} the five reported statistics determine the per-donor "
        "multiset exactly (n_donors <= 3: the sorted vector literally is (min, median, max); "
        "n_donors == 4: the two middle values sum to 2*median by definition), so both "
        "dispersion extremes coincide and median_drift is exactly 0. Neither tier analysed here "
        "exercises this case (every group holds >= 8 donors by the >= 8v8 tier filter); it is "
        "exercised by tests/test_a2_feasibility.py."
    ),
    "limitation": (
        "the committed rows carry five summary statistics per group, not the per-donor vector, so "
        "the true donor-permutation assignment-total distribution is not computable from committed "
        "data. Every number in this artifact is a bracket between two extremes consistent with "
        "those five statistics, not a measurement of the real distribution; where the two extremes "
        "disagree, the true figure lies between them at an unknown point."
    ),
}

CLASSIFICATION_METHOD = {
    "enumerated_branch": (
        f"permutation_count <= A2_ENUM_CAP ({A2_ENUM_CAP:,}): every design's assignment total is "
        "computed exactly by a subset-sum dynamic program (no sampling). The window count m is "
        f"exact and classified full (m >= {A2_FULL_MIN}) / coarse ({A2_COARSE_MIN} <= m < "
        f"{A2_FULL_MIN}) / void (m < {A2_COARSE_MIN}), matching plan sec 1.2's enumerated branch."
    ),
    "sampled_branch": (
        f"permutation_count > A2_ENUM_CAP: q is estimated from MC_DRAWS ({MC_DRAWS:,}) balanced "
        f"draws (seeded per stratum, per dispersion -- see 'seeding' below) and classified against "
        f"the acceptance-rate form of the plan's sample-budget thresholds: full at q >= "
        f"{Q_FULL_THRESHOLD:g}, coarse at {Q_COARSE_THRESHOLD:g} <= q < {Q_FULL_THRESHOLD:g}, void "
        f"below {Q_COARSE_THRESHOLD:g}. This MC_DRAWS budget is smaller than the real machinery's "
        "own A2_SAMPLE_BUDGET (1e6, escalating to 1e7) -- see the module docstring for why -- and "
        "every sampled q carries its Monte-Carlo standard error (q_mc_se) so the reader can see "
        "how much noise a projected status near a threshold carries, particularly near the coarse "
        "boundary where the expected hit count in 20,000 draws is a handful."
    ),
    "exact_moments": (
        "mean and variance of the assignment-total distribution have closed forms under balanced "
        "sampling without replacement (finite-population sampling theory) that need neither "
        "enumeration nor sampling; z is therefore exact under both dispersion extremes regardless "
        "of branch. Only the percentile and the tau-ladder q are enumerated-exact or MC-estimated. "
        "Cross-checked against the exact enumerated distribution in tests/test_a2_feasibility.py."
    ),
    "seeding": (
        "each (stratum, dispersion) Monte-Carlo draw set is seeded from "
        "numpy.random.SeedSequence([MC_SEED, sha256-derived-key, dispersion_index]), where the key "
        "is derived from (dataset_id, cell_type, disease, reference) -- independent of iteration "
        "order, so this script's own row order can change in a future revision without changing "
        "any stratum's draws."
    ),
}


def build_artifact() -> dict:
    stratum_list = load_stratum_list()
    manifest = load_manifest()

    frozen_rows = select_frozen_tier(stratum_list)
    manifest_rows = select_manifest_tier(manifest)

    declared = {
        FROZEN_TIER: {"n_strata": 150, "min_permutation_count": 24_310, "n_fully_enumerable": 32},
        MANIFEST_TIER: {"n_strata": 554, "min_permutation_count": 12_870, "n_fully_enumerable": 88},
    }
    for tier_name, rows in ((FROZEN_TIER, frozen_rows), (MANIFEST_TIER, manifest_rows)):
        exp = declared[tier_name]
        if len(rows) != exp["n_strata"]:
            raise DeclaredFigureMismatch(
                f"{tier_name} tier holds {len(rows)} strata at >= 8v8, expected {exp['n_strata']}"
            )
        min_pc = min(r["permutation_count"] for r in rows)
        if min_pc != exp["min_permutation_count"]:
            raise DeclaredFigureMismatch(
                f"{tier_name} tier minimum permutation_count is {min_pc}, expected "
                f"{exp['min_permutation_count']}"
            )
        n_enum = sum(1 for r in rows if r["permutation_count"] <= A2_ENUM_CAP)
        if n_enum != exp["n_fully_enumerable"]:
            raise DeclaredFigureMismatch(
                f"{tier_name} tier holds {n_enum} strata enumerable at the {A2_ENUM_CAP:,} cap, "
                f"expected {exp['n_fully_enumerable']}"
            )

    frozen_analysed = [analyze_stratum(r) for r in frozen_rows]
    manifest_analysed = [analyze_stratum(r) for r in manifest_rows]

    header = {
        "act": (
            "the tau-ladder VOID-projection bracket for the A2 cell-count stratification "
            "correction (plan sec 6, task T1). An evidence artifact a future, separately dated "
            "Amendment 6 may cite when choosing A2_TAU; not itself an amendment and not a claim "
            "about any real dataset's biology."
        ),
        "generated_by": "scripts/a2_feasibility.py",
        "sources": {
            "stratum_list": {
                "json": STRATUM_LIST_JSON.name, "json_sha256": STRATUM_LIST_SHA256,
                "json_bytes": STRATUM_LIST_BYTES,
                "csv": STRATUM_LIST_CSV.name, "csv_sha256": STRATUM_LIST_CSV_SHA256,
                "csv_bytes": STRATUM_LIST_CSV_BYTES,
                "frozen_date": STRATUM_LIST_FROZEN_DATE,
            },
            "census_candidates_manifest": {
                "json": MANIFEST_JSON.name, "json_sha256": MANIFEST_SHA256,
                "json_bytes": MANIFEST_BYTES,
                "csv": MANIFEST_CSV.name, "csv_sha256": MANIFEST_CSV_SHA256,
                "csv_bytes": MANIFEST_CSV_BYTES,
                "ci_run_id": CI_RUN_ID, "generated_utc": MANIFEST_GENERATED_UTC,
            },
        },
        "tier_definitions": {
            FROZEN_TIER: (
                "stratum_list rows (role = analysis_set, the 251-stratum pre-registered analysis "
                "set) filtered to min(n_donors_A, n_donors_B) >= 8."
            ),
            MANIFEST_TIER: (
                "census_candidates rows with gate_status == 'candidate' (the whole pinned "
                "Census's candidate set, 1197 rows across 68 datasets) filtered to "
                "min(n_donors_A, n_donors_B) >= 8. The frozen tier is a strict subset of this "
                "tier -- every frozen stratum is also a manifest stratum -- reported separately "
                "because the plan asks for the bracket on both the twelve-dataset analysis set "
                "and the wider manifest."
            ),
        },
        "tau_ladder": list(TAU_LADDER),
        "constants": {
            "enum_cap": A2_ENUM_CAP, "full_min": A2_FULL_MIN, "coarse_min": A2_COARSE_MIN,
            "q_full_threshold": Q_FULL_THRESHOLD, "q_coarse_threshold": Q_COARSE_THRESHOLD,
            "mc_draws": MC_DRAWS, "mc_seed": MC_SEED,
        },
        "tiers": {
            FROZEN_TIER: {
                "n_strata": len(frozen_rows),
                "min_permutation_count": min(r["permutation_count"] for r in frozen_rows),
                "n_fully_enumerable_at_cap": sum(
                    1 for r in frozen_rows if r["permutation_count"] <= A2_ENUM_CAP
                ),
            },
            MANIFEST_TIER: {
                "n_strata": len(manifest_rows),
                "min_permutation_count": min(r["permutation_count"] for r in manifest_rows),
                "n_fully_enumerable_at_cap": sum(
                    1 for r in manifest_rows if r["permutation_count"] <= A2_ENUM_CAP
                ),
            },
        },
        "test_group_convention": (
            "N_test = cells in group A (census_select's disease/test group; group B is the "
            "reference group -- see census_select.group_of). f_real = n_cells_test_real / n_cells. "
            "The tau-window condition |N(pi) - N_real| <= t* is symmetric under swapping which "
            "group is called 'test' (it is defined on a difference), so this convention affects no "
            "count in this artifact, only which real total is reported as n_cells_test_real."
        ),
        "bracket_construction": BRACKET_CONSTRUCTION,
        "classification_method": CLASSIFICATION_METHOD,
    }

    return {
        "header": header,
        "tau_ladder_summary": {
            FROZEN_TIER: summarise_tier(frozen_analysed, TAU_LADDER),
            MANIFEST_TIER: summarise_tier(manifest_analysed, TAU_LADDER),
        },
        "extreme_f_real": {
            FROZEN_TIER: extreme_f_summary(frozen_analysed),
            MANIFEST_TIER: extreme_f_summary(manifest_analysed),
        },
        "void_index": {
            FROZEN_TIER: void_index(frozen_analysed, TAU_LADDER),
            MANIFEST_TIER: void_index(manifest_analysed, TAU_LADDER),
        },
        "strata": {
            FROZEN_TIER: frozen_analysed,
            MANIFEST_TIER: manifest_analysed,
        },
    }


# ---------------------------------------------------------------------------
# Rendering — JSON and its CSV twin, byte-pinned line endings (freeze_stratum_list.py's own idiom).
# ---------------------------------------------------------------------------


def _round_floats(obj):
    """Round every float to 6 decimal places so the artifact does not carry binary float noise."""
    if isinstance(obj, float):
        return round(obj, 6)
    if isinstance(obj, dict):
        return {k: _round_floats(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_round_floats(v) for v in obj]
    return obj


def render_json(artifact: dict) -> str:
    return json.dumps(_round_floats(artifact), indent=1, ensure_ascii=False) + "\n"


CSV_FIELDS = (
    "tier", "dataset_id", "dataset_short", "cell_type", "disease", "reference", "role",
    "n_donors_A", "n_donors_B", "n_donors_total", "n_cells", "n_cells_test_real", "f_real",
    "permutation_count", "resolution", "mean_assignment_total",
    "std_low", "z_low", "percentile_low", "median_drift_A_low", "median_drift_B_low",
    "std_high", "z_high", "percentile_high", "median_drift_A_high", "median_drift_B_high",
    "tau", "t_star",
    "q_low", "n_survivors_low", "survivors_denominator_low", "q_mc_se_low", "status_low", "void_low",
    "q_high", "n_survivors_high", "survivors_denominator_high", "q_mc_se_high", "status_high",
    "void_high",
)


def _csv_rows(artifact: dict):
    for tier in (FROZEN_TIER, MANIFEST_TIER):
        for s in artifact["strata"][tier]:
            for i, tau in enumerate(TAU_LADDER):
                lo_t, hi_t = s["low"]["per_tau"][i], s["high"]["per_tau"][i]
                yield {
                    "tier": tier, "dataset_id": s["dataset_id"], "dataset_short": s["dataset_short"],
                    "cell_type": s["cell_type"], "disease": s["disease"], "reference": s["reference"],
                    "role": s["role"], "n_donors_A": s["n_donors_A"], "n_donors_B": s["n_donors_B"],
                    "n_donors_total": s["n_donors_total"], "n_cells": s["n_cells"],
                    "n_cells_test_real": s["n_cells_test_real"], "f_real": round(s["f_real"], 6),
                    "permutation_count": s["permutation_count"], "resolution": s["resolution"],
                    "mean_assignment_total": round(s["mean_assignment_total"], 6),
                    "std_low": round(s["low"]["std"], 6),
                    "z_low": None if s["low"]["z"] is None else round(s["low"]["z"], 6),
                    "percentile_low": round(s["low"]["percentile"], 6),
                    "median_drift_A_low": round(s["low"]["median_drift_A"], 6),
                    "median_drift_B_low": round(s["low"]["median_drift_B"], 6),
                    "std_high": round(s["high"]["std"], 6),
                    "z_high": None if s["high"]["z"] is None else round(s["high"]["z"], 6),
                    "percentile_high": round(s["high"]["percentile"], 6),
                    "median_drift_A_high": round(s["high"]["median_drift_A"], 6),
                    "median_drift_B_high": round(s["high"]["median_drift_B"], 6),
                    "tau": tau, "t_star": lo_t["t_star"],
                    "q_low": round(lo_t["q"], 6), "n_survivors_low": lo_t["n_survivors"],
                    "survivors_denominator_low": lo_t["survivors_denominator"],
                    "q_mc_se_low": None if lo_t["q_mc_se"] is None else round(lo_t["q_mc_se"], 6),
                    "status_low": lo_t["status"], "void_low": lo_t["void"],
                    "q_high": round(hi_t["q"], 6), "n_survivors_high": hi_t["n_survivors"],
                    "survivors_denominator_high": hi_t["survivors_denominator"],
                    "q_mc_se_high": None if hi_t["q_mc_se"] is None else round(hi_t["q_mc_se"], 6),
                    "status_high": hi_t["status"], "void_high": hi_t["void"],
                }


def _csv_value(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return ""
    return value


def render_csv(artifact: dict) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(CSV_FIELDS))
    writer.writeheader()
    for row in _csv_rows(artifact):
        writer.writerow({k: _csv_value(row[k]) for k in CSV_FIELDS})
    return buffer.getvalue()


def write(artifact: dict, out_dir: Path, *, stem: str = OUT_STEM) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{stem}.json"
    csv_path = out_dir / f"{stem}.csv"
    with open(json_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(render_json(artifact))
    with open(csv_path, "w", encoding="utf-8", newline="") as fh:
        fh.write(render_csv(artifact))
    return json_path, csv_path


def summarise(artifact: dict) -> str:
    header = artifact["header"]
    lines = [
        f"# A2 feasibility bracket -- {FREEZE_DATE}",
        f"  sources     : {header['sources']['stratum_list']['json']} "
        f"(sha256 {header['sources']['stratum_list']['json_sha256'][:12]}...)",
        f"                {header['sources']['census_candidates_manifest']['json']} "
        f"(sha256 {header['sources']['census_candidates_manifest']['json_sha256'][:12]}...)",
        "",
    ]
    for tier in (FROZEN_TIER, MANIFEST_TIER):
        t = header["tiers"][tier]
        lines.append(
            f"  {tier:<9} : {t['n_strata']} strata >= 8v8, min permutation_count "
            f"{t['min_permutation_count']:,}, {t['n_fully_enumerable_at_cap']} fully enumerable "
            f"at {A2_ENUM_CAP:,}"
        )
    lines.append("")
    lines.append("  tau     | frozen VOID low/high (%)      | manifest VOID low/high (%)")
    frozen_rows = artifact["tau_ladder_summary"][FROZEN_TIER]
    manifest_rows = artifact["tau_ladder_summary"][MANIFEST_TIER]
    for fr, mr in zip(frozen_rows, manifest_rows):
        lines.append(
            f"  {fr['tau']:<7} | {fr['void_low']:>3} / {fr['void_high']:>3} "
            f"({fr['void_low_pct']:>5.1f} / {fr['void_high_pct']:>5.1f})    | "
            f"{mr['void_low']:>3} / {mr['void_high']:>3} "
            f"({mr['void_low_pct']:>5.1f} / {mr['void_high_pct']:>5.1f})"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI.
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--output-dir", type=Path, default=PREREG_DIR,
        help=f"where {OUT_STEM}.json / .csv are written",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="verify the committed JSON and CSV are byte-identical to a fresh run; write nothing",
    )
    return parser


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)
    try:
        artifact = build_artifact()
    except (SourceArtifactMismatch, DeclaredFigureMismatch) as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2

    if args.check:
        for suffix, expected in ((".json", render_json(artifact)), (".csv", render_csv(artifact))):
            actual_path = args.output_dir / f"{OUT_STEM}{suffix}"
            if not actual_path.exists():
                print(f"MISSING: {actual_path}", file=sys.stderr)
                return 1
            with open(actual_path, encoding="utf-8", newline="") as fh:
                actual = fh.read()
            if actual != expected:
                print(f"DRIFT: {actual_path} is not what this script now produces", file=sys.stderr)
                return 1
            print(f"OK: {actual_path} matches a fresh run")
        return 0

    json_path, csv_path = write(artifact, args.output_dir)
    print(summarise(artifact))
    print(f"\n  wrote {json_path}")
    print(f"  wrote {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
