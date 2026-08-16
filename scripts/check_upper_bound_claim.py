"""Reproducer for Amendment 4 Part A's probe table: is ``sqrt(s0^2) * ln 2`` an upper bound on sigma?

WHY THIS EXISTS
---------------
``docs/AMENDMENTS.md``, Amendment 4 Part A ("The measurement behind Correction 1"), carries a
ten-row table of probe cells which demotes ``sqrt(s0²)·ln 2`` from Amendment 3's claimed **upper
bound** on ``donor_sigma`` to *an audit quantity of unknown error sign*. That table was produced
while the amendment was being written, from code that was never committed, on development seeds
that were never named, and with a per-row seed count that varied (six rows at 16 seeds, two at 8,
two at a single seed). A claim that overturns a previous amendment cannot rest on numbers nobody
can re-derive, so this script is the missing instrument: same ten cells, **one** named seed set for
all of them, an artifact that carries every per-seed realisation, and a ``header`` block that pins
the code the numbers came out of.

It is a reproducer, not a new method. Nothing here re-implements the estimator — the quantity is
read out of the pseudobulk arm's own entry point.

THE QUANTITY
------------
Per replicate, at ``n_de = 0`` (a null oracle: no condition effect, donor random effect present),

    sqrt(s0²) · ln 2

where ``s0²`` is the **untrended** empirical-Bayes prior variance the moderated arm fits on
``log2(CPM + 1)`` — ``moderation["s0_squared"]`` out of
:func:`pbcheck.methods.moderated.ebayes_from_pdata`. The ``ln 2`` converts a log2-scale standard
deviation to the natural-log scale on which ``synthetic.oracles.simulate``'s donor random effect
``re[g,d] = exp(N(−σ²/2, σ²))`` is defined, so the quantity is directly comparable with the
simulator's ``donor_sigma`` and Amendment 3's claim is exactly "this number is ≥ σ".

``trend=False`` is asserted, not assumed: the trended variant makes ``s0²`` a per-gene vector whose
``s0_squared`` metadata is a median over genes, which is a different quantity and would silently
change what the table means. Amendment 2 Change 1 measured the trended variant and did not select
it; the arm's default is untrended and this script refuses to run on anything else.

THE CODE PATH IS THE ARM'S, NOT A COPY
--------------------------------------
Each replicate goes through exactly what the pseudobulk arm does, in the arm's order, using the
same call convention as ``scripts/pb_calibration_probe.py::_aggregate``:

    synthetic.oracles.simulate(..., n_de=0)                  the pre-registered generative model
    pbcheck.methods.pseudobulk.build_pseudobulk(...)         decoupler donor aggregation + the
                                                             pre-registered thin-donor filter
    pbcheck.gene_universe.frozen_universe(...)               label-agnostic frozen universe (A3)
    pdata[:, universe].copy()                                restriction, once
    pbcheck.methods.moderated.ebayes_from_pdata(...)         log_cpm → wls_two_group → fit_f_dist

``frozen_universe`` is called with its defaults and **not** with ``min_size=MIN_UNIVERSE_SIZE``,
matching the probe; the realised universe size is recorded per replicate so a reader can check it
never came near the C5 floor of 200 (it does not: the smallest cell here simulates 1000 genes).

RESTRICT ONCE — AND WHY ``universe=`` IS NOT PASSED TO THE ARM
--------------------------------------------------------------
The matrix handed to ``ebayes_from_pdata`` is **already** restricted to ``universe``, in
``universe`` order, by the ``pdata[:, universe].copy()`` line above — the probe's structural
guarantee, restrict once and hand every consumer the same object. Passing ``universe=`` as well
would only make :func:`~pbcheck.methods.moderated.ebayes_from_pdata` redo the identical
restriction, and that line is **quadratic in the universe size**::

    keep = [g for g in universe if g in set(pdata.var_names)]     # moderated.py:336

``set(pdata.var_names)`` sits in the comprehension's condition, so it is rebuilt from the pandas
Index once per gene: G² element insertions. Measured on this machine at the 15 000-gene cells,
that single line costs **113.7 s** per call against 0.02 s for the whole rest of the fit — it, and
not the simulation, is what makes those cells expensive (simulate 1.4 s, ``build_pseudobulk``
0.4 s, ``frozen_universe`` 0.01 s, restriction 0.01 s). The same construction appears verbatim at
``pseudobulk.py:157`` and ``pb_calibration_probe.py:809``. **It is reported, not patched**: the
shipped arm's code is pinned by ``tests/test_moderated.py`` and by the amendments' numbers, and
fixing it is a separate change with its own review, not a side effect of a reproducer.

Skipping the redundant argument is therefore a saving of pure duplicated work, and it is checked
rather than asserted: ``--smoke`` runs
:func:`check_universe_argument_is_redundant`, which fits one cell both ways and requires
bit-identical p-values, log2 fold changes, gene order **and** the entire ``moderation`` dict. It
was also verified once at G = 15 000, where the amendment's heaviest cells live.

WHAT ELSE IS RECORDED, AND WHY
------------------------------
The amendment's table reads the reversal as **ordered by depth relative to a closed-form
threshold**, so the artifact has to carry the design quantities that ordering rests on, per
replicate and not merely per cell:

``L_tilde``      median over donors of the row sum of the **universe-restricted** pseudobulk
                 matrix — the library size ``moderated.log_cpm`` actually divides by, so no gene
                 outside the frozen universe leaks into it.
``median_cpm``   median over all (donor, gene) entries of ``1e6 · T / L`` on that same matrix.
``a_bar``        median over all (donor, gene) of ``CPM / (1 + CPM)``, the attenuation factor the
                 ``+1`` in ``log2(CPM + 1)`` imposes.
``prior_df_d0``  the realised ``d0``. ``null`` encodes ``+inf`` (complete pooling), which is a
                 legitimate estimate and is flagged by ``prior_is_complete_pooling`` beside it —
                 the repo's convention, because JSON has no infinity literal ``jq`` will parse.

Derived per cell: the mean / sd / SE of the mean of the quantity over seeds, ``ratio_to_sigma``,
``n_below_sigma`` (strictly below the true σ — the count that decides whether "upper bound" is even
approximately true), ``L_crit = 1e6/(2σ²)`` and ``L_tilde/L_crit``. Every aggregate is recomputable
from the per-seed block in the same artifact; nothing is reported that cannot be re-derived.

``L_crit`` is quoted for orientation only, exactly as the amendment quotes it: it is a median-gene
heuristic that drops the ``φ r C²`` term, whereas ``fit_f_dist``'s location is a log-scale central
value over the whole gene distribution. **No code reads it.**

SEEDS
-----
All ten cells run on the same set, ``--seeds 1-16`` by default. These are development seeds inside
the repository's disclosed development range [1, 999] and are disclosed as such: they are *not*
confirmatory seeds and no pre-registered claim may be read off them. What running all ten cells on
one named set buys is the thing the original table could not offer — the rows are comparable with
each other, and every row is backed by the same number of realisations.

COST
----
The four 15 000-gene cells dominate: 16 donors × 300 cells × 15 000 genes per replicate, ~1.5 GB
peak and tens of seconds each. Cells run cheapest-first so that a ``--max-minutes`` stop loses the
least informative work, and every completed replicate is appended to a ``replicates.jsonl`` ledger
in the output directory, keyed by cell id, seed **and a fingerprint of the cell's simulator
parameters**. A re-run reuses the ledger and computes only what is missing; a ledger entry whose
fingerprint disagrees is ignored rather than trusted, so an edited cell definition can never be
silently served from cache. If a run stops early, the affected cell reports ``is_complete: false``
with its realised ``n_seeds`` — a short row is visible, never quietly averaged as if it were full.

USAGE
-----
    python scripts/check_upper_bound_claim.py                       # all ten cells, seeds 1-16
    python scripts/check_upper_bound_claim.py --smoke               # 2 cheap cells, 2 seeds, ~20 s
    python scripts/check_upper_bound_claim.py --cells c06,c09 --seeds 1-8
    python scripts/check_upper_bound_claim.py --list-only

Writes ``upper_bound_check_<date>.json`` and a companion ``.csv`` of the per-cell summary into
``--out`` (default ``pilot/upper_bound_check/``), matching the ``pilot/testsel/`` convention of a
committed JSON+CSV pair beside the bulky per-run material.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata as _md
import json
import math
import platform
import subprocess
import sys
import time
import warnings
from datetime import date as _date
from pathlib import Path

import numpy as np
from anndata import ImplicitModificationWarning

# TARGETED, not blanket. This is the one warning ``pyproject.toml``'s ``filterwarnings`` block
# documents as benign and silences by message+category for the test suite: anndata fires it because
# ``synthetic/oracles.py`` builds its ``AnnData`` from an obs frame with a default integer index,
# before ``obs_names`` is set. It would otherwise print once per replicate, 160 times in a full run.
# It is NOT ``warnings.filterwarnings("ignore")`` — anything new out of a future numpy / pandas /
# scanpy stack must still surface here rather than be swallowed with it.
warnings.filterwarnings("ignore", message=r"Transforming to str index\.",
                        category=ImplicitModificationWarning)

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "synthetic"))

from oracles import simulate  # noqa: E402

from pbcheck.gene_universe import frozen_universe  # noqa: E402
from pbcheck.methods.moderated import ebayes_from_pdata  # noqa: E402
from pbcheck.methods.pseudobulk import build_pseudobulk  # noqa: E402

#: One line, carried into the artifact header so the file explains itself without this script.
QUANTITY = (
    "sqrt(s0^2) * ln 2, where s0^2 is the untrended empirical-Bayes prior variance fitted by the "
    "pseudobulk arm (pbcheck.methods.moderated.ebayes_from_pdata, moderation['s0_squared']) on "
    "log2(CPM+1) of the universe-restricted donor pseudobulk of a null oracle (n_de = 0); ln 2 "
    "puts it on the natural-log scale of the simulator's donor_sigma, so Amendment 3's claim is "
    "exactly that this number is >= donor_sigma."
)

#: Held fixed across every cell (Amendment 4 Part A: 8 v 8 donors, dispersion 0.2, simulator
#: defaults otherwise). A cell may override only what appears in :data:`CELLS`.
COMMON_SIM = {
    "n_donors_per_group": 8,
    "dispersion": 0.2,
    "n_de": 0,
    "mean_log_mu": 1.0,
    "mean_log_sigma": 1.2,
    "depth_log_sigma": 0.3,
}

#: The ten probe cells of Amendment 4 Part A's table. ``c06`` is ``gate_config.ORACLE_SIM``'s
#: geometry at its calibration sigma — the gate's own operating point, and the row the amendment
#: reads hardest.
CELLS: tuple[dict, ...] = (
    {"id": "c01", "n_genes": 1000, "n_cells_per_donor": 30, "mean_log_mu": 0.0,
     "donor_sigma": 0.35},
    {"id": "c02", "n_genes": 1000, "n_cells_per_donor": 30, "donor_sigma": 0.35},
    {"id": "c03", "n_genes": 1000, "n_cells_per_donor": 300, "mean_log_mu": 0.0,
     "donor_sigma": 0.35},
    {"id": "c04", "n_genes": 1000, "n_cells_per_donor": 300, "donor_sigma": 0.35},
    {"id": "c05", "n_genes": 1500, "n_cells_per_donor": 250, "donor_sigma": 0.35},
    {"id": "c06", "n_genes": 1500, "n_cells_per_donor": 250, "donor_sigma": 0.50,
     "note": "gate_config.ORACLE_SIM geometry at CALIBRATION_EVAL_SIGMA"},
    {"id": "c07", "n_genes": 15000, "n_cells_per_donor": 300, "donor_sigma": 0.20},
    {"id": "c08", "n_genes": 15000, "n_cells_per_donor": 100, "donor_sigma": 0.35},
    {"id": "c09", "n_genes": 15000, "n_cells_per_donor": 300, "donor_sigma": 0.35},
    {"id": "c10", "n_genes": 15000, "n_cells_per_donor": 300, "donor_sigma": 0.50},
)

#: ``--smoke``: the two cheapest cells on two seeds. Enough to exercise every code path end to end
#: (simulate → aggregate → universe → eBayes → aggregate → artifact) in well under a minute.
SMOKE_CELLS = "c01,c02"
SMOKE_SEEDS = "1-2"

LN2 = math.log(2.0)

#: Source files the measurement depends on. Hashed into the header: a commit alone does not pin the
#: numbers when the working tree may carry uncommitted edits.
HASHED_SOURCES = (
    "scripts/check_upper_bound_claim.py",
    "src/pbcheck/methods/moderated.py",
    "src/pbcheck/methods/pseudobulk.py",
    "src/pbcheck/gene_universe.py",
    "synthetic/oracles.py",
)


def sim_params(cell: dict) -> dict:
    """The full ``simulate`` keyword set for a cell: the common block with the cell's overrides."""
    p = dict(COMMON_SIM)
    p.update({k: v for k, v in cell.items() if k not in ("id", "note")})
    return p


def cell_fingerprint(cell: dict) -> str:
    """Short hash of a cell's simulator parameters, so a stale ledger entry cannot be reused."""
    blob = json.dumps(sim_params(cell), sort_keys=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:12]


def cell_cost(cell: dict) -> int:
    """Crude cost proxy (genes × cells per donor) used only to order the run cheapest-first."""
    return int(cell["n_genes"]) * int(cell["n_cells_per_donor"])


def describe(cell: dict) -> str:
    """The table's left-hand column: how the amendment names the cell."""
    p = sim_params(cell)
    bits = [f"{p['n_genes']} genes", f"{p['n_cells_per_donor']} cells"]
    if p["mean_log_mu"] != COMMON_SIM["mean_log_mu"]:
        bits.append(f"mean_log_mu = {p['mean_log_mu']:g}")
    text = ", ".join(bits)
    return f"{text} [{cell['note']}]" if "note" in cell else text


def parse_int_list(spec: str, *, what: str) -> list[int]:
    """Parse ``"1-16"``, ``"1,2,5"`` or ``"1-4,9,12-14"`` into a sorted list of unique ints."""
    out: set[int] = set()
    for chunk in (c.strip() for c in spec.split(",")):
        if not chunk:
            continue
        if "-" in chunk.lstrip("-"):
            lo_s, _, hi_s = chunk.partition("-")
            try:
                lo, hi = int(lo_s), int(hi_s)
            except ValueError:
                raise argparse.ArgumentTypeError(f"bad {what} range {chunk!r}") from None
            if hi < lo:
                raise argparse.ArgumentTypeError(f"{what} range {chunk!r} runs backwards")
            out.update(range(lo, hi + 1))
        else:
            try:
                out.add(int(chunk))
            except ValueError:
                raise argparse.ArgumentTypeError(f"bad {what} value {chunk!r}") from None
    if not out:
        raise argparse.ArgumentTypeError(f"no {what} given")
    return sorted(out)


def parse_cells(spec: str) -> list[dict]:
    """Resolve a comma list of cell ids against :data:`CELLS`, preserving the table's order."""
    known = {c["id"]: c for c in CELLS}
    wanted = [s.strip() for s in spec.split(",") if s.strip()]
    unknown = [w for w in wanted if w not in known]
    if unknown:
        raise argparse.ArgumentTypeError(
            f"unknown cell id(s) {unknown}; valid ids are {sorted(known)}"
        )
    keep = set(wanted)
    return [c for c in CELLS if c["id"] in keep]


# ===========================================================================
# One replicate: the arm's own path, and the design quantities beside it.
# ===========================================================================

def measure_replicate(cell: dict, seed: int) -> dict:
    """Run one null-oracle replicate of ``cell`` at ``seed`` and read the quantity off the arm."""
    t0 = time.perf_counter()
    params = sim_params(cell)

    oracle = simulate(seed=seed, **params)
    pdata = build_pseudobulk(oracle.adata)
    del oracle                                   # ~0.3-0.9 GB at the 15 000-gene cells

    universe = frozen_universe(pdata)            # engine defaults, exactly as the probe calls it
    pdata_u = pdata[:, universe].copy()          # restrict ONCE; see the module docstring
    del pdata

    res = ebayes_from_pdata(pdata_u, trend=False)
    mod = res.moderation
    if mod["s0_squared_is_trended"]:
        raise RuntimeError(
            "s0_squared came back trended: it would be a per-gene median, not the scalar prior "
            "variance this measurement is about. trend must stay False."
        )
    s0_squared = float(mod["s0_squared"])
    if not (s0_squared > 0 and math.isfinite(s0_squared)):
        raise RuntimeError(f"non-positive or non-finite s0_squared {s0_squared!r} at {cell['id']}"
                           f" seed {seed}")

    # The same matrix and the same library definition moderated.log_cpm uses: row sums of the
    # universe-restricted matrix, floored at 1 so an empty donor cannot divide by zero.
    X = np.asarray(pdata_u.X, dtype=float)
    lib = X.sum(axis=1)
    cpm = X / np.maximum(lib, 1.0)[:, None] * 1e6
    a = cpm / (1.0 + cpm)

    return {
        "cell_id": cell["id"],
        "seed": int(seed),
        "fingerprint": cell_fingerprint(cell),
        "s0_squared": s0_squared,
        "sqrt_s0_squared_times_ln2": math.sqrt(s0_squared) * LN2,
        "prior_df_d0": mod["prior_df_d0"],          # None encodes +inf; see the flag below
        "prior_is_complete_pooling": bool(mod["prior_is_complete_pooling"]),
        "residual_df_d": int(mod["residual_df_d"]),
        "L_tilde": float(np.median(lib)),
        "median_cpm": float(np.median(cpm)),
        "a_bar": float(np.median(a)),
        "universe_size": int(len(universe)),
        "n_donors_retained": int(pdata_u.n_obs),
        "n_genes_nonfinite_s2": int(mod["n_genes_nonfinite_s2"]),
        "seconds": float(time.perf_counter() - t0),
    }


def check_universe_argument_is_redundant(cell: dict, seed: int) -> dict:
    """Prove that not passing ``universe=`` to the arm changes nothing, on a real replicate.

    :func:`measure_replicate` restricts the pseudobulk matrix to the frozen universe itself and
    then calls ``ebayes_from_pdata`` without ``universe=``. That is only legitimate if the two
    calls are the *same* call on an already-restricted matrix, so this fits one cell both ways and
    demands bit-identical p-values, log2 fold changes, gene order and moderation metadata — no
    tolerance. Run from ``--smoke``, where the universe is small enough that the quadratic
    ``universe=`` branch costs well under a second.

    Raises
    ------
    RuntimeError
        On any difference at all, which would mean the fast path is not the arm's path.
    """
    params = sim_params(cell)
    oracle = simulate(seed=seed, **params)
    pdata = build_pseudobulk(oracle.adata)
    del oracle
    universe = frozen_universe(pdata)
    pdata_u = pdata[:, universe].copy()
    del pdata

    t0 = time.perf_counter()
    ref = ebayes_from_pdata(pdata_u, universe=universe, trend=False)
    t_with = time.perf_counter() - t0
    t0 = time.perf_counter()
    got = ebayes_from_pdata(pdata_u, trend=False)
    t_without = time.perf_counter() - t0

    problems = []
    if list(ref.table.index) != list(got.table.index):
        problems.append("gene order")
    for col in ("pval", "log2fc"):
        if not np.array_equal(ref.table[col].to_numpy(dtype=float),
                              got.table[col].to_numpy(dtype=float), equal_nan=True):
            problems.append(col)
    differing = [k for k in ref.moderation if ref.moderation[k] != got.moderation[k]]
    problems.extend(f"moderation[{k}]" for k in differing)
    if problems:
        raise RuntimeError(
            "ebayes_from_pdata(pdata_u) is NOT identical to ebayes_from_pdata(pdata_u, "
            f"universe=universe) on {cell['id']} seed {seed}: {problems}. The fast path is not "
            "the arm's path and this measurement would be invalid."
        )
    return {"cell_id": cell["id"], "seed": int(seed), "universe_size": len(universe),
            "seconds_with_universe_arg": t_with, "seconds_without": t_without}


# ===========================================================================
# Ledger: append-only, so a killed run resumes instead of restarting.
# ===========================================================================

def load_ledger(path: Path) -> dict[tuple[str, int, str], dict]:
    """Read completed replicates, keyed by (cell id, seed, fingerprint). Bad lines are skipped."""
    done: dict[tuple[str, int, str], dict] = {}
    if not path.exists():
        return done
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                key = (rec["cell_id"], int(rec["seed"]), rec["fingerprint"])
            except (json.JSONDecodeError, KeyError, TypeError, ValueError):
                continue
            done[key] = rec
    return done


def append_ledger(path: Path, rec: dict) -> None:
    """Append one replicate and flush, so the ledger survives a kill between replicates."""
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        fh.flush()


# ===========================================================================
# Per-cell aggregation. Every number here is recomputable from `replicates`.
# ===========================================================================

def summarise_cell(cell: dict, recs: list[dict], seeds_requested: list[int]) -> dict:
    """Aggregate a cell's replicates into the table's row, flagging a short row as incomplete."""
    recs = sorted(recs, key=lambda r: r["seed"])
    sigma = float(sim_params(cell)["donor_sigma"])
    v = np.asarray([r["sqrt_s0_squared_times_ln2"] for r in recs], dtype=float)
    n = int(v.size)
    l_tilde = float(np.median([r["L_tilde"] for r in recs]))
    l_crit = 1e6 / (2.0 * sigma ** 2)
    d0_finite = [r["prior_df_d0"] for r in recs if r["prior_df_d0"] is not None]

    mean = float(v.mean())
    sd = float(v.std(ddof=1)) if n > 1 else None
    sem = float(v.std(ddof=1) / math.sqrt(n)) if n > 1 else None

    return {
        "cell_id": cell["id"],
        "description": describe(cell),
        "sim_params": sim_params(cell),
        "fingerprint": cell_fingerprint(cell),
        "donor_sigma": sigma,
        "n_seeds": n,
        "n_seeds_requested": len(seeds_requested),
        "is_complete": n == len(seeds_requested),
        "seeds": [int(r["seed"]) for r in recs],

        "mean_sqrt_s0_squared_times_ln2": mean,
        "sd_sqrt_s0_squared_times_ln2": sd,
        "sem_sqrt_s0_squared_times_ln2": sem,
        "min_sqrt_s0_squared_times_ln2": float(v.min()),
        "max_sqrt_s0_squared_times_ln2": float(v.max()),
        "ratio_to_sigma": mean / sigma,
        "n_below_sigma": int(np.sum(v < sigma)),        # strictly below the truth
        "frac_below_sigma": float(np.mean(v < sigma)),

        "L_tilde_median_over_seeds": l_tilde,
        "L_crit": l_crit,
        "L_tilde_over_L_crit": l_tilde / l_crit,
        "median_cpm_median_over_seeds": float(np.median([r["median_cpm"] for r in recs])),
        "a_bar_median_over_seeds": float(np.median([r["a_bar"] for r in recs])),

        # Median over the replicates whose d0 came back FINITE. Complete pooling (d0 = +inf) is a
        # legitimate estimate and is counted beside it rather than folded in as a large number,
        # which would understate it, or dropped silently, which would hide it.
        "median_prior_df_d0_finite_only": float(np.median(d0_finite)) if d0_finite else None,
        "n_replicates_complete_pooling": int(sum(r["prior_is_complete_pooling"] for r in recs)),
        "median_universe_size": float(np.median([r["universe_size"] for r in recs])),
        "seconds_total": float(sum(r["seconds"] for r in recs)),

        "replicates": recs,
    }


# ===========================================================================
# Provenance
# ===========================================================================

def git_commit() -> str | None:
    """``git rev-parse HEAD``, read-only. ``None`` if git is unavailable or this is not a repo."""
    try:
        p = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                           text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return p.stdout.strip() if p.returncode == 0 else None


def source_hashes() -> dict:
    """SHA-256 of every file the measurement runs through — stronger than the commit alone."""
    out = {}
    for rel in HASHED_SOURCES:
        f = REPO / rel
        out[rel] = (hashlib.sha256(f.read_bytes()).hexdigest() if f.exists() else None)
    return out


def versions() -> dict:
    v = {"python": platform.python_version(), "platform": platform.platform()}
    for pkg in ("pbcheck", "numpy", "scipy", "anndata", "pandas", "decoupler", "scanpy",
                "statsmodels", "pydeseq2"):
        try:
            v[pkg] = _md.version(pkg)
        except _md.PackageNotFoundError:
            v[pkg] = "unavailable"
    return v


# ===========================================================================

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Development seeds, disclosed as such: nothing here is confirmatory evidence.",
    )
    ap.add_argument("--seeds", default="1-16",
                    help="seed set, identical for every cell. Ranges and comma lists: "
                         "'1-16' (default), '1,2,5', '1-4,9'. The repo's disclosed development "
                         "range is [1, 999].")
    ap.add_argument("--cells", default=",".join(c["id"] for c in CELLS),
                    help="comma list of cell ids (default: all ten of the amendment's table)")
    ap.add_argument("--out", type=Path, default=REPO / "pilot" / "upper_bound_check",
                    help="artifact directory; the JSON, the CSV and the resume ledger land here")
    ap.add_argument("--smoke", action="store_true",
                    help=f"fast end-to-end check: cells {SMOKE_CELLS} on seeds {SMOKE_SEEDS}, "
                         "written to a '_smoke'-suffixed artifact so it cannot be mistaken for, "
                         "or overwrite, the real run. Explicit --cells/--seeds still win.")
    ap.add_argument("--max-minutes", type=float, default=None,
                    help="stop cleanly before starting a replicate once this wall clock is "
                         "exceeded. Cells run cheapest-first, completed replicates are kept in "
                         "the ledger, and any short cell is marked is_complete: false rather "
                         "than averaged as if it were full.")
    ap.add_argument("--date-stamp", default=_date.today().isoformat(),
                    help="date in the artifact filename (default: today, ISO)")
    ap.add_argument("--list-only", action="store_true",
                    help="print the run plan (cells, seeds, replicate count) and exit")
    a = ap.parse_args(argv)

    # --smoke only supplies DEFAULTS; an explicit --cells / --seeds on the same command line wins,
    # so a smoke invocation can still be narrowed or widened without the flag silently overriding.
    given = list(sys.argv[1:] if argv is None else [str(x) for x in argv])
    if a.smoke:
        if not any(g.startswith("--cells") for g in given):
            a.cells = SMOKE_CELLS
        if not any(g.startswith("--seeds") for g in given):
            a.seeds = SMOKE_SEEDS

    seeds = parse_int_list(a.seeds, what="seed")
    cells = parse_cells(a.cells)
    stem = f"upper_bound_check{'_smoke' if a.smoke else ''}_{a.date_stamp}"

    if a.list_only:
        for c in sorted(cells, key=cell_cost):
            p = sim_params(c)
            print(f"{c['id']}  {describe(c):<52} sigma={p['donor_sigma']:.2f} "
                  f"cost={cell_cost(c):>9,d}")
        print(f"\n{len(cells)} cells x {len(seeds)} seeds = {len(cells) * len(seeds)} replicates; "
              f"seeds = {seeds}", flush=True)
        return 0

    a.out.mkdir(parents=True, exist_ok=True)
    ledger_path = a.out / "replicates.jsonl"
    ledger = load_ledger(ledger_path)

    total = len(cells) * len(seeds)
    print(f"[ubc] {len(cells)} cells x {len(seeds)} seeds = {total} replicates; "
          f"{len(ledger)} in the ledger; out={a.out}", flush=True)
    print(f"[ubc] quantity: {QUANTITY}", flush=True)

    equivalence = None
    if a.smoke:
        equivalence = check_universe_argument_is_redundant(cells[0], seeds[0])
        print(f"[ubc] universe= argument is redundant on an already-restricted matrix: "
              f"bit-identical on {equivalence['cell_id']} seed {equivalence['seed']} "
              f"(G={equivalence['universe_size']}, "
              f"{equivalence['seconds_with_universe_arg']:.2f}s with it vs "
              f"{equivalence['seconds_without']:.2f}s without)", flush=True)

    t0 = time.time()
    deadline = t0 + a.max_minutes * 60 if a.max_minutes else None
    computed = reused = 0
    stopped_early = False

    for cell in sorted(cells, key=cell_cost):
        fp = cell_fingerprint(cell)
        for seed in seeds:
            key = (cell["id"], seed, fp)
            if key in ledger:
                reused += 1
                continue
            if deadline is not None and time.time() > deadline:
                print(f"[ubc] budget of {a.max_minutes} min exhausted before {cell['id']} "
                      f"seed {seed} - stopping cleanly", flush=True)
                stopped_early = True
                break
            rec = measure_replicate(cell, seed)
            ledger[key] = rec
            append_ledger(ledger_path, rec)
            computed += 1
            print(f"[ubc] {cell['id']} seed {seed:>3}  "
                  f"sqrt(s0^2)*ln2={rec['sqrt_s0_squared_times_ln2']:.4f}  "
                  f"(sigma={sim_params(cell)['donor_sigma']:.2f})  "
                  f"d0={rec['prior_df_d0'] if rec['prior_df_d0'] is None else round(rec['prior_df_d0'], 1)}  "
                  f"G={rec['universe_size']}  {rec['seconds']:.1f}s  "
                  f"[{computed + reused}/{total}]", flush=True)
        if stopped_early:
            break

    summaries = []
    for cell in cells:
        fp = cell_fingerprint(cell)
        recs = [ledger[(cell["id"], s, fp)] for s in seeds if (cell["id"], s, fp) in ledger]
        if not recs:
            print(f"[ubc] {cell['id']}: no replicates - omitted from the artifact", flush=True)
            continue
        summaries.append(summarise_cell(cell, recs, seeds))

    elapsed = time.time() - t0
    artifact = {
        "header": {
            "script": f"scripts/{Path(__file__).name}",
            "argv": sys.argv[1:] if argv is None else [str(x) for x in argv],
            "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "git_commit": git_commit(),
            "source_sha256": source_hashes(),
            "quantity": QUANTITY,
            "seeds": seeds,
            "seed_provenance": "development seeds, disclosed; the repo's disclosed development "
                               "range is [1, 999]. NOT confirmatory seeds.",
            "reproduces": "docs/AMENDMENTS.md, Amendment 4 Part A, 'The measurement behind "
                          "Correction 1'. The amendment's own table used unnamed development "
                          "seeds and a per-row seed count of 16/8/1; this runs every cell on one "
                          "named set, so exact agreement with it is not expected.",
            "design_columns_note": "L_tilde, median CPM and a_bar are reported as MEDIANS OVER "
                                   "SEEDS, and L_tilde_over_L_crit is built from the median "
                                   "L_tilde. The amendment's table carries the single-realisation "
                                   "(first-seed) values of those three columns, which run 1-3% "
                                   "below the medians; per-seed values are in `replicates` so "
                                   "either convention can be re-derived.",
            "arm_code_path": "synthetic.oracles.simulate -> pseudobulk.build_pseudobulk -> "
                             "gene_universe.frozen_universe -> pdata[:, universe].copy() -> "
                             "moderated.ebayes_from_pdata (log_cpm -> wls_two_group -> "
                             "fit_f_dist), trend=False",
            "arm_call_note": "ebayes_from_pdata is called WITHOUT universe=, because the matrix "
                             "handed to it is already restricted to the frozen universe in "
                             "universe order. Passing it would only repeat that restriction "
                             "through moderated.py:336, whose 'set(pdata.var_names)' inside the "
                             "comprehension condition is quadratic in G (113.7 s per call at "
                             "G = 15000 against 0.02 s for the rest of the fit). The two calls "
                             "were verified bit-identical in p-values, log2fc, gene order and the "
                             "whole moderation dict; --smoke re-checks it on every run.",
            "equivalence_check": equivalence,
            "common_sim_params": dict(COMMON_SIM),
            "cell_sim_params": {c["id"]: sim_params(c) for c in cells},
            "L_crit_definition": "1e6 / (2 * donor_sigma^2); a median-gene heuristic quoted for "
                                 "orientation only. No code reads it.",
            "n_replicates_computed_this_run": computed,
            "n_replicates_reused_from_ledger": reused,
            "run_complete": (not stopped_early) and all(s["is_complete"] for s in summaries)
                            and len(summaries) == len(cells),
            "stopped_on_time_budget": stopped_early,
            "wall_clock_seconds": float(elapsed),
            "versions": versions(),
        },
        "cells": summaries,
    }

    json_path = a.out / f"{stem}.json"
    json_path.write_text(json.dumps(artifact, indent=1, ensure_ascii=False), encoding="utf-8")

    csv_fields = [
        "cell_id", "description", "donor_sigma", "n_genes", "n_cells_per_donor", "mean_log_mu",
        "n_seeds", "is_complete", "mean_sqrt_s0_squared_times_ln2",
        "sd_sqrt_s0_squared_times_ln2", "sem_sqrt_s0_squared_times_ln2", "ratio_to_sigma",
        "n_below_sigma", "L_tilde_median_over_seeds", "L_crit", "L_tilde_over_L_crit",
        "median_cpm_median_over_seeds", "a_bar_median_over_seeds",
        "median_prior_df_d0_finite_only", "n_replicates_complete_pooling",
        "median_universe_size",
    ]
    csv_path = a.out / f"{stem}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=csv_fields, extrasaction="ignore")
        w.writeheader()
        for s in summaries:
            row = dict(s)
            row.update({k: s["sim_params"][k]
                        for k in ("n_genes", "n_cells_per_donor", "mean_log_mu")})
            w.writerow(row)

    print(f"\n[ubc] {'cell':<5} {'sigma':>5} {'mean':>8} {'ratio':>7} {'below':>7} "
          f"{'L~/Lcrit':>9}", flush=True)
    for s in sorted(summaries, key=lambda r: r["L_tilde_over_L_crit"]):
        print(f"[ubc] {s['cell_id']:<5} {s['donor_sigma']:>5.2f} "
              f"{s['mean_sqrt_s0_squared_times_ln2']:>8.4f} {s['ratio_to_sigma']:>7.3f} "
              f"{s['n_below_sigma']:>3d}/{s['n_seeds']:<3d} {s['L_tilde_over_L_crit']:>9.3f}"
              + ("" if s["is_complete"] else "   INCOMPLETE"), flush=True)
    print(f"[ubc] computed={computed} reused={reused} elapsed={elapsed / 60:.1f} min", flush=True)
    print(f"[ubc] artifact -> {json_path}", flush=True)
    print(f"[ubc] summary  -> {csv_path}", flush=True)
    return 0 if artifact["header"]["run_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
