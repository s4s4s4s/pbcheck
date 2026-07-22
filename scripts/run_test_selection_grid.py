"""Phase 0 test-selection grid — the long, deterministic run.

Answers one question: **as gene-level variance heterogeneity rises, do the variance-borrowing
tests (limma-style eBayes / limma-trend / voom) become anti-conservative the way DESeq2 does,
while the unshrunken donor-level tests stay calibrated?** And, among whatever remains
calibrated, which test has the most power at the pre-registered oracle?

This is a driver, not a method: every number comes from ``pb_calibration_probe.py``, which is
the validated instrument. The driver only decides which cells to run, in what order, and
aggregates the results.

Design notes that matter for reading the output:

* **Calibration is measured on FRESH nulls only.** For donor-level tests the donor-permutation
  null is calibrated by construction (donor profiles are fixed, only labels move), so it is not
  evidence. ``--null-mode fresh`` simulates an independent dataset per replicate and uses the
  real labels.
* **The primary axis is ``sigma_het``**, not ``sigma_donor``. With every gene sharing one donor
  variance the true prior degrees of freedom are infinite — precisely the regime empirical Bayes
  was built for, in which moderation cannot fail. A calibration pass there proves nothing.
* **``directnb`` carries the high-replicate work.** It draws donor-level counts straight from the
  negative-binomial model with no cell layer, so it is far cheaper per replicate *and* it is not
  our own simulator, which makes it the stronger evidence. The cell-level arms get fewer
  replicates for breadth.
* **Power is replicated.** A single positive-oracle dataset gives a sensitivity with an
  across-dataset SD near 0.05, which is hopeless against a binding 0.60 gate, so every cell runs
  ``--n-power-reps`` and reports the mean and SD.

Resumable: a cell whose JSON already exists is skipped, so the run can be killed and restarted.
Time-budgeted: cells run in priority order and the driver stops cleanly at ``--max-hours``, so a
partial run is still a usable result rather than a loss.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PROBE = REPO / "scripts" / "pb_calibration_probe.py"

# The unshrunken baselines and the variance-borrowing candidates. wald is the known-failing
# reference and is handled separately because it costs ~2.8 s per fit rather than ~5 ms.
UNSHRUNKEN = ("ttest", "pooled_t", "wilcoxon")
MODERATED = ("ebayes", "ebayes_trend", "voom")
CHEAP = UNSHRUNKEN + MODERATED

SIGMA_HET = (0.0, 0.4, 0.8, 1.2)
N_GENES = 1500
DISPERSION = 0.2
CELLS_PER_DONOR = 250
POWER_REPS = 25
SEED = 7


def cell_key(c: dict) -> str:
    return (
        f"{c['tier']}__{c['gen_arm']}__{c['test']}__het{c['sigma_het']}"
        f"__sd{c['sigma_donor']}__n{c['n_donors']}__r{c['n_perm']}"
    )


def build_cells() -> list[dict]:
    """Cells in priority order: the decisive comparison first, breadth afterwards.

    Ordering is the whole safety net for a time-budgeted run — if the budget runs out, what is
    missing must be the least informative work, not a random half of the grid.
    """
    cells: list[dict] = []

    def add(tier, gen_arm, tests, hets, sds, ns, n_perm):
        for h in hets:
            for sd in sds:
                for n in ns:
                    for t in tests:
                        cells.append(dict(
                            tier=tier, gen_arm=gen_arm, test=t, sigma_het=h,
                            sigma_donor=sd, n_donors=n, n_perm=n_perm,
                        ))

    # Tier A — the decisive cell. One (sigma_donor, n), the full heterogeneity axis, every cheap
    # test, on the arm that is not our simulator, with enough replicates that a 0.05-vs-0.08
    # difference is actually resolvable (3000 reps => SE ~0.004).
    add("A", "directnb", CHEAP, SIGMA_HET, (0.5,), (8,), 3000)

    # Tier B — does the tier-A answer hold as donor count and donor variance move?
    add("B", "directnb", CHEAP, SIGMA_HET, (0.35, 0.7), (8,), 1500)
    add("B", "directnb", CHEAP, SIGMA_HET, (0.5,), (4, 12), 1500)

    # Tier C — does it hold on the cell-level arms, including our own pre-registered simulator?
    add("C", "lognormal", CHEAP, (0.0, 0.8), (0.5,), (8,), 400)
    add("C", "gamma", CHEAP, (0.0, 0.8), (0.5,), (8,), 400)

    # Tier D — the failing reference, so the grid contains its own positive control for
    # anti-conservativity. Few cells, few replicates: wald is ~600x more expensive per fit.
    add("D", "directnb", ("wald",), (0.0, 0.8), (0.5,), (8,), 150)

    return cells


def run_cell(c: dict, out_dir: Path, python: str, n_cpus: int) -> tuple[bool, str]:
    out = out_dir / f"{cell_key(c)}.json"
    if out.exists():
        return True, "cached"
    cmd = [
        python, str(PROBE),
        "--sigma-donor", str(c["sigma_donor"]),
        "--n-donors", str(c["n_donors"]),
        "--test", c["test"],
        "--n-genes", str(N_GENES),
        "--n-cells-per-donor", str(CELLS_PER_DONOR),
        "--dispersion", str(DISPERSION),
        "--n-perm", str(c["n_perm"]),
        "--seed", str(SEED),
        "--null-mode", "fresh",
        "--gen-arm", c["gen_arm"],
        "--sigma-het", str(c["sigma_het"]),
        "--n-power-reps", str(POWER_REPS),
        "--n-cpus", str(n_cpus),
        "--out", str(out),
    ]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except subprocess.TimeoutExpired:
        return False, "timeout after 2h"
    if p.returncode != 0 or not out.exists():
        return False, (p.stderr or p.stdout or "no output")[-400:].replace("\n", " ")
    return True, "ok"


SUMMARY_FIELDS = [
    "tier", "gen_arm", "test", "sigma_het", "sigma_donor", "n_donors", "n_perm",
    "lambda_pb", "fp_rate", "fp_binom_p", "mean_rejections",
    "power", "power_sd", "n_power_reps", "prior_df_d0", "shrinkage", "seconds",
]


def summarise(out_dir: Path, summary_dir: Path) -> int:
    rows = []
    for f in sorted(out_dir.glob("*.json")):
        try:
            o = json.load(open(f, encoding="utf-8"))
        except Exception:
            continue
        if "lambda_pb" not in o:
            continue
        d = o.get("test_diagnostics") or {}
        rows.append({
            "tier": f.name.split("__")[0],
            "gen_arm": o.get("gen_arm"),
            "test": o.get("test"),
            "sigma_het": o.get("sigma_het"),
            "sigma_donor": o.get("sigma_donor"),
            "n_donors": o.get("n_donors_per_group"),
            "n_perm": o.get("n_perm"),
            "lambda_pb": o.get("lambda_pb"),
            "fp_rate": o.get("fp_rate_perm_null"),
            "fp_binom_p": o.get("fp_rate_binomial_p_greater"),
            "mean_rejections": o.get("mean_rejections_perm_null"),
            "power": o.get("power_at_pre_registered"),
            "power_sd": o.get("power_sd"),
            "n_power_reps": o.get("n_power_reps"),
            "prior_df_d0": d.get("prior_df_d0"),
            "shrinkage": d.get("shrinkage_factor_d0_over_d0_plus_d"),
            "seconds": o.get("seconds_elapsed"),
        })
    summary_dir.mkdir(parents=True, exist_ok=True)
    with open(summary_dir / "summary.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    json.dump(rows, open(summary_dir / "summary.json", "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    return len(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="per-cell JSONs land here; existing files are treated as done")
    ap.add_argument("--summary-dir", type=Path, default=REPO / "pilot" / "testsel",
                    help="summary.csv / summary.json (small, meant to be committed)")
    ap.add_argument("--max-hours", type=float, default=10.0)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--n-cpus", type=int, default=4)
    ap.add_argument("--list-only", action="store_true")
    a = ap.parse_args()

    cells = build_cells()
    a.out_dir.mkdir(parents=True, exist_ok=True)

    if a.list_only:
        for c in cells:
            print(cell_key(c))
        print(f"\n{len(cells)} cells", flush=True)
        return 0

    t0 = time.time()
    deadline = t0 + a.max_hours * 3600
    done = failed = skipped = 0

    print(f"[grid] {len(cells)} cells, budget {a.max_hours} h, out={a.out_dir}", flush=True)

    for i, c in enumerate(cells, 1):
        if time.time() > deadline:
            print(f"[grid] budget exhausted at cell {i}/{len(cells)} — stopping cleanly", flush=True)
            break
        t = time.time()
        ok, msg = run_cell(c, a.out_dir, a.python, a.n_cpus)
        if msg == "cached":
            skipped += 1
        elif ok:
            done += 1
            print(f"[{i}/{len(cells)}] {cell_key(c)}  {time.time()-t:.1f}s", flush=True)
        else:
            failed += 1
            print(f"[{i}/{len(cells)}] FAILED {cell_key(c)}: {msg}", flush=True)
        if (done + failed) % 10 == 0:
            summarise(a.out_dir, a.summary_dir)

    n = summarise(a.out_dir, a.summary_dir)
    el = (time.time() - t0) / 3600
    print(f"[grid] done={done} cached={skipped} failed={failed} rows={n} elapsed={el:.2f}h", flush=True)
    print(f"[grid] summary -> {a.summary_dir/'summary.csv'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
