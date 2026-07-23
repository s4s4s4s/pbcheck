"""Read the test-selection grid summary and state, in numbers, which test the pseudobulk arm
should use — and whether the variance-shrinkage hypothesis held.

Consumes ``pilot/testsel/summary.json`` (written by ``run_test_selection_grid.py``). Produces no
verdict of its own beyond what the numbers force: it tabulates calibration and power, applies the
pre-registered thresholds, and prints the one comparison the whole run was built to make.

The hypothesis under test: variance-borrowing methods (eBayes, limma-trend, voom) become
anti-conservative as gene-level variance heterogeneity rises — i.e. as the realised prior degrees
of freedom d0 fall — while the unshrunken donor-level tests (Welch t, pooled Student, Wilcoxon)
stay calibrated. The decisive axis is therefore ``sigma_het``, read via the realised ``prior_df_d0``.

Nothing here is a scientific finding: every number is from synthetic data. The binding real check
is the Mathys 2019 anchor, which has not run.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ALPHA = 0.05
LAMBDA_BAND = (0.9, 1.1)
POWER_TARGET = 0.60

UNSHRUNKEN = {"ttest", "pooled_t", "wilcoxon"}
MODERATED = {"ebayes", "ebayes_trend", "voom"}


def _f(x, d="  .  "):
    return d if x is None else (f"{x:.3f}" if isinstance(x, float) else str(x))


def calibrated(row) -> bool | None:
    """A cell is calibrated if its permutation/fresh-null FP rate is not above alpha AND lambda is
    in band. ``fp_binom_p`` is P(observed >= FP | true rate = alpha); a small value means the arm
    is significantly anti-conservative. Returns None when the cell lacks the fields."""
    fp_p = row.get("fp_binom_p")
    lam = row.get("lambda_pb")
    if fp_p is None or lam is None:
        return None
    # Not-anti-conservative: we cannot reject "rate <= alpha" at 0.05 in the liberal direction.
    fp_ok = fp_p > 0.05
    lam_ok = LAMBDA_BAND[0] <= lam <= LAMBDA_BAND[1]
    return bool(fp_ok and lam_ok)


def load(summary: Path) -> list[dict]:
    rows = json.load(open(summary, encoding="utf-8"))
    # Normalise numeric fields that a CSV round-trip may have stringified.
    for r in rows:
        for k in ("sigma_het", "sigma_donor", "lambda_pb", "fp_rate", "fp_binom_p",
                  "power", "power_sd", "prior_df_d0", "shrinkage"):
            v = r.get(k)
            if isinstance(v, str) and v not in ("", "None", "inf", "Infinity"):
                try:
                    r[k] = float(v)
                except ValueError:
                    pass
    return rows


def shrinkage_table(rows: list[dict]) -> str:
    """FP rate vs sigma_het at the decisive cell (directnb, sigma_donor=0.5, n=8), per test.

    This is the table that answers the hypothesis: read each moderated test's FP rate as the
    realised d0 falls, against the unshrunken baselines on the same cells.
    """
    sel = [r for r in rows
           if r.get("gen_arm") == "directnb" and r.get("sigma_donor") == 0.5
           and r.get("n_donors") == 8]
    hets = sorted({r["sigma_het"] for r in sel if r.get("sigma_het") is not None})
    tests = [t for t in ("ttest", "pooled_t", "wilcoxon", "ebayes", "ebayes_trend", "voom")
             if any(r.get("test") == t for r in sel)]

    by = {(r["test"], r["sigma_het"]): r for r in sel}
    out = ["", "## Shrinkage hypothesis — FP rate (fresh null, directnb, sigma_donor=0.5, 8v8)",
           "  test          " + "".join(f"  het={h:<4}" for h in hets) + "   d0@max_het"]
    for t in tests:
        cells = []
        for h in hets:
            r = by.get((t, h))
            cells.append(_f(r.get("fp_rate")) if r else "  .  ")
        rmax = by.get((t, hets[-1])) if hets else None
        d0 = _f(rmax.get("prior_df_d0")) if rmax else "  .  "
        tag = "SHRINK" if t in MODERATED else "flat  "
        out.append(f"  {t:12} {tag}" + "".join(f"  {c:>7}" for c in cells) + f"   {d0:>8}")
    out.append(f"  (target FP <= {ALPHA}; a value climbing with het = shrinkage breaks under heterogeneity)")
    return "\n".join(out)


def power_frontier(rows: list[dict]) -> str:
    """Smallest n_donors for power >= 0.60 at the pre-registered oracle, per calibrated test,
    per sigma_donor — computed only on CALIBRATED cells (an uncalibrated test's power is void)."""
    out = ["", "## Power frontier — min donors/group for power >= 0.60 (calibrated cells only)",
           "  test          sigma=0.35   sigma=0.5   sigma=0.7"]
    # Use the het=0 cells for the frontier (matched to the pre-registered homoscedastic oracle);
    # calibration is judged on the full het axis separately.
    for t in ("ttest", "pooled_t", "wilcoxon", "ebayes", "ebayes_trend", "voom"):
        cols = []
        for sd in (0.35, 0.5, 0.7):
            reached = None
            cand = sorted([r for r in rows if r.get("test") == t and r.get("sigma_donor") == sd
                           and r.get("sigma_het") == 0.0 and r.get("n_donors") is not None],
                          key=lambda r: r["n_donors"])
            for r in cand:
                if calibrated(r) and (r.get("power") or 0) >= POWER_TARGET:
                    reached = r["n_donors"]
                    break
            cols.append(str(reached) if reached else ">grid")
        out.append(f"  {t:12}   " + "   ".join(f"{c:>9}" for c in cols))
    return "\n".join(out)


def verdict(rows: list[dict]) -> str:
    """State the hypothesis outcome from the directnb decisive column, with the deciding numbers."""
    sel = [r for r in rows if r.get("gen_arm") == "directnb" and r.get("sigma_donor") == 0.5
           and r.get("n_donors") == 8 and r.get("sigma_het") is not None]
    if not sel:
        return "\n## Verdict\n  (decisive directnb cells not present in the summary yet)"
    max_het = max(r["sigma_het"] for r in sel)
    hi = {r["test"]: r for r in sel if r["sigma_het"] == max_het}

    lines = ["", f"## Verdict (at the strongest heterogeneity in the grid, het={max_het})"]
    broke = [t for t in MODERATED if t in hi and calibrated(hi[t]) is False]
    held = [t for t in MODERATED if t in hi and calibrated(hi[t]) is True]
    flat_ok = [t for t in UNSHRUNKEN if t in hi and calibrated(hi[t]) is True]

    for t in ("ttest", "pooled_t", "wilcoxon", "ebayes", "ebayes_trend", "voom"):
        r = hi.get(t)
        if not r:
            continue
        c = calibrated(r)
        tag = "CALIBRATED" if c else ("ANTI-CONSERVATIVE" if c is False else "?")
        lines.append(f"  {t:12} FP={_f(r.get('fp_rate'))} (binom_p={_f(r.get('fp_binom_p'))}) "
                     f"lambda={_f(r.get('lambda_pb'))} d0={_f(r.get('prior_df_d0'))}  -> {tag}")

    if broke and flat_ok:
        lines.append(f"\n  => SHRINKAGE HYPOTHESIS SUPPORTED: {', '.join(broke)} break under heterogeneity "
                     f"while {', '.join(flat_ok)} stay calibrated.")
    elif held and not broke:
        lines.append(f"\n  => SHRINKAGE HYPOTHESIS NOT SUPPORTED: moderated tests ({', '.join(held)}) "
                     "stay calibrated even at strong heterogeneity. The moderated method is then the "
                     "better choice IF it also wins on FDR-matched power; read the power frontier.")
    else:
        lines.append("\n  => MIXED / INCONCLUSIVE at this cell — inspect the full het column above.")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", type=Path, default=REPO / "pilot" / "testsel" / "summary.json")
    a = ap.parse_args()
    if not a.summary.exists():
        print(f"no summary at {a.summary} — run the grid first")
        return 1
    rows = load(a.summary)
    print(f"# test-selection analysis  ({len(rows)} cells)")
    counts = defaultdict(int)
    for r in rows:
        counts[r.get("test")] += 1
    print("  cells per test:", dict(counts))
    print(shrinkage_table(rows))
    print(power_frontier(rows))
    print(verdict(rows))
    print("\n(synthetic only; the binding check is the Mathys 2019 real anchor, not yet run)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
