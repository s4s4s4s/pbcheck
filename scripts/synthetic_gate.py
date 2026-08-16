"""Run the full corrected pbcheck instrument on synthetic oracles and emit a GO/NO-GO-style readout.

This is the calibration of the measurement instrument on KNOWN truth, per docs/PHASE0_SPEC.md §8
oracles (b) synthetic null and (c) synthetic positive. It exercises the whole verified pipeline:
frozen label-agnostic universe (with the C5 minimum enforced) -> naive & moderated-eBayes pseudobulk
over that universe -> ONE paired BH across both arms -> donor permutation null -> genomic-inflation
lambda + permutation floor + pseudobulk calibration & power.

The real GO/NO-GO gate runs the same code over real CELLxGENE datasets (needs cellxgene-census).
Here we only prove the instrument reads correctly where we already know the answer.

Thresholds come from :mod:`pbcheck.gate_config` and nowhere else; the pre-registered ones are frozen
by the spec and changing one requires an AMENDMENTS.md entry written first. Defaults for every
command-line option are the pre-registered values, so ``python scripts/synthetic_gate.py`` with no
arguments is the protocol run.

As of Amendment 2 (2026-08-15) the pseudobulk arm's test is moderated eBayes, both arms are
corrected over one common tested set, and the run records the realised prior (the moderated analog
of DESeq2's ``fitType``, C5) plus a null p-value uniformity check (the replacement for the retired
C3 independent-filtering cross-check).

As of Amendment 3 (2026-08-15) the two halves of the validity gate are evaluated at **two different
regimes**, each labelled in the output and in the artifact:

* **calibration** (λ band, perm-null FP rate, and the naive-arm instrument checks) at the hard point
  ``donor_sigma`` = :data:`pbcheck.gate_config.CALIBRATION_EVAL_SIGMA` — unchanged, unrelaxed;
* **power** at ``donor_sigma`` = :data:`pbcheck.gate_config.POWER_EVAL_SIGMA`, the boundary of the
  declared operating envelope, at the UNCHANGED pre-registered effect size (log2FC = 1.0, K = 200).

The envelope itself is printed and written to the artifact, because a verdict of "valid" here is
valid **within a stated domain** and must not be readable as unconditional. A stratum outside the
envelope is outside what this instrument claims; see docs/AMENDMENTS.md, Amendment 3 Change 1.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
from scipy.stats import kstest  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "synthetic"))
from oracles import null_oracle, positive_oracle  # noqa: E402

from pbcheck import gate_config as gc  # noqa: E402
from pbcheck import metrics, mtc  # noqa: E402
from pbcheck.gene_universe import frozen_universe  # noqa: E402
from pbcheck.methods import naive_de, pseudobulk_de  # noqa: E402
from pbcheck.methods.pseudobulk import build_pseudobulk  # noqa: E402
from pbcheck.permutation import run_null  # noqa: E402

FDR = gc.ALPHA


def paired_real_label(oracle, universe):
    """Both arms on the real labels, BH-corrected over ONE common tested set (spec §5).

    Amendment 2 Change 2: this used the per-arm path, which gave the two arms different BH
    denominators whenever the pseudobulk arm NA'd a gene the naive arm tested.
    """
    nv = naive_de(oracle.adata, genes=universe)
    pb = pseudobulk_de(oracle.adata, universe=universe)
    paired = mtc.bh_both_arms(nv, pb, universe, alpha=FDR)
    return paired, getattr(pb, "moderation", {})


def null_pvalue_uniformity(pval_matrix):
    """Median per-permutation KS distance of the arm's null p-values from Uniform(0, 1).

    Replaces the retired C3 cross-check (Amendment 2 Change 5): DESeq2's native
    independent-filtered ``padj`` has no counterpart in a moderated arm, but C3's purpose — showing
    the arm's calls are not an artifact of one p-value-processing choice — is served directly by
    asking whether the null p-values are uniform at all.

    Reported as the KS **statistic**, not a p-value: genes within one permutation are dependent, so
    a KS p-value here would be badly anti-conservative and would read as a false alarm. D is a
    distance and stays interpretable under dependence.
    """
    ds = []
    for row in np.asarray(pval_matrix, dtype=float):
        p = row[~np.isnan(row)]
        if p.size:
            ds.append(float(kstest(p, "uniform").statistic))
    return {"ks_D_median": float(np.median(ds)) if ds else float("nan"),
            "ks_D_max": float(np.max(ds)) if ds else float("nan"),
            "n_perm": len(ds)}


def naive_floors(null_res, n_genes):
    """The naive arm's two permutation floors, each over exactly ONE BH convention.

    ``run_null`` returns the naive #DEG as two labelled series (spec §5 / Amendment 2 Change 2 for
    the paired one, the arm's own BH for the other). They are different statistics and this
    function is the only place the gate turns either into a floor, so that the choice is made once
    and visibly:

    * ``paired`` is the **headline**. Two of the gate's criteria read it, and one of them compares
      it against the pseudobulk floor — a cross-arm comparison, which only the common-tested-set
      BH supports. The other ("most of the genome") is a naive-arm statement that the solo series
      would estimate at higher resolution, but it is left on the paired series so the gate has one
      headline floor rather than two that differ by which permutations they saw.
    * ``solo`` is reported alongside, never mixed in. It is the naive arm's own floor over all
      ``n_perm`` permutations and is what the Phase 1 false-discovery map reads.

    They coincide exactly while ``gate_config.N_PERM == N_PERM_PB``, which is why the defect this
    replaces — reading a floor off an array that was paired below ``n_paired`` and solo above it —
    was invisible in the committed run and would not have been at the pre-registered sweep counts.
    """
    return {
        "paired": metrics.perm_floor(null_res["naive_ndeg_paired"], n_genes),
        "solo": metrics.perm_floor(null_res["naive_ndeg_solo"], n_genes),
    }


def envelope_lines():
    """The declared operating envelope, as console lines (Amendment 3 Change 1(c)).

    Printed on every run, next to the verdict, so that "INSTRUMENT VALID" can never be read as an
    unconditional claim: the arm is valid for strata inside this region and is not declared valid
    outside it.
    """
    out = ["  operating envelope (Amendment 3 Change 1) — the arm is declared valid ONLY inside it:",
           "      sigma_donor   min donors/group   grid support"]
    for r in gc.OPERATING_ENVELOPE:
        out.append(f"      {r['sigma_donor']:>11}   {r['min_donors_per_group']:>16}   "
                   f"{r['grid_support']}")
    out.append("      (synthetic throughout; sigma_donor is an UNANCHORED simulator knob, so whether")
    out.append("       any real stratum falls inside this envelope is still unknown — Amendment 1's")
    out.append("       closing concern, restated by Amendments 2 and 3 and still open)")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--seed", type=int, default=1,
                    help="base seed; the null oracle uses it and the positive oracle uses seed+1")
    ap.add_argument("--n-perm", type=int, default=gc.N_PERM,
                    help=f"donor-label permutations for the naive arm (default {gc.N_PERM})")
    ap.add_argument("--n-perm-pb", type=int, default=gc.N_PERM_PB,
                    help="permutations on which BOTH arms run and are paired-BH'd "
                         f"(default {gc.N_PERM_PB}); the cross-arm numbers come from these")
    ap.add_argument("--out", type=Path, default=None,
                    help="write the full readout as JSON to this path (the run artifact)")
    a = ap.parse_args(argv)

    t0 = time.time()
    sim = dict(gc.ORACLE_SIM)
    # Amendment 3 Change 1: donor_sigma is no longer one number for the whole run. The calibration
    # arms stay at the HARD regime (the conservative choice for a false-positive criterion, and
    # exactly where every previous run measured them); the power oracle moves to the envelope
    # boundary, which is the regime at which the arm actually claims ≥ POWER_TARGET. It is popped
    # here so that neither call can silently inherit the other's sigma.
    sim.pop("donor_sigma")
    calib_sigma = gc.CALIBRATION_EVAL_SIGMA
    power_sigma = gc.POWER_EVAL_SIGMA
    n_d = sim["n_donors_per_group"]

    print(f"Building synthetic NULL oracle (donor structure, truth = 0 DE, "
          f"donor_sigma={calib_sigma} — the hard calibration regime)...")
    null = null_oracle(seed=a.seed, donor_sigma=calib_sigma, **sim)
    pdata = build_pseudobulk(null.adata)          # thin-donor filter runs here (Change 7)
    # C5, restored (Amendment 2 Change 4): the minimum-universe gate is ENFORCED, not accepted and
    # ignored. A stratum below it must SKIP rather than contribute a meaningless lambda and floor.
    universe = frozen_universe(pdata, min_size=gc.MIN_UNIVERSE_SIZE)
    G = len(universe)
    thin = pdata.uns.get("thin_donor_filter", {})
    print(f"  frozen label-agnostic universe: {G} genes (C5 minimum {gc.MIN_UNIVERSE_SIZE})")
    print(f"  thin-donor filter: {thin.get('n_dropped', 0)} of {thin.get('n_profiles_before', '?')}"
          f" profiles dropped (min_cells={gc.MIN_CELLS}, min_counts={gc.MIN_COUNTS})")

    paired, moderation = paired_real_label(null, universe)
    nv, pb = paired.naive, paired.pseudobulk
    print(f"  real-label #DEG    naive={nv.n_significant(fdr=FDR):4d}   "
          f"pseudobulk={pb.n_significant(fdr=FDR):3d}   "
          f"(one common tested set: {paired.n_tested_common}/{G} genes)")

    print(f"Running donor-permutation null (both arms, {a.n_perm} perms, "
          f"{a.n_perm_pb} paired)...")
    null_res = run_null(null.adata, universe, n_perm=a.n_perm, n_perm_pb=a.n_perm_pb, fdr=FDR,
                        seed=a.seed)

    lam_naive = metrics.lambda_over_permutations(null_res["naive_pvals"])
    lam_pb = metrics.lambda_over_permutations(null_res["pb_pvals"])
    naive_floor = naive_floors(null_res, G)
    floor, floor_solo = naive_floor["paired"], naive_floor["solo"]
    pb_floor = metrics.perm_floor(null_res["pb_ndeg"], G)
    pb_fp_rate = null_res["monte_carlo"]["pb_fp_rate"]
    pb_uniformity = null_pvalue_uniformity(null_res["pb_pvals"])

    # B5's construction (Amendment 2 Change 3), reported as what it is: a check that the permutation
    # set is exchangeable with the real labeling. It is uniform BY CONSTRUCTION under the sharp null
    # and therefore cannot be, and is not, a calibration criterion — see metrics.empirical_perm_pvalues.
    lam_emp = {
        arm: metrics.genomic_inflation(
            metrics.empirical_perm_pvalues(null_res[f"{arm}_pvals_real"], null_res[f"{arm}_pvals"]))
        for arm in ("naive", "pb")
    }

    # Spec §8(c) pre-registers K=200 at |log2FC|=1.0. Running an easier oracle here would be
    # exactly the goalpost-move pre-registration exists to prevent. See docs/AMENDMENTS.md.
    # The EFFECT SIZE is untouched by Amendment 3 and stays exactly as §8(c) froze it; what
    # Amendment 3 Change 1(b) changes is only WHERE the threshold binds — donor_sigma moves from
    # the arbitrary hard point to the boundary of the declared envelope. Power 0.60 at
    # donor_sigma=0.5 with 8v8 remains unmet and unclaimed; it is not re-tested here and is not
    # asserted anywhere.
    print(f"\nBuilding synthetic POSITIVE oracle ({gc.ORACLE_K} true DE genes, "
          f"|log2FC|={gc.ORACLE_LOG2FC}, shared; donor_sigma={power_sigma} — the envelope "
          f"boundary)...")
    pos = positive_oracle(n_de=gc.ORACLE_K, log2fc=gc.ORACLE_LOG2FC, seed=a.seed + 1,
                          donor_sigma=power_sigma, **sim)
    pos_uni = frozen_universe(build_pseudobulk(pos.adata), min_size=gc.MIN_UNIVERSE_SIZE)
    pos_paired, pos_moderation = paired_real_label(pos, pos_uni)
    pos_pb = pos_paired.pseudobulk
    truth = set(pos.de_genes) & set(pos_uni)
    hits = pos_pb.significant(fdr=FDR)
    sensitivity = len(hits & truth) / max(len(truth), 1)
    pos_fp = len(hits - set(pos.de_genes))
    emp_fdr = pos_fp / max(len(hits), 1)  # observed false-discovery rate on the positive oracle

    # Calibration yardstick (spec §8(a)/B1). NOTE: an earlier version of this file argued that
    # "fraction of permutations with >=1 rejection <= alpha" is a family-wise bound that BH does not
    # control. That was WRONG and is retracted: under the COMPLETE null every rejection is false, so
    # V = R, FDP = 1{R>=1}, and FDR = E[FDP] = P(R>=1). FWER and FDR coincide there, so BH does bound
    # it and the spec's original criterion stands. See docs/AMENDMENTS.md, Amendment 1.
    pb_mean_ndeg = float(np.mean(null_res["pb_ndeg"].counts))
    lo, hi = gc.LAMBDA_BAND

    def band(x):
        return "calibrated" if lo <= x <= hi else ("INFLATED" if x > hi else "under")

    print("\n" + "=" * 72)
    print(f"  SYNTHETIC GATE READOUT (instrument calibration on known truth, {n_d}v{n_d} donors)")
    print("=" * 72)
    print(f"  pseudobulk arm                      = {pb.method}")
    print(f"  calibration regime                  = donor_sigma {calib_sigma} (hard point)")
    print(f"  power regime                        = donor_sigma {power_sigma} (envelope boundary)")
    print(f"  lambda_naive (perm-null p-values)   = {lam_naive['lambda']:6.2f}  -> {band(lam_naive['lambda'])}")
    print(f"  lambda_pseudobulk                   = {lam_pb['lambda']:6.2f}  -> {band(lam_pb['lambda'])}")
    print(f"  pseudobulk perm-null FP rate        = {pb_fp_rate:6.2f}  (target <= {FDR}, "
          f"MC SE {null_res['monte_carlo']['pb_fp_rate_mc_se']:.3f})")
    print(f"  naive perm-floor  (median #DEG)     = {floor['median_count']:6.0f}  ({100*floor['median_frac']:.1f}% of {G} genes)")
    print(f"    BH convention = {floor['bh_mode']} over {floor['n_perm']} perms (cross-arm comparable); "
          f"naive-only BH = {floor_solo['median_count']:.0f} over {floor_solo['n_perm']} perms")
    print(f"  pseudobulk perm-floor (median #DEG) = {pb_floor['median_count']:6.0f}  (mean {pb_mean_ndeg:.2f})")
    print(f"  pseudobulk power (positive oracle)  = {sensitivity:6.2f}  (target >= {gc.POWER_TARGET}"
          f", at donor_sigma {power_sigma})")
    print(f"  pseudobulk empirical FDR (positive) = {emp_fdr:6.2f}  (incl. chance donor imbalance; "
          f"true FDR calib. is the perm-null above)")
    print("-" * 72)
    print("  moderation (the fitType analog, spec C5):")
    print(f"    prior df d0 = {moderation.get('prior_df_d0')}   complete pooling = "
          f"{moderation.get('prior_is_complete_pooling')}   residual df = {moderation.get('residual_df_d')}")
    print(f"    shrinkage d0/(d0+d) = {moderation.get('shrinkage_factor_d0_over_d0_plus_d')}   "
          f"trended prior = {moderation.get('s0_squared_is_trended')}   "
          f"untestable genes = {moderation.get('n_pval_nan')}")
    print(f"    null p-value uniformity (retired-C3 replacement): median KS D = "
          f"{pb_uniformity['ks_D_median']:.4f}, max {pb_uniformity['ks_D_max']:.4f}")
    print("  fairness (spec §5, one common tested set — Amendment 2 Change 2):")
    print(f"    real labels: {paired.n_tested_common}/{G} genes common, "
          f"{paired.n_na_naive} naive NA, {paired.n_na_pseudobulk} pseudobulk NA")
    print(f"    permutations: {null_res['n_perm_paired']} paired of {null_res['n_perm_naive']}")
    print("  B5 empirical permutation p-values (machinery check, NOT a calibration criterion):")
    print(f"    lambda_naive_empirical = {lam_emp['naive']:.2f}   "
          f"lambda_pb_empirical = {lam_emp['pb']:.2f}   "
          f"(near 1 by construction; resolution floor p >= {1/(1+null_res['n_perm_paired']):.3f}, "
          f"so few permutations drive these BELOW 1 by discreteness alone)")
    print("-" * 72)

    # The pseudobulk validity gate is BINDING (spec decision rule item 1): if the denominator arm is
    # not both calibrated and powered, every inflation number above is uninterpretable and the study
    # is VOID. The two pre-registered calibration criteria are therefore wired in as pass/fail here,
    # at the spec's own thresholds — not narrated afterwards in prose.
    #
    # Each criterion carries the regime it was evaluated at and the document that pins it, because
    # after Amendment 3 the two halves of the gate no longer run at the same donor_sigma and a bare
    # PASS/FAIL list would hide that.
    sanity = "instrument sanity — NOT pre-registered, not the decision rule"
    criteria = [
        ("naive flagged inflated (lambda >> 1)",
         lam_naive["lambda"] > gc.INSTRUMENT_LAMBDA_NAIVE_MIN,
         f"donor_sigma={calib_sigma}, {n_d}v{n_d}", sanity),
        ("naive floor is most of the genome",
         floor["median_frac"] > gc.INSTRUMENT_NAIVE_FLOOR_FRAC_MIN,
         f"donor_sigma={calib_sigma}, {n_d}v{n_d}", sanity),
        (f"pseudobulk lambda in pre-registered band [{lo}, {hi}]",
         lo <= lam_pb["lambda"] <= hi,
         f"donor_sigma={calib_sigma} (hard regime), {n_d}v{n_d}",
         "spec §8(a) / decision rule item 1; regime pinned by Amendment 3 Change 1(a), "
         "threshold unchanged"),
        ("pseudobulk perm-null FP rate <= alpha (spec §8a)",
         pb_fp_rate <= FDR,
         f"donor_sigma={calib_sigma} (hard regime), {null_res['n_perm_paired']} paired perms",
         "spec §8(a) / decision rule item 1 (restored binding by Amendment 1 Change 1); "
         "resolution raised to 200 permutations by Amendment 3 Change 2"),
        ("naive floor >> pseudobulk floor",
         floor["median_count"] > gc.INSTRUMENT_FLOOR_RATIO_MIN * max(pb_floor["median_count"], 1),
         f"donor_sigma={calib_sigma}, {n_d}v{n_d}", sanity),
        (f"pseudobulk powered (sens >= {gc.POWER_TARGET} at pre-registered "
         f"log2FC={gc.ORACLE_LOG2FC}, K={gc.ORACLE_K})",
         sensitivity >= gc.POWER_TARGET,
         f"donor_sigma={power_sigma} (envelope boundary), {n_d}v{n_d}",
         "spec §8(c) / decision rule item 1 — effect size and threshold UNCHANGED; evaluation "
         "regime re-scoped from a point to the envelope boundary by Amendment 3 Change 1(b)"),
    ]
    checks = {name: ok for name, ok, _, _ in criteria}
    for name, ok, at, _ in criteria:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        print(f"         evaluated at {at}")
    print("-" * 72)
    for line in envelope_lines():
        print(line)
    verdict = ("INSTRUMENT VALID WITHIN THE STATED OPERATING ENVELOPE "
               "-> ready for real CELLxGENE sweep, envelope-gated"
               if all(checks.values()) else "INSTRUMENT NEEDS ATTENTION")
    elapsed = time.time() - t0
    print("=" * 72)
    print(f"  {verdict}")
    print("=" * 72)
    print(f"  ({elapsed:.1f} s)")

    if a.out:
        artifact = {
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "machine": {"platform": platform.platform(), "python": platform.python_version()},
            "runtime_seconds": round(elapsed, 1),
            "params": {"seed": a.seed, "n_perm": a.n_perm, "n_perm_pb": a.n_perm_pb,
                       "calibration_eval_sigma": calib_sigma, "power_eval_sigma": power_sigma,
                       "n_donors_per_group": n_d},
            "config": gc.manifest(),
            # Amendment 3 Change 1(c): the validity domain is a first-class output, not a footnote.
            "operating_envelope": [dict(r) for r in gc.OPERATING_ENVELOPE],
            "operating_envelope_source":
                "docs/AMENDMENTS.md Amendment 3 Change 1 (frontier from Amendment 1, corroborated "
                "by pilot/testsel/summary.json @ 72dec7b). SYNTHETIC: sigma_donor is an unanchored "
                "simulator knob; whether any real stratum falls inside is unknown.",
            "validity_scope":
                "The pseudobulk arm is declared valid ONLY for strata inside operating_envelope. "
                "Power 0.60 at donor_sigma=0.5 with 8 donors/group remains UNMET (measured 0.35 on "
                "2026-08-15) and is NOT claimed.",
            "arm": {"pseudobulk": pb.method, "naive": nv.method},
            # Which naive engine produced the permutation numbers. The fast engine is a
            # bitwise-identical reformulation of the same statistic (one ranking pass per stratum
            # instead of one scanpy re-run per permutation), but an artifact that does not record
            # which code path ran cannot be audited against a later one that changed it.
            "naive_engine": null_res["naive_engine"],
            "naive_engine_setup": null_res["naive_engine_setup"],
            "universe_size": G,
            "thin_donor_filter": thin,
            "real_label_ndeg": {"naive": nv.n_significant(fdr=FDR),
                                "pseudobulk": pb.n_significant(fdr=FDR)},
            "lambda_naive": lam_naive["lambda"],
            "lambda_pseudobulk": lam_pb["lambda"],
            "lambda_empirical_perm_b5": lam_emp,
            "pb_perm_null_fp_rate": pb_fp_rate,
            # Two floors, each tagged with the BH convention it was taken over; never one array
            # holding both (see permutation.run_null / metrics.NDegSeries).
            "naive_perm_floor": floor,
            "naive_perm_floor_solo": floor_solo,
            "pb_perm_floor": pb_floor,
            "pb_mean_ndeg": pb_mean_ndeg,
            "pb_null_pvalue_uniformity": pb_uniformity,
            "power_sensitivity": sensitivity,
            "power_empirical_fdr": emp_fdr,
            "moderation_null_oracle": moderation,
            "moderation_positive_oracle": pos_moderation,
            "paired_bh_real": null_res["paired_bh_real"],
            "monte_carlo": null_res["monte_carlo"],
            "real_split_inside_perm_range": null_res["real_split_inside_perm_range"],
            "checks": {k: bool(v) for k, v in checks.items()},
            "criteria": [{"name": n, "passed": bool(ok), "evaluated_at": at, "source": src}
                         for n, ok, at, src in criteria],
            "verdict": verdict,
            "amendments_applied": ["Amendment 1 (2026-07-22)", "Amendment 2 (2026-08-15)",
                                   "Amendment 3 (2026-08-15)"],
        }
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(artifact, indent=2, default=str), encoding="utf-8")
        print(f"  artifact -> {a.out}")

    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
