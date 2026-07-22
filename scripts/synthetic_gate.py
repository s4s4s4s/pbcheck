"""Run the full corrected pbcheck instrument on synthetic oracles and emit a GO/NO-GO-style readout.

This is the calibration of the measurement instrument on KNOWN truth, per docs/PHASE0_SPEC.md §8
oracles (b) synthetic null and (c) synthetic positive. It exercises the whole verified pipeline:
frozen label-agnostic universe -> naive & pseudobulk over that universe -> identical BH -> donor
permutation null -> genomic-inflation lambda + permutation floor + pseudobulk calibration & power.

The real GO/NO-GO gate runs the same code over real CELLxGENE datasets (needs cellxgene-census).
Here we only prove the instrument reads correctly where we already know the answer.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "synthetic"))
from oracles import null_oracle, positive_oracle  # noqa: E402

from pbcheck.methods import naive_de, pseudobulk_de  # noqa: E402
from pbcheck.methods.pseudobulk import build_pseudobulk  # noqa: E402
from pbcheck.gene_universe import frozen_universe  # noqa: E402
from pbcheck import mtc, metrics  # noqa: E402
from pbcheck.permutation import run_null  # noqa: E402

FDR = 0.05


def real_label(oracle, universe):
    pb = mtc.bh_over_universe(pseudobulk_de(oracle.adata, universe=universe), universe, alpha=FDR)
    nv = mtc.bh_over_universe(naive_de(oracle.adata, genes=universe), universe, alpha=FDR)
    return nv, pb


def main() -> None:
    # 8v8 donors: the regime the spec says carries the headline (permutations ~orthogonal to truth,
    # pseudobulk both calibrated AND powered). 4v4 is deliberately under-powered for pseudobulk.
    kw = dict(n_genes=1500, n_donors_per_group=8, n_cells_per_donor=250, dispersion=0.2)

    print("Building synthetic NULL oracle (donor structure, truth = 0 DE)...")
    null = null_oracle(seed=1, donor_sigma=0.5, **kw)
    pdata = build_pseudobulk(null.adata)
    universe = frozen_universe(pdata)
    G = len(universe)
    print(f"  frozen label-agnostic universe: {G} genes")

    nv, pb = real_label(null, universe)
    print(f"  real-label #DEG    naive={nv.n_significant(fdr=FDR):4d}   pseudobulk={pb.n_significant(fdr=FDR):3d}")

    print("Running donor-permutation null (both arms)...")
    null_res = run_null(null.adata, universe, n_perm=40, n_perm_pb=20, fdr=FDR)

    lam_naive = metrics.lambda_over_permutations(null_res["naive_pvals"])
    lam_pb = metrics.lambda_over_permutations(null_res["pb_pvals"])
    floor = metrics.perm_floor(null_res["naive_ndeg"], G)
    pb_floor = metrics.perm_floor(null_res["pb_ndeg"], G)
    pb_fp_rate = float(np.mean(null_res["pb_ndeg"] >= 1))

    # Spec §8(c) pre-registers K=200 at |log2FC|=1.0. Running an easier oracle here would be
    # exactly the goalpost-move pre-registration exists to prevent. See docs/AMENDMENTS.md.
    print("\nBuilding synthetic POSITIVE oracle (200 true DE genes, |log2FC|=1.0, shared)...")
    pos = positive_oracle(n_de=200, log2fc=1.0, seed=2, donor_sigma=0.5, **kw)
    pos_uni = frozen_universe(build_pseudobulk(pos.adata))
    pos_pb = mtc.bh_over_universe(pseudobulk_de(pos.adata, universe=pos_uni), pos_uni, alpha=FDR)
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
    pb_mean_ndeg = float(np.mean(null_res["pb_ndeg"]))

    def band(x, lo, hi):
        return "calibrated" if lo <= x <= hi else ("INFLATED" if x > hi else "under")

    print("\n" + "=" * 72)
    print("  SYNTHETIC GATE READOUT (instrument calibration on known truth, 8v8 donors)")
    print("=" * 72)
    print(f"  lambda_naive (perm-null p-values)   = {lam_naive['lambda']:6.2f}  -> {band(lam_naive['lambda'],0.9,1.1)}")
    print(f"  lambda_pseudobulk                   = {lam_pb['lambda']:6.2f}  -> {band(lam_pb['lambda'],0.9,1.1)}")
    print(f"  pseudobulk perm-null FP rate        = {pb_fp_rate:6.2f}  (target <= {FDR})")
    print(f"  naive perm-floor  (median #DEG)     = {floor['median_count']:6.0f}  ({100*floor['median_frac']:.1f}% of {G} genes)")
    print(f"  pseudobulk perm-floor (median #DEG) = {pb_floor['median_count']:6.0f}  (mean {pb_mean_ndeg:.2f})")
    print(f"  pseudobulk power (positive oracle)  = {sensitivity:6.2f}  (target >= 0.60)")
    print(f"  pseudobulk empirical FDR (positive) = {emp_fdr:6.2f}  (incl. chance donor imbalance; "
          f"true FDR calib. is the perm-null above)")
    print("-" * 72)

    # The pseudobulk validity gate is BINDING (spec decision rule item 1): if the denominator arm is
    # not both calibrated and powered, every inflation number above is uninterpretable and the study
    # is VOID. The two pre-registered calibration criteria are therefore wired in as pass/fail here,
    # at the spec's own thresholds — not narrated afterwards in prose.
    checks = {
        "naive flagged inflated (lambda >> 1)": lam_naive["lambda"] > 1.5,
        "naive floor is most of the genome": floor["median_frac"] > 0.3,
        "pseudobulk lambda in pre-registered band [0.9, 1.1]": 0.9 <= lam_pb["lambda"] <= 1.1,
        "pseudobulk perm-null FP rate <= alpha (spec §8a)": pb_fp_rate <= FDR,
        "naive floor >> pseudobulk floor": floor["median_count"] > 10 * max(pb_floor["median_count"], 1),
        "pseudobulk powered (sens >= 0.6 at pre-registered log2FC=1.0, K=200)": sensitivity >= 0.60,
    }
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    verdict = ("INSTRUMENT VALID -> ready for real CELLxGENE sweep" if all(checks.values())
               else "INSTRUMENT NEEDS ATTENTION")
    print("=" * 72)
    print(f"  {verdict}")
    print("=" * 72)


if __name__ == "__main__":
    main()
