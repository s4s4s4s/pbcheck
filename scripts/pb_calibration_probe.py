"""Pseudobulk-arm calibration probe (DIAGNOSTIC — not part of the frozen Phase 0 pipeline).

WHY THIS EXISTS
---------------
``scripts/synthetic_gate.py`` reports the pseudobulk arm as valid while, at 8v8 donors, it
actually shows ``lambda_pseudobulk ~ 1.25`` (pre-registered band [0.9, 1.1]) and a
permutation-null false-positive rate far above alpha. The spec (docs/PHASE0_SPEC.md, decision
rule item 1 and section 8a/8c) makes pseudobulk calibration+power a BINDING gate, so that
failure has to be diagnosed before the real sweep's stratum list is pre-registered.

Two rival explanations with opposite remedies:

  H1  DESeq2-Wald is anti-conservative at n ~ 8 (dispersion estimated from 8 samples).
      -> simpler donor-level tests on the SAME aggregated matrix would be calibrated,
         and the remedy is to change the pseudobulk arm's TEST.
  H2  the simulator in ``synthetic/oracles.py`` is misspecified (a log-normal donor random
      effect on top of an NB does not give a donor-level NB), so DESeq2's dispersion model
      is wrong and NO test would be calibrated.
      -> the remedy is to recalibrate the SIMULATOR (spec section 8(b) already demands a real
         mean-dispersion trend), not to swap the test.

This script is the measurement instrument for that question. It holds EVERYTHING fixed
(generative model, seed, donor aggregation, frozen gene universe, BH over that universe,
alpha = 0.05) and varies exactly one thing: ``--test``.

STRUCTURAL GUARANTEE
--------------------
The donor pseudobulk matrix is aggregated ONCE per oracle via the engine's
``build_pseudobulk``; the frozen universe is built ONCE from it via the engine's
``frozen_universe``; the matrix is then restricted to that universe ONCE. Every test
receives the identical ``(pdata_u, universe, condition_vector)`` triple through the single
dispatch function :func:`run_test`. There is no code path by which two tests can see
different data. All four go through the engine's ``pbcheck.mtc.bh_over_universe`` and the
engine's ``pbcheck.metrics`` for lambda and the floor.

WHAT IS REUSED VS BYPASSED
--------------------------
Reused (the real pipeline, not a reimplementation):
  * ``synthetic.oracles.null_oracle`` / ``positive_oracle``  — the generative model.
  * ``pbcheck.methods.pseudobulk.build_pseudobulk``          — decoupler donor aggregation.
  * ``pbcheck.gene_universe.frozen_universe``                — label-agnostic frozen universe.
  * ``pbcheck.methods.pseudobulk.deseq_from_pdata``          — the Wald arm, verbatim.
  * ``pbcheck.permutation.build_perms`` / ``_labels_for``    — the donor-permutation null.
  * ``pbcheck.mtc.bh_over_universe``                         — identical BH for every arm.
  * ``pbcheck.metrics.genomic_inflation`` / ``lambda_over_permutations`` / ``perm_floor``.

BYPASSED, deliberately, each flagged again at the call site:
  1. ``pbcheck.permutation.run_null`` is NOT used. It hard-codes the Wald test
     (``deseq_from_pdata``) and always runs the expensive naive per-cell arm, which this
     probe does not measure. The permutation loop here is a copy of ``run_null``'s
     pseudobulk half with the test made pluggable; the permutation SET is still produced by
     the engine's ``build_perms``, so the null itself is not reimplemented.
  2. ``--test lrt`` is hand-rolled. PyDESeq2 0.5.4 exposes NO likelihood-ratio test: the
     installed ``pydeseq2/ds.py`` has only ``run_wald_test``, ``DeseqStats.__init__`` has no
     ``test``/``reduced`` parameter, and grepping the installed package for
     ``lrt|likelihood_ratio|reduced`` returns nothing. See :func:`_deseq_lrt` for exactly
     what is computed instead and how it maps onto DESeq2's ``nbinomLRT``.
  3. ``--test ttest`` / ``--test wilcoxon`` have no engine counterpart at all; they are
     defined here on the aggregated donor matrix.

CAVEAT THE CALLER MUST NOT FORGET
---------------------------------
Under a DONOR-LABEL PERMUTATION null the donor profiles are fixed and only labels move, so
any permutation-invariant test is calibrated BY CONSTRUCTION. The Mann-Whitney U null IS the
donor-label permutation distribution, so ``--test wilcoxon --null-mode permutation`` is
guaranteed to give lambda ~ 1 and FP rate ~ 0 and is therefore NOT evidence about H1 vs H2.
The Welch t-test is nearly so. Distinguishing H1 from H2 requires calibration measured
against the GENERATIVE model, i.e. independent fresh simulations — ``--null-mode fresh``,
which draws ``--n-perm`` independent datasets with truth = 0 DE and the real (unpermuted)
labels. Both modes emit the same JSON keys; read ``null_mode`` before interpreting them.

Second structural trap, reported as ``bh_min_genes_at_min_p``. With D donors split n/n the
smallest attainable two-sided exact Mann-Whitney p-value is 2/C(D, n) — at 8v8 that is
2/12870 = 1.554e-4, which is ABOVE the Bonferroni line alpha/G = 3.33e-5 at G = 1500. That
does NOT mean the arm cannot reject: BH is a step-up procedure, so it rejects as soon as
p_(i) <= alpha*i/G for some i, and k genes tied at the floor value need only
k >= ceil(min_p * G / alpha) = 5 to clear it. (This was checked numerically before it was
written down: the wilcoxon arm did in fact reject 11 true genes at 8v8/G=1500.) What the
number does mean is that the wilcoxon arm has a hard granularity floor: it cannot report
fewer than ~5 discoveries, and at small D its p-value support is so coarse that lambda is
dominated by discreteness (at 4v4, min_p = 2/70 = 0.0286 and lambda collapses to ~0.36).
Treat wilcoxon lambda as uninterpretable below ~6v6 and read ``min_attainable_pvalue``
before using any wilcoxon number.

USAGE
-----
    python scripts/pb_calibration_probe.py \
        --sigma-donor 0.5 --n-donors 8 --test wald \
        --n-genes 1500 --n-cells-per-donor 250 --dispersion 0.2 \
        --n-perm 20 --seed 7 --out out.json

All results belong OUTSIDE the repo. Nothing here writes into the repo.
"""

from __future__ import annotations

import argparse
import importlib.metadata as _md
import json
import platform
import sys
import time
import warnings
from math import comb, log
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.stats import chi2, mannwhitneyu, ttest_ind

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))
sys.path.insert(0, str(_REPO / "synthetic"))

from oracles import null_oracle, positive_oracle  # noqa: E402

from pbcheck import metrics, mtc  # noqa: E402
from pbcheck.gene_universe import frozen_universe  # noqa: E402
from pbcheck.methods.de import DEResult  # noqa: E402
from pbcheck.methods.pseudobulk import build_pseudobulk, deseq_from_pdata  # noqa: E402
from pbcheck.permutation import _labels_for, build_perms  # noqa: E402

# ---------------------------------------------------------------------------
# Frozen constants. These are the pre-registered values (docs/PHASE0_SPEC.md).
# They are NOT exposed as flags: the whole point of the audit is that the gate
# script silently substituted easier ones (log2FC 1.5 / K 150).
# ---------------------------------------------------------------------------
ALPHA = 0.05                    # spec section 5: alpha = 0.05 everywhere
ORACLE_LOG2FC = 1.0             # spec section 8(c): fixed log2FC = 1.0
ORACLE_K = 200                  # spec section 8(c): K = 200 injected genes
ORACLE_K_REFERENCE_GENES = 1500  # universe size K = 200 was written against
LAMBDA_BAND = (0.9, 1.1)        # spec decision rule item 1
POWER_TARGET = 0.60             # spec section 8(c)
TESTS = ("wald", "lrt", "ttest", "wilcoxon")


# ===========================================================================
# The four tests. Every one takes the SAME (pdata_u, condition_values) pair.
# ===========================================================================

def _counts_frame(pdata_u) -> pd.DataFrame:
    """Integer donor x gene count matrix, exactly as ``deseq_from_pdata`` builds it."""
    return pd.DataFrame(
        np.asarray(pdata_u.X).round().astype(int),
        index=pdata_u.obs_names,
        columns=pdata_u.var_names,
    )


def _log_cpm(pdata_u) -> np.ndarray:
    """log2(CPM + 1) of the SAME aggregated donor matrix. Shape (n_donors, n_genes).

    Library size is the row sum of the universe-restricted matrix, so the transform is a
    function of exactly the data the DESeq2 arms see (no genes outside the frozen universe
    leak into the normalisation).

    Note: Welch's t is invariant to the log base (it is a constant rescaling of every value)
    and Mann-Whitney is rank-based, so log2 vs ln changes no p-value. log2 is used only so
    the reported ``log2fc`` column is directly interpretable.
    """
    X = np.asarray(pdata_u.X, dtype=float)
    lib = X.sum(axis=1, keepdims=True)
    cpm = X / np.maximum(lib, 1.0) * 1e6
    return np.log2(cpm + 1.0)


def _donor_level_result(pdata_u, cond, test_level, ref_level, method_name, kind) -> DEResult:
    """Welch t-test or Mann-Whitney U across donors on log2(CPM+1). BYPASS: no engine equivalent."""
    Y = _log_cpm(pdata_u)
    cond = np.asarray(cond, dtype=object)
    a = Y[cond == test_level, :]
    b = Y[cond == ref_level, :]
    if a.shape[0] < 2 or b.shape[0] < 2:
        raise ValueError(f"need >= 2 donors per group, got {a.shape[0]} vs {b.shape[0]}")

    with np.errstate(invalid="ignore", divide="ignore"):
        if kind == "ttest":
            pval = ttest_ind(a, b, axis=0, equal_var=False).pvalue
        elif kind == "wilcoxon":
            # method="exact" is refused when ties are present; "auto" falls back to the
            # normal approximation with tie correction, matching the naive arm's convention.
            pval = mannwhitneyu(a, b, axis=0, alternative="two-sided", method="auto").pvalue
        else:  # pragma: no cover
            raise ValueError(kind)
    pval = np.asarray(pval, dtype=float)
    log2fc = a.mean(axis=0) - b.mean(axis=0)

    table = pd.DataFrame(
        {"pval": pval, "padj": np.nan, "log2fc": log2fc},
        index=pd.Index(list(pdata_u.var_names), name="gene"),
    )
    return DEResult(method=method_name, table=table, contrast=("condition", test_level, ref_level))


def _deseq_lrt(pdata_u, cond, test_level, ref_level, *, n_cpus=4) -> DEResult:
    """DESeq2 likelihood-ratio test, ``~condition`` vs ``~1``, hand-rolled.

    BYPASS + API SURPRISE. PyDESeq2 0.5.4 has no LRT: ``DeseqStats`` exposes only
    ``run_wald_test`` and its signature is
    ``(dds, contrast, alpha, cooks_filter, independent_filter, prior_LFC_var, lfc_null,
    alt_hypothesis, inference, quiet, n_cpus)`` — no ``test=``, no ``reduced=``. So the LRT is
    assembled from the installed internals, mirroring DESeq2's ``nbinomLRT``:

      1. Build the DeseqDataSet EXACTLY as ``pbcheck.methods.pseudobulk.deseq_from_pdata``
         does (poscounts size factors C1, refit_cooks=True, design ~condition) and run
         ``dds.deseq2()``. Dispersions are therefore the same MAP dispersions the Wald arm
         uses, estimated under the FULL design — which is what DESeq2's LRT also does.
      2. Refit the full model and fit the reduced model ``~1`` with ``dds.inference.irls``,
         both on ``dds.X`` with those same dispersions, size factors and ``min_mu``. Fitting
         the full model here too (rather than reading ``obsm['_mu_LFC']``) makes the two
         deviances internally consistent by construction; it was verified to reproduce
         ``_mu_LFC`` to 0.0 absolute error at 8v8/400 genes.
      3. stat = 2 * (nll_reduced - nll_full) using ``pydeseq2.utils.nb_nll``, p = chi2.sf(stat, 1).

    Known departure from R DESeq2: Cook's outlier count REPLACEMENT (refit_cooks) affects the
    dispersion estimates but the deviances here are computed on the unreplaced ``dds.X``. At
    8v8 on this simulator ``dds.var['replaced'].sum()`` was 0, so the departure was inert; the
    JSON reports ``n_cooks_replaced_genes`` so a caller can see when it is not.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.utils import nb_nll

    counts = _counts_frame(pdata_u)
    metadata = pd.DataFrame(
        {"condition": np.asarray(cond, dtype=object).astype(str)}, index=pdata_u.obs_names
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design="~condition",
            ref_level=["condition", ref_level],
            size_factors_fit_type="poscounts",  # C1, same as the Wald arm
            refit_cooks=True,
            quiet=True,
            n_cpus=n_cpus,
        )
        dds.deseq2()

        nz = dds.var["non_zero"].to_numpy()
        disp = dds.var.loc[dds.var["non_zero"], "dispersions"].to_numpy()
        sf = dds.obs["size_factors"].to_numpy()
        Y = dds.X[:, dds.non_zero_idx]
        X_full = dds.obsm["design_matrix"].values
        X_red = np.ones((dds.n_obs, 1))

        beta_full, mu_full, _, conv_full = dds.inference.irls(
            counts=Y, size_factors=sf, design_matrix=X_full, disp=disp,
            min_mu=dds.min_mu, beta_tol=dds.beta_tol,
        )
        _, mu_red, _, conv_red = dds.inference.irls(
            counts=Y, size_factors=sf, design_matrix=X_red, disp=disp,
            min_mu=dds.min_mu, beta_tol=dds.beta_tol,
        )

    stat = 2.0 * (nb_nll(Y, mu_red, disp) - nb_nll(Y, mu_full, disp))
    stat = np.maximum(np.asarray(stat, dtype=float), 0.0)
    pval_nz = chi2.sf(stat, df=1)

    cols = list(dds.obsm["design_matrix"].columns)
    coef = cols.index([c for c in cols if c != "Intercept"][0])

    pval = np.full(dds.n_vars, np.nan)
    log2fc = np.full(dds.n_vars, np.nan)
    pval[nz] = pval_nz
    log2fc[nz] = beta_full[:, coef] / log(2.0)  # pydeseq2 betas are natural log

    table = pd.DataFrame(
        {"pval": pval, "padj": np.nan, "log2fc": log2fc},
        index=pd.Index(list(dds.var_names), name="gene"),
    )
    res = DEResult(
        method="pseudobulk[DESeq2:poscounts,LRT ~condition vs ~1]",
        table=table,
        contrast=("condition", test_level, ref_level),
    )
    res._diag = {  # attached for provenance, not part of DEResult's contract
        "n_cooks_replaced_genes": int(dds.var["replaced"].sum()) if "replaced" in dds.var else -1,
        "irls_full_converged_frac": float(np.mean(conv_full)),
        "irls_reduced_converged_frac": float(np.mean(conv_red)),
        "disp_function_type": str(dds.uns.get("disp_function_type", "?")),
    }
    return res


def run_test(pdata_u, universe, cond, test, *, test_level="disease", ref_level="ctrl", n_cpus=4):
    """THE single dispatch point. Every arm gets the identical aggregated matrix + universe.

    ``pdata_u`` is already restricted to ``universe`` in ``universe`` order, so no test can
    silently see a different gene set. Returns a BH-corrected DEResult over exactly ``universe``.
    """
    if test == "wald":
        # Engine function, verbatim — this IS the pseudobulk arm as the repo runs it today.
        res = deseq_from_pdata(
            pdata_u, condition_col="condition", test_level=test_level, ref_level=ref_level,
            universe=universe, condition_values=cond, fdr=ALPHA, n_cpus=n_cpus,
        )
    elif test == "lrt":
        res = _deseq_lrt(pdata_u, cond, test_level, ref_level, n_cpus=n_cpus)
    elif test in ("ttest", "wilcoxon"):
        name = {"ttest": "donor[Welch t on log2CPM]", "wilcoxon": "donor[MWU on log2CPM]"}[test]
        res = _donor_level_result(pdata_u, cond, test_level, ref_level, name, test)
    else:
        raise ValueError(f"unknown test {test!r}; expected one of {TESTS}")
    out = mtc.bh_over_universe(res, universe, alpha=ALPHA)
    if hasattr(res, "_diag"):
        out._diag = res._diag
    return out


# ===========================================================================
# Oracle plumbing
# ===========================================================================

def _aggregate(adata):
    """Aggregate donors ONCE, freeze the universe ONCE, restrict ONCE."""
    pdata = build_pseudobulk(adata)                  # engine
    universe = frozen_universe(pdata)                # engine, label-agnostic (A3)
    keep = [g for g in universe if g in set(pdata.var_names)]
    assert keep == universe, "frozen universe must be a subset of the pseudobulk var_names"
    pdata_u = pdata[:, universe].copy()
    donor_names = pd.Index(
        pdata_u.obs["donor"].astype(str) if "donor" in pdata_u.obs else pdata_u.obs_names
    )
    return pdata_u, universe, donor_names


def _effective_k(n_genes: int) -> tuple[int, str]:
    """Pre-registered K = 200; scaled proportionally only if it cannot fit in n_genes."""
    if ORACLE_K <= n_genes:
        return ORACLE_K, "none (pre-registered K=200 used as-is)"
    k = max(1, int(round(ORACLE_K * n_genes / ORACLE_K_REFERENCE_GENES)))
    k = min(k, n_genes)
    return k, (f"K=200 exceeds n_genes={n_genes}; scaled proportionally against the "
               f"reference universe of {ORACLE_K_REFERENCE_GENES} genes -> K={k}")


# ===========================================================================
# Null: calibration
# ===========================================================================

def _null_permutation(sim_kw, test, n_perm, seed, n_cpus):
    """Donor-label permutation null on ONE simulated dataset (spec section 4 / oracle 8a).

    BYPASS: ``pbcheck.permutation.run_null`` is not called — it hard-codes the Wald test and
    always runs the naive per-cell arm, which this probe does not measure. The permutation SET
    still comes from the engine's ``build_perms`` (balanced, identity and complement removed).
    """
    oracle = null_oracle(seed=seed, **sim_kw)
    pdata_u, universe, donor_names = _aggregate(oracle.adata)
    G = len(universe)

    tmap = (oracle.adata.obs[["donor", "condition"]].astype(str).drop_duplicates()
            .set_index("donor")["condition"])
    donors = list(tmap.index)
    true_test = set(tmap.index[tmap == "disease"])
    perms = build_perms(donors, true_test, n_perm=n_perm, seed=seed)[:n_perm]

    pvals = np.full((len(perms), G), np.nan)
    ndeg = np.zeros(len(perms), dtype=int)
    diag = {}
    for i, tset in enumerate(perms):
        cond = _labels_for(donor_names, tset, "disease", "ctrl")   # engine
        res = run_test(pdata_u, universe, cond, test, n_cpus=n_cpus)
        pvals[i] = res.table["pval"].reindex(pd.Index(universe, name="gene")).to_numpy()
        ndeg[i] = res.n_significant(fdr=ALPHA)
        diag = getattr(res, "_diag", diag)
    return pvals, ndeg, G, len(perms), diag, len(donors)


def _null_fresh(sim_kw, test, n_rep, seed, n_cpus):
    """Independent-realisation null: ``n_rep`` fresh datasets, truth = 0 DE, REAL labels.

    This is the mode that can actually separate H1 from H2, because it measures each test
    against the GENERATIVE model rather than against a permutation distribution that any
    permutation-invariant test satisfies by construction.
    """
    pvals, ndeg, G, diag, n_donors = None, [], None, {}, None
    for i in range(n_rep):
        rep_seed = int((seed * 1_000_003 + i) % (2**31 - 1))
        oracle = null_oracle(seed=rep_seed, **sim_kw)
        pdata_u, universe, _ = _aggregate(oracle.adata)
        cond = np.asarray(pdata_u.obs["condition"].astype(str))
        res = run_test(pdata_u, universe, cond, test, n_cpus=n_cpus)
        p = res.table["pval"].reindex(pd.Index(universe, name="gene")).to_numpy()
        if pvals is None:
            G = len(universe)
            pvals = np.full((n_rep, G), np.nan)
            n_donors = pdata_u.n_obs
        # Universes are re-derived per replicate and can differ in size by a few genes;
        # pad/truncate to the first replicate's G only for the lambda matrix (lambda is a
        # median over genes, so a few genes at the margin cannot move it).
        m = min(G, p.size)
        pvals[i, :m] = p[:m]
        ndeg.append(res.n_significant(fdr=ALPHA))
        diag = getattr(res, "_diag", diag)
    return pvals, np.asarray(ndeg, dtype=int), G, n_rep, diag, n_donors


# ===========================================================================
# Positive oracle: power at the PRE-REGISTERED parameters
# ===========================================================================

def _power(sim_kw, test, seed, n_cpus):
    n_genes = sim_kw["n_genes"]
    k, k_note = _effective_k(n_genes)
    pos_seed = int((seed + 1000) % (2**31 - 1))
    pos = positive_oracle(n_de=k, log2fc=ORACLE_LOG2FC, seed=pos_seed, **sim_kw)
    pdata_u, universe, _ = _aggregate(pos.adata)
    cond = np.asarray(pdata_u.obs["condition"].astype(str))
    res = run_test(pdata_u, universe, cond, test, n_cpus=n_cpus)

    truth = set(pos.de_genes) & set(universe)
    hits = res.significant(fdr=ALPHA)
    sensitivity = len(hits & truth) / max(len(truth), 1)
    fp = len(hits - set(pos.de_genes))
    return {
        "power_at_pre_registered": float(sensitivity),
        "fdr_observed_on_positive_oracle": float(fp / max(len(hits), 1)),
        "positive_oracle": {
            "log2fc": ORACLE_LOG2FC,
            "k_requested": ORACLE_K,
            "k_injected": int(k),
            "k_scaling_note": k_note,
            "n_true_de_in_universe": int(len(truth)),
            "n_hits": int(len(hits)),
            "n_true_hits": int(len(hits & truth)),
            "n_false_hits": int(fp),
            "universe_size": int(len(universe)),
            "seed": pos_seed,
        },
    }


# ===========================================================================

def _versions() -> dict:
    v = {"python": platform.python_version(), "platform": platform.platform()}
    for pkg in ("pydeseq2", "decoupler", "scanpy", "anndata", "numpy", "scipy",
                "statsmodels", "pandas"):
        try:
            v[pkg] = _md.version(pkg)
        except Exception:
            v[pkg] = "unavailable"
    return v


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sigma-donor", type=float, required=True)
    ap.add_argument("--n-donors", type=int, required=True, help="donors PER GROUP")
    ap.add_argument("--test", choices=TESTS, required=True)
    ap.add_argument("--n-genes", type=int, required=True)
    ap.add_argument("--n-cells-per-donor", type=int, required=True)
    ap.add_argument("--dispersion", type=float, required=True)
    ap.add_argument("--n-perm", type=int, required=True,
                    help="permutations (null-mode=permutation) or independent replicates (fresh)")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    # Optional; defaults preserve the contracted behaviour exactly.
    ap.add_argument("--null-mode", choices=("permutation", "fresh"), default="permutation")
    ap.add_argument("--n-cpus", type=int, default=4)
    ap.add_argument("--skip-power", action="store_true",
                    help="benchmark aid: skip the positive oracle (power fields become null)")
    a = ap.parse_args(argv)

    t0 = time.perf_counter()
    sim_kw = dict(
        n_genes=a.n_genes,
        n_donors_per_group=a.n_donors,
        n_cells_per_donor=a.n_cells_per_donor,
        dispersion=a.dispersion,
        donor_sigma=a.sigma_donor,
    )

    t_null = time.perf_counter()
    runner = _null_permutation if a.null_mode == "permutation" else _null_fresh
    pvals, ndeg, G, n_used, diag, n_donors_total = runner(
        sim_kw, a.test, a.n_perm, a.seed, a.n_cpus
    )
    sec_null = time.perf_counter() - t_null

    lam = metrics.lambda_over_permutations(pvals)        # engine
    floor = metrics.perm_floor(ndeg, G)                  # engine
    per_perm = np.asarray(lam.get("per_perm", []), dtype=float)
    fp_rate = float(np.mean(ndeg >= 1))
    n_hit = int(np.sum(ndeg >= 1))

    # Exact binomial tail P(X >= n_hit | n, alpha) — how surprising the FP rate is if the arm
    # were calibrated. Computed with scipy (not asserted by hand).
    from scipy.stats import binomtest
    p_binom = float(binomtest(n_hit, n_used, ALPHA, alternative="greater").pvalue) if n_used else float("nan")

    # p-value tail profile: anti-conservativity should be monotone if it is systematic.
    flat = pvals[~np.isnan(pvals)]
    tails = {f"P(p<={t})": float(np.mean(flat <= t)) for t in (0.05, 0.01, 0.001)}

    # Discreteness floor of the exact Mann-Whitney null. BH is a STEP-UP procedure, so the
    # relevant quantity is not "min_p <= alpha/G" (that is the Bonferroni line) but how many
    # genes must sit at the floor value before p_(k) <= alpha*k/G can hold.
    if a.test == "wilcoxon":
        min_p = 2.0 / comb(2 * a.n_donors, a.n_donors)
        k_needed = int(np.ceil(min_p * G / ALPHA))
        attainable = bool(k_needed <= G)
    else:
        min_p, k_needed, attainable = None, None, True

    out = {
        "sigma_donor": a.sigma_donor,
        "n_donors_per_group": a.n_donors,
        "test": a.test,
        "n_genes": a.n_genes,
        "n_cells_per_donor": a.n_cells_per_donor,
        "dispersion": a.dispersion,
        "n_perm": n_used,
        "n_perm_requested": a.n_perm,
        "seed": a.seed,
        "null_mode": a.null_mode,
        "alpha": ALPHA,
        "universe_size_G": int(G),
        "n_donors_total": int(n_donors_total),

        "lambda_pb": float(lam["lambda"]),
        "lambda_pb_iqr": float(lam["lambda_iqr"]),
        "lambda_pb_min": float(per_perm.min()) if per_perm.size else float("nan"),
        "lambda_pb_max": float(per_perm.max()) if per_perm.size else float("nan"),
        "lambda_pb_per_perm": [float(x) for x in per_perm],
        "lambda_in_preregistered_band": bool(LAMBDA_BAND[0] <= lam["lambda"] <= LAMBDA_BAND[1]),

        "fp_rate_perm_null": fp_rate,
        "n_perm_with_rejection": n_hit,
        "fp_rate_binomial_p_greater": p_binom,
        "fp_rate_mc_se": float(np.sqrt(fp_rate * (1 - fp_rate) / n_used)) if n_used else float("nan"),
        "mean_rejections_perm_null": float(np.mean(ndeg)),
        "median_rejections_perm_null": float(np.median(ndeg)),
        "max_rejections_perm_null": float(floor["max_count"]),
        "rejections_per_perm": [int(x) for x in ndeg],
        "pvalue_tail_profile": tails,

        "bh_rejection_attainable": attainable,
        "min_attainable_pvalue": min_p,
        "bh_min_genes_at_min_p": k_needed,

        "power_at_pre_registered": None,
        "fdr_observed_on_positive_oracle": None,
        "positive_oracle": None,

        "power_target": POWER_TARGET,
        "lambda_band": list(LAMBDA_BAND),
        "test_diagnostics": diag,
        "cli": " ".join(sys.argv[1:]) if argv is None else " ".join(map(str, argv)),
        "versions": _versions(),
    }

    sec_power = 0.0
    if not a.skip_power:
        tp = time.perf_counter()
        out.update(_power(sim_kw, a.test, a.seed, a.n_cpus))
        sec_power = time.perf_counter() - tp

    out["seconds_elapsed"] = float(time.perf_counter() - t0)
    out["seconds_null"] = float(sec_null)
    out["seconds_power"] = float(sec_power)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2), encoding="utf-8")

    pw = out["power_at_pre_registered"]
    print(f"[{a.test}] mode={a.null_mode} G={G} n_perm={n_used}  "
          f"lambda_pb={out['lambda_pb']:.3f} (band {LAMBDA_BAND})  "
          f"FP={fp_rate:.3f} ({n_hit}/{n_used}, binom p={p_binom:.2g})  "
          f"mean_rej={out['mean_rejections_perm_null']:.2f}  "
          f"power={'n/a' if pw is None else f'{pw:.3f}'}  "
          f"{out['seconds_elapsed']:.1f}s -> {a.out}")
    if min_p is not None:
        print(f"  !  exact-MWU discreteness floor: min attainable p = {min_p:.3e}; BH needs "
              f">= {k_needed} genes at that floor before it can reject anything"
              + ("" if attainable else "  -- IMPOSSIBLE, this arm cannot reject at all"))
    if a.null_mode == "permutation" and a.test in ("wilcoxon", "ttest"):
        print("  !! permutation null is (near-)exact for this test BY CONSTRUCTION; "
              "use --null-mode fresh for H1-vs-H2 evidence.")
    return out


if __name__ == "__main__":
    main()
