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
  4. ``--test pooled_t`` is the EQUAL-VARIANCE two-sample Student t on the donor matrix, i.e.
     the exact d0 = 0 limit of the moderated t. It exists because ``--test ttest`` is WELCH
     and every eBayes arm POOLS the residual variance, so "ebayes beats ttest" bundles two
     changes. The grid must decompose Welch -> pooled -> moderated and credit moderation only
     with the second step.
  5. ``--test ebayes`` / ``ebayes_trend`` / ``voom`` are pure-Python re-implementations of
     limma's empirical-Bayes machinery (Smyth 2004) and of voom (Law et al. 2014). There is
     no R, no rpy2, no limma and no edgeR on this machine and no system C compiler, and
     pbcheck is a Python-native tool for the scanpy/scverse ecosystem, so an rpy2+R
     dependency is not an available option. The estimators are therefore implemented from
     the published formulas and validated numerically (four checks: d0 -> 0 reduces to the
     ordinary t; homoscedastic-null calibration by KS against t_{d+d0}; a power gain over
     the ordinary t in the homoscedastic case; recovery of a KNOWN d0). See
     :func:`_fit_f_dist`, :func:`_moderated_t`, :func:`_voom_precision_weights`.

GENERATIVE ARMS
---------------
``--gen-arm {lognormal,gamma,directnb}`` selects the donor random effect (see the block above
:func:`_make_dataset`). ``lognormal`` is the repo simulator called verbatim and is the
default, so the probe's previous behaviour is unchanged. ``gamma`` is variance-matched to it;
``directnb`` draws donor counts straight from DESeq2's own NB model with no cells and no
aggregation. Every ``--test`` works with every ``--gen-arm`` and with both ``--null-mode``
values; calibration claims must be read off ``--null-mode fresh``.

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

READING THE DIAGNOSTICS
-----------------------
``test_diagnostics`` is the LAST replicate's diagnostics only and is NOT stable across
replicates (``prior_df_d0`` has been observed to flip between +inf and ~8e2 at one fixed
configuration). ``test_diagnostics_is_last_replicate_only`` is emitted next to it as a
permanent warning. Any statement about whether moderation was actually ACTIVE must be read
off ``test_diagnostics_summary``, which carries median/min/max/mean/n for every numeric
diagnostic over all replicates and a true-fraction for every boolean. For ``prior_df_d0``,
``n_null`` counts the complete-pooling (d0 = +inf) replicates.

WHAT --sigma-het IS FOR (and what it is not)
--------------------------------------------
The pre-registered simulator gives EVERY gene the same ``donor_sigma``: it is homoscedastic
across genes, so the true prior degrees of freedom d0 are INFINITE. That is exactly the
regime empirical Bayes was designed for and the regime in which shrinkage cannot hurt, so a
calibration pass at ``--sigma-het 0`` is NOT evidence that moderation is safe -- it is a
measurement taken where the hypothesis could not have been refuted. ``--sigma-het h > 0``
makes the per-gene donor variance log-normal, which is what drives the realised d0 down to
single digits and finally stresses the shrinkage. It is DIAGNOSTIC and NOT pre-registered.
Caveat for power comparisons, MEASURED, not asserted. The knob preserves E[sigma_g^2] on the
LATENT parameter to within 3e-4 relative for h <= 1.2, but that is not the variance a
donor-level test sees on log2(CPM+1). Measured per-gene residual variance at 8v8 /
sigma_donor 0.5, h = 0 / 0.6 / 1.2 / 2.0:
    lognormal   mean 0.523 / 0.522 / 0.512 / 0.396   median 0.497 / 0.409 / 0.243 / 0.075
    gamma       mean 0.680 / 0.804 / 0.980 / 0.684   median 0.632 / 0.500 / 0.263 / 0.074
    directnb    mean 0.683 / 0.795 / 0.973 / 0.689   median 0.635 / 0.506 / 0.268 / 0.073
So the TYPICAL gene gets 2-9x quieter as h rises (and on gamma/directnb the mean rises ~43%
first). Calibration is scale-free, so FP rates and lambda remain comparable across h; POWER
is not. Compare tests against each other WITHIN one h, never across h.

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

import anndata as ad
import numpy as np
import pandas as pd
from scipy.special import digamma, polygamma
from scipy.stats import chi2, mannwhitneyu, t as t_dist, ttest_ind

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))
sys.path.insert(0, str(_REPO / "synthetic"))

import oracles  # noqa: E402  (module handle: _nb_counts is reused by the gamma generator arm)
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
TESTS = ("wald", "lrt", "ttest", "wilcoxon", "pooled_t", "ebayes", "ebayes_trend", "voom")
GEN_ARMS = ("lognormal", "gamma", "directnb")

VOOM_SPAN = 0.5                 # Law et al. 2014 voom() default span
EBAYES_TREND_DF = 4             # limma's fitFDist(covariate=) default: splinedf = min(4, ...)


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


# ---------------------------------------------------------------------------
# limma-style empirical Bayes (Smyth 2004), pure Python.
#
# There is no R / rpy2 / limma / edgeR on this machine and pbcheck is a Python-native
# tool for the scanpy ecosystem, so an rpy2 dependency is not an option. Everything
# below is a from-scratch implementation of the published estimators, validated
# numerically (see the four self-validation checks reported alongside this file).
# ---------------------------------------------------------------------------

def _trigamma_inverse(x):
    """Solve trigamma(y) = x for y > 0. Newton iteration on limma's parameterisation.

    limma::trigammaInverse iterates on y with step
        dy = trigamma(y) * (1 - trigamma(y)/x) / psigamma(y, deriv=2)
    which is Newton's method applied to the (nearly linear, convex) function
    1/trigamma(y); psigamma(y, 2) = polygamma(2, y) is the tetragamma function.

    Closed-form branches at the extremes (limma's, and both are the correct asymptotics):
      trigamma(y) ~ 1/y^2 as y -> 0   =>  x very large  =>  y = 1/sqrt(x)
      trigamma(y) ~ 1/y   as y -> inf =>  x very small  =>  y = 1/x
    Accuracy is checked numerically: max |trigamma(f(x)) / x - 1| over x in [1e-8, 1e10].
    """
    x = np.asarray(x, dtype=float)
    scalar = x.ndim == 0
    x = np.atleast_1d(x).astype(float)
    y = np.full_like(x, np.nan)
    big = x > 1e7          # trigamma(y) ~ 1/y for tiny y
    small = x < 1e-6       # trigamma(y) ~ 1/y for huge y  => y ~ 1/x
    mid = ~(big | small) & np.isfinite(x) & (x > 0)
    y[big] = 1.0 / np.sqrt(x[big])
    y[small] = 1.0 / x[small]
    if mid.any():
        yy = 0.5 + 1.0 / x[mid]
        xm = x[mid]
        for _ in range(60):
            tri = polygamma(1, yy)
            dif = tri * (1.0 - tri / xm) / polygamma(2, yy)
            yy = yy + dif
            if np.max(-dif / yy) < 1e-8:
                break
        y[mid] = yy
    return float(y[0]) if scalar else y


def _natural_spline_basis(x, df=EBAYES_TREND_DF):
    """Natural cubic regression spline basis, ``df`` columns INCLUDING the intercept.

    Truncated-power construction (Hastie/Tibshirani/Friedman eq. 5.4-5.5): with knots
    xi_1 < ... < xi_K,
        N_1 = 1, N_2 = x, N_{k+2} = d_k(x) - d_{K-1}(x),  k = 1..K-2
        d_k(x) = [ (x - xi_k)_+^3 - (x - xi_K)_+^3 ] / (xi_K - xi_k)
    which spans the natural cubic splines (linear beyond the boundary knots) in K columns.
    Knots are placed at K = df equally spaced quantiles of x.

    DEPARTURE FROM R, stated explicitly: limma calls ``splines::ns(covariate, df=4,
    intercept=TRUE)``, an orthonormalised B-spline basis. The SPAN of the two bases is the
    same space of natural cubic splines when the knot sets agree, so the OLS *fitted values*
    -- the only thing fitFDist uses -- are identical up to conditioning; R's ns() places
    interior knots at quantiles and boundary knots at the data range, which is not
    bit-identical to the quantile grid used here. This is a smoother-choice difference, not a
    different estimator.
    """
    x = np.asarray(x, dtype=float)
    # Affine standardisation: the space of natural cubic splines is invariant under
    # x -> (x - m)/s together with the same map on the knots, so this only helps conditioning.
    m, s = float(np.mean(x)), float(np.std(x))
    z = (x - m) / (s if s > 0 else 1.0)
    K = max(int(df), 2)
    knots = np.unique(np.quantile(z, np.linspace(0.0, 1.0, K)))
    if knots.size < 3:
        return np.column_stack([np.ones_like(z), z])
    K = knots.size

    def d(k):
        num = np.maximum(z - knots[k], 0.0) ** 3 - np.maximum(z - knots[K - 1], 0.0) ** 3
        return num / (knots[K - 1] - knots[k])

    cols = [np.ones_like(z), z] + [d(k) - d(K - 2) for k in range(K - 2)]
    return np.column_stack(cols)


def _fit_f_dist(s2, df, covariate=None, spline_df=EBAYES_TREND_DF):
    """Smyth (2004) fitFDist: method of moments on log(s2) via digamma/trigamma.

    MODEL.  s2_g | sigma2_g ~ sigma2_g * chisq_d / d  and  d0 * s02 / sigma2_g ~ chisq_d0,
    so marginally s2_g ~ s02 * F(d, d0).

    ESTIMATOR (derived, not copied).  Write z_g = log s2_g. Because
    log(chisq_d / d) has mean digamma(d/2) - log(d/2) and variance trigamma(d/2), and
    log(sigma2_g) = log(d0 s02) - log(chisq_d0) has mean log(s02) + log(d0/2) - digamma(d0/2)
    and variance trigamma(d0/2),

        E[z_g]   = log(s02) + digamma(d/2) - log(d/2) + log(d0/2) - digamma(d0/2)
        Var[z_g] = trigamma(d/2) + trigamma(d0/2)

    Define e_g = z_g - digamma(d/2) + log(d/2). Then, with n usable genes,

        evar = sum_g (e_g - ebar_g)^2 / (n - rank)        [rank = 1, or the spline rank]
        trigamma(d0/2) = evar - trigamma(d/2)
        d0  = 2 * trigammaInverse(evar - trigamma(d/2))
        s02 = exp(ebar + digamma(d0/2) - log(d0/2))

    If evar - trigamma(d/2) <= 0 the between-gene variance is not distinguishable from the
    within-gene sampling variance, and the correct limit is d0 = +inf (complete pooling),
    s02 = exp(ebar). This is limma's behaviour and is reported, never silently clipped.

    ``covariate`` (limma-trend): ebar becomes a per-gene FITTED value from an OLS regression
    of e on a natural cubic spline in the covariate, so s02 becomes a per-gene TREND. d0 stays
    scalar, exactly as in limma.

    Returns ``(d0, s02, info)`` where s02 is a scalar (no covariate) or a length-G vector.
    """
    s2 = np.asarray(s2, dtype=float)
    G = s2.size
    ok = np.isfinite(s2) & (s2 > 0)
    n = int(ok.sum())
    if n < 3:
        raise ValueError(f"fitFDist needs >= 3 genes with s2 > 0, got {n}")
    z = np.log(s2[ok])
    e_ok = z - digamma(df / 2.0) + log(df / 2.0)

    if covariate is None:
        rank = 1
        ebar_ok = np.full(n, float(np.mean(e_ok)))
        ebar_all = np.full(G, float(np.mean(e_ok)))
    else:
        cov = np.asarray(covariate, dtype=float)
        B_all = _natural_spline_basis(cov, df=spline_df)
        B_ok = B_all[ok, :]
        coef, _, rank, _ = np.linalg.lstsq(B_ok, e_ok, rcond=None)
        rank = int(rank)
        ebar_ok = B_ok @ coef
        ebar_all = B_all @ coef

    resid = e_ok - ebar_ok
    denom = max(n - rank, 1)
    evar = float(np.sum(resid ** 2) / denom)
    evar_between = evar - float(polygamma(1, df / 2.0))

    if evar_between > 0:
        d0 = 2.0 * _trigamma_inverse(evar_between)
        s02 = np.exp(ebar_all + digamma(d0 / 2.0) - log(d0 / 2.0))
    else:
        d0 = np.inf
        s02 = np.exp(ebar_all)

    if covariate is None:
        s02 = float(s02[0])
    info = {
        "d0": float(d0),
        "evar_total": evar,
        "evar_between_genes": float(evar_between),
        "trigamma_d_over_2": float(polygamma(1, df / 2.0)),
        "n_genes_used_in_prior": n,
        "n_genes_dropped_zero_var": int(G - n),
        "prior_rank": int(rank if covariate is not None else 1),
    }
    return d0, s02, info


def _wls_two_group(Y, x, weights=None, return_fitted=False):
    """Per-gene (weighted) least squares for the design ``~ 1 + x`` with x a 0/1 indicator.

    ``Y`` is (n_donors, n_genes); ``weights`` is None (ordinary LS) or (n_donors, n_genes).
    Returns ``(beta, s2, d, unscaled_var)``:
      beta          contrast estimate (group x=1 minus group x=0), length G
      s2            residual variance sum_i w_i e_i^2 / d, length G
      d             residual degrees of freedom = n - 2 (scalar)
      unscaled_var  (X' W X)^{-1}[1,1], length G  -> se(beta) = sqrt(s2 * unscaled_var)

    Sanity anchor: with weights = 1 this reduces to beta = mean_A - mean_B,
    unscaled_var = 1/n_A + 1/n_B, and s2 the pooled variance -- i.e. Student's equal-variance
    two-sample t. That identity is checked numerically in the validation script.
    """
    Y = np.asarray(Y, dtype=float)
    x = np.asarray(x, dtype=float).reshape(-1, 1)
    n = Y.shape[0]
    W = np.ones_like(Y) if weights is None else np.asarray(weights, dtype=float)

    S0 = W.sum(axis=0)                 # sum w
    S1 = (W * x).sum(axis=0)           # sum w x  (= sum w x^2, x is 0/1)
    Sy = (W * Y).sum(axis=0)           # sum w y
    Sxy = (W * x * Y).sum(axis=0)      # sum w x y
    det = S0 * S1 - S1 ** 2            # = S1 * (S0 - S1)
    with np.errstate(divide="ignore", invalid="ignore"):
        beta0 = (S1 * Sy - S1 * Sxy) / det          # intercept  (weighted mean of x=0 group)
        beta1 = (S0 * Sxy - S1 * Sy) / det          # contrast
        unscaled_var = S0 / det                     # (X'WX)^{-1}[1,1]
        fitted = beta0[None, :] + beta1[None, :] * x
        resid = Y - fitted
        d = n - 2
        s2 = (W * resid ** 2).sum(axis=0) / d
    if return_fitted:
        return beta1, s2, d, unscaled_var, fitted
    return beta1, s2, d, unscaled_var


def _moderated_t(beta, s2, d, unscaled_var, d0, s02):
    """Smyth's moderated t and its two-sided p-value.

        s~_g^2 = (d0 * s02_g + d * s2_g) / (d0 + d)
        t_g    = beta_g / sqrt(s~_g^2 * unscaled_var_g)
        p_g    = 2 * P(T_{d + d0} > |t_g|)

    d0 = +inf is handled in the limit: s~^2 -> s02 and T_{inf} -> N(0, 1). Passing d0 = 0
    reproduces the ordinary (pooled, equal-variance) t exactly -- that is `--test pooled_t`.

    A gene with a non-finite s2 is UNTESTABLE and must come out as p = NaN in every branch, so
    that ``mtc.bh`` drops it instead of publishing a number derived entirely from the prior.
    See the F1 comment in the d0 = +inf branch: that was a live silent-p-value bug.
    """
    s2 = np.asarray(s2, dtype=float)
    s02_arr = np.broadcast_to(np.asarray(s02, dtype=float), s2.shape)
    finite_s2 = np.isfinite(s2)
    if np.isinf(d0):
        # FIX F1 (silent p-value bug). The posterior variance in the complete-pooling limit is
        # the prior, which DISCARDS the gene's own s2 entirely. Before this fix a gene whose s2
        # was NaN -- reachable through the voom arm, which writes NaN precision weights whenever
        # its lowess mean-variance curve is non-positive -- silently received a finite,
        # valid-looking p-value computed purely from the prior, and n_pval_nan stayed 0. NaN must
        # propagate here exactly as it already does in the finite-d0 branch below, so that an
        # untestable gene is dropped from BH rather than given a fabricated p-value.
        s2_post = np.where(finite_s2, s02_arr.astype(float), np.nan)
        df_total = np.inf
    elif d0 == 0:
        # Unmoderated limit, written out exactly rather than as (0*s02 + d*s2)/(0+d): this is
        # the `--test pooled_t` baseline (equal-variance two-sample Student t) and it must be
        # bit-identical to it, not merely equal to rounding. s02 is unused here by construction.
        s2_post = s2.astype(float, copy=True)
        df_total = float(d)
    else:
        s2_post = (d0 * s02_arr + d * s2) / (d0 + d)
        df_total = d + d0
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(s2_post * np.asarray(unscaled_var, dtype=float))
        tstat = np.asarray(beta, dtype=float) / se
        if np.isinf(df_total):
            from scipy.stats import norm
            pval = 2.0 * norm.sf(np.abs(tstat))
        else:
            pval = 2.0 * t_dist.sf(np.abs(tstat), df_total)
    return tstat, np.asarray(pval, dtype=float), s2_post, df_total


def _voom_precision_weights(counts, x):
    """voom (Law, Chen, Shi & Smyth 2014, Genome Biology 15:R29) precision weights.

    Reproduced from the published algorithm, on the DONOR pseudobulk count matrix (which is
    exactly where a limma-voom pseudobulk pipeline applies it):

      1. lib_d   = sum_g counts[d, g]
         y[d, g] = log2( (counts[d, g] + 0.5) / (lib_d + 1) * 1e6 )        (voom's offset)
      2. OLS of y on the full design ~ 1 + x  ->  fitted values, residual sd sigma_g
      3. sx_g = mean_d y[d, g] + mean_d log2(lib_d + 1) - log2(1e6)        (mean log2-count)
         sy_g = sqrt(sigma_g)                                              (sqrt of the sd)
      4. lowess(sy ~ sx, f = 0.5) -> a monotone-ish mean-variance curve, then linear
         interpolation with constant extrapolation (R's approxfun(rule = 2))
      5. per-observation predicted log2-count
             lam[d, g] = fitted[d, g] + log2(lib_d + 1) - log2(1e6)
         and weights w[d, g] = 1 / curve(lam[d, g])^4
         (^4 because the curve models sqrt(sd), so curve^4 is a variance)

    Returns ``(y, w, info)``. ``y`` is voom's log-CPM (with the +0.5/+1 offsets), NOT the
    log2(CPM+1) used by the other donor-level tests -- that difference is voom's own definition.

    R's lowess and statsmodels' lowess are both Cleveland's LOWESS with the same default 3
    robustness iterations; the tricube local-linear fit is the same algorithm. This is the one
    place where the numbers can differ from R in the last digits.
    """
    from statsmodels.nonparametric.smoothers_lowess import lowess

    C = np.asarray(counts, dtype=float)
    lib = C.sum(axis=1)
    y = np.log2((C + 0.5) / (lib[:, None] + 1.0) * 1e6)

    # Step 2: ordinary LS on the FULL design; `fitted` comes from the same solve, so the
    # weights cannot drift away from the fit they are derived from.
    _, s2, _, _, fitted = _wls_two_group(y, x, weights=None, return_fitted=True)
    sigma = np.sqrt(np.maximum(s2, 0.0))

    amean = y.mean(axis=0)
    sx = amean + float(np.mean(np.log2(lib + 1.0))) - log(1e6) / log(2.0)
    sy = np.sqrt(sigma)
    good = np.isfinite(sx) & np.isfinite(sy) & (C.sum(axis=0) > 0) & (sigma > 0)
    if good.sum() < 10:
        raise ValueError("voom: fewer than 10 usable genes for the mean-variance trend")

    lo = lowess(sy[good], sx[good], frac=VOOM_SPAN, it=3, return_sorted=True)
    lx, ly = lo[:, 0], lo[:, 1]
    # R's approxfun(..., rule = 2): linear interpolation, constant outside the range.
    lam = fitted + np.log2(lib + 1.0)[:, None] - log(1e6) / log(2.0)
    curve = np.interp(lam, lx, ly, left=ly[0], right=ly[-1])
    n_extrap = int(np.sum((lam < lx[0]) | (lam > lx[-1])))
    with np.errstate(divide="ignore", invalid="ignore"):
        w = 1.0 / curve ** 4
    w[~np.isfinite(w)] = np.nan
    info = {
        "voom_span": VOOM_SPAN,
        "voom_n_genes_in_trend": int(good.sum()),
        "voom_n_obs_extrapolated": n_extrap,
        "voom_curve_min": float(ly.min()),
        "voom_curve_max": float(ly.max()),
        "voom_weight_median": float(np.nanmedian(w)),
        "voom_weight_iqr": float(np.nanpercentile(w, 75) - np.nanpercentile(w, 25)),
    }
    return y, w, info


def _ebayes_result(pdata_u, cond, test_level, ref_level, kind) -> DEResult:
    """``pooled_t`` / ``ebayes`` / ``ebayes_trend`` / ``voom`` on the SAME aggregated donor matrix.

    All four consume ``pdata_u`` -- already restricted to the frozen universe by
    :func:`_aggregate` -- so they see byte-identical data to wald / lrt / ttest / wilcoxon.

      pooled_t      OLS on log2(CPM + 1), EQUAL-VARIANCE (pooled) two-sample Student t, i.e.
                    exactly the d0 = 0 limit of the moderated t. NO moderation at all.
      ebayes        OLS on log2(CPM + 1) [the same _log_cpm the ttest arm uses]  + eBayes,
                    single global prior (no mean-variance trend).
      ebayes_trend  identical, but the prior variance is a spline trend in average log-CPM.
      voom          voom's own log-CPM with precision weights -> WLS -> eBayes (flat prior),
                    which is the standard limma-voom pipeline.

    FIX F2 (power attribution). ``--test ttest`` is WELCH (unequal-variance, Satterthwaite df)
    while every eBayes arm pools the residual variance across the two groups, so an
    "ebayes vs ttest" gap silently bundles two changes: Welch -> pooled, and unmoderated ->
    moderated. ``pooled_t`` is the missing middle term; the grid must read the decomposition
    Welch (ttest) -> pooled Student (pooled_t) -> moderated (ebayes) and attribute only the
    second step to moderation.
    """
    cond = np.asarray(cond, dtype=object)
    x = (cond == test_level).astype(float)
    nA, nB = int(x.sum()), int((cond == ref_level).sum())
    if nA < 2 or nB < 2:
        raise ValueError(f"need >= 2 donors per group, got {nA} vs {nB}")
    if nA + nB != len(cond):
        raise ValueError("condition vector contains levels other than test/ref")

    diag = {}
    if kind == "voom":
        counts = np.asarray(pdata_u.X, dtype=float).round()
        Y, W, vinfo = _voom_precision_weights(counts, x)
        diag.update(vinfo)
        beta, s2, d, uvar = _wls_two_group(Y, x, weights=W)
        covariate = None
        name = "donor[voom + limma eBayes]"
    else:
        Y = _log_cpm(pdata_u)                       # identical transform to --test ttest
        beta, s2, d, uvar = _wls_two_group(Y, x, weights=None)
        covariate = Y.mean(axis=0) if kind == "ebayes_trend" else None
        name = {"pooled_t": "donor[pooled Student t on log2CPM, UNMODERATED]",
                "ebayes": "donor[limma eBayes on log2CPM, no trend]",
                "ebayes_trend": "donor[limma eBayes on log2CPM, mean-variance trend]"}[kind]

    # The prior is fitted for EVERY arm, but it is APPLIED only by the moderated ones. For
    # pooled_t the fit is diagnostic-only: it is what tells the grid how heterogeneous the
    # gene-level variances actually were in the very same dataset the unmoderated test saw.
    d0_fit, s02_fit, finfo = _fit_f_dist(s2, d, covariate=covariate, spline_df=EBAYES_TREND_DF)
    if kind == "pooled_t":
        d0, s02 = 0.0, 0.0          # d0 = 0 => s~^2 = s2 exactly; s02 is unused (see _moderated_t)
    else:
        d0, s02 = d0_fit, s02_fit
    tstat, pval, s2_post, df_total = _moderated_t(beta, s2, d, uvar, d0, s02)

    # d0 = +inf is a legitimate, and here a COMMON, estimate (complete pooling: the between-gene
    # variance is not separable from sampling noise). JSON has no infinity literal that `jq` will
    # parse, so it is emitted as null with an explicit boolean beside it. A runner must read
    # `prior_is_complete_pooling`: when it is true EVERY gene was given the same variance, which
    # is the maximum-shrinkage regime and the one the anti-conservativity hypothesis is about.
    d0_finite = bool(np.isfinite(d0))
    d0_fit_finite = bool(np.isfinite(d0_fit))
    diag.update(finfo)
    # finfo["d0"] is the RAW fit and may be +inf, which json.dumps writes as the bare literal
    # `Infinity` (not parseable by jq). Overwrite both with the null-encoded form.
    diag["d0"] = float(d0) if d0_finite else None
    diag.update({
        "prior_df_d0": float(d0) if d0_finite else None,
        "prior_is_complete_pooling": (not d0_finite),
        "residual_df_d": int(d),
        "total_df": float(df_total) if np.isfinite(df_total) else None,
        "s0_squared": (float(s02) if np.isscalar(s02) or np.ndim(s02) == 0
                       else float(np.median(s02))),
        "s0_squared_is_trended": covariate is not None,
        "median_s2_raw": float(np.nanmedian(s2)),
        "median_s2_posterior": float(np.nanmedian(s2_post)),
        # 1.0 = the posterior variance is entirely the prior (no gene-specific information left);
        # 0.0 = no shrinkage at all (the unmoderated pooled_t baseline).
        "shrinkage_factor_d0_over_d0_plus_d": (1.0 if not d0_finite else float(d0 / (d0 + d))),
        "moderation_applied": (kind != "pooled_t"),
        # For pooled_t: the prior that WOULD have been used, so the grid can read the realised
        # gene-level variance heterogeneity off the unmoderated arm too. None = +inf.
        # None here means EITHER "not applicable" (a moderated arm, which applied its d0) OR
        # "+inf"; the boolean beside it disambiguates.
        "d0_fitted_but_not_applied": (None if kind != "pooled_t"
                                      else (float(d0_fit) if d0_fit_finite else None)),
        "d0_fitted_but_not_applied_is_infinite": (kind == "pooled_t" and not d0_fit_finite),
        # FIX F1 monitoring: how many genes ended up untestable. mtc.bh drops NaN p-values, so a
        # nonzero count here means this arm's BH denominator is smaller than the universe.
        "n_genes_nonfinite_s2": int(np.sum(~np.isfinite(np.asarray(s2, dtype=float)))),
        "n_pval_nan": int(np.sum(~np.isfinite(np.asarray(pval, dtype=float)))),
    })

    # beta is the contrast on whichever log2 scale the arm used, so it IS a log2 fold change.
    table = pd.DataFrame(
        {"pval": pval, "padj": np.nan, "log2fc": beta},
        index=pd.Index(list(pdata_u.var_names), name="gene"),
    )
    res = DEResult(method=name, table=table, contrast=("condition", test_level, ref_level))
    res._diag = diag
    return res


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
    elif test in ("pooled_t", "ebayes", "ebayes_trend", "voom"):
        # STRUCTURAL GUARANTEE, enforced HERE: the moderated arms receive the very same
        # `pdata_u` object (already universe-restricted by _aggregate) and the very same
        # `cond` vector as wald/lrt/ttest/wilcoxon. No arm re-aggregates, re-filters or
        # re-simulates; the only thing that varies across --test is the estimator.
        res = _ebayes_result(pdata_u, cond, test_level, ref_level, test)
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


# ---------------------------------------------------------------------------
# GENERATIVE ARMS (--gen-arm). Previously these lived only in a scratch script
# (scratchpad/pbprobe/gen_lens.py); they are wired in here so that EVERY --test,
# including the new moderated ones, can be judged on more than our own simulator.
#
#   lognormal  the repo simulator, called through synthetic.oracles VERBATIM. This is the
#              default and is byte-identical to the probe's previous behaviour.
#   gamma      identical cell-level model except the donor random effect is Gamma with the
#              SAME mean (1) and SAME variance (e^{s^2} - 1); the donor pseudobulk is then
#              essentially exactly NB, i.e. DESeq2 is well specified.
#   directnb   donor counts drawn straight from the NB model with no cells and no
#              aggregation at all: Y_gd ~ NB(mu_g * S_d, alpha_d) with
#              alpha_d = (e^{s^2} - 1) + dispersion * e^{s^2} * Q_d / S_d^2,
#              S_d = sum_c depth_c, Q_d = sum_c depth_c^2 -- the moment-matched donor-level
#              law implied by the cell-level arms.
# ---------------------------------------------------------------------------

def _per_gene_sigma(donor_sigma, sigma_het, n_genes, rng):
    """Per-gene donor SD. ``sigma_het = 0`` (the default) returns the constant pre-registered
    value AND CONSUMES NO RNG, so every default result is bit-identical to the previous probe.

    For ``sigma_het = h > 0`` the per-gene donor VARIANCE is log-normal about the pre-registered
    value, mean-preserving on the variance scale:

        sigma_g^2 = donor_sigma^2 * exp(N(-h^2/2, h^2))   =>   E[sigma_g^2] = donor_sigma^2

    This exists for exactly one reason. The pre-registered simulator gives EVERY gene the same
    donor-level variance, i.e. it is homoscedastic across genes on the log scale -- which is
    precisely the regime empirical-Bayes moderation was designed for and in which shrinkage
    cannot hurt. A calibration pass there is therefore NOT evidence that a moderated method is
    safe. ``--sigma-het`` is the knob that makes the true donor variance heterogeneous, which is
    the condition under which the shrinkage hypothesis predicts anti-conservative tails.
    It is a DIAGNOSTIC knob and is NOT part of the pre-registered generative model (spec 8b).
    """
    if sigma_het <= 0:
        return np.full(n_genes, float(donor_sigma))
    var = (donor_sigma ** 2) * np.exp(
        rng.normal(-0.5 * sigma_het ** 2, sigma_het, size=n_genes)
    )
    return np.sqrt(var)


def _sim_cell_level(*, n_genes, n_donors_per_group, n_cells_per_donor, dispersion,
                    donor_sigma, re_dist, n_de=0, log2fc=1.0, seed=0, sigma_het=0.0,
                    mean_log_mu=1.0, mean_log_sigma=1.2, depth_log_sigma=0.3):
    """oracles.simulate() with the donor-RE distribution swappable. Returns (adata, de_genes)."""
    rng = np.random.default_rng(seed)
    n_donors = 2 * n_donors_per_group
    donors = [f"d{i:02d}" for i in range(n_donors)]
    groups = np.array(["ctrl"] * n_donors_per_group + ["disease"] * n_donors_per_group)
    genes = [f"g{j:04d}" for j in range(n_genes)]

    mu_g = np.exp(rng.normal(mean_log_mu, mean_log_sigma, size=n_genes))
    cond_fc = np.ones(n_genes)
    de_idx = rng.choice(n_genes, size=n_de, replace=False) if n_de > 0 else np.array([], dtype=int)
    signs = rng.choice([-1.0, 1.0], size=len(de_idx))
    cond_fc[de_idx] = 2.0 ** (signs * log2fc)
    de_genes = [genes[i] for i in de_idx]

    sig = _per_gene_sigma(donor_sigma, sigma_het, n_genes, rng)   # no RNG use when sigma_het = 0
    v = np.exp(sig ** 2) - 1.0                                    # per-gene RE variance
    if donor_sigma > 0:
        if re_dist == "lognormal":
            donor_re = np.exp(rng.normal(-0.5 * sig[:, None] ** 2, sig[:, None],
                                         size=(n_genes, n_donors)))
        elif re_dist == "gamma":
            donor_re = rng.gamma(shape=1.0 / v[:, None],
                                 scale=np.broadcast_to(v[:, None], (n_genes, n_donors)))
        else:
            raise ValueError(re_dist)
    else:
        donor_re = np.ones((n_genes, n_donors))

    blocks, obs_donor, obs_group = [], [], []
    for d_idx, donor in enumerate(donors):
        fc = cond_fc if groups[d_idx] == "disease" else np.ones(n_genes)
        donor_mean = mu_g * donor_re[:, d_idx] * fc
        depth = np.exp(rng.normal(0.0, depth_log_sigma, size=n_cells_per_donor))
        blocks.append(oracles._nb_counts(depth[:, None] * donor_mean[None, :], dispersion, rng))
        obs_donor.extend([donor] * n_cells_per_donor)
        obs_group.extend([groups[d_idx]] * n_cells_per_donor)

    X = np.vstack(blocks).astype(np.float32)
    obs = pd.DataFrame({
        "donor": pd.Categorical(obs_donor),
        "condition": pd.Categorical(obs_group, categories=["ctrl", "disease"]),
        "cell_type": pd.Categorical(["T_cell"] * X.shape[0]),
    })
    adata = ad.AnnData(X=X, obs=obs)
    adata.var_names = genes
    adata.obs_names = [f"c{i:06d}" for i in range(X.shape[0])]
    return adata, de_genes


def _sim_directnb(*, n_genes, n_donors_per_group, n_cells_per_donor, dispersion,
                  donor_sigma, n_de=0, log2fc=1.0, seed=0, sigma_het=0.0,
                  mean_log_mu=1.0, mean_log_sigma=1.2, depth_log_sigma=0.3):
    """Donor-level counts from DESeq2's own exact-NB model. Returns (pdata, de_genes)."""
    rng = np.random.default_rng(seed)
    n_donors = 2 * n_donors_per_group
    donors = [f"d{i:02d}" for i in range(n_donors)]
    groups = ["ctrl"] * n_donors_per_group + ["disease"] * n_donors_per_group
    genes = [f"g{j:04d}" for j in range(n_genes)]
    mu_g = np.exp(rng.normal(mean_log_mu, mean_log_sigma, size=n_genes))
    cond_fc = np.ones(n_genes)
    de_idx = rng.choice(n_genes, size=n_de, replace=False) if n_de > 0 else np.array([], dtype=int)
    signs = rng.choice([-1.0, 1.0], size=len(de_idx))
    cond_fc[de_idx] = 2.0 ** (signs * log2fc)
    de_genes = [genes[i] for i in de_idx]
    sig = _per_gene_sigma(donor_sigma, sigma_het, n_genes, rng)   # no RNG use when sigma_het = 0
    w = np.exp(sig ** 2)                                          # per-gene, scalar when het = 0

    rows = []
    for d_idx in range(n_donors):
        fc = cond_fc if groups[d_idx] == "disease" else np.ones(n_genes)
        depth = np.exp(rng.normal(0.0, depth_log_sigma, size=n_cells_per_donor))
        S, Q = depth.sum(), (depth ** 2).sum()
        alpha = (w - 1.0) + dispersion * w * Q / S ** 2           # per-gene NB dispersion
        mean = mu_g * fc * S
        if np.all(alpha <= 0):
            rows.append(rng.poisson(mean))
        else:
            rows.append(rng.poisson(rng.gamma(shape=1.0 / alpha, scale=mean * alpha)))
    X = np.vstack(rows).astype(np.float32)
    obs = pd.DataFrame(
        {"donor": donors,
         "condition": pd.Categorical(groups, categories=["ctrl", "disease"]),
         "cell_type": ["T_cell"] * n_donors},
        index=donors,
    )
    pdata = ad.AnnData(X=X, obs=obs)
    pdata.var_names = genes
    return pdata, de_genes


def _make_dataset(sim_kw, gen_arm, seed, *, n_de=0, log2fc=ORACLE_LOG2FC, sigma_het=0.0):
    """One simulated dataset -> (pdata_u, universe, donor_names, de_genes).

    The default arm ('lognormal') goes through synthetic.oracles VERBATIM, so the probe's
    pre-existing behaviour and RNG stream are unchanged. Aggregation/universe freezing happen
    exactly ONCE per dataset here, for every arm and every test alike.
    """
    if gen_arm == "lognormal" and sigma_het <= 0:
        # Default path: the engine's own simulator, called VERBATIM. Unchanged.
        oracle = (positive_oracle(n_de=n_de, log2fc=log2fc, seed=seed, **sim_kw) if n_de > 0
                  else null_oracle(seed=seed, **sim_kw))
        pdata_u, universe, donor_names = _aggregate(oracle.adata)
        return pdata_u, universe, donor_names, list(oracle.de_genes)
    if gen_arm in ("lognormal", "gamma"):
        # _sim_cell_level(re_dist='lognormal', sigma_het=0) is a bit-for-bit clone of
        # oracles.simulate (verified: identical X and identical de_genes at n_de = 0 and 200).
        adata, de_genes = _sim_cell_level(re_dist=gen_arm, n_de=n_de, log2fc=log2fc,
                                          seed=seed, sigma_het=sigma_het, **sim_kw)
        pdata_u, universe, donor_names = _aggregate(adata)
        return pdata_u, universe, donor_names, de_genes
    if gen_arm == "directnb":
        # No cells exist in this arm, so there is nothing to aggregate: the donor matrix IS
        # the simulated object. Universe freezing and restriction are still done once, by the
        # engine's frozen_universe, exactly as in _aggregate.
        pdata, de_genes = _sim_directnb(n_de=n_de, log2fc=log2fc, seed=seed,
                                        sigma_het=sigma_het, **sim_kw)
        universe = frozen_universe(pdata)
        pdata_u = pdata[:, universe].copy()
        donor_names = pd.Index(pdata_u.obs["donor"].astype(str))
        return pdata_u, universe, donor_names, de_genes
    raise ValueError(f"unknown gen arm {gen_arm!r}; expected one of {GEN_ARMS}")


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

def _summarise_diags(diags):
    """FIX F3. ``test_diagnostics`` in the JSON is the LAST replicate's diagnostics only, and it
    is not stable: at the pre-registered 8v8 / sigma 0.5 log-normal setting ``prior_df_d0`` came
    out +inf on some replicates and finite (order 1e2-1e3) on others, so whether a reader
    concludes "moderation was fully active" or "moderation was inert" depended on which replicate
    happened to be last. The grid's conclusion turns on that, so every numeric diagnostic is
    summarised over ALL replicates here (median, min, max, mean, n), booleans as a true-fraction.

    ``prior_df_d0`` is emitted as ``null`` when d0 = +inf, so for that key ``n_null`` is exactly
    the number of complete-pooling replicates (and equals ``prior_is_complete_pooling.n_true``).
    """
    diags = [d for d in diags if isinstance(d, dict) and d]
    if not diags:
        return {}
    keys = []
    for dd in diags:
        for k in dd:
            if k not in keys:
                keys.append(k)
    out = {"n_replicates_summarised": len(diags)}
    for k in keys:
        vals = [dd.get(k, None) for dd in diags]
        present = [v for v in vals if v is not None]
        n_null = len(vals) - len(present)
        if present and all(isinstance(v, bool) for v in present):
            n_true = int(sum(bool(v) for v in present))
            out[k] = {"n_true": n_true, "n_false": len(present) - n_true,
                      "frac_true": float(n_true / len(present)), "n_null": n_null}
        elif present and all(isinstance(v, (int, float, np.integer, np.floating))
                             and not isinstance(v, bool) for v in present):
            arr = np.asarray([float(v) for v in present], dtype=float)
            fin = arr[np.isfinite(arr)]
            out[k] = {
                "median": float(np.median(fin)) if fin.size else None,
                "min": float(fin.min()) if fin.size else None,
                "max": float(fin.max()) if fin.size else None,
                "mean": float(fin.mean()) if fin.size else None,
                "n": int(fin.size),
                "n_null": n_null,
                "n_nonfinite": int(arr.size - fin.size),
            }
        else:
            out[k] = {"distinct_values": sorted({str(v) for v in present})[:10],
                      "n": len(present), "n_null": n_null}
    return out


def _null_permutation(sim_kw, test, n_perm, seed, n_cpus, gen_arm="lognormal", sigma_het=0.0):
    """Donor-label permutation null on ONE simulated dataset (spec section 4 / oracle 8a).

    BYPASS: ``pbcheck.permutation.run_null`` is not called — it hard-codes the Wald test and
    always runs the naive per-cell arm, which this probe does not measure. The permutation SET
    still comes from the engine's ``build_perms`` (balanced, identity and complement removed).
    """
    pdata_u, universe, donor_names, _ = _make_dataset(sim_kw, gen_arm, seed, sigma_het=sigma_het)
    G = len(universe)

    tmap = (pdata_u.obs[["donor", "condition"]].astype(str).drop_duplicates()
            .set_index("donor")["condition"])
    donors = list(tmap.index)
    true_test = set(tmap.index[tmap == "disease"])
    perms = build_perms(donors, true_test, n_perm=n_perm, seed=seed)[:n_perm]

    pvals = np.full((len(perms), G), np.nan)
    ndeg = np.zeros(len(perms), dtype=int)
    diag, diags, nan_counts = {}, [], []
    for i, tset in enumerate(perms):
        cond = _labels_for(donor_names, tset, "disease", "ctrl")   # engine
        res = run_test(pdata_u, universe, cond, test, n_cpus=n_cpus)
        p_i = res.table["pval"].reindex(pd.Index(universe, name="gene")).to_numpy()
        pvals[i] = p_i
        ndeg[i] = res.n_significant(fdr=ALPHA)
        diag = getattr(res, "_diag", diag)
        diags.append(getattr(res, "_diag", None))         # FIX F3: keep EVERY replicate
        nan_counts.append(int(np.sum(~np.isfinite(np.asarray(p_i, dtype=float)))))
    return pvals, ndeg, G, len(perms), diag, len(donors), diags, nan_counts


def _null_fresh(sim_kw, test, n_rep, seed, n_cpus, gen_arm="lognormal", sigma_het=0.0):
    """Independent-realisation null: ``n_rep`` fresh datasets, truth = 0 DE, REAL labels.

    This is the mode that can actually separate H1 from H2, because it measures each test
    against the GENERATIVE model rather than against a permutation distribution that any
    permutation-invariant test satisfies by construction.
    """
    pvals, ndeg, G, diag, n_donors = None, [], None, {}, None
    diags, nan_counts = [], []
    for i in range(n_rep):
        rep_seed = int((seed * 1_000_003 + i) % (2**31 - 1))
        pdata_u, universe, _, _ = _make_dataset(sim_kw, gen_arm, rep_seed, sigma_het=sigma_het)
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
        diags.append(getattr(res, "_diag", None))         # FIX F3: keep EVERY replicate
        # counted on THIS replicate's own universe, before the pad/truncate above
        nan_counts.append(int(np.sum(~np.isfinite(np.asarray(p, dtype=float)))))
    return pvals, np.asarray(ndeg, dtype=int), G, n_rep, diag, n_donors, diags, nan_counts


# ===========================================================================
# Positive oracle: power at the PRE-REGISTERED parameters
# ===========================================================================

def _power(sim_kw, test, seed, n_cpus, gen_arm="lognormal", sigma_het=0.0, n_reps=1):
    """Sensitivity at the PRE-REGISTERED oracle (log2FC = 1.0, K = 200, BH < 0.05).

    ``n_reps = 1`` (the default) is the historical behaviour EXACTLY: one positive-oracle
    dataset at ``seed + 1000``, and every emitted key is unchanged. That single draw is,
    however, one realisation of a random variable whose across-dataset SD was measured at
    0.046-0.062 at 8v8 -- roughly +-0.10 at 95% against a BINDING 0.60 gate -- so any run whose
    power number is going to be compared with POWER_TARGET must pass ``--n-power-reps`` > 1.
    Replicate j uses seed ``seed + 1000 + j``; ``power_at_pre_registered`` is then the MEAN over
    replicates and ``power_se`` is its Monte-Carlo standard error.
    """
    n_genes = sim_kw["n_genes"]
    k, k_note = _effective_k(n_genes)
    base = int((seed + 1000) % (2**31 - 1))
    sens_all, fdr_all, reps = [], [], []
    first = None
    for j in range(max(int(n_reps), 1)):
        pos_seed = int((base + j) % (2**31 - 1))
        pdata_u, universe, _, de_genes = _make_dataset(
            sim_kw, gen_arm, pos_seed, n_de=k, log2fc=ORACLE_LOG2FC, sigma_het=sigma_het
        )
        cond = np.asarray(pdata_u.obs["condition"].astype(str))
        res = run_test(pdata_u, universe, cond, test, n_cpus=n_cpus)

        truth = set(de_genes) & set(universe)
        hits = res.significant(fdr=ALPHA)
        sensitivity = len(hits & truth) / max(len(truth), 1)
        fp = len(hits - set(de_genes))
        detail = {
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
        }
        if first is None:
            first = detail
        sens_all.append(float(sensitivity))
        fdr_all.append(float(fp / max(len(hits), 1)))
        reps.append({"seed": pos_seed, "sensitivity": float(sensitivity),
                     "true_fdp": float(fp / max(len(hits), 1)), "n_hits": int(len(hits))})

    s = np.asarray(sens_all, dtype=float)
    n = s.size
    return {
        "power_at_pre_registered": float(s.mean()),
        "fdr_observed_on_positive_oracle": float(np.mean(fdr_all)),
        "positive_oracle": first,
        # ---- added by the F2/power-MC fix; identical-valued and single-element when n_reps = 1
        "power_n_reps": int(n),
        "power_mean": float(s.mean()),
        "power_sd_across_datasets": float(s.std(ddof=1)) if n > 1 else None,
        "power_se": float(s.std(ddof=1) / np.sqrt(n)) if n > 1 else None,
        "power_min": float(s.min()),
        "power_max": float(s.max()),
        "power_replicates": reps,
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
    ap.add_argument("--gen-arm", choices=GEN_ARMS, default="lognormal",
                    help="generative donor random effect; 'lognormal' = the repo simulator "
                         "verbatim (default, unchanged behaviour)")
    ap.add_argument("--sigma-het", type=float, default=0.0,
                    help="DIAGNOSTIC, NOT pre-registered. 0 (default) = the pre-registered "
                         "homoscedastic donor RE, bit-identical to before. h > 0 makes the "
                         "per-gene donor VARIANCE log-normal about sigma_donor^2 with "
                         "log-sd h, mean-preserving: sigma_g^2 = sigma_donor^2 * "
                         "exp(N(-h^2/2, h^2)). This is the condition under which variance "
                         "shrinkage (eBayes / edgeR-QLF) is predicted to become "
                         "anti-conservative; with h = 0 the simulator is homoscedastic and "
                         "moderation cannot be stressed at all.")
    ap.add_argument("--n-power-reps", type=int, default=1,
                    help="independent positive-oracle datasets averaged into "
                         "power_at_pre_registered. 1 (default) = the historical single draw, "
                         "bit-identical. The across-dataset SD of sensitivity was measured at "
                         "0.046-0.062 at 8v8, i.e. ~+-0.10 at 95%% on ONE draw, against a "
                         "BINDING 0.60 gate -- so any run whose power is compared with the gate "
                         "must use > 1 and read power_se.")
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
    pvals, ndeg, G, n_used, diag, n_donors_total, diags, nan_counts = runner(
        sim_kw, a.test, a.n_perm, a.seed, a.n_cpus, a.gen_arm, a.sigma_het
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
        "gen_arm": a.gen_arm,
        "sigma_het": a.sigma_het,
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
        "power_n_reps": None,
        "power_sd_across_datasets": None,
        "power_se": None,

        "power_target": POWER_TARGET,
        "lambda_band": list(LAMBDA_BAND),
        "test_diagnostics": diag,
        # ---- keys added by the F1/F3 fixes; every pre-existing key above is untouched ----
        # `test_diagnostics` is the LAST replicate only and is NOT stable across replicates.
        # Read `test_diagnostics_summary` instead for anything that affects a conclusion.
        "test_diagnostics_is_last_replicate_only": True,
        "test_diagnostics_summary": _summarise_diags(diags),
        # FIX F1 monitoring, computed for EVERY --test on each replicate's own universe.
        # mtc.bh drops NaN p-values, so nonzero here means the BH denominator was < G.
        "n_pval_nan_total": int(np.sum(nan_counts)),
        "n_pval_nan_max_per_replicate": int(np.max(nan_counts)) if nan_counts else 0,
        "cli": " ".join(sys.argv[1:]) if argv is None else " ".join(map(str, argv)),
        "versions": _versions(),
    }

    sec_power = 0.0
    if not a.skip_power:
        tp = time.perf_counter()
        out.update(_power(sim_kw, a.test, a.seed, a.n_cpus, a.gen_arm, a.sigma_het,
                          n_reps=a.n_power_reps))
        sec_power = time.perf_counter() - tp

    out["seconds_elapsed"] = float(time.perf_counter() - t0)
    out["seconds_null"] = float(sec_null)
    out["seconds_power"] = float(sec_power)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2), encoding="utf-8")

    pw = out["power_at_pre_registered"]
    if pw is None:
        pw_str = "n/a"
    elif out.get("power_se") is None:
        pw_str = f"{pw:.3f} (n=1 dataset, MC error NOT estimated)"
    else:
        pw_str = f"{pw:.3f} +- {out['power_se']:.3f} (n={out['power_n_reps']})"
    print(f"[{a.test}] gen={a.gen_arm} mode={a.null_mode} G={G} n_perm={n_used}  "
          f"lambda_pb={out['lambda_pb']:.3f} (band {LAMBDA_BAND})  "
          f"FP={fp_rate:.3f} ({n_hit}/{n_used}, binom p={p_binom:.2g})  "
          f"mean_rej={out['mean_rejections_perm_null']:.2f}  "
          f"power={pw_str}  "
          f"{out['seconds_elapsed']:.1f}s -> {a.out}")
    if min_p is not None:
        print(f"  !  exact-MWU discreteness floor: min attainable p = {min_p:.3e}; BH needs "
              f">= {k_needed} genes at that floor before it can reject anything"
              + ("" if attainable else "  -- IMPOSSIBLE, this arm cannot reject at all"))
    if out["n_pval_nan_total"]:
        print(f"  !  {out['n_pval_nan_total']} NaN p-values across the null "
              f"(max {out['n_pval_nan_max_per_replicate']} in one replicate): mtc.bh drops them, "
              f"so this arm's BH denominator was smaller than G = {G}.")
    if a.null_mode == "permutation" and a.test in (
        "wilcoxon", "ttest", "pooled_t", "ebayes", "ebayes_trend", "voom"
    ):
        print("  !! permutation null is (near-)exact for this test BY CONSTRUCTION; "
              "use --null-mode fresh for H1-vs-H2 evidence.")
    return out


if __name__ == "__main__":
    main()
