"""The pseudobulk arm's test: moderated t with empirical-Bayes shrunk variances (Smyth 2004).

Selected by **Amendment 2 (2026-08-15), Change 1**, which closes the replacement-test question
Amendment 1 deliberately left open. The deciding evidence is the 146-cell test-selection grid
(``pilot/testsel/summary.{csv,json}``, committed at ``72dec7b``): the shrinkage hypothesis is not
supported — ``ebayes`` stays calibrated at the strongest gene-variance heterogeneity in the grid
(FP = 0.046 against alpha = 0.05, lambda = 0.998, realised ``d0`` = 2.13 against d = 14) — and only
the moderated tests reach the pre-registered power target (spec §8(c): sensitivity >= 0.60 at
log2FC = 1.0, K = 200) anywhere in the grid. It replaces DESeq2-Wald, which Amendment 1 showed
fails the binding validity gate of decision rule item 1 / §8(a).

Spec sections this implements: §3 (pseudobulk DE model — the DESeq2 block is superseded; C1
poscounts and C4 covariate-df lapse with DESeq2, C2's no-asymmetric-NA requirement survives and is
carried by :func:`pbcheck.mtc.bh_both_arms`), §5 (raw p-values only; BH is applied outside), §9
item 5. The realised prior reported in :func:`ebayes_from_pdata`'s ``moderation`` metadata is the
moderated analog of DESeq2's ``fitType``, per Amendment 2 Change 4 (spec C5).

PROVENANCE. The statistics below are **not a rewrite**. They are lifted verbatim from
``scripts/pb_calibration_probe.py`` — the code that produced the grid's ``prior_df_d0`` and
``shrinkage`` columns — so the arm this package ships is numerically the arm the selection evidence
was measured on. ``tests/test_moderated.py`` pins that equivalence against the probe's own
functions. Changes made in the lift are confined to naming, docstrings, and the public entry point;
no formula, constant, branch or guard is altered.

WHY PURE PYTHON. There is no R, rpy2, limma or edgeR available here and pbcheck targets the
scanpy/scverse ecosystem, so wrapping limma is not an option. These are from-scratch
implementations of the published estimators (Smyth, G. K. 2004, *Linear models and empirical Bayes
methods for assessing differential expression in microarray experiments*, Stat. Appl. Genet. Mol.
Biol. 3, Article 3), validated numerically by the probe's four self-checks: d0 -> 0 reduces to the
ordinary t; homoscedastic-null calibration by KS against ``t_{d+d0}``; a power gain over the
ordinary t in the homoscedastic case; recovery of a known d0.
"""

from __future__ import annotations

from math import log

import numpy as np
import pandas as pd
from scipy.special import digamma, polygamma
from scipy.stats import norm, t as t_dist

from pbcheck.methods.de import DEResult

#: limma's ``fitFDist(covariate=)`` default: ``splinedf = min(4, ...)``. Only used when
#: ``trend=True``, which Amendment 2 Change 1 measured and did **not** select.
EBAYES_TREND_DF = 4


def log_cpm(pdata) -> np.ndarray:
    """log2(CPM + 1) of the aggregated donor matrix. Shape (n_donors, n_genes).

    Library size is the row sum of the **universe-restricted** matrix, so the transform is a
    function of exactly the data the test sees and no gene outside the frozen universe leaks into
    the normalisation. This is the moderated arm's counterpart to a DESeq2 size factor (Amendment 2
    Change 1); log2 rather than ln so the reported ``log2fc`` column is directly interpretable.
    """
    X = np.asarray(pdata.X, dtype=float)
    lib = X.sum(axis=1, keepdims=True)
    cpm = X / np.maximum(lib, 1.0) * 1e6
    return np.log2(cpm + 1.0)


def trigamma_inverse(x):
    """Solve trigamma(y) = x for y > 0. Newton iteration on limma's parameterisation.

    ``limma::trigammaInverse`` iterates on y with step
    ``dy = trigamma(y) * (1 - trigamma(y)/x) / psigamma(y, deriv=2)``, which is Newton's method
    applied to the (nearly linear, convex) function ``1/trigamma(y)``; ``psigamma(y, 2)`` is
    ``polygamma(2, y)``, the tetragamma function.

    Closed-form branches at the extremes (limma's, and both the correct asymptotics):
    ``trigamma(y) ~ 1/y^2`` as ``y -> 0`` (x very large => ``y = 1/sqrt(x)``), and
    ``trigamma(y) ~ 1/y`` as ``y -> inf`` (x very small => ``y = 1/x``).
    """
    x = np.asarray(x, dtype=float)
    scalar = x.ndim == 0
    x = np.atleast_1d(x).astype(float)
    y = np.full_like(x, np.nan)
    big = x > 1e7          # trigamma(y) ~ 1/y^2 for tiny y
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


def natural_spline_basis(x, df: int = EBAYES_TREND_DF):
    """Natural cubic regression spline basis, ``df`` columns INCLUDING the intercept.

    Truncated-power construction (Hastie/Tibshirani/Friedman eq. 5.4-5.5): with knots
    ``xi_1 < ... < xi_K``, ``N_1 = 1``, ``N_2 = x``, ``N_{k+2} = d_k(x) - d_{K-1}(x)`` for
    ``k = 1..K-2``, where ``d_k(x) = [(x - xi_k)_+^3 - (x - xi_K)_+^3] / (xi_K - xi_k)``. This spans
    the natural cubic splines (linear beyond the boundary knots) in K columns. Knots are placed at
    K = ``df`` equally spaced quantiles of x.

    DEPARTURE FROM R, stated explicitly: limma calls ``splines::ns(covariate, df=4,
    intercept=TRUE)``, an orthonormalised B-spline basis. The SPAN of the two bases is the same
    space of natural cubic splines when the knot sets agree, so the OLS *fitted values* — the only
    thing ``fit_f_dist`` uses — are identical up to conditioning; R's ``ns()`` places interior knots
    at quantiles and boundary knots at the data range, which is not bit-identical to the quantile
    grid used here. A smoother-choice difference, not a different estimator.

    Only reachable via ``trend=True``, which Amendment 2 Change 1 measured and did not select.
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


def fit_f_dist(s2, df, covariate=None, spline_df: int = EBAYES_TREND_DF):
    """Smyth (2004) ``fitFDist``: method of moments on log(s2) via digamma/trigamma.

    MODEL. ``s2_g | sigma2_g ~ sigma2_g * chisq_d / d`` and ``d0 * s02 / sigma2_g ~ chisq_d0``, so
    marginally ``s2_g ~ s02 * F(d, d0)``.

    ESTIMATOR. Write ``z_g = log s2_g``. Because ``log(chisq_d / d)`` has mean
    ``digamma(d/2) - log(d/2)`` and variance ``trigamma(d/2)``, and ``log(sigma2_g)`` has mean
    ``log(s02) + log(d0/2) - digamma(d0/2)`` and variance ``trigamma(d0/2)``::

        E[z_g]   = log(s02) + digamma(d/2) - log(d/2) + log(d0/2) - digamma(d0/2)
        Var[z_g] = trigamma(d/2) + trigamma(d0/2)

    Define ``e_g = z_g - digamma(d/2) + log(d/2)``. Then, with n usable genes::

        evar = sum_g (e_g - ebar_g)^2 / (n - rank)        [rank = 1, or the spline rank]
        d0   = 2 * trigammaInverse(evar - trigamma(d/2))
        s02  = exp(ebar + digamma(d0/2) - log(d0/2))

    If ``evar - trigamma(d/2) <= 0`` the between-gene variance is not distinguishable from the
    within-gene sampling variance, and the correct limit is ``d0 = +inf`` (complete pooling),
    ``s02 = exp(ebar)``. This is limma's behaviour and is reported, never silently clipped.

    ``covariate`` (limma-trend): ``ebar`` becomes a per-gene FITTED value from an OLS regression of
    e on a natural cubic spline in the covariate, so ``s02`` becomes a per-gene TREND. ``d0`` stays
    scalar, exactly as in limma.

    Returns ``(d0, s02, info)``; ``s02`` is a scalar (no covariate) or a length-G vector.
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
        B_all = natural_spline_basis(cov, df=spline_df)
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
        d0 = 2.0 * trigamma_inverse(evar_between)
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


def wls_two_group(Y, x, weights=None, return_fitted: bool = False):
    """Per-gene (weighted) least squares for the design ``~ 1 + x``, x a 0/1 indicator.

    ``Y`` is (n_donors, n_genes); ``weights`` is None (ordinary LS) or (n_donors, n_genes).
    Returns ``(beta, s2, d, unscaled_var)``:

    - ``beta`` contrast estimate (group x=1 minus group x=0), length G
    - ``s2`` residual variance ``sum_i w_i e_i^2 / d``, length G
    - ``d`` residual degrees of freedom = n - 2 (scalar)
    - ``unscaled_var`` ``(X' W X)^{-1}[1,1]``, length G, so ``se(beta) = sqrt(s2 * unscaled_var)``

    Sanity anchor: with weights = 1 this reduces to ``beta = mean_A - mean_B``,
    ``unscaled_var = 1/n_A + 1/n_B``, and ``s2`` the pooled variance — i.e. Student's
    equal-variance two-sample t.
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
        beta0 = (S1 * Sy - S1 * Sxy) / det          # intercept (weighted mean of the x=0 group)
        beta1 = (S0 * Sxy - S1 * Sy) / det          # contrast
        unscaled_var = S0 / det                     # (X'WX)^{-1}[1,1]
        fitted = beta0[None, :] + beta1[None, :] * x
        resid = Y - fitted
        d = n - 2
        s2 = (W * resid ** 2).sum(axis=0) / d
    if return_fitted:
        return beta1, s2, d, unscaled_var, fitted
    return beta1, s2, d, unscaled_var


def moderated_t(beta, s2, d, unscaled_var, d0, s02):
    """Smyth's moderated t and its two-sided p-value::

        s~_g^2 = (d0 * s02_g + d * s2_g) / (d0 + d)
        t_g    = beta_g / sqrt(s~_g^2 * unscaled_var_g)
        p_g    = 2 * P(T_{d + d0} > |t_g|)

    ``d0 = +inf`` is handled in the limit: ``s~^2 -> s02`` and ``T_inf -> N(0, 1)``. Passing
    ``d0 = 0`` reproduces the ordinary (pooled, equal-variance) t exactly.

    A gene with a non-finite ``s2`` is UNTESTABLE and comes out as ``p = NaN`` in every branch, so
    that :func:`pbcheck.mtc.bh` drops it instead of publishing a number derived entirely from the
    prior. See the comment in the ``d0 = +inf`` branch: that was a live silent-p-value bug.
    """
    s2 = np.asarray(s2, dtype=float)
    s02_arr = np.broadcast_to(np.asarray(s02, dtype=float), s2.shape)
    finite_s2 = np.isfinite(s2)
    if np.isinf(d0):
        # The posterior variance in the complete-pooling limit is the prior, which DISCARDS the
        # gene's own s2 entirely. Before this guard a gene whose s2 was NaN silently received a
        # finite, valid-looking p-value computed purely from the prior. NaN must propagate here
        # exactly as it does in the finite-d0 branch, so an untestable gene is dropped from BH
        # rather than given a fabricated p-value.
        s2_post = np.where(finite_s2, s02_arr.astype(float), np.nan)
        df_total = np.inf
    elif d0 == 0:
        # Unmoderated limit, written out exactly rather than as (0*s02 + d*s2)/(0+d), so that it is
        # bit-identical to the equal-variance two-sample Student t rather than merely equal to
        # rounding. s02 is unused here by construction.
        s2_post = s2.astype(float, copy=True)
        df_total = float(d)
    else:
        s2_post = (d0 * s02_arr + d * s2) / (d0 + d)
        df_total = d + d0
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(s2_post * np.asarray(unscaled_var, dtype=float))
        tstat = np.asarray(beta, dtype=float) / se
        if np.isinf(df_total):
            pval = 2.0 * norm.sf(np.abs(tstat))
        else:
            pval = 2.0 * t_dist.sf(np.abs(tstat), df_total)
    return tstat, np.asarray(pval, dtype=float), s2_post, df_total


def ebayes_from_pdata(
    pdata,
    *,
    condition_col: str = "condition",
    test_level: str = "disease",
    ref_level: str = "ctrl",
    universe: list[str] | None = None,
    condition_values=None,
    trend: bool = False,
) -> DEResult:
    """The pseudobulk arm: moderated eBayes across donors (Amendment 2 Change 1).

    Same interface and role as the superseded
    :func:`pbcheck.methods.pseudobulk.deseq_from_pdata`, so a permutation harness can aggregate
    donors **once** and only relabel donor->condition per permutation (donor profiles are invariant
    under relabeling).

    Parameters
    ----------
    pdata
        Donor-pseudobulk ``AnnData`` (one row per donor), raw summed counts in ``.X``.
    universe
        Frozen label-agnostic gene universe (spec §3/A3); the matrix is restricted to it before the
        log-CPM transform, so library sizes are computed over exactly the tested genes.
    condition_values
        Optional per-donor condition vector overriding ``pdata.obs[condition_col]``, aligned with
        ``pdata.obs_names`` (used by the permutation null).
    trend
        Fit the prior variance as a spline trend in mean log-CPM (limma-trend). **Default False**:
        Amendment 2 Change 1 measured the trended variant and did not select it (mean power
        difference −0.0002 over 24 matched grid cells). Exposed so the choice stays inspectable and
        the ``s0_squared_is_trended`` metadata means something, not so it is toggled casually —
        changing the arm's default is a protocol change requiring an amendment.

    Returns
    -------
    :class:`~pbcheck.methods.de.DEResult` with **raw** p-values (``padj`` is NaN: BH is applied
    outside, over the frozen universe, identically for both arms — spec §5). The realised prior and
    the rest of the moderation metadata (the ``fitType`` analog, Amendment 2 Change 4 / spec C5) is
    attached as the ``.moderation`` attribute.
    """
    if universe is not None:
        # The set is built ONCE, outside the comprehension. It used to sit inside the condition, so
        # the pandas Index was rebuilt into a set per candidate gene and the restriction was O(G^2):
        # profiled at G = 15 000 that single line took 113.7 s of a 115.5 s call, 98 % of the cost,
        # against 0.02 s hoisted. Invisible at the 1500-gene ORACLE_SIM the gate has always run at,
        # and paid on every stratum, every permutation and every jackknife replicate of a real sweep
        # against a realistic frozen universe. The result is unchanged and the ordering is still the
        # universe's, not the matrix's. Amendment 4 Part A, Correction 1's reproducer is what
        # surfaced it. NOTE: the same construction survives in ``methods/pseudobulk.py`` (the
        # superseded DESeq2 arm) and in ``scripts/pb_calibration_probe.py`` (the frozen reference
        # instrument, deliberately not edited); neither is on this arm's path.
        present = set(pdata.var_names)
        keep = [g for g in universe if g in present]
        pdata = pdata[:, keep].copy()
    if pdata.n_vars == 0:
        raise ValueError("No genes to test after universe restriction.")

    cond = (pd.Series(condition_values, index=pdata.obs_names) if condition_values is not None
            else pdata.obs[condition_col])
    cond = np.asarray(cond.astype(str).to_numpy(), dtype=object)

    x = (cond == test_level).astype(float)
    nA, nB = int(x.sum()), int((cond == ref_level).sum())
    if nA < 2 or nB < 2:
        raise ValueError(f"need >= 2 donors per group, got {nA} vs {nB}")
    if nA + nB != len(cond):
        raise ValueError("condition vector contains levels other than test/ref")

    Y = log_cpm(pdata)
    beta, s2, d, uvar = wls_two_group(Y, x, weights=None)
    covariate = Y.mean(axis=0) if trend else None
    d0, s02, finfo = fit_f_dist(s2, d, covariate=covariate, spline_df=EBAYES_TREND_DF)
    tstat, pval, s2_post, df_total = moderated_t(beta, s2, d, uvar, d0, s02)

    # d0 = +inf is a legitimate estimate (complete pooling: the between-gene variance is not
    # separable from sampling noise). JSON has no infinity literal, so it is emitted as null with an
    # explicit boolean beside it. A reader must check `prior_is_complete_pooling`: when true, EVERY
    # gene was given the same variance, i.e. the maximum-shrinkage regime.
    d0_finite = bool(np.isfinite(d0))
    moderation = dict(finfo)
    moderation["d0"] = float(d0) if d0_finite else None
    moderation.update({
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
        # 0.0 = no shrinkage at all.
        "shrinkage_factor_d0_over_d0_plus_d": (1.0 if not d0_finite else float(d0 / (d0 + d))),
        "moderation_applied": True,
        # mtc drops NaN p-values, so a nonzero count here means this arm's BH denominator is
        # smaller than the universe — which is exactly what mtc.bh_both_arms equalises (Change 2).
        "n_genes_nonfinite_s2": int(np.sum(~np.isfinite(np.asarray(s2, dtype=float)))),
        "n_pval_nan": int(np.sum(~np.isfinite(np.asarray(pval, dtype=float)))),
        "n_donors_test": nA,
        "n_donors_ref": nB,
    })

    # beta is the contrast on the log2 scale the arm used, so it IS a log2 fold change.
    table = pd.DataFrame(
        {"pval": pval, "padj": np.nan, "log2fc": beta},
        index=pd.Index(list(pdata.var_names), name="gene"),
    )
    res = DEResult(
        method=("pseudobulk[eBayes:log2CPM,trend]" if trend else "pseudobulk[eBayes:log2CPM]"),
        table=table,
        contrast=(condition_col, test_level, ref_level),
    )
    res.moderation = moderation
    res.tstat = tstat
    return res
