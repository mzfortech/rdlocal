"""Sharp and fuzzy RD point estimation with robust bias-corrected inference (CCT 2014)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

from ._base import Result, clean
from ._kernels import kernel_name
from ._locpoly import Fit
from .bandwidth import rd_bandwidth, sigma2_nn


@dataclass(frozen=True)
class RDResult(Result):
    coef: float          # conventional estimate, order p, bandwidth h
    coef_bc: float       # bias-corrected estimate
    se_conv: float
    se_robust: float
    z_conv: float
    z_robust: float
    pv_conv: float
    pv_robust: float
    ci_conv: tuple
    ci_robust: tuple     # centred at coef_bc, robust SE
    n_left: int          # effective observations (inside h)
    n_right: int
    n_total_left: int
    n_total_right: int
    h_left: float
    h_right: float
    b_left: float
    b_right: float
    kernel: str
    p: int
    q: int
    level: float
    cutoff: float
    vce: str
    bwselect: str
    fuzzy: bool
    first_stage: float | None  # fuzzy only: jump in treatment take-up


def _hc_adj(vce, lev, n_s, k):
    if vce == "hc0":
        return np.ones_like(lev)
    if vce == "hc1":
        return np.full_like(lev, n_s / max(n_s - k, 1))
    if vce == "hc2":
        return 1 / (1 - lev)
    if vce == "hc3":
        return 1 / (1 - lev) ** 2
    raise ValueError("vce must be hc0, hc1, hc2, hc3 or nn")


def _variance(z, W, fh, hs, p, vce, cl, xc, fb=None, q=None):
    """Variance of the linear statistic sum_i W_i z_i, residuals of z from the order-p fit at h.

    CCT (2014) supplement S.2.4: plug-in residuals eps_i = Y_i - mu_hat_p(X_i; h) for the conventional
    variance.  For the bias-corrected statistic the residuals "are of the same form, with only the
    bandwidth and polynomial order changed" (CCF 2020 supplement S.4): observations in the b-window use the
    order-q fit at b (`fb`), the rest the order-p fit at h.  DEVIATION: the HC0-HC3 / cluster refinements are
    the textbook ones (MacKinnon 2013; CCF 2020 S.5.2 only cites them), with the weighted hat matrix of the
    fit that produced the residual.
    """
    n = len(z)
    if vce == "nn" and cl is None:
        return float(np.sum(W**2 * sigma2_nn(xc, z)))
    r = np.zeros(n)
    adj = np.ones(n)
    for s, f in fh.items():
        beta = f.Ainv @ (f.X.T @ (f.w * z[f.idx]))
        side = (xc >= 0) if s > 0 else (xc < 0)
        pw = (xc[side] / f.h)[:, None] ** np.arange(p + 1)
        r[side] = z[side] - pw @ beta
        if cl is None:
            _, lev = f.residuals(z)
            adj[f.idx] = _hc_adj(vce, lev, f.n, p + 1)
    if fb is not None:
        for f in fb.values():
            r_, lev = f.residuals(z)
            r[f.idx] = r_
            if cl is None:
                adj[f.idx] = _hc_adj(vce, lev, f.n, q + 1)
    if cl is None:
        return float(np.sum(W**2 * adj * r**2))
    # cluster-robust: sum over clusters of squared cluster totals of the score W_i r_i
    G = len(np.unique(cl))
    tot = np.bincount(np.unique(cl, return_inverse=True)[1], weights=W * r)
    nn_ = sum(f.n for f in fh.values())
    corr = G / max(G - 1, 1) * (nn_ - 1) / max(nn_ - 2 * (p + 1), 1)
    return float(np.sum(tot**2) * corr)


def rd_estimate(y, x, cutoff=0.0, *, fuzzy=None, p=1, q=None, kernel="triangular", h=None, b=None,
                bwselect="mserd", vce="hc1", cluster=None, level=0.95, scaleregul=1.0) -> RDResult:
    """Local polynomial RD estimate with conventional and robust bias-corrected inference.

    Treatment is x >= cutoff.  Conventional: order-p fit at bandwidth h.  Bias-corrected: the
    order-p estimate minus the estimated leading bias from an order-q fit at bandwidth b
    (CCT 2014; finite-sample form as in CCF 2020 S.5.3).  Robust SE: variance of the
    bias-corrected linear statistic.  With `cluster`, `vce` is ignored (cluster-robust).
    """
    kernel = kernel_name(kernel)
    q = p + 1 if q is None else q
    if q <= p:
        raise ValueError("q must exceed p")
    (y, x, d, cl), _ = clean(y, x, fuzzy, cluster)
    xc = x - cutoff
    bwsel = bwselect
    if h is None:
        bw = rd_bandwidth(y, x, cutoff, fuzzy=d, p=p, q=q, kernel=kernel, bwselect=bwselect, cluster=cl, scaleregul=scaleregul)
        hl, hr = bw.h_left, bw.h_right
        bl, br = (bw.b_left, bw.b_right) if b is None else (b, b)
    else:
        hl, hr = (h, h) if np.isscalar(h) else tuple(h)
        bwsel = "manual"
        bl, br = (hl, hr) if b is None else ((b, b) if np.isscalar(b) else tuple(b))
    hs = {1: hr, -1: hl}
    bs = {1: br, -1: bl}

    n = len(xc)
    Wc = np.zeros(n)
    Wb = np.zeros(n)
    fh, fbd = {}, {}
    for s in (1, -1):
        f = Fit(xc, hs[s], p, kernel, s)
        fb = Fit(xc, bs[s], q, kernel, s)
        fh[s] = f
        fbd[s] = fb
        wo = f.coef_weights(0)
        kappa = float(wo @ f.u ** (p + 1))  # bias of the order-p intercept per unit of m_{p+1} h^(p+1)
        Wc[f.idx] += s * wo
        Wb[f.idx] += s * wo
        # CCT 2014 / CCF 2020 S.5.3: tau_bc = tau - h^(p+1) * kappa * m_hat_{p+1}(b), m = coef of x^(p+1)
        Wb[fb.idx] -= s * hs[s] ** (p + 1) * kappa * fb.coef_weights(p + 1) / bs[s] ** (p + 1)

    if d is None:
        tau, tau_bc = Wc @ y, Wb @ y
        se_c = np.sqrt(_variance(y, Wc, fh, hs, p, vce, cl, xc))
        se_r = np.sqrt(_variance(y, Wb, fh, hs, p, vce, cl, xc, fbd, q))
        fs = None
    else:
        ty, td = Wc @ y, Wc @ d
        tau = ty / td
        # CCT supplement (proof of Thm A.2): zeta_bc = zeta - h^(p+1) B_F with the *linearised* bias
        # B_F = bias_Y/tau_T - tau_Y bias_T/tau_T^2, i.e. (bias of y - tau d)/tau_T.
        tau_bc = tau - ((Wc - Wb) @ (y - tau * d)) / td
        # delta method == variance of the linear statistic in the pseudo-outcome y - tau*d, scaled by 1/td
        # (CCT supplement Lemma 2: V_F with sigma^2_YY - 2 tau sigma_YT + tau^2 sigma_TT).  The bias-corrected
        # linear statistic has weights Wb/td; the pseudo-outcome keeps the conventional ratio tau (the
        # linearisation point of the supplement's expansion around zeta_hat_p(h)).
        se_c = np.sqrt(_variance(y - tau * d, Wc / td, fh, hs, p, vce, cl, xc))
        se_r = np.sqrt(_variance(y - tau * d, Wb / td, fh, hs, p, vce, cl, xc, fbd, q))
        fs = float(td)
    zc = stats.norm.ppf(1 - (1 - level) / 2)
    z1, z2 = tau / se_c, tau_bc / se_r
    return RDResult(
        coef=float(tau), coef_bc=float(tau_bc), se_conv=float(se_c), se_robust=float(se_r),
        z_conv=float(z1), z_robust=float(z2),
        pv_conv=float(2 * stats.norm.sf(abs(z1))), pv_robust=float(2 * stats.norm.sf(abs(z2))),
        ci_conv=(float(tau - zc * se_c), float(tau + zc * se_c)),
        ci_robust=(float(tau_bc - zc * se_r), float(tau_bc + zc * se_r)),
        n_left=fh[-1].n, n_right=fh[1].n,
        n_total_left=int((xc < 0).sum()), n_total_right=int((xc >= 0).sum()),
        h_left=float(hl), h_right=float(hr), b_left=float(bl), b_right=float(br),
        kernel=kernel, p=p, q=q, level=level, cutoff=float(cutoff), vce="cluster" if cl is not None else vce,
        bwselect=bwsel, fuzzy=d is not None, first_stage=fs,
    )
