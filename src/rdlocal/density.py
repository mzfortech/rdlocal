"""Local polynomial density estimator and manipulation (density discontinuity) test.

Cattaneo, Jansson & Ma (2020, JASA), "Simple Local Polynomial Density Estimators":
f_hat(x) = slope of a kernel-weighted order-p regression of the empirical CDF on (x_i - x)
(Section 2); variance estimator of Theorem 2; test statistic of Section 4 / Corollary 1; robust
statistic T_{p+1}(h_p) (Section 4).  The "unrestricted" model is used: each side of the cutoff has its
own polynomial, with the whole-sample empirical CDF.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

from ._base import Result, clean
from ._kernels import kernel_name, kvals, min_gamma, mom
from .bandwidth import _global_coef, _pre_pilot


@dataclass(frozen=True)
class DensityResult(Result):
    f_left: float
    f_right: float
    se_left: float
    se_right: float
    f_left_bc: float
    f_right_bc: float
    se_left_bc: float
    se_right_bc: float
    statistic: float      # robust bias-corrected T_{p+1}(h), the recommended test
    p_value: float
    statistic_conv: float  # conventional T_p(h)
    p_value_conv: float
    h_left: float
    h_right: float
    n_left: int
    n_right: int
    n_eff_left: int
    n_eff_right: int
    p: int
    q: int
    kernel: str
    cutoff: float


def _side(xs, F, n, c, h, order, kernel, side):
    """Return (f_hat, variance, n_eff) for one side.  xs sorted, F = empirical CDF at xs."""
    xc = xs - c
    m = ((xc >= 0) if side > 0 else (xc < 0)) & (np.abs(xc) <= h)
    u = xc[m] / h
    w = kvals(kernel, u)
    pos = w > 0
    u, w, Fm = u[pos], w[pos], F[m][pos]
    if len(u) < order + 2 or len(np.unique(u)) < order + 1:
        raise ValueError(f"too few observations within bandwidth {h:.6g} on side {side:+d}")
    X = u[:, None] ** np.arange(order + 1)
    Ainv = np.linalg.inv(X.T @ (w[:, None] * X))
    cj = (Ainv[1] @ (X * w[:, None]).T) / h  # f_hat = sum_j c_j F(x_j)  (CJM 2020 sec. 2)
    f = float(cj @ Fm)
    # CJM 2020 Thm 2: Var = (1/n) sum_{j,k} c_j c_k [F(x_j ^ x_k) - F(x_j) F(x_k)]  (algebraically the
    # printed e1' A^-1 B A^-1 e1 / (n h) with the triple-sum B).  Suffix sums give the min() term in O(m).
    suff = np.cumsum((cj)[::-1])[::-1] - cj  # sum_{k>j} c_k
    T = np.sum(cj * Fm * (cj + 2 * suff))
    var = float(max(T - f**2, 0.0) / n)
    return f, var, len(u)


def _pilot_density(xs, F, n, c, kernel, p):
    v = _pre_pilot(xs - c, None, p, kernel)
    out = {}
    for s in (1, -1):
        h = v
        for _ in range(60):
            try:
                out[s] = max(_side(xs, F, n, c, h, p, kernel, s)[0], 1e-12)
                break
            except ValueError:
                h *= 1.15
    return out


def _h_mse(xs, F, n, c, p, kernel):
    """MSE-optimal common bandwidth for f_hat(+) - f_hat(-).

    CJM 2020 Sec. 3.1: h = (V / (2 p B^2))^(1/(1+2p)) n^(-1/(1+2p)), B = F^{(p+1)}/(p+1)! e1' S^-1 c,
    V = f e1' S^-1 Gamma S^-1 e1.  DEVIATION: the paper's data-driven plug-in for the *difference* is in
    the omitted SA text; here F^{(p+1)}/(p+1)! comes from a global order-(p+3) polynomial fit of the
    empirical CDF on each side (Fan-Gijbels style, as CCF 2020 do for derivatives), f from a Silverman
    pilot window, and the one-sided biases/variances are combined (sign (-1)^p for the left bias).
    """
    xc = xs - c
    S = np.array([[mom(kernel, j + k) for k in range(p + 1)] for j in range(p + 1)])
    Si = np.linalg.inv(S)
    cvec = np.array([mom(kernel, p + 1 + j) for j in range(p + 1)])
    Gm = min_gamma(kernel, p)
    Bc = float(Si[1] @ cvec)
    Vc = float(Si[1] @ Gm @ Si[1])
    fp = _pilot_density(xs, F, n, c, kernel, p)
    mR = _global_coef(xc, F, p + 3, p + 1, 1)
    mL = _global_coef(xc, F, p + 3, p + 1, -1)
    B = Bc * (mR - (-1) ** p * mL)
    V = Vc * (fp[1] + fp[-1])
    cap = min(np.abs(xc[xc >= 0]).max(), np.abs(xc[xc < 0]).max())
    if B == 0 or not np.isfinite(B):
        return cap
    return float(min((V / (2 * p * B**2)) ** (1 / (1 + 2 * p)) * n ** (-1 / (1 + 2 * p)), cap))


def rd_density_test(x, cutoff=0.0, *, p=2, h=None, kernel="triangular") -> DensityResult:
    """Manipulation test: H0 is continuity of the running-variable density at the cutoff.

    h: None (data-driven common bandwidth), a float, or (h_left, h_right).
    """
    kernel = kernel_name(kernel)
    if p < 1:
        raise ValueError("p must be >= 1")
    (x,), _ = clean(x)
    n = len(x)
    xs = np.sort(x)
    # F(x_i) = (1/(n-1)) sum_{j != i} 1(x_j <= x_i): the leave-in-free form of CJM 2020 SA sec. 5.2
    # (the self term only adds O(1/n) "leave-in bias").  DEVIATION from the sec. 2 definition, which
    # keeps j = i; the difference is O(1/n) in f_hat.
    F = (np.searchsorted(xs, xs, side="right") - 1) / (n - 1)
    if h is None:
        hl = hr = _h_mse(xs, F, n, cutoff, p, kernel)
    else:
        hl, hr = (h, h) if np.isscalar(h) else tuple(h)
    q = p + 1
    res = {}
    for order in (p, q):
        fl, vl, nl = _side(xs, F, n, cutoff, hl, order, kernel, -1)
        fr, vr, nr = _side(xs, F, n, cutoff, hr, order, kernel, 1)
        res[order] = (fl, vl, nl, fr, vr, nr)
    t = {o: (r[3] - r[0]) / np.sqrt(r[1] + r[4]) for o, r in res.items()}
    pv = {o: float(2 * stats.norm.sf(abs(t[o]))) for o in t}
    a, b = res[p], res[q]
    return DensityResult(
        f_left=a[0], f_right=a[3], se_left=float(np.sqrt(a[1])), se_right=float(np.sqrt(a[4])),
        f_left_bc=b[0], f_right_bc=b[3], se_left_bc=float(np.sqrt(b[1])), se_right_bc=float(np.sqrt(b[4])),
        statistic=float(t[q]), p_value=pv[q], statistic_conv=float(t[p]), p_value_conv=pv[p],
        h_left=float(hl), h_right=float(hr), n_left=int((x < cutoff).sum()), n_right=int((x >= cutoff).sum()),
        n_eff_left=a[2], n_eff_right=a[5], p=p, q=q, kernel=kernel, cutoff=float(cutoff),
    )
