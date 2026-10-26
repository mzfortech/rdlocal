"""Bandwidth selectors: CCT (2014) plug-in MSE (mserd / msetwo), CCF (2020) rule-of-thumb CER, IK (2012)."""
from __future__ import annotations

from dataclasses import dataclass
import re
from math import factorial

import numpy as np

from ._base import Result, clean
from ._kernels import kernel_name, lead_const, silverman_const
from ._locpoly import Fit


@dataclass(frozen=True)
class BandwidthResult(Result):
    h_left: float
    h_right: float
    b_left: float
    b_right: float
    bwselect: str
    kernel: str
    p: int
    q: int
    n_left: int
    n_right: int

    @property
    def h(self):
        return self.h_left if self.h_left == self.h_right else (self.h_left, self.h_right)


# ---------------------------------------------------------------- helpers

def _nn_sigma2(x, z, J=3):
    """Fixed-J nearest-neighbour variance, per observation of ONE side.

    CCT (2014) supplement S.2.4: sigma2_hat(X_i) = J/(J+1) (Y_i - mean_j Y_{j(i)})^2.
    """
    from scipy.spatial import cKDTree

    m = len(x)
    if m <= J:
        raise ValueError("too few observations for nearest-neighbour variance")
    d, ind = cKDTree(x[:, None]).query(x[:, None], k=J + 1)
    selfm = ind == np.arange(m)[:, None]
    selfm[~selfm.any(axis=1), -1] = True  # tie case: self not among neighbours -> drop farthest
    nb = ind[~selfm].reshape(m, J)
    return J / (J + 1) * (z - z[nb].mean(axis=1)) ** 2


def sigma2_nn(xc, z, J=3):
    """NN variance for every observation (each side matched within its own side)."""
    s2 = np.zeros(len(xc))
    for side in (xc >= 0, xc < 0):
        ix = np.flatnonzero(side)
        s2[ix] = _nn_sigma2(xc[ix], z[ix], J)
    return s2


def fit_enough(xc, h, order, kernel, side, nmin=None):
    """Fit, widening h by 10% steps until the window supports the fit (and nmin points)."""
    nmin = nmin or 2 * (order + 1)
    xs = np.abs(xc[(xc >= 0) if side > 0 else (xc < 0)])
    if len(xs) < nmin:
        raise ValueError(f"too few observations on side {side:+d}")
    h = max(h, np.sort(xs)[nmin - 1] * 1.0001)
    for _ in range(80):
        try:
            f = Fit(xc, h, order, kernel, side)
            if f.n >= nmin:
                return f
        except ValueError:
            pass
        h *= 1.1
    raise ValueError("could not form a local fit")


def deriv_var(fit: Fit, sig2, nu):
    """V-hat_{nu,p}(h): nu!^2 e_nu' A^-1 (X' W^2 diag(sig2) X) A^-1 e_nu / h^(2 nu), one side.

    CCT (2014) supplement S.2.6: V_hat_{+,nu,p}(h) = nu!^2 e' Gamma^-1 Psi_hat Gamma^-1 e / (n h^(2 nu)),
    which algebraically equals this sandwich in scaled units.
    """
    e = fit.Ainv[nu]
    M = fit.X.T @ ((fit.w**2 * sig2[fit.idx])[:, None] * fit.X)
    return float(factorial(nu) ** 2 * (e @ M @ e) / fit.h ** (2 * nu))


def coef(fit: Fit, y, k):
    """k-th Taylor coefficient (unscaled x-units) of the local fit."""
    return float(fit.coef_weights(k) @ y[fit.idx] / fit.h**k)


def _global_coef(xc, y, deg, k, side):
    m = (xc >= 0) if side > 0 else (xc < 0)
    s = np.abs(xc[m]).max()
    X = (xc[m] / s)[:, None] ** np.arange(deg + 1)
    beta = np.linalg.lstsq(X, y[m], rcond=None)[0]
    return float(beta[k] / s**k)


def _pre_pilot(xc, y, p, kernel):
    n = len(xc)  # y unused (kept for call symmetry)
    sx = min(np.std(xc, ddof=1), (np.percentile(xc, 75) - np.percentile(xc, 25)) / 1.349)
    return silverman_const(kernel) * sx * n ** (-0.2)


# ---------------------------------------------------------------- CCT plug-in

def _cct(xc, y, p, q, kernel):
    """CCT (2014) supplement S.2.6 plug-in, for the three estimands named in CCFT (2017, sec. 2.3):
    rd (difference of one-sided estimators), two (each side alone), sum (sum of one-sided estimators).

    Returns {variant: (h_left, h_right, b_left, b_right)}.  The variants differ only in which one-sided
    bias terms are combined (difference, one side, sum) and in whose variances are added.
    """
    n = len(xc)
    v = _pre_pilot(xc, y, p, kernel)
    s2 = sigma2_nn(xc, y)
    S = (1, -1)

    # Step 0 (ii): pilot c_n via global (q+2)-th order polynomial fits.  Estimand-independent.
    fq1 = {s: fit_enough(xc, v, q + 1, kernel, s) for s in S}
    Vq1 = sum(deriv_var(fq1[s], s2, q + 1) for s in S)
    Bq1 = lead_const(kernel, q + 1, q + 1)
    g = _global_coef(xc, y, q + 2, q + 2, 1) - (-1) ** q * _global_coef(xc, y, q + 2, q + 2, -1)
    Cc = (2 * q + 3) * n * v ** (2 * q + 3) * Vq1 / (2 * Bq1**2 * g**2)
    c_n = Cc ** (1 / (2 * q + 5)) * n ** (-1 / (2 * q + 5)) if np.isfinite(Cc) and Cc > 0 else v
    fc = {s: fit_enough(xc, c_n, q + 1, kernel, s) for s in S}
    cq = {s: coef(fc[s], y, q + 1) for s in S}
    Vcc = {s: deriv_var(fc[s], s2, q + 1) / factorial(q + 1) ** 2 for s in S}  # DEVIATION 1 (coefficient scale)
    fv = {s: fit_enough(xc, v, q, kernel, s) for s in S}
    Vpq = {s: deriv_var(fv[s], s2, p + 1) for s in S}
    Bpq = lead_const(kernel, q, p + 1)
    fp = {s: fit_enough(xc, v, p, kernel, s) for s in S}
    V0 = {s: deriv_var(fp[s], s2, 0) for s in S}
    B0 = lead_const(kernel, p, 0)

    # sign with which the LEFT one-sided bias enters: difference -> -(-1)^k, sum -> +(-1)^k, k = order+1
    def comb(a, bm, k, kind):
        if kind == "rd":
            return a[1] - (-1) ** k * a[-1]
        if kind == "sum":
            return a[1] + (-1) ** k * a[-1]
        raise ValueError

    out = {}
    for kind in ("rd", "sum", "two"):
        sides = [(1,), (-1,)] if kind == "two" else [S]
        hh, bb = {}, {}
        for grp in sides:
            def tot(d):
                return sum(d[s] for s in grp)
            # Step 1: b
            if kind == "two":
                d1 = cq[grp[0]]
            else:
                d1 = comb(cq, q + 1, q + 1, kind)
            Cb = (2 * p + 3) * n * v ** (2 * p + 3) * tot(Vpq) / (
                2 * (q - p) * Bpq**2 * (d1**2 + 3 * tot(Vcc)))
            b = Cb ** (1 / (2 * q + 3)) * n ** (-1 / (2 * q + 3))
            fb = {s: fit_enough(xc, b, q, kernel, s) for s in grp}
            mb = {s: coef(fb[s], y, p + 1) for s in grp}
            Vb = sum(deriv_var(fb[s], s2, p + 1) for s in grp) / factorial(p + 1) ** 2
            d2 = mb[grp[0]] if kind == "two" else comb(mb, p + 1, p + 1, kind)
            Ch = n * v * tot(V0) / (2 * (p + 1) * B0**2 * (d2**2 + 3 * Vb))
            h = Ch ** (1 / (2 * p + 3)) * n ** (-1 / (2 * p + 3))
            for s in grp:
                hh[s], bb[s] = h, b
        out[kind] = (hh[-1], hh[1], bb[-1], bb[1])
    return out


# ---------------------------------------------------------------- IK

def _ik(xc, y):
    """Imbens-Kalyanaraman (2012) Section 4.2, eqs. (4.8)-(4.13); edge-kernel constant C_K = 3.4375."""
    N = len(xc)
    sx = np.std(xc, ddof=1)
    h1 = 1.84 * sx * N ** (-0.2)
    L = (xc >= -h1) & (xc < 0)
    R = (xc >= 0) & (xc <= h1)
    nL, nR = L.sum(), R.sum()
    f = (nL + nR) / (N * h1)  # (4.8)
    s2 = (((y[L] - y[L].mean()) ** 2).sum() + ((y[R] - y[R].mean()) ** 2).sum()) / (nL + nR)  # (4.9)
    medL, medR = np.median(xc[xc < 0]), np.median(xc[xc >= 0])
    keep = (xc >= medL) & (xc <= medR)
    xk = xc[keep]
    X = np.column_stack([np.ones(keep.sum()), (xk >= 0).astype(float), xk, xk**2, xk**3])  # (4.10)
    m3 = 6 * np.linalg.lstsq(X, y[keep], rcond=None)[0][4]
    base = (s2 / (f * max(m3**2, 0.01))) ** (1 / 7)  # (4.11)
    Nside = {1: int((xc >= 0).sum()), -1: int((xc < 0).sum())}
    m2, r = {}, {}
    for s in (1, -1):
        h2 = 3.56 * base * Nside[s] ** (-1 / 7)
        w = (xc >= 0) & (xc <= h2) if s > 0 else (xc < 0) & (xc >= -h2)
        Xq = np.column_stack([np.ones(w.sum()), xc[w], xc[w] ** 2])
        m2[s] = 2 * np.linalg.lstsq(Xq, y[w], rcond=None)[0][2]
        r[s] = 720 * s2 / (w.sum() * h2**4)  # (4.12)
    return 3.4375 * (2 * s2 / (f * ((m2[1] - m2[-1]) ** 2 + r[1] + r[-1]))) ** 0.2 * N ** (-0.2)  # (4.13)


# ---------------------------------------------------------------- public

def rd_bandwidth(y, x, cutoff=0.0, *, fuzzy=None, p=1, q=None, kernel="triangular",
                 bwselect="mserd", cluster=None) -> BandwidthResult:
    """Data-driven bandwidths for RD estimation.

    bwselect: "mse"/"cer" + "rd" (common, default "mserd"), "two" (side-specific), "sum" (MSE of the sum of
    one-sided estimators), "comb1" (min of rd and sum), "comb2" (per-side median of two, rd, sum), as listed in
    CCFT (2017) sec. 2.3.  CER = n^(-p/((3+p)(3+2p))) times the MSE bandwidth.  "ik": IK (2012).
    """
    kernel = kernel_name(kernel)
    q = p + 1 if q is None else q
    if q <= p:
        raise ValueError("q must exceed p")
    (y, x, fz, cl), _ = clean(y, x, fuzzy, cluster)
    xc = x - cutoff
    nL, nR = int((xc < 0).sum()), int((xc >= 0).sum())
    if bwselect == "ik":
        if p != 1 or kernel != "triangular" or fz is not None:
            raise ValueError("bwselect='ik' is implemented for sharp designs, p=1, triangular kernel only")
        h = float(_ik(xc, y))
        return BandwidthResult(h, h, h, h, "ik", kernel, p, q, nL, nR)
    z = y
    if fz is not None:
        # pilot fuzzy effect from local fits at the Silverman pilot bandwidth; the MSE problem for
        # the ratio is the sharp problem for the pseudo-outcome y - tau*d (CCT supp. Lemma 2).
        v = _pre_pilot(xc, y, p, kernel)
        fits = {s: fit_enough(xc, v, p, kernel, s) for s in (1, -1)}
        ty = sum(s * (fits[s].coef_weights(0) @ y[fits[s].idx]) for s in (1, -1))
        td = sum(s * (fits[s].coef_weights(0) @ fz[fits[s].idx]) for s in (1, -1))
        z = y - ty / td * fz
    m = re.fullmatch(r"(mse|cer)(rd|two|sum|comb1|comb2)", bwselect)
    if not m:
        raise ValueError("bwselect must be mse|cer + rd|two|sum|comb1|comb2, or ik")
    kind, variant = m.groups()
    allv = _cct(xc, z, p, q, kernel)
    rd, two, sm = allv["rd"], allv["two"], allv["sum"]   # (hl, hr, bl, br)
    if variant in ("rd", "two", "sum"):
        hl, hr, bl, br = allv[variant]
    elif variant == "comb1":  # CCFT 2017 sec. 2.3: h_comb = min(h_rd, h_sum); b follows the chosen one
        hl, hr, bl, br = rd if rd[0] <= sm[0] else sm
    else:  # comb2: per side the median of (two, rd, sum); b of the selected estimator
        pick = {}
        for i_, nm in ((0, "l"), (1, "r")):
            cands = sorted([(two[i_], two[i_ + 2]), (rd[i_], rd[i_ + 2]), (sm[i_], sm[i_ + 2])])
            pick[nm] = cands[1]
        hl, bl = pick["l"]
        hr, br = pick["r"]
    scale = 1.0
    if kind == "cer":  # CCFT 2017 sec. 2.3 / CCF 2020 sec. 5: n^(-p/((3+p)(3+2p))) h_mse
        N = len(np.unique(cl)) if cl is not None else len(xc)
        scale = N ** (-p / ((3 + p) * (3 + 2 * p)))
    cap = {1: np.abs(xc[xc >= 0]).max(), -1: np.abs(xc[xc < 0]).max()}
    hh = {1: float(min(hr * scale, cap[1])), -1: float(min(hl * scale, cap[-1]))}
    bb = {1: float(min(br, cap[1])), -1: float(min(bl, cap[-1]))}
    return BandwidthResult(hh[-1], hh[1], bb[-1], bb[1], bwselect, kernel, p, q, nL, nR)
