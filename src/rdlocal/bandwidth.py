"""Bandwidth selectors: CCT (2014) plug-in MSE (mserd / msetwo), CCF (2020) rule-of-thumb CER, IK (2012)."""
from __future__ import annotations

from dataclasses import dataclass
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

def _cct(xc, y, p, q, kernel, two_sided):
    n = len(xc)
    sgn = (-1) ** (p + 1)
    v = _pre_pilot(xc, y, p, kernel)

    # NN variances are residuals of the data themselves, independent of the fit.
    s2 = sigma2_nn(xc, y)

    # Step 0 (ii): pilot c_n via global (q+2)-th order polynomial fits (CCT supp. S.2.6).
    fq1 = {s: fit_enough(xc, v, q + 1, kernel, s) for s in (1, -1)}
    Vq1 = sum(deriv_var(fq1[s], s2, q + 1) for s in (1, -1))
    Bq1 = lead_const(kernel, q + 1, q + 1)
    g = _global_coef(xc, y, q + 2, q + 2, 1) - (-1) ** q * _global_coef(xc, y, q + 2, q + 2, -1)
    Cc = (2 * q + 3) * n * v ** (2 * q + 3) * Vq1 / (2 * Bq1**2 * g**2)
    c_n = (Cc * 1.0) ** (1 / (2 * q + 5)) * n ** (-1 / (2 * q + 5)) if np.isfinite(Cc) and Cc > 0 else v

    # Step 1: b_n.
    fc = {s: fit_enough(xc, c_n, q + 1, kernel, s) for s in (1, -1)}
    dq = coef(fc[1], y, q + 1) - (-1) ** (q + 1) * coef(fc[-1], y, q + 1)
    Vcc = sum(deriv_var(fc[s], s2, q + 1) for s in (1, -1))
    fv = {s: fit_enough(xc, v, q, kernel, s) for s in (1, -1)}
    Vpq = sum(deriv_var(fv[s], s2, p + 1) for s in (1, -1))
    Bpq = lead_const(kernel, q, p + 1)
    # DEVIATION: the printed bracket adds 3*V_hat (derivative-scale variance, carries nu!^2) to a squared
    # Taylor coefficient; the regularisation is applied in coefficient scale (V_hat / nu!^2) so both terms
    # share units.  Verbatim form would give materially smaller bandwidths (see CLEANROOM.md).
    Cb = (2 * p + 3) * n * v ** (2 * p + 3) * Vpq / (2 * (q - p) * Bpq**2 * (dq**2 + 3 * Vcc / factorial(q + 1) ** 2))
    b = Cb ** (1 / (2 * q + 3)) * n ** (-1 / (2 * q + 3))

    # Step 2: h_n.
    fb = {s: fit_enough(xc, b, q, kernel, s) for s in (1, -1)}
    mb = {s: coef(fb[s], y, p + 1) for s in (1, -1)}
    Vb = {s: deriv_var(fb[s], s2, p + 1) for s in (1, -1)}
    fp = {s: fit_enough(xc, v, p, kernel, s) for s in (1, -1)}
    V0 = {s: deriv_var(fp[s], s2, 0) for s in (1, -1)}
    B0 = lead_const(kernel, p, 0)
    if not two_sided:
        d = mb[1] - sgn * mb[-1]
        Ch = n * v * (V0[1] + V0[-1]) / (2 * (p + 1) * B0**2 * (d**2 + 3 * (Vb[1] + Vb[-1]) / factorial(p + 1) ** 2))
        hh = Ch ** (1 / (2 * p + 3)) * n ** (-1 / (2 * p + 3))
        return {1: hh, -1: hh}, {1: b, -1: b}
    # DEVIATION (msetwo): the papers obtained do not state the two-bandwidth selector; the
    # same Step-2 formula is applied to each side separately (own variance, own bias term).
    out = {}
    for s in (1, -1):
        Ch = n * v * V0[s] / (2 * (p + 1) * B0**2 * (mb[s] ** 2 + 3 * Vb[s] / factorial(p + 1) ** 2))
        out[s] = Ch ** (1 / (2 * p + 3)) * n ** (-1 / (2 * p + 3))
    return out, {1: b, -1: b}


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

    bwselect: "mserd" (common MSE), "msetwo" (side-specific MSE), "cerrd" (rule-of-thumb CER,
    CCF 2020 sec. 5: n^(-p/((2p+3)(p+3))) times the MSE bandwidth), "ik" (IK 2012).
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
    base = {"mserd": "mserd", "msetwo": "msetwo", "cerrd": "mserd", "certwo": "msetwo"}.get(bwselect)
    if base is None:
        raise ValueError("bwselect must be one of mserd, msetwo, cerrd, certwo, ik")
    h, b = _cct(xc, z, p, q, kernel, base == "msetwo")
    scale = 1.0
    if bwselect.startswith("cer"):
        N = len(np.unique(cl)) if cl is not None else len(xc)
        scale = N ** (-p / ((2 * p + 3) * (p + 3)))
    cap = {1: np.abs(xc[xc >= 0]).max(), -1: np.abs(xc[xc < 0]).max()}
    hh = {s: float(min(h[s] * scale, cap[s])) for s in (1, -1)}
    bb = {s: float(min(b[s], cap[s])) for s in (1, -1)}
    return BandwidthResult(hh[-1], hh[1], bb[-1], bb[1], bwselect, kernel, p, q, nL, nR)
