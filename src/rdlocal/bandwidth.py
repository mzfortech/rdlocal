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
    """CCFT supplement S12.1: v_hat = C_K min(s, IQR/1.349) n^(-1/5); C_K = 2.576 (triangular), 1.843 (uniform)."""
    n = len(xc)  # y unused (kept for call symmetry)
    sx = min(np.std(xc, ddof=1), (np.percentile(xc, 75) - np.percentile(xc, 25)) / 1.349)
    return silverman_const(kernel) * sx * n ** (-0.2)


# ---------------------------------------------------------------- CCFT (2019) supplement, Section 12

def _fac(fit: Fit, nu: int, order: int) -> float:
    """O_{nu,order}(v): nu! * sum_i w_i u_i^(order+1), the finite-sample bias constant of the nu-th
    coefficient of an order-`order` local fit (CCFT supplement S12, "Bias Constants"; weights w = row nu of
    (X'WX)^-1 X'W in scaled units, so no bandwidth power appears here)."""
    return float(factorial(nu) * (fit.coef_weights(nu) @ fit.u ** (order + 1)))


def _coef_var(fit: Fit, sig2, nu: int) -> float:
    """Variance of the nu-th Taylor coefficient of the fit (deriv_var without the nu!^2 factor)."""
    return deriv_var(fit, sig2, nu) / factorial(nu) ** 2


def h_given_b(xc, y, p, q, kernel, b, scaleregul=1.0, kind="rd"):
    """Section 12.3 (Step 3) alone: h for a given bias bandwidth b (common b for rd / sum)."""
    lam = 3.0 * scaleregul
    v = _pre_pilot(xc, y, p, kernel)
    s2 = sigma2_nn(xc, y)
    S = (1, -1)
    fb = {s: fit_enough(xc, b, q, kernel, s) for s in S}
    fp = {s: fit_enough(xc, v, p, kernel, s) for s in S}
    V0 = sum(v * deriv_var(fp[s], s2, 0) for s in S)
    kap = {s: _fac(fp[s], 0, p) for s in S}
    sg = -1 if kind == "rd" else 1
    Bh = kap[1] * coef(fb[1], y, p + 1) + sg * kap[-1] * coef(fb[-1], y, p + 1)
    Rh = sum(kap[s] ** 2 * _coef_var(fb[s], s2, p + 1) for s in S)
    return float((V0 / (2 * (p + 1) * (Bh**2 + lam * Rh))) ** (1 / (3 + 2 * p)))


def _cct(xc, y, p, q, kernel, scaleregul=1.0):
    """CCFT (2019) supplement Section 12 (pp. 65-70), sharp designs, no covariates.

    Step 1: v_hat = C_K C_sd n^(-1/5).
    Step 2: b_hat = [ (3+2p)/(2(q-p)) * (V/n) / B^2 ]^(1/(3+2q)), V/n and O evaluated at v_hat, the
            (1+q)-th Taylor coefficient estimated by an order-(1+q) fit at the pilot d_hat.
    Step 3: h_hat = [ (1/(2(1+p))) * (V/n) / B^2 ]^(1/(3+2p)), B = O_{0,p}(v_hat) * m_{p+1}(b_hat).
    Regularisation: B^2 -> B^2 + 3*scaleregul*Var(B); the factor 3 and the coefficient scale are not printed
    in Section 12 and were identified from black-box agreement (see CLEANROOM.md).
    DEVIATION (d_hat): Section 12.2 defines d_hat through B_ring built from mu_hat^(1+q)(x_pm) at the
    third/first quartiles x_pm, whose meaning is not recoverable from the text; here d_hat uses the CCT (2014)
    supplement S.2.6 Step 0 rule (global polynomial of degree q+2) with the finite-sample O of Section 12.

    Returns {variant: (h_left, h_right, b_left, b_right)} for variants rd / sum / two.
    """
    n = len(xc)
    lam = 3.0 * scaleregul
    v = _pre_pilot(xc, y, p, kernel)
    s2 = sigma2_nn(xc, y)
    S = (1, -1)

    fq1 = {s: fit_enough(xc, v, q + 1, kernel, s) for s in S}
    Vd = {s: v ** (2 * q + 3) * deriv_var(fq1[s], s2, q + 1) for s in S}
    Od = {s: _fac(fq1[s], q + 1, q + 1) for s in S}
    Bd = {s: Od[s] * _global_coef(xc, y, q + 2, q + 2, s) for s in S}
    fv = {s: fit_enough(xc, v, q, kernel, s) for s in S}
    Vb = {s: v ** (2 * p + 3) * deriv_var(fv[s], s2, p + 1) for s in S}
    Ob = {s: _fac(fv[s], p + 1, q) for s in S}
    fp = {s: fit_enough(xc, v, p, kernel, s) for s in S}
    V0 = {s: v * deriv_var(fp[s], s2, 0) for s in S}
    kap = {s: _fac(fp[s], 0, p) for s in S}

    def comb(a, kind, grp):
        if kind == "two":
            return a[grp[0]]
        return a[1] - a[-1] if kind == "rd" else a[1] + a[-1]

    out = {}
    for kind in ("rd", "sum", "two"):
        hh, bb = {}, {}
        for grp in ([(1,), (-1,)] if kind == "two" else [S]):
            tot = lambda d: sum(d[s] for s in grp)
            dh = ((3 + 2 * q) / 2 * tot(Vd) / comb(Bd, kind, grp) ** 2) ** (1 / (5 + 2 * q))
            fd = {s: fit_enough(xc, dh, q + 1, kernel, s) for s in grp}
            Bb = comb({s: Ob[s] * coef(fd[s], y, q + 1) for s in grp}, kind, grp)
            Rb = sum(Ob[s] ** 2 * _coef_var(fd[s], s2, q + 1) for s in grp)
            b = ((3 + 2 * p) / (2 * (q - p)) * tot(Vb) / (Bb**2 + lam * Rb)) ** (1 / (3 + 2 * q))
            fb = {s: fit_enough(xc, b, q, kernel, s) for s in grp}
            Bh = comb({s: kap[s] * coef(fb[s], y, p + 1) for s in grp}, kind, grp)
            Rh = sum(kap[s] ** 2 * _coef_var(fb[s], s2, p + 1) for s in grp)
            h = (tot(V0) / (2 * (p + 1) * (Bh**2 + lam * Rh))) ** (1 / (3 + 2 * p))
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
                 bwselect="mserd", cluster=None, scaleregul=1.0) -> BandwidthResult:
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
    allv = _cct(xc, z, p, q, kernel, scaleregul)
    rd, two, sm = allv["rd"], allv["two"], allv["sum"]   # (hl, hr, bl, br)
    if variant in ("rd", "two", "sum"):
        hl, hr, bl, br = allv[variant]
    elif variant == "comb1":  # CCFT 2017 sec. 2.3: min(rd, sum), applied to h and to b separately
        hl, hr, bl, br = (min(rd[i_], sm[i_]) for i_ in range(4))
    else:  # comb2: per side the median of (two, rd, sum), for h and for b separately
        med = lambda i_: sorted([two[i_], rd[i_], sm[i_]])[1]
        hl, hr, bl, br = (med(i_) for i_ in range(4))
    scale = 1.0
    if kind == "cer":  # CCFT 2017 sec. 2.3 / CCF 2020 sec. 5: n^(-p/((3+p)(3+2p))) h_mse
        N = len(np.unique(cl)) if cl is not None else len(xc)
        scale = N ** (-p / ((3 + p) * (3 + 2 * p)))
    cap = {1: np.abs(xc[xc >= 0]).max(), -1: np.abs(xc[xc < 0]).max()}
    hh = {1: float(min(hr * scale, cap[1])), -1: float(min(hl * scale, cap[-1]))}
    bb = {1: float(min(br, cap[1])), -1: float(min(bl, cap[-1]))}
    return BandwidthResult(hh[-1], hh[1], bb[-1], bb[1], bwselect, kernel, p, q, nL, nR)
