"""Kernels (support [-1, 1]) and exact one-sided kernel moments.

Constants such as the bias and variance factors of local polynomial fits are
integrals of polynomials times K on [0, 1]; kernels are piecewise polynomials so
the integrals are exact.  Kernel scale cancels in every ratio used downstream.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial import polynomial as P

_POLY = {"triangular": [1.0, -1.0], "epanechnikov": [1.0, 0.0, -1.0], "uniform": [1.0]}
_ALIAS = {"tri": "triangular", "epa": "epanechnikov", "uni": "uniform"}


def kernel_name(kernel: str) -> str:
    k = _ALIAS.get(kernel, kernel)
    if k not in _POLY:
        raise ValueError(f"kernel must be one of {sorted(_POLY)}, got {kernel!r}")
    return k


def kvals(kernel: str, u):
    """K(u) for |u| <= 1, else 0."""
    a = np.abs(np.asarray(u, float))
    return np.where(a <= 1, P.polyval(a, _POLY[kernel_name(kernel)]), 0.0)


def mom(kernel: str, j: int, sq: bool = False) -> float:
    """Integral over [0, 1] of u**j * K(u) (or K(u)**2 when sq)."""
    c = np.array(_POLY[kernel_name(kernel)], float)
    if sq:
        c = P.polymul(c, c)
    return float(P.polyval(1.0, P.polyint(P.polymul(c, [0.0] * j + [1.0]))))


def gamma_inv(kernel: str, r: int) -> np.ndarray:
    G = np.array([[mom(kernel, j + k) for k in range(r + 1)] for j in range(r + 1)])
    return np.linalg.inv(G)


def sandwich(kernel: str, r: int, nu: int) -> tuple[float, float]:
    """(bias const, variance const) for the nu-th coefficient of an order-r fit.

    bias  = e_nu' G^-1 L_r         (multiplies m_{r+1} h^(r+1-nu))
    var   = e_nu' G^-1 Psi G^-1 e_nu  (multiplies sigma^2 / (n f h^(1+2nu)))
    """
    Gi = gamma_inv(kernel, r)
    L = np.array([mom(kernel, r + 1 + j) for j in range(r + 1)])
    Psi = np.array([[mom(kernel, j + k, True) for k in range(r + 1)] for j in range(r + 1)])
    e = Gi[nu]
    return float(e @ L), float(e @ Psi @ e)


def lead_const(kernel: str, r: int, nu: int) -> float:
    """B_{nu,r} = nu! e_nu' Gamma_r^{-1} vartheta_{r,r+1}, one-sided kernel moments.

    Source: CCT (2014) supplement, Lemma 1 notation (B_{nu,p,r,s}) and Section S.2.6
    ("Let B_{nu,p} = nu! e_nu' Gamma_p^{-1} vartheta_{p,p+1}").
    """
    from math import factorial

    Gi = gamma_inv(kernel, r)
    th = np.array([mom(kernel, r + 1 + j) for j in range(r + 1)])
    return float(factorial(nu) * Gi[nu] @ th)


def silverman_const(kernel: str) -> float:
    """Constant c in c * min(sd, IQR/1.349) * n^(-1/5).

    Source: CCT (2014) supplement S.2.6, Step 0 (i): the printed value for the triangular
    kernel is 2.58; other kernels use the printed Silverman formula
    (8 sqrt(pi) int K^2 / (3 (int u^2 K)^2))^(1/5).
    """
    k = kernel_name(kernel)
    if k == "triangular":
        return 2.58
    k2 = 2 * mom(k, 0, True)
    m2 = 2 * mom(k, 2)
    return float((8 * np.sqrt(np.pi) * k2 / (3 * m2**2)) ** 0.2)


def min_gamma(kernel: str, r: int) -> np.ndarray:
    """Gamma_{jk} = int int min(u,v) u^j v^k K(u) K(v) du dv over [0,1]^2 (CJM 2020, Thm 1, matrix B)."""
    K = np.array(_POLY[kernel_name(kernel)], float)

    def I(a):  # polynomial of v -> int_0^v u^a K(u) du
        return P.polyint(P.polymul(K, [0.0] * a + [1.0]))

    def tot(j, k):  # int_0^1 K(v) v^k I_{j+1}(v) dv
        poly = P.polymul(P.polymul(K, [0.0] * k + [1.0]), I(j + 1))
        return float(P.polyval(1.0, P.polyint(poly)))

    return np.array([[tot(j, k) + tot(k, j) for k in range(r + 1)] for j in range(r + 1)])
