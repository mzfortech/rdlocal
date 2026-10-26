"""Weighted one-sided local polynomial fits."""
from __future__ import annotations

import numpy as np

from ._kernels import kvals


class Fit:
    """Kernel-weighted order-p fit on one side of the cutoff, in scaled units u = x / h."""

    def __init__(self, xc, h, p, kernel, side):
        m = (xc >= 0) if side > 0 else (xc < 0)
        idx = np.flatnonzero(m & (np.abs(xc) <= h))
        w = kvals(kernel, xc[idx] / h)
        pos = w > 0
        self.idx, self.w, self.u = idx[pos], w[pos], xc[idx[pos]] / h
        self.h, self.p, self.n = h, p, int(pos.sum())
        if self.n < p + 1 or len(np.unique(self.u)) < p + 1:
            raise ValueError(f"too few observations ({self.n}) within bandwidth {h:.6g} for order {p}")
        self.X = self.u[:, None] ** np.arange(p + 1)
        self.Ainv = np.linalg.inv(self.X.T @ (self.w[:, None] * self.X))

    def coef_weights(self, k: int) -> np.ndarray:
        """Weights w_i with coef_k (in u-units) = sum_i w_i y_i."""
        return self.Ainv[k] @ (self.X * self.w[:, None]).T

    def residuals(self, z):
        beta = self.Ainv @ (self.X.T @ (self.w * z[self.idx]))
        r = z[self.idx] - self.X @ beta
        lev = self.w * np.einsum("ij,jk,ik->i", self.X, self.Ainv, self.X)
        return r, lev
