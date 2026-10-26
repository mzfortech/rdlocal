"""Binned means for RD plots."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._base import Result, clean


@dataclass(frozen=True)
class BinnedResult(Result):
    side: list       # -1 left of the cutoff, +1 right
    x_mean: list
    y_mean: list
    n: list
    lo: list
    hi: list
    spacing: str
    cutoff: float


def rd_binned(y, x, cutoff=0.0, nbins=20, *, spacing="evenly") -> BinnedResult:
    """Means of y in bins of x, separately on each side of the cutoff (x >= cutoff is the right side).

    nbins: int (per side) or (n_left, n_right).  spacing: "evenly" (equal-width bins) or "quantile"
    (equal-count bins).  Empty bins are dropped.
    """
    if spacing not in ("evenly", "quantile"):
        raise ValueError("spacing must be 'evenly' or 'quantile'")
    (y, x), _ = clean(y, x)
    nb = {-1: nbins, 1: nbins} if np.isscalar(nbins) else {-1: nbins[0], 1: nbins[1]}
    out = {k: [] for k in ("side", "x_mean", "y_mean", "n", "lo", "hi")}
    for s in (-1, 1):
        m = (x >= cutoff) if s > 0 else (x < cutoff)
        xs, ys = x[m], y[m]
        if len(xs) == 0:
            continue
        k = int(nb[s])
        if k < 1:
            raise ValueError("nbins must be >= 1")
        if spacing == "evenly":
            edges = np.linspace(cutoff, xs.max(), k + 1) if s > 0 else np.linspace(xs.min(), cutoff, k + 1)
        else:
            edges = np.quantile(xs, np.linspace(0, 1, k + 1))
            edges = np.unique(edges)
        # right side bins are [lo, hi); the last bin is closed so the extreme point is kept
        idx = np.clip(np.searchsorted(edges, xs, side="right") - 1, 0, len(edges) - 2)
        for j in range(len(edges) - 1):
            mj = idx == j
            if mj.any():
                for key, v in zip(out, (s, xs[mj].mean(), ys[mj].mean(), int(mj.sum()), edges[j], edges[j + 1])):
                    out[key].append(v)
    return BinnedResult(**out, spacing=spacing, cutoff=float(cutoff))
