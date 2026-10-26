import json

import numpy as np
import pytest

from rdlocal import rd_density_test
from rdlocal.density import _side


def brute_variance(xs, F, n, h, order, kernel_w):
    """Triple-sum variance estimator of CJM (2020) Theorem 2, written naively, one (right) side."""
    u = xs / h
    m = (xs >= 0) & (np.abs(u) <= 1)
    um, Fm = u[m], F[m]
    w = kernel_w(um)
    R = um[:, None] ** np.arange(order + 1)
    A = (R * w[:, None]).T @ R / (n * h)
    B = np.zeros((order + 1, order + 1))
    for j in range(len(um)):
        for k in range(len(um)):
            ind = np.minimum(xs[m][j], xs[m][k])
            cov = ((xs <= xs[m][j]).astype(float) - F[m][j]) @ ((xs <= xs[m][k]).astype(float) - F[m][k])
            B += np.outer(R[j] * w[j], R[k] * w[k]) * cov
    B /= n**3 * h**3
    Ai = np.linalg.inv(A)
    return (Ai @ B @ Ai)[1, 1] / (n * h)


def test_variance_matches_paper_triple_sum():
    rng = np.random.default_rng(3)
    n = 120
    x = np.sort(rng.normal(size=n))
    F = np.searchsorted(x, x, side="right") / n  # Thm 2 uses the plain empirical CDF
    h = 0.9
    f, var, _ = _side(x, F, n, 0.0, h, 2, "triangular", 1)
    v2 = brute_variance(x, F, n, h, 2, lambda u: np.maximum(1 - np.abs(u), 0))
    # B-hat of Thm 2 centres by F-hat inside the sum; ours is the algebraic equivalent up to the
    # (F-hat vs F-hat(min)) bookkeeping of the discrete sum: they must agree closely.
    assert var == pytest.approx(v2, rel=0.05)


def test_no_discontinuity_not_rejected_and_discontinuity_rejected():
    rng = np.random.default_rng(5)
    x = rng.normal(size=4000)
    r = rd_density_test(x, 0.3)
    assert r.p_value > 0.01
    json.dumps(r.to_dict())
    keep = (x < 0) | (rng.uniform(size=4000) < 0.5)  # halve the density right of the cutoff
    r2 = rd_density_test(x[keep], 0.0)
    assert r2.p_value < 1e-3 and r2.f_right < r2.f_left


def test_density_estimator_recovers_normal_density():
    rng = np.random.default_rng(8)
    x = rng.normal(size=20000)
    r = rd_density_test(x, 0.0, h=0.5)
    assert r.f_left == pytest.approx(0.3989, abs=0.03) and r.f_right == pytest.approx(0.3989, abs=0.03)


@pytest.mark.blackbox
def test_density_agrees_with_rddensity_fixed_h(designs, reference):
    """Bias-corrected density estimates agree to ~1e-9; SEs within 2% (their jackknife SE, Sec. 5.2, is
    a different consistent estimator from the Thm 2 one used here); statistic within 2%."""
    for name, rec in reference["rddensity"].items():
        x = designs[name]["x"]
        for f in rec["fixed"] + [rec["default"]]:
            r = rd_density_test(x, h=f["h"])
            assert r.f_left_bc == pytest.approx(f["hat"][0], rel=1e-9)
            assert r.f_right_bc == pytest.approx(f["hat"][1], rel=1e-9)
            assert r.se_left_bc == pytest.approx(f["sd_jk"][0], rel=2e-2)
            assert r.se_right_bc == pytest.approx(f["sd_jk"][1], rel=2e-2)
            assert abs(r.statistic - f["t_jk"]) < 0.02 * max(1, abs(f["t_jk"]))
