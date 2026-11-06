"""Nearest-neighbour variance on a running variable with ties: the neighbours must not depend on row order or the search tree."""

import numpy as np
import pytest

from rdlocal import rd_bandwidth, rd_estimate
from rdlocal.bandwidth import _nn_sigma2, sigma2_nn


def brute_force(x, z, J=3, tol=1e-14):
    """The definition, O(n^2): every other point within the distance of the J-th nearest other point."""
    out = np.zeros(len(x))
    for i in range(len(x)):
        others = np.delete(np.arange(len(x)), i)
        dist = np.abs(x[others] - x[i])
        chosen = others[dist <= np.sort(dist)[J - 1] + tol]
        M = len(chosen)
        out[i] = M / (M + 1) * (z[i] - z[chosen].mean()) ** 2
    return out


def tied_design(n=2000, decimals=2, seed=7):
    rng = np.random.default_rng(seed)
    x = rng.normal(0, 1, n).round(decimals)
    y = 5 + 1.5 * x + 2 * (x >= 0) + rng.normal(0, 1, n)
    return x, y


def test_a_case_worked_by_hand():
    # J = 2. x = 1: the 2nd nearest is at distance 1, and all 4 others are within it.
    x, z = np.array([0.0, 0.0, 0.0, 1.0, 2.0]), np.array([1.0, 2.0, 3.0, 4.0, 6.0])
    got = _nn_sigma2(x, z, J=2)
    # x = 0, z = 1: others at distance <= 0 are the two other zeros, mean 2.5, M = 2.
    assert got[0] == pytest.approx(2 / 3 * (1 - 2.5) ** 2)
    assert got[3] == pytest.approx(4 / 5 * (4 - 3) ** 2)
    assert got[4] == pytest.approx(4 / 5 * (6 - 2.5) ** 2)  # x = 2: 2nd nearest is at distance 2, all 4 others qualify


@pytest.mark.parametrize("decimals", [0, 1, 2, 3])
def test_matches_the_definition_on_tied_data(decimals):
    x, y = tied_design(300, decimals)
    side = x >= 0
    assert _nn_sigma2(x[side], y[side]) == pytest.approx(brute_force(x[side], y[side]), rel=1e-9)


def test_untied_data_keep_exactly_j_neighbours():
    rng = np.random.default_rng(1)
    x, z = rng.uniform(0.001, 1, 200), rng.normal(0, 1, 200)
    assert _nn_sigma2(x, z) == pytest.approx(brute_force(x, z), rel=1e-9)
    near = np.argsort(np.abs(x - x[0]))[1:4]
    assert _nn_sigma2(x, z)[0] == pytest.approx(3 / 4 * (z[0] - z[near].mean()) ** 2)


def test_points_equidistant_on_a_decimal_grid_are_ties_despite_rounding():
    x = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])  # 0.3 - 0.2 != 0.2 - 0.1 in floating point
    z = np.array([0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0])
    got = _nn_sigma2(x, z, J=1)
    assert got[1] == pytest.approx(2 / 3 * (1 - 2.0) ** 2)  # neighbours of 0.2 are 0.1 and 0.3


@pytest.mark.parametrize("decimals", [1, 2, 3])
def test_row_order_does_not_change_any_result(decimals):
    x, y = tied_design(decimals=decimals)
    shuffle = np.random.default_rng(3).permutation(len(x))
    assert sigma2_nn(x[shuffle], y[shuffle]) == pytest.approx(sigma2_nn(x, y)[shuffle], rel=1e-12)
    a, b = rd_bandwidth(y, x), rd_bandwidth(y[shuffle], x[shuffle])
    assert (a.h_left, a.b_left) == pytest.approx((b.h_left, b.b_left), rel=1e-9)  # summation order leaves about 1e-12
    ra, rb = rd_estimate(y, x), rd_estimate(y[shuffle], x[shuffle])
    assert (ra.coef, ra.se_robust) == pytest.approx((rb.coef, rb.se_robust), rel=1e-9)


def test_the_estimator_stays_unbiased_at_mass_points():
    rng = np.random.default_rng(11)
    x = rng.integers(0, 12, 20000).astype(float)  # like a test score: every value heavily tied
    z = 2.0 * rng.normal(size=20000)  # constant conditional mean, variance 4
    assert _nn_sigma2(x, z).mean() == pytest.approx(4.0, rel=0.03)


@pytest.mark.parametrize("outlier", [1e3, 1e6, 1e9])
def test_an_outlier_in_x_does_not_merge_distinct_points(outlier):
    rng = np.random.default_rng(5)
    x, z = rng.uniform(0, 1, 2000), rng.normal(size=2000)
    x[0] = outlier
    assert _nn_sigma2(x, z) == pytest.approx(brute_force(x, z, tol=0.0), rel=1e-9)


def test_a_large_outcome_mean_loses_no_accuracy():
    rng = np.random.default_rng(6)
    x, z = rng.normal(size=5000), rng.normal(size=5000)
    assert _nn_sigma2(x, z + 1e8) == pytest.approx(_nn_sigma2(x, z), rel=1e-5, abs=1e-5)


def test_heavy_ties_use_the_same_definition():
    rng = np.random.default_rng(8)
    x, z = rng.integers(0, 3, 400).astype(float), rng.normal(size=400)
    assert _nn_sigma2(x, z) == pytest.approx(brute_force(x, z), rel=1e-9)


def test_a_non_finite_running_variable_is_refused():
    with pytest.raises(ValueError, match="finite"):
        _nn_sigma2(np.array([0.0, 1.0, np.nan, 3.0, 4.0]), np.zeros(5))


def test_integer_inputs_work():
    x, z = np.array([0, 0, 1, 2, 3, 5]), np.array([1, 2, 3, 4, 5, 6])
    assert _nn_sigma2(x, z) == pytest.approx(_nn_sigma2(x.astype(float), z.astype(float)))


def test_too_few_observations_still_raise():
    with pytest.raises(ValueError, match="too few"):
        _nn_sigma2(np.array([0.0, 1.0, 2.0]), np.zeros(3), J=3)
