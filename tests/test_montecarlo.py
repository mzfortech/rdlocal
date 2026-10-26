import numpy as np
import pytest

from rdlocal import rd_estimate


@pytest.mark.slow
def test_robust_ci_coverage_lee_model1():
    """Model 1 of CCT (2014) supplement S.3.1, n=1000, 500 replications; true jump = 0.52 - 0.48 = 0.04."""
    cl = [0.48, 1.27, 7.18, 20.21, 21.54, 7.33]
    cr = [0.52, 0.84, -3.00, 7.99, -9.01, 3.56]
    rng = np.random.default_rng(2024)
    hit = hit_c = 0
    reps = 500
    for _ in range(reps):
        x = 2 * rng.beta(2, 4, 1000) - 1
        y = np.where(x < 0, np.polyval(cl[::-1], x), np.polyval(cr[::-1], x)) + rng.normal(0, 0.1295, 1000)
        r = rd_estimate(y, x, vce="hc1")
        hit += r.ci_robust[0] <= 0.04 <= r.ci_robust[1]
        hit_c += r.ci_conv[0] <= 0.04 <= r.ci_conv[1]
    cov, cov_c = hit / reps, hit_c / reps
    print(f"coverage robust={cov:.3f} conventional={cov_c:.3f}")
    assert 0.90 <= cov <= 0.98  # observed 0.922 (MC s.e. 0.010)
