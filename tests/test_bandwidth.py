import json

import numpy as np
import pytest

from rdlocal import rd_bandwidth, rd_estimate


def test_mserd_common_and_b_larger(designs):
    d = designs["lee1"]
    bw = rd_bandwidth(d["y"], d["x"])
    assert bw.h_left == bw.h_right and bw.b_left > bw.h_left > 0
    json.dumps(bw.to_dict())


def test_msetwo_sides_differ(designs):
    d = designs["lee1"]
    bw = rd_bandwidth(d["y"], d["x"], bwselect="msetwo")
    assert bw.h_left != bw.h_right
    r = rd_estimate(d["y"], d["x"], bwselect="msetwo")
    assert (r.h_left, r.h_right) == (bw.h_left, bw.h_right)


def test_cer_is_rule_of_thumb_shrinkage(designs):
    d = designs["lee1"]
    n, p = len(d["x"]), 1
    mse = rd_bandwidth(d["y"], d["x"])
    cer = rd_bandwidth(d["y"], d["x"], bwselect="cerrd")
    # CCF (2020) sec. 5: h_cer = n^(-p/((2p+3)(p+3))) h_mse  (p=1: n^(-1/20))
    assert cer.h_right == pytest.approx(mse.h_right * n ** (-1 / 20), rel=1e-12)


def test_ik_runs_and_validates(designs):
    d = designs["lee1"]
    bw = rd_bandwidth(d["y"], d["x"], bwselect="ik")
    assert 0.05 < bw.h_left < 1 and bw.h_left == bw.b_left
    with pytest.raises(ValueError):
        rd_bandwidth(d["y"], d["x"], bwselect="ik", p=2)
    with pytest.raises(ValueError):
        rd_bandwidth(d["y"], d["x"], bwselect="nope")


def test_ik_on_smooth_curvature_reference():
    # With a pure quadratic regression function and no noise-free degeneracy, IK's h must shrink as curvature grows.
    rng = np.random.default_rng(0)
    x = rng.uniform(-1, 1, 3000)
    h = []
    for c in (0.5, 4.0):
        y = c * x**2 * np.sign(x) + rng.normal(0, 0.2, 3000)
        h.append(rd_bandwidth(y, x, bwselect="ik").h_left)
    assert h[1] < h[0]


def test_bandwidth_scale_invariance(designs):
    d = designs["lee1"]
    a = rd_bandwidth(d["y"], d["x"])
    b = rd_bandwidth(3.0 * d["y"] + 7, 5.0 * d["x"], 0.0)
    assert b.h_left == pytest.approx(5.0 * a.h_left, rel=1e-6)
