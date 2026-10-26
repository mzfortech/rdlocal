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
    # CCF (2020) sec. 5: h_cer = n^(-p/((3+p)(3+2p))) h_mse  (p=1: n^(-1/20))
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


def test_published_cer_scaling_senate_table():
    """CCFT (2017, Stata Journal 17(2)) p.400, rdbwselect output for the U.S. Senate data (n = 1297, p = 1):
    mserd h = 17.708 and cerrd h = 12.374, i.e. h_cer = n^(-p/((3+p)(3+2p))) h_mse = n^(-1/20) h_mse.
    The same table lists cersum = 12.806 for msesum = 18.326, certwo = 11.288/12.585 for msetwo = 16.154/18.009.
    """
    n, p = 1297, 1
    f = n ** (-p / ((3 + p) * (3 + 2 * p)))
    for mse, cer in [(17.708, 12.374), (18.326, 12.806), (16.154, 11.288), (18.009, 12.585)]:
        assert mse * f == pytest.approx(cer, rel=5e-4)


def test_combination_rules_and_variants(designs):
    d = designs["lee1"]
    r = {k: rd_bandwidth(d["y"], d["x"], bwselect=k) for k in
         ("mserd", "msetwo", "msesum", "msecomb1", "msecomb2", "cerrd", "certwo", "cersum", "cercomb1", "cercomb2")}
    # CCFT (2017) sec. 2.3: comb1 = min(rd, sum), comb2 = per-side median(two, rd, sum)
    assert r["msecomb1"].h_left == min(r["mserd"].h_left, r["msesum"].h_left)
    for side in ("left", "right"):
        hs = sorted(getattr(r[k], f"h_{side}") for k in ("msetwo", "mserd", "msesum"))
        assert getattr(r["msecomb2"], f"h_{side}") == hs[1]
    n = len(d["x"])
    assert r["cersum"].h_left == pytest.approx(r["msesum"].h_left * n ** (-1 / 20), rel=1e-12)
    assert r["msetwo"].b_left != r["msetwo"].b_right  # side-specific b, as in the published table
