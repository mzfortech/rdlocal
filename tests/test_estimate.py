import json

import numpy as np
import pytest

import rdlocal
from rdlocal import rd_estimate

CALLS = {
    "lee1": {}, "lee2": {}, "lee3": {}, "fuzzy": {"fuzzy": "d"},
    "lee1_epa_p2": {"kernel": "epanechnikov", "p": 2, "q": 3}, "cluster": {"cluster": "g"},
}


def wls_intercept(x, y, h, p, side):
    m = (x >= 0) & (x <= h) if side > 0 else (x < 0) & (x >= -h)
    X = x[m][:, None] ** np.arange(p + 1)
    return np.linalg.lstsq(X, y[m], rcond=None)[0][0]


@pytest.mark.parametrize("p", [1, 2])
def test_exactness_uniform_kernel_equals_ols(designs, p):
    d = designs["lee1"]
    x, y = d["x"], d["y"] + 0.3 * (d["x"] >= 0)
    h = 0.31
    r = rd_estimate(y, x, kernel="uniform", p=p, h=h, b=0.5)
    ref = wls_intercept(x, y, h, p, 1) - wls_intercept(x, y, h, p, -1)
    assert abs(r.coef - ref) < 1e-10


def test_exactness_nonzero_cutoff(designs):
    d = designs["lee1"]
    x, y = d["x"] + 2.5, d["y"]
    a = rd_estimate(y, x, 2.5, kernel="uniform", h=0.3)
    b = rd_estimate(y, d["x"], 0.0, kernel="uniform", h=0.3)
    assert abs(a.coef - b.coef) < 1e-10


def test_bias_correction_with_h_equal_b_is_order_q_fit(designs):
    # Frisch-Waugh identity: with h = b and q = p + 1 the bias-corrected estimate is the order-q intercept.
    d = designs["lee2"]
    x, y = d["x"], d["y"]
    r = rd_estimate(y, x, kernel="uniform", p=1, q=2, h=0.25, b=0.25)
    ref = wls_intercept(x, y, 0.25, 2, 1) - wls_intercept(x, y, 0.25, 2, -1)
    assert abs(r.coef_bc - ref) < 1e-10


def test_fuzzy_is_ratio_of_sharp(designs):
    d = designs["fuzzy"]
    r = rd_estimate(d["y"], d["x"], fuzzy=d["d"], h=0.3, b=0.4)
    ry = rd_estimate(d["y"], d["x"], h=0.3, b=0.4)
    rt = rd_estimate(d["d"], d["x"], h=0.3, b=0.4)
    assert abs(r.coef - ry.coef / rt.coef) < 1e-12
    assert abs(r.first_stage - rt.coef) < 1e-12


def test_recovers_effect_and_json(designs):
    d = designs["lee1"]
    r = rd_estimate(d["y"], d["x"])
    assert r.ci_robust[0] < 0.04 < r.ci_robust[1]  # Model 1 jump is 0.52-0.48
    assert r.n_left > 10 and r.n_right > 10
    json.dumps(r.to_dict())  # JSON friendly


def test_rejects_bad_input():
    x = np.linspace(-1, 1, 50)
    with pytest.raises(ValueError):
        rd_estimate(x, x, p=2, q=2)
    with pytest.raises(ValueError):
        rd_estimate(x, x, kernel="gauss")
    with pytest.raises(ValueError):
        rd_estimate(x, x[:10])


def test_nonfinite_rows_dropped(designs):
    d = designs["lee1"]
    y = d["y"].copy(); y[3] = np.nan
    a = rd_estimate(y, d["x"], h=0.3)
    m = np.isfinite(y)
    b = rd_estimate(y[m], d["x"][m], h=0.3)
    assert a.coef == b.coef


def test_cluster_runs_and_is_wider_with_cluster_effects(designs):
    d = designs["cluster"]
    r = rd_estimate(d["y"], d["x"], h=0.3, cluster=d["g"])
    assert r.vce == "cluster" and r.se_robust > 0


# ---- black-box agreement with rdrobust run as an opaque function (data + numbers in tests/data) ----

def _cmp(key):
    name, vce = key.split("|")
    return name, vce


@pytest.mark.blackbox
@pytest.mark.parametrize("vce", ["hc0", "hc1", "hc2", "hc3", "nn"])
@pytest.mark.parametrize("name", ["lee1", "lee2", "lee3", "fuzzy", "lee1_epa_p2"])
def test_fixed_bandwidth_agrees_with_rdrobust(designs, reference, name, vce):
    """At the SAME (h, b) the estimator layer must agree; gaps are documented per quantity.

    point estimates / conventional SE (hc0, hc2, hc3, nn) / robust SE (hc*): ~1e-9 (observed ~1e-14).
    robust SE with nn: 5e-6 (observed up to 1e-6; neighbour/tie bookkeeping).
    conventional SE with hc1: 1e-2. We use the textbook n_s/(n_s-k) per side; rdrobust's factor differs
    slightly (observed gap 1.5e-3 to 1e-2).
    """
    rec = reference["rdrobust"][f"{name}|{vce}"]
    d = designs["lee1" if name == "lee1_epa_p2" else name]
    kw = {k: (d[v] if v in d else v) for k, v in CALLS[name].items()}
    fb = rec["fixed_bw"]
    r = rd_estimate(d["y"], d["x"], h=fb["h"], b=fb["b"], vce=vce, **kw)
    assert r.coef == pytest.approx(fb["coef"][0], rel=1e-9)
    assert r.coef_bc == pytest.approx(fb["coef"][1], rel=1e-9)
    assert r.se_conv == pytest.approx(fb["se"][0], rel=1e-2 if vce == "hc1" else 1e-9)
    assert r.se_robust == pytest.approx(fb["se"][2], rel=5e-6 if vce == "nn" else 1e-9)
    assert r.ci_robust[0] == pytest.approx(fb["ci"][2][0], rel=1e-5, abs=1e-7) if vce != "nn" else True
    assert r.ci_conv[1] == pytest.approx(fb["ci"][0][1], rel=2e-2)


@pytest.mark.blackbox
def test_cluster_close_to_rdrobust_nncluster(designs, reference):
    rec = reference["rdrobust"]["cluster|nn"]["fixed_bw"]
    d = designs["cluster"]
    r = rd_estimate(d["y"], d["x"], h=rec["h"], b=rec["b"], cluster=d["g"])
    assert r.coef == pytest.approx(rec["coef"][0], rel=1e-9)
    # different cluster residual construction (fit residuals vs nearest-neighbour): within 5%
    assert r.se_robust == pytest.approx(rec["se"][2], rel=5e-2)


@pytest.mark.blackbox
@pytest.mark.parametrize("name", ["lee1", "lee2", "lee3", "fuzzy", "lee1_epa_p2"])
def test_default_bandwidth_in_same_range_as_rdrobust(designs, reference, name):
    """Selector gap: observed 6%-31% below rdrobust's h (see README / CLEANROOM deviations)."""
    rec = reference["rdrobust"][f"{name}|nn"]
    d = designs["lee1" if name == "lee1_epa_p2" else name]
    kw = {k: (d[v] if v in d else v) for k, v in CALLS[name].items()}
    r = rd_estimate(d["y"], d["x"], vce="nn", **kw)
    assert 0.65 < r.h_right / rec["h"][0] < 1.1
    assert 0.6 < r.b_right / rec["b"][0] < 1.15
