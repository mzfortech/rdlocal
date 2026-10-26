"""Regenerate tests/data/ (simulated designs + reference numbers).

rdrobust / rddensity are run strictly as opaque black boxes; only INPUTS (the simulated data and the call
arguments) and OUTPUT NUMBERS are recorded.  Run with an interpreter that has them installed:
    python make_reference.py
"""
import json
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).parent / "data"
HERE.mkdir(exist_ok=True)
COEF = {  # CCT (2014) supplement S.3.1: Models 1-3, polynomial coefficients (x^0..x^5)
    1: ([0.48, 1.27, 7.18, 20.21, 21.54, 7.33], [0.52, 0.84, -3.00, 7.99, -9.01, 3.56]),
    2: ([3.71, 2.30, 3.28, 1.45, 0.23, 0.03], [0.26, 18.49, -54.81, 74.30, -45.02, 9.83]),
    3: ([0.48, 1.27, -0.5 * 7.18, 0.7 * 20.21, 1.1 * 21.54, 1.5 * 7.33],
        [0.52, 0.84, -3.01 * 3.00, -0.3 * 7.99, -0.1 * 9.01, 3.56]),
}


def dgp(model, n, seed, sigma=0.1295):
    rng = np.random.default_rng(seed)
    x = 2 * rng.beta(2, 4, n) - 1
    cl, cr = COEF[model]
    mu = np.where(x < 0, np.polyval(cl[::-1], x), np.polyval(cr[::-1], x))
    return rng, x, mu + rng.normal(0, sigma, n)


def designs():
    D = {}
    for name, model, n, seed in [("lee1", 1, 1000, 11), ("lee2", 2, 1000, 12), ("lee3", 3, 1000, 13)]:
        _, x, y = dgp(model, n, seed)
        D[name] = dict(x=x, y=y)
    rng, x, y = dgp(1, 1500, 14)
    d = (rng.uniform(size=1500) < 0.15 + 0.7 * (x >= 0)).astype(float)
    D["fuzzy"] = dict(x=x, y=y + 0.4 * d, d=d)
    D["lee1_epa_p2"] = dict(x=D["lee1"]["x"], y=D["lee1"]["y"])
    D["lee1_uni"] = dict(x=D["lee1"]["x"], y=D["lee1"]["y"])
    D["lee2_p2"] = dict(x=D["lee2"]["x"], y=D["lee2"]["y"])
    D["lee3_p2"] = dict(x=D["lee3"]["x"], y=D["lee3"]["y"])
    rng, x, y = dgp(1, 1200, 16)
    D["cluster"] = dict(x=x, y=y, g=rng.integers(0, 60, 1200).astype(float))
    return D


def main():
    from rddensity import rddensity
    from rdrobust import rdrobust

    D = designs()
    ref = {"rdrobust": {}, "rddensity": {}}
    calls = {
        "lee1": dict(), "lee2": dict(), "lee3": dict(), "fuzzy": dict(fuzzy="d"),
        "lee1_epa_p2": dict(kernel="epanechnikov", p=2, q=3), "cluster": dict(cluster="g"),
        "lee1_uni": dict(kernel="uniform"), "lee2_p2": dict(p=2, q=3), "lee3_p2": dict(p=2, q=3),
    }
    for name, extra in calls.items():
        d = D[name]
        kw = {k: (d[v] if isinstance(v, str) and v in d else v) for k, v in extra.items()}
        for vce in ("nn", "hc0", "hc1", "hc2", "hc3"):
            if (name == "cluster" or name in ("lee1_uni", "lee2_p2", "lee3_p2")) and vce != "nn":
                continue
            r = rdrobust(d["y"], d["x"], c=0, vce=vce, **kw)
            rec = dict(args={k: (v if not isinstance(v, np.ndarray) else "<array>") for k, v in {**extra, "vce": vce}.items()},
                       coef=r.coef.values.ravel().tolist(), se=r.se.values.ravel().tolist(),
                       ci=r.ci.values.tolist(), pv=r.pv.values.ravel().tolist(),
                       h=r.bws.values[0].tolist(), b=r.bws.values[1].tolist(),
                       N_h=[int(v) for v in r.N_h])
            # same call with their bandwidths fixed: isolates the estimator from the selector
            r2 = rdrobust(d["y"], d["x"], c=0, vce=vce, h=rec["h"][0], b=rec["b"][0], **kw)
            rec["fixed_bw"] = dict(h=rec["h"][0], b=rec["b"][0], coef=r2.coef.values.ravel().tolist(),
                                   se=r2.se.values.ravel().tolist(), ci=r2.ci.values.tolist())
            if vce == "nn":  # bandwidth selectors: no regularisation, and every selector variant
                bws = {}
                for reg in (0, 1):
                    for sel in ("mserd", "msetwo", "msesum", "msecomb1", "msecomb2", "cerrd", "certwo", "cersum"):
                        if sel not in ("mserd", "cerrd") and name not in ("lee1", "lee2"):
                            continue
                        rr = rdrobust(d["y"], d["x"], c=0, vce="nn", bwselect=sel, scaleregul=reg, **kw)
                        bws[f"{sel}|reg{reg}"] = dict(h=rr.bws.values[0].tolist(), b=rr.bws.values[1].tolist())
                rec["bws"] = bws
            ref["rdrobust"][f"{name}|{vce}"] = rec
    for name in ("lee1", "lee2", "lee3"):
        x = D[name]["x"]
        r = rddensity(x, c=0.0)
        rec = dict(default=dict(h=r.h.values.tolist(), hat=r.hat.values.tolist(), sd_jk=r.sd_jk.values.tolist(),
                                t_jk=float(r.test.values[1]), p_jk=float(r.test.values[3])))
        rec["fixed"] = []
        for h in ([0.3, 0.3], [0.25, 0.45]):
            r = rddensity(x, c=0.0, h=h)
            rec["fixed"].append(dict(h=h, hat=r.hat.values.tolist(), sd_jk=r.sd_jk.values.tolist(),
                                     t_jk=float(r.test.values[1]), p_jk=float(r.test.values[3])))
        ref["rddensity"][name] = rec
    np.savez(HERE / "designs.npz", **{f"{k}__{kk}": vv for k, v in D.items() for kk, vv in v.items()})
    (HERE / "reference.json").write_text(json.dumps(ref, indent=1))


if __name__ == "__main__":
    main()
