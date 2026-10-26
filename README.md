# rdlocal

Regression discontinuity (RD) estimation in plain NumPy/SciPy/pandas, BSD-3-Clause.

* `rd_estimate`: sharp and fuzzy RD, local polynomial of order p, kernels triangular / epanechnikov / uniform,
  conventional and bias-corrected estimates, conventional and robust SEs, CIs and p-values
  (Calonico-Cattaneo-Titiunik 2014), HC0-HC3, nearest-neighbour and cluster-robust variances.
* `rd_bandwidth`: MSE-optimal common (`mserd`) and side-specific (`msetwo`) bandwidths (CCFT 2019 supplement, Section 12), `msesum`, the
  `comb1`/`comb2` combinations, their CER versions (`cerrd`, `certwo`, `cersum`, `cercomb1/2`; CCFT 2017 sec. 2.3),
  and Imbens-Kalyanaraman (`ik`, 2012).
* `rd_density_test`: local polynomial density manipulation test (Cattaneo-Jansson-Ma 2020).
* `rd_binned`: evenly spaced or quantile-spaced binned means for RD plots.

All functions return frozen dataclasses with a JSON-friendly `to_dict()`.

```python
import numpy as np
from rdlocal import rd_estimate, rd_density_test, rd_binned

rng = np.random.default_rng(0)
x = rng.uniform(-1, 1, 2000)
y = 0.5 * x + 0.4 * (x >= 0) + rng.normal(0, 0.3, 2000)

r = rd_estimate(y, x, cutoff=0.0)            # p=1, triangular, mserd, hc1
print(r.coef, r.coef_bc, r.ci_robust, r.h_left, r.n_left, r.n_right)
fz = rd_estimate(y, x, fuzzy=(x >= 0).astype(float))
print(rd_density_test(x).p_value)
print(rd_binned(y, x, nbins=20).to_dict()["y_mean"][:3])
```

Treatment is `x >= cutoff`.  Non-finite rows are dropped.  With `cluster=...` the `vce` argument is ignored.

## Validation

`pytest` runs the fast suite; `pytest -m slow` runs the Monte Carlo coverage check.  At a fixed bandwidth with
a uniform kernel the sharp estimate equals plain OLS on each side (1e-10).  Against `rdrobust`/`rddensity`,
run only as black boxes (inputs and numbers in `tests/data`, see `tests/make_reference.py`), estimates, SEs and
bias-corrected density values agree to ~1e-9 or better at the same bandwidths; bandwidth h given b agrees to <=1%; the full selector agrees to 1% on 14 of 56 tested cells and ~1% on the
Senate default (the d-pilot of Section 12.2 is not recoverable, see `CLEANROOM.md`)
(see `CLEANROOM.md`, deviations, and the measured gaps in the test docstrings).

Published examples: the CCFT (2017) Senate output (pp.392-393, 400) is encoded in `tests/test_published_senate.py`,
which skips with an explanatory message unless you supply the Senate vote-margin data (not shipped).  With the data
it reproduces p.393 (7.4253, SE 1.4954, robust CI [4.06975, 10.9833], N 343/310).  The CER scaling and combination structure of the published table are
checked without the data.  Simulation designs of CCT (2014) supplement S.3.1 provide the other checks.

## Roadmap (out of scope for 0.1.0)

Covariate adjustment (CCFT 2019), multiple cutoffs/scores, kink designs, McCrary (2008) test, restricted
density model and jackknife SE, CER arg-min bandwidth, local randomisation / power tools.

## Methodology

[`docs/METHODOLOGY.md`](docs/METHODOLOGY.md) explains the statistics (local polynomial RD, bandwidth selection, robust bias-corrected inference, manipulation testing), how the package was built, and the validation evidence, including the reproduced US Senate results and the one open problem.

## Provenance

An independent clean-room implementation written from the papers listed in `CLEANROOM.md`.  It is not
affiliated with, endorsed by, or derived from the authors' `rdpackages` software.

## Citations

* Calonico, Cattaneo, Titiunik (2014), Econometrica 82(6): 2295-2326.
* Calonico, Cattaneo, Farrell (2020), Econometrics Journal 23(2): 192-210.
* Calonico, Cattaneo, Farrell, Titiunik (2019), Review of Economics and Statistics 101(3): 442-451.
* Cattaneo, Jansson, Ma (2020), Journal of the American Statistical Association 115(531): 1449-1455.
* Cattaneo, Jansson, Ma (2018), Stata Journal 18(1): 234-261.
* Imbens, Kalyanaraman (2012), Review of Economic Studies 79(3): 933-959.
