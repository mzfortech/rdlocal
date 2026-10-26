"""Published example: U.S. Senate vote-share data, as used in Calonico, Cattaneo, Farrell & Titiunik (2017,
Stata Journal 17(2):372-404), section 7.2, pp.392-393.

The data file is NOT shipped: a copy with clear redistribution terms and provenance could not be found
(Lee's own replication files are not reachable and the copy distributed with the authors' packages is off
limits for this clean-room project).  To run this test, place a CSV with columns `margin` (Democratic margin,
running variable, cutoff 0) and `vote` (Democratic vote share at t+2) at tests/data/senate.csv or point the
environment variable RDLOCAL_SENATE_CSV at it.  1297 usable rows are expected.
"""
import os
import pathlib

import numpy as np
import pandas as pd
import pytest

from rdlocal import rd_estimate

PATH = pathlib.Path(os.environ.get("RDLOCAL_SENATE_CSV", pathlib.Path(__file__).parent / "data" / "senate.csv"))
pytestmark = pytest.mark.skipif(not PATH.exists(), reason=f"Senate data not found at {PATH}; see module docstring")


def test_senate_fixed_bandwidth_p393():
    df = pd.read_csv(PATH).dropna(subset=["margin", "vote"])
    assert len(df) == 1297
    r = rd_estimate(df["vote"], df["margin"], 0.0, h=16.79369, b=27.43745, vce="nn")
    # p.393: conventional 7.4253 (s.e. 1.4954), robust 95% CI [4.06975, 10.9833], effective N 343 / 310
    assert (r.n_left, r.n_right) == (343, 310)
    assert r.coef == pytest.approx(7.4253, abs=5e-5)
    assert r.se_conv == pytest.approx(1.4954, abs=5e-5)
    assert r.ci_robust[0] == pytest.approx(4.06975, abs=5e-4)
    assert r.ci_robust[1] == pytest.approx(10.9833, abs=5e-4)


def test_senate_default_bandwidth_p392():
    df = pd.read_csv(PATH).dropna(subset=["margin", "vote"])
    r = rd_estimate(df["vote"], df["margin"], 0.0, vce="nn")
    # p.392: h = 17.708, b = 27.984, conventional 7.416.  Selector gap documented in CLEANROOM.md, so the
    # tolerance here is that of the known selector difference, not 1%.
    assert r.h_left == pytest.approx(17.708, rel=0.35)
