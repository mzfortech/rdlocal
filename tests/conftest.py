import json
import pathlib

import numpy as np
import pytest

DATA = pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def designs():
    z = np.load(DATA / "designs.npz")
    out = {}
    for k in z.files:
        name, var = k.split("__")
        out.setdefault(name, {})[var] = z[k]
    return out


@pytest.fixture(scope="session")
def reference():
    return json.loads((DATA / "reference.json").read_text())
