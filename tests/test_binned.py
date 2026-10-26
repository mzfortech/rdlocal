import json

import numpy as np
import pytest

from rdlocal import rd_binned


def test_binned_even_and_quantile():
    rng = np.random.default_rng(2)
    x = rng.uniform(-1, 1, 1000)
    y = x + rng.normal(0, 0.1, 1000)
    for sp in ("evenly", "quantile"):
        b = rd_binned(y, x, 0.0, nbins=(10, 8), spacing=sp)
        assert sum(b.n) == 1000
        assert set(b.side) == {-1, 1}
        assert max(b.hi[i] for i in range(len(b.n)) if b.side[i] < 0) <= 0.0 + 1e-12
        json.dumps(b.to_dict())
    q = rd_binned(y, x, 0.0, nbins=10, spacing="quantile")
    left_n = [n for n, s in zip(q.n, q.side) if s < 0]
    assert max(left_n) - min(left_n) <= 2
    with pytest.raises(ValueError):
        rd_binned(y, x, spacing="x")
