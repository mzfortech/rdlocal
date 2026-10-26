"""Shared helpers: JSON-friendly serialisation and input cleaning."""
from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np


def jsonable(v: Any) -> Any:
    if isinstance(v, np.ndarray):
        return [jsonable(a) for a in v.tolist()]
    if isinstance(v, (np.floating, float)):
        return float(v) if np.isfinite(v) else None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (list, tuple)):
        return [jsonable(a) for a in v]
    if isinstance(v, dict):
        return {k: jsonable(a) for k, a in v.items()}
    return v


class Result:
    """Mixin for frozen dataclass results."""

    def to_dict(self) -> dict:
        return {f.name: jsonable(getattr(self, f.name)) for f in dataclasses.fields(self)}


def clean(*arrays, cutoff=0.0):
    """Float arrays with rows containing any non-finite value dropped (None entries pass through)."""
    cols = [None if a is None else np.asarray(a, dtype=float).ravel() for a in arrays]
    n = next(len(c) for c in cols if c is not None)
    if any(c is not None and len(c) != n for c in cols):
        raise ValueError("inputs must have the same length")
    ok = np.ones(n, bool)
    for c in cols:
        if c is not None:
            ok &= np.isfinite(c)
    return [None if c is None else c[ok] for c in cols], ok
