"""rdlocal: clean-room regression discontinuity tools (BSD-3-Clause)."""
from .bandwidth import BandwidthResult, rd_bandwidth
from .estimate import RDResult, rd_estimate

__version__ = "0.1.0"
__all__ = ["rd_estimate", "rd_bandwidth", "RDResult", "BandwidthResult"]
