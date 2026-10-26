"""rdlocal: clean-room regression discontinuity tools (BSD-3-Clause)."""
from .bandwidth import BandwidthResult, rd_bandwidth
from .density import DensityResult, rd_density_test
from .estimate import RDResult, rd_estimate

__version__ = "0.1.0"
__all__ = ["rd_estimate", "rd_bandwidth", "rd_density_test", "RDResult", "BandwidthResult", "DensityResult"]
