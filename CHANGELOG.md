# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.1.1] - 2026-11-06

### Fixed
- Nearest-neighbour variance on a running variable with ties. Among tied or equidistant points the
  neighbours were whichever `scipy.spatial.cKDTree.query` returned, so results depended on row order and
  on the scipy version. Simulated, n = 2000, `x` rounded to one decimal, 30 seeds: shuffling the rows moved
  the bandwidth by 6% on average and by up to 22%. All points within the distance of the J-th nearest are
  now neighbours, with the J/(J+1) factor following their count, so row order no longer matters and results
  on scipy 1.18 and 1.14 agree to about 1e-12. Neighbours come from the sorted data with no search tree,
  which is also far faster on tied data, at about 1.9 times the peak memory on very large inputs.

### Changed
- Results with `vce="nn"` (the default) and data-driven bandwidths change on data with ties. In the same
  simulation, against 0.1.0: bandwidth by 4.2% on average (up to 11.5%) and robust standard error by 1.8%
  (up to 5.2%) at one decimal; 0.4% (up to 1.4%) and 0.2% (up to 0.5%) at three decimals. Results on data
  without ties do not change beyond rounding error.

## [0.1.0] - 2026-10-26

Initial public release.

### Added
- `rd_estimate`: sharp and fuzzy regression discontinuity estimation with robust bias-corrected
  inference and HC, nearest-neighbour and cluster-robust variance estimators.
- `rd_bandwidth`: MSE-optimal and CER bandwidth selectors, and the Imbens-Kalyanaraman (IK) selector.
- `rd_density_test`: local-polynomial density (manipulation) test at the cutoff.
- `rd_binned`: binned scatter summaries for RD plots.
- Typed package (`py.typed`).

### Known limitations
- The data-driven bandwidths do not reproduce every reference implementation to 1%: the d-hat pilot
  (CCFT 2019 supplement, Section 12.2) is not fully specified in the permitted sources. See the
  deviations in [CLEANROOM.md](CLEANROOM.md).
