# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

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
