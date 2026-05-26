
# Changelog

All notable changes to this project will be documented in this file.

## [v1.1.0] - 2026-05-25

### Added
- R package documentation and help pages generated with roxygen2
- Documentation for `generate_sde_ellipses()`
- Probabilistic SDE implementation with MVN coverage targets
- Repeated MVN simulation validation workflow
- Synthetic rotation validation against CrimeStat
- Expanded validation framework:
  - ArcGIS validation
  - CrimeStat validation
  - Synthetic validation datasets
- Validation documentation and figures
- Discussion of MVN assumptions and covariance ellipse interpretation

### Changed
- Reorganized validation structure into:
  - `validation/data/ArcGISValid/`
  - `validation/data/CrimeStatValid/`
  - `validation/data/Synthetic/`
  - `validation/scripts/`
  - `validation/figures/`
- Improved package examples and workflow documentation
- Updated README installation and usage guidance

### Fixed
- CRS handling improvements
- Geometry construction updates to avoid duplicate ring vertices
- Relative path handling in validation workflows
- Rotation validation workflow cleanup****

## [v1.0.0] - 2025-07-26
### Added
- Initial release of the Standard Deviational Ellipse (SDE) Tool
- Functions to generate 1, 2, and 3 SD ellipses
- Weighted point support
- Yuill + √2 correction and degrees of freedom adjustment
- Simulated examples with visualization
- Validation against CrimeStat ellipses

