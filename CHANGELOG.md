## SDEtool 1.1.0

### Added
- R package documentation and help pages generated via roxygen2
- Documentation for `generate_sde_ellipses()`
- Expanded validation workflow including repeated MVN simulation validation
- Added `validation/probabilistic_validation.R`
- Discussion of MVN assumptions and probabilistic interpretation limits
- README installation and workflow improvements

### Changed
- Improved package examples and workflow documentation

### Fixed
- CRS handling improvements
- Geometry construction updates to avoid duplicate ring vertices


# Changelog

All notable changes to this project will be documented in this file.

## [v1.0.0] - 2025-07-26
### Added
- Initial release of the Standard Deviational Ellipse (SDE) Tool
- Functions to generate 1, 2, and 3 SD ellipses
- Weighted point support
- Yuill + √2 correction and degrees of freedom adjustment
- Simulated examples with visualization
- Validation against CrimeStat ellipses

