# Validation workflow

This directory contains validation data, figures, and scripts used to evaluate SDEtool against external software and probabilistic expectations.

## Components

### ArcGIS parity validation

Compares:

- ArcGIS output
- SDEtool (`mode="arcgis"`)

Metrics:

- IoU overlap
- angle agreement
- centroid distance
- axis differences

Reference figures:

- `ArcGIS_R_vs_ArcRef_1x.png`
- `ArcGIS_R_vs_ArcRef_2x.png`

---

### CrimeStat parity validation

Compares:

- CrimeStat output
- SDEtool (`mode="crimestat"`)

Metrics:

- IoU (Intersection over Union) overlap
- angle agreement
- metric parity in projected CRS

Reference figures:

- `CrimeStat_R_vs_CSRef_1x.png`
- `CrimeStat_R_vs_CSRef_2x.png`

---

### Synthetic rotation validation

Stress test using rotated synthetic data.

Script:

- `cs_rotate_validation.R`

Output figure:

- `CrimeStat_RotateValidation.png`

---

### Probabilistic validation

Evaluates long-run probabilistic coverage behavior under MVN assumptions.

Script:

- `probabilistic_validation.R`

Output:

- repeated simulation statistics
- probabilistic coverage figures

---

Main validation summary:

`SDE_validation.md`
