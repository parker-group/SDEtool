# Validation workflow

This directory contains validation data, figures, and scripts used to evaluate SDEtool against external software and probabilistic expectations.

## Validation metrics

**IoU (Intersection over Union)** measures geometric overlap between two polygons:

IoU = Area(A ∩ B) / Area(A ∪ B)

Interpretation:

- IoU = 1.0 → perfect overlap
- IoU ≈ 0.99 → near-identical geometry
- IoU = 0 → no overlap

SDEtool uses IoU throughout validation to compare ellipse geometry against ArcGIS and CrimeStat outputs.

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
