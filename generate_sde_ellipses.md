# `generate_sde_ellipses()` — Full Reference

Dictionary for every argument, return field, and common recipes for **SDEtool**.

---

## Usage

```r
sde_sf <- generate_sde_ellipses(
  sf_data        = sf_pts_proj,  # sf POINTS
  group_vars     = "Region",     # or NULL for all points together
  sd_levels      = c(1, 2, 3),
  min_points     = 5,
  sqrt2_scaling  = TRUE,         # legacy; ignored when 'mode' is used
  dof_correction = TRUE,         # legacy; ignored when 'mode' is used
  weight_col     = NULL,         # e.g., "w" for weights
  mode           = "arcgis",     # "arcgis" | "crimestat" | "prob"
  compute_in     = "working",    # "input" | "working"
  working_crs    = "auto_utm",   # e.g., 32648 or "auto_utm" (when compute_in = "working")
  output_crs     = "working",    # "input" | "working"
  return_metric  = FALSE,        # adds 'geom_metric' when useful
  coverage       = c(0.6827, 0.95, 0.9973)  # for mode = "prob"
)
```

---

## Parameters

- **`sf_data`** *(sf POINTS; required)*  
  Input point layer. Must have a valid CRS.

- **`group_vars`** *(character vector \| NULL; default in script)*  
  Columns to group by before fitting ellipses. Set `NULL` to treat **all points as one group**.

- **`sd_levels`** *(numeric; default `c(1,2,3)`)*  
  Which SD multipliers to draw **when `mode` is not `"prob"`**.

- **`min_points`** *(integer; default `5`)*  
  Skip groups smaller than this.

- **`sqrt2_scaling`, `dof_correction`** *(logical; legacy)*  
  Back-compat flags. **Ignored when `mode` is supplied.**

- **`weight_col`** *(character \| NULL; default `NULL`)*  
  Optional column name with **nonnegative weights**. If present, centers/covariances are weighted and the **Kish** effective sample size `n_eff = (Σw)^2 / Σ(w^2)` is used for df:
  - `arcgis`: df = `n_eff`
  - `crimestat`: df = `n_eff − 2`
  - `prob`: df = `n_eff − 1`

- **`mode`** *(character; default `"arcgis"`)*  
  Controls df, scaling, and the reported angle basis:

  | mode          | df used              | scale factor                               | angle basis (reported) |
  |---------------|----------------------|--------------------------------------------|------------------------|
  | `"arcgis"`    | `n` or `n_eff`       | `k · √2` (for `sd_levels`)                 | `north_cw`             |
  | `"crimestat"` | `n − 2` or `n_eff−2` | `k · √2`                                   | `north_cw`             |
  | `"prob"`      | `n − 1` or `n_eff−1` | `√(qchisq(p, df = 2))` from `coverage`     | `north_cw` (internally `east_ccw`) |

- **`compute_in`** *(character; `"input"` or `"working"`)*  
  Where calculations happen:  
  - `"input"`: compute in the **original CRS** (e.g., EPSG:4326 degrees). Best for **byte-match parity** with ArcGIS/CrimeStat shapefiles.  
  - `"working"`: compute in a **projected CRS** (e.g., UTM) for **metric** axes/areas.

- **`working_crs`** *(EPSG integer \| `"auto_utm"` \| NULL)*  
  Used only when `compute_in = "working"`. `"auto_utm"` picks a WGS84 UTM zone from the data, or pass a specific EPSG (e.g., `32648`).

- **`output_crs`** *(character; `"input"` or `"working"`)*  
  CRS of the **returned ellipse geometry**: `"input"` returns in original CRS; `"working"` returns in projected CRS.

- **`return_metric`** *(logical; default `FALSE`)*  
  If `TRUE` **and** `compute_in = "working"` **and** `output_crs = "input"`, adds a second geometry column `geom_metric` (UTM) so you can measure in meters while keeping returned ellipses in degrees.

- **`coverage`** *(numeric vector in (0,1); default `c(0.6827, 0.95, 0.9973)`)*  
  **Only used when `mode = "prob"`**; drives the scale via `√(qchisq(p, df = 2))`. `sd_levels` is ignored in this mode.

---

### CRS settings — exact behavior

- `compute_in="input"`: compute in the CRS of `sf_data`. `working_crs` is ignored; no transform.
- `compute_in="working"`: compute in `working_crs`. Transform happens **only** if `working_crs` ≠ `st_crs(sf_data)`.

**Examples**

| `st_crs(sf_data)` | `compute_in` | `working_crs`   | What happens                    |
|-------------------|--------------|------------------|----------------------------------|
| EPSG:4326 (deg)   | `"input"`    | (ignored)        | Compute in degrees (no transform) |
| EPSG:4326 (deg)   | `"working"`  | `"auto_utm"`     | Transform to UTM, compute in meters |
| EPSG:4326 (deg)   | `"working"`  | `4326`           | Compute in degrees (no transform) |
| EPSG:32611 (UTM)  | `"input"`    | (ignored)        | Compute in meters (no transform)  |
| EPSG:32611 (UTM)  | `"working"`  | `"auto_utm"`     | Stays in 32611 (identity)         |
| EPSG:32611 (UTM)  | `"working"`  | `4326`           | Transform to degrees, compute there |

- `output_crs="input"` returns geometry in the original CRS.
- `output_crs="working"` returns geometry in the working CRS.
- `return_metric=TRUE` is useful when you **compute in UTM** but want to **return degrees** (adds a `geom_metric` column in meters).
---
## Output columns (key fields)

- `sd_level` (for non-`prob`) / `target_coverage` (for `prob`)  
- `n_points`, `count_inside`, `percent_inside` (weights respected if provided)  
- `area` (in the **returned geometry’s CRS**)  
- `major_axis`, `minor_axis` (semi-axes; in **compute** CRS units)  
- `orientation_deg` (east-CCW), `angle_north_cw` (north-CW), `angle_basis`  
- `sde_mode`, `df_used`, `scale_rule`, `implied_coverage`  
- `center_x`, `center_y`, `lambda1`, `lambda2`  
- `geometry` (returned in `output_crs`), optional `geom_metric` (if requested)

---

## Quick recipes

**ArcGIS parity in degrees (recommended for comparisons):**
```r
generate_sde_ellipses(
  sf_pts, group_vars = "Region",
  mode = "arcgis",
  compute_in = "input", output_crs = "input"
)
```

**CrimeStat parity in degrees:**
```r
generate_sde_ellipses(
  sf_pts, group_vars = "Region",
  mode = "crimestat",
  compute_in = "input", output_crs = "input"
)
```

**Metric axes/areas (UTM), keep output in meters:**
```r
generate_sde_ellipses(
  sf_pts, group_vars = "Region",
  mode = "arcgis",
  compute_in = "working", working_crs = "auto_utm",
  output_crs = "working"
)
```

**Probabilistic coverage (e.g., 68/95/99.73%):**
```r
generate_sde_ellipses(
  sf_pts, group_vars = NULL,
  mode = "prob", compute_in = "input",
  coverage = c(0.6827, 0.95, 0.9973)
)
```

**Weighted example (df uses Kish n_eff):**
```r
generate_sde_ellipses(
  sf_pts_w, group_vars = "Region",
  weight_col = "w",
  mode = "arcgis", compute_in = "input"
)
```

---

## Notes & edge cases

- Groups with `< min_points` are skipped (message printed). Default is 5 min_points; it seems silly to work with less than 5. 
- Zero-variance groups (coincident points) do not produce ellipses.  
- For shapefile export, consider mapping attribute names to ≤10 chars (DBF-safe).  
- Angles are reported as `north_cw` in summaries; internal math uses `east_ccw` (converted for reporting).

---

## See also

- Validation results and methodology: [`validation/SDE_validation.md`](https://github.com/parker-group/SDEtool/blob/main/validation/SDE_validation.md)  
- Minimal [README](https://github.com/parker-group/SDEtool/blob/main/README.md) usage example that links to this page.
