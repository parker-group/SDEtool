# 🧪 Standard Deviational Ellipse (SDE) Validation

This page validates the **R-based SDEtool** against **ArcGIS** and **CrimeStat** using the same dataset (**Lepto**, *n = 24*, WGS84). We report geometric overlap (IoU), angles, centers, and % point coverage. All angles reported here use **clockwise-from-north** (`north_cw`) to match ArcGIS/CrimeStat UI.

> **Files used (WGS84 degrees):**  
> • Points: `validation/data/Lepto.csv`  
> • ArcGIS SDEs: `validation/data/1sd.shp`, `validation/data/2sd.shp`  
> • CrimeStat SDEs: `validation/data/SDECS_Lepto.shp`, `validation/data/2SDECS_Lepto.shp`  
> • R outputs written by the code below: `validation/data/LeptoSDE_RtoolArcGIS.shp`, `validation/data/LeptoSDE_RtoolCrimeStat.shp`

---

## 🔍 Summary of Results

### ArcGIS preset (R) vs ArcGIS shapefiles (1×, 2×)

| Metric | 1 SD | 2 SD |
|---|---:|---:|
| **IoU (overlap)** | **0.999873** | **0.999873** |
| **Angle Δ (north_cw)** | **0.00027°** | **0.00027°** |
| **Centroid distance** | ~**9×10⁻⁶°** (~**1 m**) | ~**9×10⁻⁶°** (~**1 m**) |
| **% points inside (unweighted)** | **75.0000%** | **91.6667%** |
| **Major axis Δ (%)** | **−0.00002%** | **−0.00002%** |
| **Minor axis Δ (%)** | **−0.00074%** | **−0.00074%** |

> **Config:** `mode="arcgis"`, `df = n`, `scale = k·√2`, `angle = north_cw`, computed **in WGS84 degrees**.

---

### CrimeStat preset (R) vs CrimeStat shapefiles (1×, 2×)

| Metric | 1 SD | 2 SD |
|---|---:|---:|
| **IoU (overlap)** | **0.996356** | **0.993754** |
| **Angle Δ (north_cw)** | **0.062°** | **0.069°** |
| **Centroid distance** | **0.003182°** (~**0.35 km**) | **0.012736°** (~**1.4 km**) |
| **% points inside (unweighted)** | **79.1667%** | **91.6667%** |

> **Config:** `mode="crimestat"`, `df = n−2`, `scale = k·√2`, `angle = north_cw`.  
> Geometry parity with CrimeStat’s shapefiles is best when computed **in degrees**.  
> **Meters parity** (SDs/areas) is best in a projected CRS (UTM); remaining diffs are within **~0.03–0.3%**.

---

## 📐 Methodology Notes

- **Angle basis:** `north_cw`. Convert to mathematically conventional `east_ccw` with `east_ccw = (90 − north_cw) mod 360`.
- **Containment:** points on the boundary count as inside.
- **Shapefile fields:** DBF truncates names (e.g., `sd_level` → `sd_levl`, `angle_north_cw` → `angl_n_`).

**Presets implemented in SDEtool**

| Mode        | df used | Scale | Angle basis | Typical use |
|-------------|---------|-------|-------------|-------------|
| `arcgis`    | `n`     | `k·√2`| `north_cw`  | Byte-match ArcGIS |
| `crimestat` | `n−2`   | `k·√2`| `north_cw`  | Match CrimeStat |
| `prob`      | `n−1`   | `√(qchisq(p, df=2))` | `east_ccw` (convert/report as needed) | Target exact coverage \(p\) |

---

## 🖼️ Visual Comparison (reproducible)

The code below builds overlays in **WGS84** and writes PNGs to `validation/figures/`.

    # --- setup ---
    library(sf)
    library(dplyr)
    library(ggplot2)
    dir.create("validation/figures", recursive = TRUE, showWarnings = FALSE)

    # Load points (WGS84)
    pts <- read.csv("validation/data/Lepto.csv")
    pts_sf <- st_as_sf(pts, coords = c("longitude", "latitude"), crs = 4326)

    # Load ArcGIS reference ellipses (1x, 2x)
    arc_1 <- st_read("validation/data/1sd.shp", quiet=TRUE)
    arc_2 <- st_read("validation/data/2sd.shp", quiet=TRUE)

    # Generate R-tool ArcGIS preset in degrees (exact parity)
    ell_arc <- generate_sde_ellipses(
      pts_sf, group_vars = character(0), sd_levels = c(1,2),
      mode = "arcgis", compute_in = "input", output_crs = "input", return_metric = FALSE
    )
    st_write(ell_arc, "validation/data/LeptoSDE_RtoolArcGIS.shp", delete_dsn = TRUE, quiet=TRUE)

    # Load CrimeStat reference ellipses (1x, 2x)
    cs_1 <- st_read("validation/data/SDECS_Lepto.shp", quiet=TRUE)
    cs_2 <- st_read("validation/data/2SDECS_Lepto.shp", quiet=TRUE)

    # Generate R-tool CrimeStat preset in degrees (geometry parity)
    ell_cs <- generate_sde_ellipses(
      pts_sf, group_vars = character(0), sd_levels = c(1,2,3),
      mode = "crimestat", compute_in = "input", output_crs = "input", return_metric = FALSE
    )
    st_write(ell_cs, "validation/data/LeptoSDE_RtoolCrimeStat.shp", delete_dsn = TRUE, quiet=TRUE)

    # Helper: overlay plot function
    plot_overlay <- function(ref_sf, rt_sf, title, outfile) {
      gg <- ggplot() +
        geom_sf(data = ref_sf, fill = NA, color = "black", linewidth = 0.7) +
        geom_sf(data = rt_sf, aes(fill = as.factor(sd_level)), alpha = 0.35, color = NA) +
        geom_sf(data = pts_sf, color = "red", size = 1.1) +
        scale_fill_brewer(palette = "Set2", name = "R sd_level") +
        coord_sf(expand = FALSE) +
        theme_minimal() +
        labs(title = title, subtitle = "Reference (black outline) vs R-tool (filled)",
             x = "Longitude (°)", y = "Latitude (°)")
      ggsave(outfile, gg, width = 7.2, height = 5.2, dpi = 300)
    }

    # Build ArcGIS overlays (1x and 2x)
    plot_overlay(arc_1, ell_arc %>% filter(sd_level==1), 
                 "ArcGIS 1× SDE: Reference vs R-tool", 
                 "validation/figures/ArcGIS_1x_overlay.png")
    plot_overlay(arc_2, ell_arc %>% filter(sd_level==2), 
                 "ArcGIS 2× SDE: Reference vs R-tool", 
                 "validation/figures/ArcGIS_2x_overlay.png")

    # Build CrimeStat overlays (1x and 2x)
    plot_overlay(cs_1, ell_cs %>% filter(sd_level==1), 
                 "CrimeStat 1× SDE: Reference vs R-tool", 
                 "validation/figures/CrimeStat_1x_overlay.png")
    plot_overlay(cs_2, ell_cs %>% filter(sd_level==2), 
                 "CrimeStat 2× SDE: Reference vs R-tool", 
                 "validation/figures/CrimeStat_2x_overlay.png")

**Suggested figure layout in this page:**

- `validation/figures/ArcGIS_1x_overlay.png`  
- `validation/figures/ArcGIS_2x_overlay.png`  
- `validation/figures/CrimeStat_1x_overlay.png`  
- `validation/figures/CrimeStat_2x_overlay.png`

---

## 📊 Reproduce the key metrics

    # IoU helper
    iou <- function(a, b) {
      inter <- suppressMessages(st_area(st_intersection(a, b)))
      uni   <- suppressMessages(st_area(st_union(a, b)))
      as.numeric(inter / uni)
    }

    # Metrics: ArcGIS (expect ≈ 0.999873 for both)
    r1 <- ell_arc %>% dplyr::filter(sd_level==1)
    r2 <- ell_arc %>% dplyr::filter(sd_level==2)
    m_arc <- tibble::tibble(
      metric = c("IoU", "IoU"),
      level  = c(1, 2),
      value  = c(iou(r1, arc_1), iou(r2, arc_2))
    )

    # Metrics: CrimeStat (expect ≈ 0.996356 / 0.993754)
    c1 <- ell_cs %>% dplyr::filter(sd_level==1)
    c2 <- ell_cs %>% dplyr::filter(sd_level==2)
    m_cs <- tibble::tibble(
      metric = c("IoU", "IoU"),
      level  = c(1, 2),
      value  = c(iou(c1, cs_1), iou(c2, cs_2))
    )

    m_arc
    m_cs

> Minor numerical jitter across platforms is expected (tolerance ~1e-6–1e-4).

---

## 🎲 (Optional) Probabilistic Ellipses

This documents an **inference-oriented** option that targets coverage levels \(p\) under a **bivariate normal** assumption:

- **Scale:** `sqrt(qchisq(p, df=2))`  
- **df:** `n−1` (sample covariance)  
- **Coverage:** e.g., `p = c(0.6827, 0.95, 0.9973)` (≈ “1σ/2σ/3σ” under MVN)

    # Generate probabilistic ellipses (coverage targets p)
    ell_prob <- generate_sde_ellipses(
      pts_sf, group_vars = character(0),
      mode = "prob", coverage = c(0.6827, 0.95, 0.9973),
      compute_in = "input", output_crs = "input", return_metric = FALSE
    )

    # Empirical coverage check (what % of points fall inside each p-ellipse?)
    emp_cov <- ell_prob %>%
      rowwise() %>%
      mutate(empirical_pct = {
        inside <- sf::st_within(pts_sf, geometry, sparse = FALSE)[,1]
        mean(inside) * 100
      }) %>%
      ungroup() %>%
      dplyr::select(target_coverage, empirical_pct) %>%
      dplyr::arrange(target_coverage)

    emp_cov

    # Optional overlay plot
    ggplot() +
      geom_sf(data = ell_prob, aes(fill = as.factor(target_coverage)), alpha = 0.35, color = NA) +
      geom_sf(data = pts_sf, color = "red", size = 1.1) +
      scale_fill_brewer(palette = "Set2", name = "Coverage") +
      theme_minimal() + labs(title = "Probabilistic Ellipses (MVN)", x = "Lon (°)", y = "Lat (°)")

> **Interpretation:** In a single dataset, empirical % inside will vary around the target \(p\). Across repeated samples from an MVN process, the long-run average tends to \(p\).

---

## ✅ Takeaways

- **ArcGIS parity** is effectively exact when computed **in degrees** with `df=n`, `k·√2`, `north_cw`.  
- **CrimeStat parity** in degrees is **very close**; meter SDs/areas match CrimeStat’s text within **~0.03–0.3%** when computed in a projected CRS.  
- A **probabilistic** mode is available to target exact coverages under MVN assumptions, documented here for transparency.
