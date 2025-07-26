# 📍 Standard Deviational Ellipse (SDE) Tool

This R tool computes **Standard Deviational Ellipses (SDEs)** for spatial point data grouped by user-defined variables.  
It supports:
- Multiple standard deviation levels (e.g., 1, 2, 3 SD)
- Weighted points (optional)
- Yuill + √2 correction (default)
- Degrees of freedom correction (default)
- Summary of ellipse shape + % of points enclosed

For background, see:  
📚 Reference: Yuill, R. S. (1971). *The Standard Deviational Ellipse: An Updated Tool for Spatial Description*. Geografiska Annaler: Series B, Human Geography, 53(1), 28–39. https://doi.org/10.2307/490885  
📖 [ArcGIS documentation on Standard Deviational Ellipses](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-directional-distribution-standard-deviationa.htm)

---

## 🔧 Workflow Steps

### 1. Load the required functions

You must first load the R functions before any other step will work:

```r
# Option 1: If running locally after cloning this repo
source("SDE_functions.r")

# Option 2: Run directly from GitHub (raw link)
source("https://raw.githubusercontent.com/parker-group/SDEtool/main/SDE_functions.r")
```

➡️ [View the SDE_functions.r script on GitHub](https://github.com/parker-group/SDEtool/blob/main/SDE_functions.r)

---

### 2. Load your data

Your dataset must:
- Include **geographic coordinates**
  - Either longitude/latitude (in degrees), or
  - Projected X/Y values (e.g., UTM in meters)
- Include a **grouping variable** (e.g., Region, Year, etc.) if you want to compare SDEs between subsets.

Example:

```r
df <- data.frame(
  Longitude = c(36.8, 36.9, 36.7),
  Latitude = c(-1.3, -1.4, -1.2),
  Region = c("A", "A", "B")
)
```

---

### 3. Convert to sf object and UTM projection

Use `convert_to_sf_utm()` to convert the data to a spatial object and project it to UTM automatically:

```r
sf_pts_proj <- convert_to_sf_utm(df)
```

If your data is already projected, you can manually specify the EPSG code:

```r
sf_pts_proj <- convert_to_sf_utm(df, input_crs = 32636, target_epsg = 32636)
```

---

### 4. Generate the SDEs

Use the main function to create ellipses for each group:

```r
sde_sf <- generate_sde_ellipses(
  sf_pts_proj,
  group_vars = "Region",
  sd_levels = c(1, 2, 3),
  min_points = 5,
  sqrt2_scaling = TRUE,
  dof_correction = TRUE,
  weight_col = NULL
)
```

---

### 5. Inspect and summarize results

Print the results or summarize how many points fall within each SD level:

```r
print(sde_sf)
aggregate(percent_inside ~ sd_level, data = sde_sf, summary)
```

---

### 6. Export to shapefile (optional)

If desired, export your SDEs to a shapefile:

```r
sf::st_write(sde_sf, "SDE_ellipses.shp", delete_dsn = TRUE)
```

---

## 🧪 Simulated Example 1: Latitude/Longitude Data

```r
set.seed(123)
n <- 100
df <- data.frame(
  Longitude = runif(n, min = 36.6, max = 37.0),
  Latitude = runif(n, min = -1.5, max = -1.0),
  Region = sample(c("East", "West"), n, replace = TRUE)
)

sf_pts_proj <- convert_to_sf_utm(df)
sde_sf <- generate_sde_ellipses(sf_pts_proj, group_vars = "Region")
print(sde_sf)
```

---

## 🧪 Simulated Example 2: UTM Projected X/Y Data

```r
set.seed(456)
n <- 100
df <- data.frame(
  X = rnorm(n, mean = 500000, sd = 10000),
  Y = rnorm(n, mean = 990000, sd = 12000),
  Region = sample(c("North", "South"), n, replace = TRUE)
)

sf_pts_proj <- convert_to_sf_utm(df, input_crs = 32636, target_epsg = 32636)
sde_sf <- generate_sde_ellipses(sf_pts_proj, group_vars = "Region")
print(sde_sf)
```

---

## 🛍 Coordinate System Tips

| Your Data Looks Like…                          | Coordinate Type              | What You Should Do                                       | Example Call                                                    |
|------------------------------------------------|------------------------------|----------------------------------------------------------|------------------------------------------------------------------|
| Values like `-1.3`, `36.8`                     | Latitude/Longitude (degrees) | Nothing special — default settings will work             | `convert_to_sf_utm(df)`                                         |
| GPS data from phone/app                        | Latitude/Longitude (degrees) | Default is fine — UTM zone will be auto-detected         | `convert_to_sf_utm(my_data)`                                   |
| X/Y values like `500000`, `1000000` (meters)   | Projected (e.g., UTM)        | You **must** specify the CRS (EPSG code)                 | `convert_to_sf_utm(df, input_crs = 32632, target_epsg = 32632)` |
| You're unsure what system your data is in      | 🤷 Unknown                   | Ask the data provider or check in GIS software           | —                                                                |
| You want to override auto-detect UTM           | Lat/lon or Projected         | Manually set `target_epsg` to force your own zone        | `convert_to_sf_utm(df, target_epsg = 32633)`                    |

**💡 How to Find EPSG Codes:**
- Visit [epsg.io](https://epsg.io)
- Use `sf::st_crs()` on known spatial data
- In QGIS: Right-click layer → Properties → CRS

---

## 🔬 What This Calculates

The Standard Deviational Ellipse (SDE) summarizes the spatial distribution of points by showing the directional trend and spread.  
Each ellipse covers approximately:
- **~63%** of points at 1 standard deviation
- **~98%** at 2 standard deviations
- **~99.9%** at 3 standard deviations  
Assumes approximately normal distribution in 2D space.

Ellipse orientation is defined by the **eigenvector** of the covariance matrix of X and Y — this shows the direction of greatest spread (i.e., the major axis of the ellipse).

---

🧪 Created for internal spatial analysis. Feel free to fork or adapt!
