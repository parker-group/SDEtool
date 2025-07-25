# 📍 Standard Deviational Ellipse (SDE) Tool

This R tool computes **Standard Deviational Ellipses (SDEs)** for spatial point data grouped by user-defined variables.  
It supports:
- Multiple standard deviation levels (e.g., 1, 2, 3 SD)
- Weighted points (optional)
- Yuill + √2 correction (default)
- Degrees of freedom correction (default)
- Summary of ellipse shape + % of points enclosed

For background, see:  
📚 Reference:Yuill, R. S. (1971). The Standard Deviational Ellipse: An Updated Tool for Spatial Description. Geografiska Annaler: Series B, Human Geography, 53(1), 28–39. https://doi.org/10.2307/490885

📖 [ArcGIS documentation on Standard Deviational Ellipses](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-directional-distribution-standard-deviationa.htm)

---

## 🔧 Usage

### 1. Source the functions

```r
# Option 1: If running locally after cloning this repo
source("SDE_functions.r")

# Option 2: Run directly from GitHub (raw link)
source("https://raw.githubusercontent.com/parker-group/SDEtool/main/SDE_functions.r")
```

➡️ [View the SDE_functions.r script on GitHub](https://github.com/parker-group/SDEtool/blob/main/SDE_functions.r)

---

### 2. Generate synthetic test data

```r

### creating 2 'groups' of points. Group 1 will have Region = SimRegion1; Group 2 will have Region = SimRegion2
### we will later create SDEs by Region
set.seed(42)
n <- 100

group1 <- data.frame(
  X = rnorm(n, mean = 0, sd = 5),
  Y = rnorm(n, mean = 0, sd = 2),
  Region = "SimRegion1"
)

group2 <- data.frame(
  X = rnorm(n, mean = 20, sd = 3),
  Y = rnorm(n, mean = 15, sd = 6),
  Region = "SimRegion2"
)

df <- rbind(group1, group2)
```

---

### 3. Convert to spatial object and auto-detect UTM

```r
sf_pts_proj <- convert_to_sf_utm(df)
```

---

### 4. Generate SDEs

```r
#note that you can modify different components of the function here
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

### 5. View or summarize output

```r
print(sde_sf)

# Summarize % of points within each ellipse level
aggregate(percent_inside ~ sd_level, data = sde_sf, summary)
```

---

### 6. Export as shapefile (optional)

```r
sf::st_write(sde_sf, "SDE_ellipses.shp", delete_dsn = TRUE)
```

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

## ✅ To Do

- [ ] Add Shiny app version
- [ ] Integrate real shapefile upload
- [ ] Validate against other GIS packages

---

🧪 Created for internal spatial analysis. Feel free to fork or adapt!
