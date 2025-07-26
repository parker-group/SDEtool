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
# note that you can modify different components of the function here
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

### 6. Plot the ellipses and points on a map

```r
# load necessary packages
library(ggplot2)
library(sf)

# ggplot function to make the map
## this map will be in UTMs. It would be possible to convert back to WGS84
ggplot() +
  geom_sf(data = sde_sf, aes(fill = as.factor(sd_level)),
          alpha = 0.3, color = "black", linetype = "solid") +
  geom_sf(data = sf_pts_proj, aes(color = Region),
          size = 1.2, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", name = "SD Level") +
  scale_color_brewer(palette = "Dark2", name = "Region") +
  theme_minimal() +
  labs(
    title = "Standard Deviational Ellipses and Point Distributions",
    subtitle = "Grouped by Region",
    x = "Easting (meters)", y = "Northing (meters)"
  )
```

---

### 7. Export as shapefile (optional)

```r
sf::st_write(sde_sf, "SDE_ellipses.shp", delete_dsn = TRUE)
```

---

## 🛍 Coordinate System Tips

| Your Data Look Like…                          | Coordinate Type              | What You Should Do                                       | Example Call                                                    |
|------------------------------------------------|------------------------------|----------------------------------------------------------|------------------------------------------------------------------|
| Values like `-1.3`, `36.8`                      | Latitude/Longitude (degrees) | Nothing special — default settings will work             | `convert_to_sf_utm(df)`                                         |
| GPS data from phone/app                         | Latitude/Longitude (degrees) | Default is fine — UTM zone will be auto-detected         | `convert_to_sf_utm(my_data)`                                   |
| X/Y values like `500000`, `1000000` (meters)    | Projected (e.g., UTM)        | You **must** specify the CRS (EPSG code)                 | `convert_to_sf_utm(df, input_crs = 32632, target_epsg = 32632)` |
| You're unsure what system your data is in       | 🤷 Unknown                   | Ask the data provider or check in GIS software           | —                                                                |
| You want to override auto-detect UTM            | Lat/lon or Projected         | Manually set `target_epsg` to force your own zone        | `convert_to_sf_utm(df, target_epsg = 32633)`                    |

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


---

🧪 Created for internal spatial analysis. Feel free to fork or adapt!
