# 📍 Standard Deviational Ellipse (SDE) Tool: SDEtool

[![](https://img.shields.io/badge/version-v1.0.0-blue.svg)](https://github.com/parker-group/SDEtool/releases/tag/v1.0.0)

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/parker-group/SDEtool/blob/main/LICENSE)


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

## 📂 Table of Contents

- [Workflow Steps](#-workflow-steps)
- [Simulated Example 1: Latitude/Longitude Data](#-simulated-example-1-latitudelongitude-data)
- [Simulated Example 2: UTM Projected X/Y Data](#-simulated-example-2-utm-projected-xy-data)
- [Simulated Example 3: Latitude/Longitude Data with Weights](#-simulated-example-3-latitudelongitude-data-with-a-count-of-peoplesamples-from-each-location-to-be-used-as-a-weight)
- [Coordinate System Tips](#-coordinate-system-tips)
- [What This Calculates](#-what-this-calculates)
- [Validation](validation/SDE_validation.md)
- [Reference PDF](The%20Standard%20Deviational%20Ellipse%20%20An%20Updated%20Tool%20for%20Spatial%20Description.pdf)


---

## 🔧 Workflow Steps

### 1. Load the required functions

You must first load the R functions before any other step will work.

There are **three ways** to do this:

- 🖱️ **Option 1: Copy + paste directly into R**  
  You can literally scroll down to the `SDE_functions.r` script, copy the functions, and paste them into your R console.

- 💻 **Option 2: Source the file locally** (if you've cloned or downloaded this repo)  
  ```r
  source("SDE_functions.r")
  ```

- 🌐 **Option 3: Source directly from GitHub**  
  ```r
  source("https://raw.githubusercontent.com/parker-group/SDEtool/main/SDE_functions.r")
  ```

➡️ [**View the full `SDE_functions.r` script on GitHub**](https://github.com/parker-group/SDEtool/blob/main/SDE_functions.r)

---

### 2. Load your data

Your dataset must:
- Include **geographic coordinates**
  - Either longitude/latitude (in degrees), or
  - Projected X/Y values (e.g., UTM in meters)
- Optionally include a **grouping variable** (e.g., Region, Year, or group_var) if you want to compute SDEs for different subsets.

The tool will automatically try to detect latitude and longitude columns using common names. Specifically, it searches for:

```r
lat_candidates <- c("latitude", "Latitude", "lat", "Lat", "y", "Y")
lon_candidates <- c("longitude", "Longitude", "lon", "Lon", "x", "X")
```

If your coordinate columns don't match these, rename them before running the functions.

If you're grouping by a variable (e.g., region, year, or category), name it or them clearly — e.g., "Region" or "genus" (`group_vars = "Region"`). You can have multiple (`group_vars = c("Region", "genus")`). 

If your data have repeats per location, you can either run the tool with multiple rows having the same location - or - you could generate a dataset that has one row per location and a count of people or samples from each location. Make sure you clearly name that "count" variable as well, and you can use that as a `weight` in the SDE function (`weight_col = count`). 

Example:

```r
df <- data.frame(
  longitude = c(36.8, 36.9, 36.7),
  latitude = c(-1.3, -1.4, -1.2),
  group_var = c("A", "A", "B")
)
```

---

### 3. Convert to sf object and UTM projection

We need to convert to an `sf` (simple feature), which tells R that the data are spatial (and special too) — meaning that the coordinates represent geometry. The `sf` object includes CRS metadata. We'll also convert to UTM so that distances and areas are straightforward in their calculations (UTMs are metric).

Use `convert_to_sf_utm()` to convert the data to a spatial object and project it to UTM automatically:

```r
sf_pts_proj <- convert_to_sf_utm(df)
```

If your data are **already projected** (i.e., not in latitude/longitude), you must **explicitly specify both the input and target EPSG codes**:

```r
sf_pts_proj <- convert_to_sf_utm(df, input_crs = 32636, target_epsg = 32636)
```

> ⚠️ **Note:** The function cannot auto-detect the UTM zone from already projected coordinates. You must supply `target_epsg` when `input_crs` is not lat/lon.

---

### 4. Generate the SDEs

Use the main function to create ellipses for each group. Note that you can set the group vars to "NULL" if you want SDEs for all points in the data.

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

##plot a map of the data
library(ggplot2)
library(sf)

ggplot() +
  geom_sf(data = sde_sf, aes(fill = as.factor(sd_level)), alpha = 0.3, color = "black") +
  geom_sf(data = sf_pts_proj, color = "red", size = 1.5) +
  scale_fill_brewer(palette = "Set2", name = "SD Level") +
  theme_minimal() +
  labs(
    title = "Standard Deviational Ellipses",
    subtitle = "With Input Points Overlaid",
    x = "Easting (m)", y = "Northing (m)"
  )
```

---

### 6. Export to shapefile (optional)

If desired, export your SDEs to a shapefile:

```r
#this will export a .shp file in the UTMs we are using here. You will need to update the folder location
sf::st_write(sde_sf, "SDE_ellipses.shp", delete_dsn = TRUE)

#this will convert back to WGS84 (normal lat/lon) before you export the .shp file. You will need to update the folder location
sde_wgs84 <- sf::st_transform(sde_sf, crs = 4326)
sf::st_write(sde_wgs84, "SDE_ellipses_WGS84.shp", delete_dsn = TRUE)
```

---

## 🧪 Simulated Example 1: Latitude/Longitude Data

```r
set.seed(123)
n <- 100
df1 <- rbind(
  data.frame(
    Longitude = runif(n / 2, min = 36.6, max = 36.7),
    Latitude = runif(n / 2, min = -1.5, max = -1.4),
    Region = "East"
  ),
  data.frame(
    Longitude = runif(n / 2, min = 36.9, max = 37.0),
    Latitude = runif(n / 2, min = -1.1, max = -1.0),
    Region = "West"
  )
)

#what do my data look like?
head(df1)

#now use the SDE tools
sf_pts_proj <- convert_to_sf_utm(df1)
sde_sf <- generate_sde_ellipses(sf_pts_proj, group_vars = "Region")
print(sde_sf)

##generate a simple map in R
# Load packages
library(ggplot2)
library(sf)

# Plot
ggplot() +
  geom_sf(data = sde_sf, aes(fill = as.factor(sd_level)), alpha = 0.3, color = "black") +
  geom_sf(data = sf_pts_proj, aes(color = Region), size = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "SD Level") +
  scale_color_brewer(palette = "Dark2", name = "Region") +
  theme_minimal() +
  labs(
    title = "Simulated Ellipses from Latitude/Longitude Data",
    subtitle = "Projected to UTM automatically",
    x = "Easting (m)", y = "Northing (m)"
  )
```

---

## 🧪 Simulated Example 2: UTM Projected X/Y Data

```r
set.seed(456)
n <- 100
df2 <- rbind(
  data.frame(
    X = rnorm(n / 2, mean = 490000, sd = 5000),
    Y = rnorm(n / 2, mean = 980000, sd = 5000),
    Region = "North"
  ),
  data.frame(
    X = rnorm(n / 2, mean = 510000, sd = 5000),
    Y = rnorm(n / 2, mean = 1000000, sd = 5000),
    Region = "South"
  )
)

#what do my data look like?
head(df2)

#now use the SDE tools
sf_pts_proj <- convert_to_sf_utm(df2, input_crs = 32636, target_epsg = 32636)
sde_sf <- generate_sde_ellipses(sf_pts_proj, group_vars = "Region")
print(sde_sf)

#generate simple map in R
# Load packages
library(ggplot2)
library(sf)
# Plot
ggplot() +
  geom_sf(data = sde_sf, aes(fill = as.factor(sd_level)), alpha = 0.3, color = "black") +
  geom_sf(data = sf_pts_proj, aes(color = Region), size = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "SD Level") +
  scale_color_brewer(palette = "Dark2", name = "Region") +
  theme_minimal() +
  labs(
    title = "Simulated Ellipses from Projected UTM Coordinates",
    subtitle = "EPSG:32636",
    x = "Easting (m)", y = "Northing (m)"
  )

```

---

## 🧪 Simulated Example 3: Latitude/Longitude Data with a count of people/samples from each location, to be used as a weight
```r
# Simulate spatial points in two regions with distinct geographic clusters
set.seed(789)
n <- 100

# Group "East"
east <- data.frame(
  Longitude = runif(n / 2, min = 36.6, max = 36.75),
  Latitude = runif(n / 2, min = -1.5, max = -1.35),
  Region = "East",
  count = sample(1:10, n / 2, replace = TRUE)
)

# Group "West" farther away
west <- data.frame(
  Longitude = runif(n / 2, min = 36.85, max = 37.0),
  Latitude = runif(n / 2, min = -1.2, max = -1.05),
  Region = "West",
  count = sample(1:10, n / 2, replace = TRUE)
)

# Combine the two
df3 <- rbind(east, west)

# what do my data look like?
head(df3)

# Convert to spatial object with UTM projection
sf_pts_proj <- convert_to_sf_utm(df3)

# Generate weighted SDEs
sde_sf <- generate_sde_ellipses(
  sf_pts_proj,
  group_vars = "Region",
  sd_levels = c(1, 2, 3),
  min_points = 5,
  sqrt2_scaling = TRUE,
  dof_correction = TRUE,
  weight_col = "count"
)

# Inspect output
print(sde_sf)

# Plot SDEs and weighted points
library(ggplot2)
library(sf)

ggplot() +
  geom_sf(data = sde_sf, aes(fill = as.factor(sd_level)),
          alpha = 0.3, color = "black") +
  geom_sf(data = sf_pts_proj, aes(size = count, color = Region),
          alpha = 0.8) +
  scale_size_continuous(range = c(1, 6)) +
  scale_fill_brewer(palette = "Set2", name = "SD Level") +
  scale_color_brewer(palette = "Dark2", name = "Region") +
  theme_minimal() +
  labs(
    title = "Weighted SDEs by Region",
    subtitle = "Point size represents 'count' used as weight",
    x = "Easting", y = "Northing"
  )


```
---

## 🛍 Coordinate System Tips

| Your Data Look Like…                           | Coordinate Type              | What You Should Do                                       | Example Call                                                    |
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

### 💡 Motivation

There are several other tools that support Standard Deviational Ellipse (SDE) calculations, including:

- 🧭 **[CrimeStat](https://nij.ojp.gov/library/publications/crimestat-40-user-manual)** – a comprehensive spatial analysis program used widely in criminology and public safety  
- 🛰️ **[ArcGIS SDE Tool](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-directional-distribution-standard-deviationa.htm)** – a built-in function in ArcGIS for directional distribution analysis  
- 🧩 **[QGIS SDE Plugin](https://github.com/havatv/qgisstandarddeviationalellipseplugin)** – a community-developed plugin for generating ellipses in QGIS  

These tools are well established and widely used within desktop GIS environments.

This R-based tool was created to provide a fully **open-source**, **script-based**, and **reproducible** alternative tailored for R users. It integrates smoothly into analytical workflows, supports **weighted points**, **custom grouping**, and offers flexible control over ellipse generation — making it especially suitable for transparent and research-grade spatial analyses.  

I developed this tool after encountering limitations with plugin-based approaches and needing a workflow that could be easily validated, customized, and shared in R. Feel free to fork or adapt!

---
