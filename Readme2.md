# 📍 Standard Deviational Ellipse (SDE) Tool

# This R tool computes Standard Deviational Ellipses (SDEs) for spatial point data grouped by user-defined variables.
# It supports:
# - Multiple standard deviation levels (e.g., 1, 2, 3 SD)
# - Weighted points (optional)
# - Yuill + √2 correction (default)
# - Degrees of freedom correction (default)
# - Summary of ellipse shape + % of points enclosed

# 📚 Reference: Yuill, R. S. (1971). The Standard Deviational Ellipse: An Updated Tool for Spatial Description.
# Geografiska Annaler: Series B, Human Geography, 53(1), 28–39. https://doi.org/10.2307/490885

# 📖 ArcGIS documentation: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-directional-distribution-standard-deviationa.htm

# --------------------------------------------------------------------------------
# 🔧 Step-by-Step Instructions

# 1. Load required functions (must be done before anything else!)
# Option 1: If running locally after cloning this repo
source("SDE_functions.r")

# Option 2: Run directly from GitHub (raw link)
source("https://raw.githubusercontent.com/parker-group/SDEtool/main/SDE_functions.r")

# --------------------------------------------------------------------------------
# 🛍 Coordinate System Tips

# | Your Data Looks Like…                          | Coordinate Type              | What You Should Do                                       | Example Call                                                    |
# |------------------------------------------------|------------------------------|----------------------------------------------------------|-----------------------------------------------------------------|
# | Values like -1.3, 36.8                         | Latitude/Longitude (degrees) | Nothing special — default settings will work             | convert_to_sf_utm(df)                                           |
# | GPS data from phone/app                        | Latitude/Longitude (degrees) | Default is fine — UTM zone will be auto-detected         | convert_to_sf_utm(my_data)                                     |
# | X/Y values like 500000, 1000000 (meters)       | Projected (e.g., UTM)        | You must specify the CRS (EPSG code)                     | convert_to_sf_utm(df, input_crs = 32632, target_epsg = 32632)   |
# | You’re unsure what system your data is in      | 🤷 Unknown                   | Ask the data provider or check in GIS software           | —                                                               |
# | You want to override auto-detect UTM           | Lat/lon or Projected         | Manually set target_epsg to force your own zone          | convert_to_sf_utm(df, target_epsg = 32633)                      |

# 💡 How to Find EPSG Codes:
# - Visit https://epsg.io
# - Use sf::st_crs() on known spatial data
# - In QGIS: Right-click layer → Properties → CRS

# --------------------------------------------------------------------------------
# 🧪 Simulated Use Case 1: Latitude / Longitude Input

# Simulate points with lat/lon
set.seed(123)
latlon_df <- data.frame(
  lon = rnorm(200, mean = 36.8, sd = 0.1),
  lat = rnorm(200, mean = -1.3, sd = 0.1),
  Group = rep(c("ZoneA", "ZoneB"), each = 100)
)

# Convert to sf + auto-detect UTM
sf_latlon <- convert_to_sf_utm(latlon_df)

# Run SDE function
sde1 <- generate_sde_ellipses(
  sf_latlon,
  group_vars = "Group",
  sd_levels = c(1, 2, 3),
  min_points = 5
)

# Plot
library(ggplot2)
library(sf)

ggplot() +
  geom_sf(data = sde1, aes(fill = as.factor(sd_level)), alpha = 0.3, color = "black") +
  geom_sf(data = sf_latlon, aes(color = Group), size = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "SD Level") +
  scale_color_brewer(palette = "Dark2", name = "Group") +
  theme_minimal() +
  labs(title = "Lat/Lon Input Example", subtitle = "Auto-detected UTM projection")

# --------------------------------------------------------------------------------
# 🧪 Simulated Use Case 2: UTM Input

# Simulate projected points
utm_df <- data.frame(
  X = rnorm(150, mean = 500000, sd = 800),
  Y = rnorm(150, mean = 9850000, sd = 1200),
  Type = rep(c("West", "East", "Central"), each = 50)
)

# Convert to sf using known EPSG (UTM Zone 36N = EPSG:32636)
sf_utm <- convert_to_sf_utm(utm_df, input_crs = 32636, target_epsg = 32636)

# Run SDE
sde2 <- generate_sde_ellipses(
  sf_utm,
  group_vars = "Type",
  sd_levels = c(1, 2),
  sqrt2_scaling = TRUE,
  dof_correction = TRUE
)

# Plot
ggplot() +
  geom_sf(data = sde2, aes(fill = as.factor(sd_level)), alpha = 0.3, color = "black") +
  geom_sf(data = sf_utm, aes(color = Type), size = 1.2) +
  scale_fill_brewer(palette = "Pastel2", name = "SD Level") +
  scale_color_brewer(palette = "Set1", name = "Type") +
  theme_minimal() +
  labs(title = "UTM Input Example", subtitle = "User-defined projection (EPSG 32636)")

# --------------------------------------------------------------------------------
# 🔬 What This Calculates

# The Standard Deviational Ellipse (SDE) summarizes spatial distribution of points by showing:
# - Directional trend (via ellipse orientation)
# - Dispersion along major and minor axes
# - % of points enclosed at each SD level (~63% at 1 SD, ~98% at 2 SD, ~99.9% at 3 SD)

# The ellipse is calculated using the eigenvectors of the covariance matrix of X and Y.

# --------------------------------------------------------------------------------
# ✅ To Do
# - [ ] Add Shiny app version
# - [ ] Integrate shapefile upload option
# - [ ] Validate against ArcGIS or QGIS SDE tools

# --------------------------------------------------------------------------------
# 🧪 Created for internal spatial analysis. Feel free to fork or adapt!
