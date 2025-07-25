################################################################################
####### step one of the process is to run these functions. #####################
#######   then either use the simulation (in readme) or use your own data ######
################################################################################


# Load required libraries
library(sf)
library(dplyr)
library(purrr)

# ----------------------------
# Helper Functions
# ----------------------------

# Detect latitude and longitude columns
# I've written in a few different ways you could label your lat/lon
# if you've not used any of these labels, you'll likely get an error message here
detect_latlon <- function(df) {
  lat_candidates <- c("latitude", "Latitude", "lat", "Lat", "y", "Y")
  lon_candidates <- c("longitude", "Longitude", "lon", "Lon", "x", "X")
  lat_col <- intersect(names(df), lat_candidates)
  lon_col <- intersect(names(df), lon_candidates)
  if (length(lat_col) == 0 || length(lon_col) == 0) {
    stop("❌ Latitude or longitude column not detected.")
  }
  list(lat = lat_col[1], lon = lon_col[1])
}

# Auto-detect UTM zone from coordinates
## we want to move to UTMs so that we can later have some straightforward 
##   calculations on distances, area, etc. 
auto_utm <- function(lon, lat) {
  zone <- floor((mean(lon, na.rm = TRUE) + 180) / 6) + 1
  hemisphere <- ifelse(mean(lat, na.rm = TRUE) >= 0, "north", "south")
  epsg <- ifelse(hemisphere == "north", 32600, 32700) + zone
  return(epsg)
}

# Convert to sf and project to detected UTM zone
## now that we know what the UTM zone should be, we project the lat/lon to that
convert_to_sf_utm <- function(df) {
  coords <- detect_latlon(df)
  sf_obj <- st_as_sf(df, coords = c(coords$lon, coords$lat), crs = 4326, remove = FALSE)
  epsg <- auto_utm(df[[coords$lon]], df[[coords$lat]])
  st_transform(sf_obj, epsg)
}

# Build an ellipse polygon
## now we move towards generating the SDEs
build_ellipse <- function(x, y, sd = 1, n_points = 100, sqrt2_scaling = TRUE) {
  if (length(x) < 2 || length(y) < 2 || anyNA(x) || anyNA(y)) return(NULL)
  cov_mat <- stats::cov(cbind(x, y))
  if (any(!is.finite(cov_mat)) || any(diag(cov_mat) == 0)) return(NULL)
  eig <- eigen(cov_mat)
  center <- c(mean(x), mean(y))
  circle <- rbind(cos(seq(0, 2 * pi, length.out = n_points)), sin(seq(0, 2 * pi, length.out = n_points)))
  axes <- diag(sqrt(eig$values))
  scale_factor <- if (sqrt2_scaling) sqrt(2) else 1
  shape <- scale_factor * sd * eig$vectors %*% axes %*% circle
  coords <- sweep(t(shape), 2, center, "+")
  coords <- rbind(coords, coords[1, , drop = FALSE])
  st_polygon(list(coords))
}

# Main ellipse generator
## note that we're automatically generating 1, 2, and 3 standard deviations
## I've set the minimum number of points to 5, less than that seems silly for stats
#### note that you could also use weighted data here, if you want: weight_col = ??
generate_sde_ellipses <- function(sf_data,
                                   group_vars = c("Location", "org1_genus"),
                                   sd_levels = c(1, 2, 3),
                                   min_points = 5,
                                   sqrt2_scaling = TRUE,
                                   dof_correction = TRUE,
                                   weight_col = NULL) {

  if (!inherits(sf_data, "sf")) stop("❌ Input must be an sf object.")
  coords <- st_coordinates(sf_data)
  df <- sf_data %>% st_drop_geometry()
  df$x <- coords[, 1]
  df$y <- coords[, 2]
  group_data <- df %>% group_by(across(all_of(group_vars))) %>% group_split()
  results <- list()
  skipped <- list()

  for (g in group_data) {
    key <- paste(g[[group_vars[1]]][1], g[[group_vars[2]]][1], sep = " | ")
    n <- nrow(g)
    if (n < min_points) {
      skipped[[key]] <- paste0("Skipped (", n, " points)")
      next
    }
    x <- g$x
    y <- g$y
    w <- if (!is.null(weight_col) && weight_col %in% names(g)) g[[weight_col]] else rep(1, n)

    for (sd in sd_levels) {
      poly <- try(build_ellipse(x, y, sd = sd, sqrt2_scaling = sqrt2_scaling), silent = TRUE)
      if (inherits(poly, "try-error") || is.null(poly)) next
      sfc_poly <- st_sfc(poly, crs = st_crs(sf_data))
      area <- as.numeric(st_area(sfc_poly))

      # Count points inside
      inside_logical <- st_within(st_as_sf(g, coords = c("x", "y"), crs = st_crs(sf_data)), sfc_poly, sparse = FALSE)[,1]
      count_inside <- sum(w[inside_logical])
	  percent_inside <- round(100 * count_inside / sum(w), 1)
	  
      # Ellipse stats
      cov_mat <- stats::cov(cbind(x, y))
      eig <- eigen(cov_mat)
      scale_factor <- if (sqrt2_scaling) sqrt(2) else 1
      major <- sd * scale_factor * sqrt(eig$values[1])
      minor <- sd * scale_factor * sqrt(eig$values[2])
      orient_rad <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
      orient_deg <- (orient_rad * 180 / pi) %% 180

      result_row <- tibble(
        sd_level = sd,
        n_points = n,
        count_inside = count_inside,
		percent_inside = percent_inside,
		area_m2 = area,
        major_axis = major,
        minor_axis = minor,
        orientation_rad = orient_rad,
        orientation_deg = orient_deg
      )
      for (v in group_vars) result_row[[v]] <- g[[v]][1]
      result_row$geometry <- sfc_poly
      results[[length(results) + 1]] <- st_as_sf(result_row)
    }
  }

  if (length(skipped) > 0) {
    cat("ℹ️ Skipped groups:\n")
    print(skipped)
  }

  bind_rows(results)
}

