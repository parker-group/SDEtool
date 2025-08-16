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
convert_to_sf_utm <- function(df, input_crs = 4326, target_epsg = NULL) {
  coords <- detect_latlon(df)

  # Create sf object with defined input CRS
  sf_obj <- st_as_sf(df, coords = c(coords$lon, coords$lat), crs = input_crs, remove = FALSE)

  # Decide target EPSG
  if (is.null(target_epsg)) {
    if (input_crs == 4326) {
      target_epsg <- auto_utm(df[[coords$lon]], df[[coords$lat]])
    } else {
      stop("❌ Cannot auto-detect UTM from non-lat/lon input. Please specify `target_epsg` explicitly.")
    }
  }

  # Message if input and target CRS are the same
  if (input_crs == target_epsg) {
    message("🤔 Are you nuts? Input and target CRS are the same — skipping transformation.")
    return(sf_obj)
  }

  # Otherwise, transform to target CRS
  st_transform(sf_obj, target_epsg)
}

# ---------- NEW HELPERS ----------
sde_scale_factor <- function(levels, mode = "arcgis", coverage = c(0.6827, 0.95, 0.9973)) {
  if (mode %in% c("arcgis","crimestat")) return(levels * sqrt(2))
  if (mode == "qgis")                    return(rep(1, length(levels)))
  if (mode == "prob")                    return(sqrt(stats::qchisq(coverage, df = 2)))
  stop("Unknown mode in sde_scale_factor().")
}

df_used_text <- function(mode) {
  if (mode %in% c("arcgis","qgis")) return("n")       # population
  if (mode == "crimestat")          return("n-2")
  if (mode == "prob")               return("n-1")     # conventional sample
  "n"
}

# multiply an (n-1) sample covariance to the desired divisor
cov_rescale_for_df <- function(cov_mat, n, mode) {
  if (!is.finite(n) || n < 2) return(cov_mat)
  txt <- df_used_text(mode)
  if (txt == "n")    return(cov_mat * ((n - 1) / n))
  if (txt == "n-2")  return(cov_mat * ((n - 1) / (n - 2)))
  if (txt == "n-1")  return(cov_mat)  # already sample cov
  cov_mat
}

angle_eastccw_to_northcw <- function(deg) (90 - deg) %% 360


#########################################################
##### now we move towards generating the SDEs    ########
#########################################################

# Build an ellipse polygon
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
generate_sde_ellipses <- function(
  sf_data,
  group_vars = c("Location", "org1_genus"),
  sd_levels = c(1, 2, 3),
  min_points = 5,
  # legacy params kept; ignored when mode is used
  sqrt2_scaling = TRUE,
  dof_correction = TRUE,
  weight_col = NULL,
  # NEW:
  mode = "arcgis",                  # "arcgis","crimestat","qgis","prob"
  compute_in = "input",             # "input","working"
  working_crs = NULL,               # e.g., 32648 or "auto_utm"
  output_crs = "input",             # "input","working"
  return_metric = FALSE,
  coverage = c(0.6827, 0.95, 0.9973)  # used in mode="prob"
) {
  stopifnot(inherits(sf_data, "sf"))
  if (!all(group_vars %in% names(sf_data))) {
    stop("❌ group_vars not found in data.")
  }

  # Choose the working layer
  work <- sf_data
  if (compute_in == "working") {
    if (is.null(working_crs) || identical(working_crs, "auto_utm")) {
      # reuse your helper to infer UTM when input is WGS84
      coords <- st_coordinates(sf_data)
      # find lon/lat columns (from sf geometry)
      lon <- coords[,1]; lat <- coords[,2]
      epsg <- auto_utm(lon, lat)
      working_crs <- epsg
    }
    work <- st_transform(sf_data, working_crs)
  }

  # After computation, decide which CRS to return as the active geometry
  to_return_crs <- if (output_crs == "working") st_crs(work) else st_crs(sf_data)

  coords <- st_coordinates(work)
  df <- work |> st_drop_geometry()
  df$x <- coords[,1]; df$y <- coords[,2]

  groups <- df |> dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |> dplyr::group_split()
  results <- list(); skipped <- list()

  for (g in groups) {
    key_vals <- sapply(group_vars, function(v) as.character(g[[v]][1]))
    key <- paste(key_vals, collapse = " | ")
    n <- nrow(g)
    if (n < min_points) { skipped[[key]] <- paste0("Skipped (", n, " points)"); next }

    x <- g$x; y <- g$y
    w <- if (!is.null(weight_col) && weight_col %in% names(g)) g[[weight_col]] else rep(1, n)

    # center & covariance
    center <- c(mean(x), mean(y))
    cov_samp <- stats::cov(cbind(x, y))           # (n-1) by default
    cov_use  <- cov_rescale_for_df(cov_samp, n, mode)

    eig <- eigen(cov_use)
    # ensure descending order
    ord <- order(eig$values, decreasing = TRUE)
    eig$values <- eig$values[ord]; eig$vectors <- eig$vectors[, ord, drop = FALSE]

    # angles
    orient_rad_east <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
    orient_deg_east <- (orient_rad_east * 180 / pi) %% 360
    orient_deg_north <- angle_eastccw_to_northcw(orient_deg_east)

    # scale factors for this mode
    if (mode == "prob") {
      # ignore sd_levels; replicate by coverage length
      levels_vec <- rep(1, length(coverage))
      scales <- sde_scale_factor(levels_vec, mode = "prob", coverage = coverage)
      levels_out <- coverage
      scale_rule <- "sqrt(qchisq)"
      implied_cov <- coverage
    } else {
      levels_vec <- sd_levels
      scales <- sde_scale_factor(levels_vec, mode = mode, coverage = coverage)
      levels_out <- sd_levels
      scale_rule <- "k*sqrt(2)"
      # implied coverage for k·√2 in 2D = 1 - exp(-k^2)
      implied_cov <- 1 - exp(-(sd_levels^2))
    }

    for (i in seq_along(levels_vec)) {
      scl <- scales[i]
      # semi-axes
      major <- scl * sqrt(eig$values[1])
      minor <- scl * sqrt(eig$values[2])

      # build ellipse polygon in the working CRS
      t <- seq(0, 2*pi, length.out = 200)
      circ <- rbind(cos(t)*major, sin(t)*minor)
      R <- matrix(c(cos(orient_rad_east), -sin(orient_rad_east),
                    sin(orient_rad_east),  cos(orient_rad_east)), nrow = 2, byrow = TRUE)
      shape <- R %*% circ
      coords_poly <- sweep(t(shape), 2, center, "+")
      coords_poly <- rbind(coords_poly, coords_poly[1, , drop = FALSE])
      sfc_poly <- st_sfc(st_polygon(list(coords_poly)), crs = st_crs(work))

      # reproject active geometry if needed
      geom_main <- if (output_crs == "working") sfc_poly else st_transform(sfc_poly, st_crs(sf_data))

      # optional metric geometry (for measuring in meters even if output in input CRS)
      geom_metric <- if (return_metric && output_crs == "input" && compute_in == "working") sfc_poly else NULL

      # area (of active geometry)
      area_val <- as.numeric(st_area(geom_main))

      # percent inside using active geometry CRS
      pts_sf <- st_as_sf(g, coords = c("x","y"), crs = st_crs(work))
      pts_for_within <- if (output_crs == "working") pts_sf else st_transform(pts_sf, st_crs(sf_data))
      inside <- st_within(pts_for_within, geom_main, sparse = FALSE)[,1]
      count_inside <- sum(w[inside]); percent_inside <- round(100 * count_inside / sum(w), 1)

      # row (keep your original fields; add new)
      row <- dplyr::tibble(
        sd_level = if (mode == "prob") NA_real_ else levels_out[i],
        target_coverage = if (mode == "prob") levels_out[i] else NA_real_,
        n_points = n,
        count_inside = count_inside,
        percent_inside = percent_inside,
        area = area_val,
        major_axis = major,
        minor_axis = minor,
        orientation_rad = orient_rad_east,
        orientation_deg = orient_deg_east,              # east_ccw (kept)
        angle_north_cw = orient_deg_north,              # new
        sde_mode = mode,
        df_used = df_used_text(mode),
        scale_rule = scale_rule,
        implied_coverage = if (mode == "prob") levels_out[i] else implied_cov[i],
        angle_basis = if (mode %in% c("arcgis","crimestat")) "north_cw" else "east_ccw",
        center_x = center[1], center_y = center[2],
        lambda1 = eig$values[1], lambda2 = eig$values[2]
      )
      for (v in group_vars) row[[v]] <- g[[v]][1]
      row$geometry <- geom_main
      if (!is.null(geom_metric)) row$geom_metric <- geom_metric
      results[[length(results) + 1]] <- st_as_sf(row)
    }
  }

  if (length(skipped) > 0) {
    cat("ℹ️ Skipped groups:\n"); print(skipped)
  }

  out <- dplyr::bind_rows(results)
  st_crs(out) <- to_return_crs
  out
}

