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

#' Detect latitude and longitude columns in a data.frame
#'
#' Tries a small set of common names for latitude and longitude.
#' If none are found, stops with a clear error.
#'
#' @param df A data.frame that contains latitude/longitude columns.
#' @return A list with elements `lat` and `lon` giving the detected column names.
#' @examples
#' detect_latlon(data.frame(lat = 1:3, lon = 1:3))
#' @keywords internal
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

#' Choose a WGS84 UTM EPSG code from lon/lat vectors
#'
#' Selects UTM zone from the **mean** lon/lat and returns the appropriate
#' WGS84 UTM EPSG (326xx north / 327xx south).
#'
#' @param lon Numeric vector of longitudes (degrees).
#' @param lat Numeric vector of latitudes (degrees).
#' @return Integer EPSG code for the chosen UTM zone.
#' @examples
#' auto_utm(c(-118,-117), c(34,35))
#' @keywords internal
auto_utm <- function(lon, lat) {
  zone <- floor((mean(lon, na.rm = TRUE) + 180) / 6) + 1
  hemisphere <- ifelse(mean(lat, na.rm = TRUE) >= 0, "north", "south")
  epsg <- ifelse(hemisphere == "north", 32600, 32700) + zone
  return(epsg)
}

#' Convert a data.frame (lon/lat) to sf and project to UTM
#'
#' Creates an `sf` object from detected lon/lat columns (or leaves coordinates
#' in place with `remove = FALSE`) and transforms to a UTM CRS. If `target_epsg`
#' is `NULL` and `input_crs` is 4326, `auto_utm()` is used.
#'
#' @param df A data.frame with lon/lat columns, or an `sf` object.
#' @param input_crs Input CRS EPSG (default 4326).
#' @param target_epsg Target EPSG (optional). If `NULL` with lon/lat input, auto-select UTM.
#' @return An `sf` object in the target CRS.
#' @examples
#' # convert_to_sf_utm(my_df_with_lon_lat)
#' @export
convert_to_sf_utm <- function(df, input_crs = 4326, target_epsg = NULL) {
  coords <- detect_latlon(df)
  sf_obj <- sf::st_as_sf(df, coords = c(coords$lon, coords$lat), crs = input_crs, remove = FALSE)
  if (is.null(target_epsg)) {
    if (input_crs == 4326) {
      target_epsg <- auto_utm(df[[coords$lon]], df[[coords$lat]])
    } else {
      stop("❌ Cannot auto-detect UTM from non-lat/lon input. Please specify `target_epsg` explicitly.")
    }
  }
  if (input_crs == target_epsg) {
    message("🤔 Are you nuts? Input and target CRS are the same — skipping transformation.")
    return(sf_obj)
  }
  sf::st_transform(sf_obj, target_epsg)
}

#' Scale factor for SDE semi-axes by mode
#'
#' Implements the preset scale rules used by desktop tools and a
#' probabilistic option:
#' - `arcgis` / `crimestat`: \(k \cdot \sqrt{2}\)
#' - `qgis`: 1 (placeholder)
#' - `prob`: \(\sqrt{\chi^2_2(p)}\) from coverage `p`
#'
#' @param levels Numeric SD levels (e.g., `c(1,2,3)`) for non-prob modes.
#' @param mode One of `"arcgis"`, `"crimestat"`, `"qgis"`, `"prob"`.
#' @param coverage Probabilities used when `mode="prob"`.
#' @return Numeric vector of scale factors.
#' @examples
#' sde_scale_factor(c(1,2,3), mode="arcgis")
#' sde_scale_factor(levels = 1, mode="prob", coverage=c(0.95))
#' @keywords internal
sde_scale_factor <- function(levels, mode = "arcgis", coverage = c(0.6827, 0.95, 0.9973)) {
  if (mode %in% c("arcgis","crimestat")) return(levels * sqrt(2))
  if (mode == "qgis")                    return(rep(1, length(levels)))
  if (mode == "prob")                    return(sqrt(stats::qchisq(coverage, df = 2)))
  stop("Unknown mode in sde_scale_factor().")
}

#' Text description of df used by preset
#'
#' Returns the divisor choice implied by a preset. Used downstream to
#' rescale a sample covariance from `(n-1)` to the desired divisor.
#'
#' @param mode One of `"arcgis"`, `"qgis"`, `"crimestat"`, `"prob"`.
#' @return Character scalar: `"n"`, `"n-2"`, or `"n-1"`.
#' @examples
#' df_used_text("crimestat")
#' @keywords internal
df_used_text <- function(mode) {
  if (mode %in% c("arcgis","qgis")) return("n")       # population
  if (mode == "crimestat")          return("n-2")
  if (mode == "prob")               return("n-1")     # conventional sample
  "n"
}

#' Rescale a sample covariance (n-1) to an alternate divisor
#'
#' Accepts a sample covariance matrix (computed with `(n-1)` denominator)
#' and rescales it to behave as if it were computed with divisor `n`, `n-2`,
#' or left as `n-1` depending on preset.
#'
#' @param cov_mat 2x2 sample covariance (from `stats::cov`).
#' @param n Sample size used for the covariance.
#' @param mode Preset mode (see `df_used_text()`).
#' @return 2x2 covariance on the requested divisor.
#' @examples
#' cov_rescale_for_df(diag(2), n = 10, mode = "arcgis")
#' @keywords internal
cov_rescale_for_df <- function(cov_mat, n, mode) {
  if (!is.finite(n) || n < 2) return(cov_mat)
  txt <- df_used_text(mode)
  if (txt == "n")    return(cov_mat * ((n - 1) / n))
  if (txt == "n-2")  return(cov_mat * ((n - 1) / (n - 2)))
  if (txt == "n-1")  return(cov_mat)  # already sample cov
  cov_mat
}

#' Convert east-CCW (math) to north-CW (GIS UI) degrees
#'
#' @param deg Numeric vector of angles (degrees) measured east and counter-clockwise.
#' @return Numeric vector of angles (degrees) measured clockwise from north.
#' @examples
#' angle_eastccw_to_northcw(0)   # 90
#' angle_eastccw_to_northcw(90)  # 0
#' @keywords internal
angle_eastccw_to_northcw <- function(deg) (90 - deg) %% 360

#' Build a valid ellipse polygon from sample points
#'
#' Computes a covariance from `x,y`, takes its eigenvectors/values, and
#' returns a single-ring polygon. The ring is closed exactly once and
#' duplicate vertices are removed to avoid GEOS/s2 issues.
#'
#' @param x Numeric vector of x coordinates.
#' @param y Numeric vector of y coordinates.
#' @param sd Scale to apply to semi-axes (already includes any preset-specific factor).
#' @param n_points Number of vertices on the ellipse boundary (default 100).
#' @param sqrt2_scaling Back-compat flag; if `TRUE` multiplies `sd` by `sqrt(2)`.
#' @return An `sfg` polygon.
#' @details
#' This helper is convenient when you want to go from raw points directly
#' to a ring. In the main generator we typically compute eigenpairs once
#' and construct rings manually for multiple scales.
#' @examples
#' # build_ellipse(x, y, sd = 1)
#' @keywords internal
build_ellipse <- function(x, y, sd = 1, n_points = 100, sqrt2_scaling = TRUE) {
  if (length(x) < 2 || length(y) < 2 || anyNA(x) || anyNA(y)) return(NULL)
  cov_mat <- stats::cov(cbind(x, y))
  if (any(!is.finite(cov_mat)) || any(diag(cov_mat) == 0)) return(NULL)
  eig <- eigen(cov_mat)
  center <- c(mean(x), mean(y))

  # angle sequence WITHOUT 2*pi to avoid duplicating the first point
  t <- seq(0, 2 * pi, length.out = n_points + 1)
  t <- t[-length(t)]
  circle <- rbind(cos(t), sin(t))
  axes <- diag(sqrt(eig$values))
  scale_factor <- if (sqrt2_scaling) sqrt(2) else 1
  shape <- scale_factor * sd * eig$vectors %*% axes %*% circle

  coords <- sweep(t(shape), 2, center, "+")
  coords <- rbind(coords, coords[1, , drop = FALSE]) # close once

  # remove accidental duplicates before closing (numeric tolerance)
  ring_no_close <- coords[-nrow(coords), , drop = FALSE]
  ring_no_dup   <- ring_no_close[!duplicated(round(ring_no_close, 12)), , drop = FALSE]
  coords        <- rbind(ring_no_dup, ring_no_dup[1, , drop = FALSE])

  sf::st_polygon(list(coords))
}


#' Kish effective sample size for nonnegative weights
#' @param w Numeric weights (NAs allowed; negatives coerced to 0).
#' @return Scalar n_eff in [0, length(w)].
#' @keywords internal
kish_neff <- function(w) {
  stopifnot(is.numeric(w))
  w <- w[is.finite(w)]
  if (!length(w)) return(0)
  w[w < 0] <- 0
  s1 <- sum(w); s2 <- sum(w * w)
  if (s1 == 0 || s2 == 0) return(0)
  min((s1 * s1) / s2, length(w))
}

#' Weighted column means
#' @keywords internal
wmean_mat <- function(X, w) {
  X <- as.matrix(X)
  stopifnot(nrow(X) == length(w))
  w <- ifelse(is.finite(w) & w > 0, w, 0)
  ws <- sum(w)
  if (ws == 0) return(colMeans(X))
  w <- w / ws
  drop(colSums(X * w))
}

#' Weighted 2x2 scatter matrix (probability weights)
#' @keywords internal
wcov_2d <- function(X, w) {
  X <- as.matrix(X)
  stopifnot(ncol(X) == 2, nrow(X) == length(w))
  w <- ifelse(is.finite(w) & w > 0, w, 0)
  ws <- sum(w)
  if (ws == 0) return(stats::cov(X))
  w <- w / ws
  mu <- wmean_mat(X, w)
  Z  <- sweep(X, 2, mu, `-`)
  t(Z) %*% (Z * w)  # population-like (Σ w = 1)
}



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
# angle sequence WITHOUT 2*pi to avoid duplicating the first point
t <- seq(0, 2 * pi, length.out = n_points + 1)
t <- t[-length(t)]

circle <- rbind(cos(t), sin(t))
axes <- diag(sqrt(eig$values))
scale_factor <- if (sqrt2_scaling) sqrt(2) else 1
shape <- scale_factor * sd * eig$vectors %*% axes %*% circle

coords <- sweep(t(shape), 2, center, "+")
# close ring exactly once
coords <- rbind(coords, coords[1, , drop = FALSE])

# optional: remove any accidental duplicates (numerical) before closing
ring_no_close <- coords[-nrow(coords), , drop = FALSE]
ring_no_dup   <- ring_no_close[!duplicated(round(ring_no_close, 12)), , drop = FALSE]
coords        <- rbind(ring_no_dup, ring_no_dup[1, , drop = FALSE])

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

  # Allow no-group usage
  if (is.null(group_vars) || length(group_vars) == 0) {
    sf_data$.__grp__ <- 1L
    group_vars <- ".__grp__"
  }

  stopifnot(inherits(sf_data, "sf"))
  if (!all(group_vars %in% names(sf_data))) {
    stop("❌ group_vars not found in data.")
  }

  # Optional: validate flags
  compute_in <- match.arg(compute_in, c("input","working"))
  output_crs <- match.arg(output_crs, c("input","working"))

  # Choose the working layer
  work <- sf_data
  if (compute_in == "working") {
    if (is.null(working_crs) || identical(working_crs, "auto_utm")) {
      coords <- sf::st_coordinates(sf_data)
      lon <- coords[,1]; lat <- coords[,2]
      working_crs <- auto_utm(lon, lat)
    }
    work <- sf::st_transform(sf_data, working_crs)
  }

  # Decide which CRS the returned active geometry will use
  to_return_crs <- if (output_crs == "working") sf::st_crs(work) else sf::st_crs(sf_data)


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

# center & covariance (weighted if weight_col present and usable)
XY <- cbind(x, y)

if (!is.null(weight_col) && weight_col %in% names(g) &&
    any(is.finite(w)) && sum(w, na.rm = TRUE) > 0 && sum(w > 0, na.rm = TRUE) >= 2) {

  # 1) weighted center & population-like scatter (Σw normalized to 1)
  center <- wmean_mat(XY, w)
  S_pop  <- wcov_2d(XY, w)

  # 2) convert to a "sample-like" covariance on (n_eff - 1) so your rescaler works
  n_eff <- kish_neff(w)
  if (n_eff > 1) {
    cov_samp <- as.matrix(S_pop) * (n_eff / (n_eff - 1))
  } else {
    cov_samp <- stats::cov(XY)  # fallback
  }

  # 3) rescale to the divisor implied by the preset using n_eff
  cov_use <- cov_rescale_for_df(cov_samp, n = n_eff, mode = mode)

} else {
  # unweighted path (exactly your current behavior)
  center   <- c(mean(x), mean(y))
  cov_samp <- stats::cov(XY)                 # (n-1) divisor
  cov_use  <- cov_rescale_for_df(cov_samp, n = n, mode = mode)
}

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
  t <- seq(0, 2*pi, length.out = 200 + 1)
  t <- t[-length(t)]  # drop 2*pi to avoid duplicating the first point

  circ <- rbind(cos(t)*major, sin(t)*minor)
  R <- matrix(c(cos(orient_rad_east), -sin(orient_rad_east),
              sin(orient_rad_east),  cos(orient_rad_east)), nrow = 2, byrow = TRUE)
  shape <- R %*% circ
  coords_poly <- sweep(t(shape), 2, center, "+")

  # close ring exactly once
  coords_poly <- rbind(coords_poly, coords_poly[1, , drop = FALSE])

  # (optional) numeric de-dup to satisfy s2 on some platforms
  ring_nc   <- coords_poly[-nrow(coords_poly), , drop = FALSE]
  ring_nodu <- ring_nc[!duplicated(round(ring_nc, 12)), , drop = FALSE]
  coords_poly <- rbind(ring_nodu, ring_nodu[1, , drop = FALSE])

  sfc_poly <- sf::st_sfc(sf::st_polygon(list(coords_poly)), crs = sf::st_crs(work))

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



