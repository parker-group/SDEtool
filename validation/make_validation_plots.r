# --- SDEtool validation figures: ArcGIS, CrimeStat, Probabilistic ---------------
# Regenerates all figures to validation/figures/ AND shows them on screen.
# Assumes your SDE functions (e.g., generate_sde_ellipses) are available.

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
base_dir <- "C:/Users/salem/OneDrive/Documents/SDEtool/validation"
data_dir <- file.path(base_dir, "data")
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Helpers (robust lat/lon detection; bbox; plotting)
# ------------------------------------------------------------------------------
detect_latlon <- function(df) {
  lat_candidates <- c("latitude","Latitude","lat","Lat","y","Y")
  lon_candidates <- c("longitude","Longitude","lon","Lon","x","X")
  lat_col <- intersect(names(df), lat_candidates)
  lon_col <- intersect(names(df), lon_candidates)
  if (length(lat_col) == 0 || length(lon_col) == 0) {
    stop("Latitude/longitude columns not found.")
  }
  list(lat = lat_col[1], lon = lon_col[1])
}

tight_bbox <- function(..., pad = 0.05) {
  bbs <- lapply(list(...), sf::st_bbox)
  xmin <- min(sapply(bbs, `[[`, "xmin"), na.rm = TRUE)
  xmax <- max(sapply(bbs, `[[`, "xmax"), na.rm = TRUE)
  ymin <- min(sapply(bbs, `[[`, "ymin"), na.rm = TRUE)
  ymax <- max(sapply(bbs, `[[`, "ymax"), na.rm = TRUE)
  dx <- xmax - xmin; dy <- ymax - ymin
  c(xmin = xmin - pad*dx, xmax = xmax + pad*dx,
    ymin = ymin - pad*dy, ymax = ymax + pad*dy)
}

# R-tool vs reference (legend included)
plot_overlay_ref <- function(ref_sf, rt_sf, pts_sf, title, outfile,
                             pad = 0.05, legend_position = "top") {
  bb <- tight_bbox(ref_sf, rt_sf, pts_sf, pad = pad)
  gg <- ggplot() +
    geom_sf(data = rt_sf,  aes(color = "R-tool",     linetype = "R-tool"),
            fill = NA, linewidth = 1.1) +
    geom_sf(data = ref_sf, aes(color = "Reference",  linetype = "Reference"),
            fill = NA, linewidth = 1.1) +
    geom_sf(data = pts_sf, aes(shape = "Points"), color = "black", size = 1.4) +
    scale_color_manual(values = c("R-tool" = "#1f77b4", "Reference" = "#d62728"),
                       name = NULL) +
    scale_linetype_manual(values = c("R-tool" = "solid", "Reference" = "dashed"),
                          name = NULL) +
    scale_shape_manual(values = c("Points" = 16), name = NULL) +
    guides(color = guide_legend(order = 1, override.aes = list(fill = NA)),
           linetype = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]),
             expand = FALSE) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_blank(),
          legend.position = legend_position) +
    labs(title = title, x = "Longitude (°)", y = "Latitude (°)")
  # show + save
  if (Sys.getenv("RSTUDIO") == "" && .Platform$OS.type == "windows") windows(10, 7.5)
  print(gg)
  ggsave(outfile, gg, width = 10, height = 7.5, dpi = 300)
  message("Wrote: ", outfile)
}

# R-only comparison (ArcGIS vs CrimeStat vs Prob)
plot_compare <- function(combo_sf, pts_sf, title, outfile, pad = 0.05) {
  bb <- tight_bbox(combo_sf, pts_sf, pad = pad)
  gg <- ggplot() +
    geom_sf(data = combo_sf, aes(color = model, linetype = model),
            fill = NA, linewidth = 1.1) +
    geom_sf(data = pts_sf, color = "black", size = 1.4) +
    scale_color_manual(values = c("ArcGIS" = "#1f77b4",
                                  "CrimeStat" = "#d62728",
                                  "Prob (68.27%)" = "#2ca02c",
                                  "Prob (95%)"    = "#2ca02c"),
                       name = NULL) +
    scale_linetype_manual(values = c("ArcGIS" = "solid",
                                     "CrimeStat" = "dashed",
                                     "Prob (68.27%)" = "dotdash",
                                     "Prob (95%)"    = "dotdash"),
                          name = NULL) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]),
             expand = FALSE) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_blank(),
          legend.position = "right") +
    labs(title = title, x = "Longitude (°)", y = "Latitude (°)")
  if (Sys.getenv("RSTUDIO") == "" && .Platform$OS.type == "windows") windows(10, 7.5)
  print(gg)
  ggsave(outfile, gg, width = 10, height = 7.5, dpi = 300)
  message("Wrote: ", outfile)
}

# Probabilistic outlines only
plot_prob <- function(prob_sf, pts_sf, outfile, pad = 0.05) {
  prob_sf <- prob_sf |> arrange(target_coverage)
  bb <- tight_bbox(prob_sf, pts_sf, pad = pad)
  gg <- ggplot() +
    geom_sf(data = prob_sf, aes(color = factor(target_coverage)),
            fill = NA, linewidth = 1.2) +
    geom_sf(data = pts_sf, color = "black", size = 1.4) +
    scale_color_brewer(palette = "Set1", name = "Coverage") +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_blank()) +
    labs(title = "Probabilistic Ellipses (MVN coverage targets)",
         x = "Longitude (°)", y = "Latitude (°)")
  if (Sys.getenv("RSTUDIO") == "" && .Platform$OS.type == "windows") windows(10, 7.5)
  print(gg)
  ggsave(outfile, gg, width = 10, height = 7.5, dpi = 300)
  message("Wrote: ", outfile)
}

# ------------------------------------------------------------------------------
# Load points (WGS84). Use your helper if present; otherwise detect columns.
# ------------------------------------------------------------------------------
Lepto <- read.csv(file.path(data_dir, "Lepto.csv"))
if (exists("convert_to_sf_utm")) {
  pts <- convert_to_sf_utm(Lepto, input_crs = 4326, target_epsg = 4326)
} else {
  cols <- detect_latlon(Lepto)
  pts  <- st_as_sf(Lepto, coords = c(cols$lon, cols$lat), crs = 4326)
}

# ------------------------------------------------------------------------------
# Generate ellipses in degrees (EPSG:4326) for all modes
# ------------------------------------------------------------------------------
ell_arc <- generate_sde_ellipses(
  pts, group_vars = character(0), sd_levels = c(1,2),
  mode = "arcgis", compute_in = "input", output_crs = "input",
  return_metric = FALSE
)

ell_cs <- generate_sde_ellipses(
  pts, group_vars = character(0), sd_levels = c(1,2,3),
  mode = "crimestat", compute_in = "input", output_crs = "input",
  return_metric = FALSE
)

ell_prob <- generate_sde_ellipses(
  pts, group_vars = character(0),
  mode = "prob", coverage = c(0.6827, 0.95, 0.9973),
  compute_in = "input", output_crs = "input",
  return_metric = FALSE
)

# Convenience pickers for 68% / 95% probabilistic rings
pick_prob <- function(x, target) {
  x |>
    mutate(diff = abs(target_coverage - target)) |>
    arrange(diff) |>
    slice(1) |>
    select(-diff)
}
prob_68 <- pick_prob(ell_prob, 0.6827)
prob_95 <- pick_prob(ell_prob, 0.95)

# ------------------------------------------------------------------------------
# Read reference shapefiles (WGS84). If CS .prj is missing, assign EPSG:4326.
# ------------------------------------------------------------------------------
arc_1_path <- file.path(data_dir, "1sd.shp")
arc_2_path <- file.path(data_dir, "2sd.shp")
cs_1_path  <- file.path(data_dir, "SDECS_Lepto.shp")
cs_2_path  <- file.path(data_dir, "2SDECS_Lepto.shp")

have_arc <- file.exists(arc_1_path) && file.exists(arc_2_path)
have_cs  <- file.exists(cs_1_path)  && file.exists(cs_2_path)

if (have_arc) {
  arc_1 <- st_read(arc_1_path, quiet = TRUE); if (is.na(st_crs(arc_1))) st_crs(arc_1) <- 4326
  arc_2 <- st_read(arc_2_path, quiet = TRUE); if (is.na(st_crs(arc_2))) st_crs(arc_2) <- 4326
} else {
  message("ArcGIS reference shapefiles not found; skipping ArcGIS overlays.")
}

if (have_cs) {
  cs_1  <- st_read(cs_1_path,  quiet = TRUE); if (is.na(st_crs(cs_1)))  st_crs(cs_1)  <- 4326
  cs_2  <- st_read(cs_2_path,  quiet = TRUE); if (is.na(st_crs(cs_2)))  st_crs(cs_2)  <- 4326
} else {
  message("CrimeStat reference shapefiles not found; skipping CrimeStat overlays.")
}

# ------------------------------------------------------------------------------
# Build combined layers for 1× and 2× R-only comparisons
# ------------------------------------------------------------------------------
mk_combo <- function(sd_level) {
  rbind(
    ell_arc |> filter(sd_level == sd_level) |> mutate(model = "ArcGIS"),
    ell_cs  |> filter(sd_level == sd_level) |> mutate(model = "CrimeStat"),
    if (sd_level == 1) prob_68 |> mutate(model = "Prob (68.27%)")
    else                prob_95 |> mutate(model = "Prob (95%)")
  )
}
combo1 <- mk_combo(1)
combo2 <- mk_combo(2)

# ------------------------------------------------------------------------------
# Generate & SHOW all figures (and save PNGs)
# ------------------------------------------------------------------------------
plot_compare(combo1, pts, "1× SDE — ArcGIS vs CrimeStat vs Probabilistic",
             file.path(fig_dir, "Compare_1x_Ronly.png"))

plot_compare(combo2, pts, "2× SDE — ArcGIS vs CrimeStat vs Probabilistic",
             file.path(fig_dir, "Compare_2x_Ronly.png"))

if (have_arc) {
  plot_overlay_ref(arc_1, dplyr::filter(ell_arc, sd_level == 1), pts,
                   "ArcGIS — R-tool vs ArcGIS reference (1×)",
                   file.path(fig_dir, "ArcGIS_R_vs_ArcRef_1x.png"))
  plot_overlay_ref(arc_2, dplyr::filter(ell_arc, sd_level == 2), pts,
                   "ArcGIS — R-tool vs ArcGIS reference (2×)",
                   file.path(fig_dir, "ArcGIS_R_vs_ArcRef_2x.png"))
}

if (have_cs) {
  plot_overlay_ref(cs_1,  dplyr::filter(ell_cs, sd_level == 1), pts,
                   "CrimeStat — R-tool vs CrimeStat reference (1×)",
                   file.path(fig_dir, "CrimeStat_R_vs_CSRef_1x.png"))
  plot_overlay_ref(cs_2,  dplyr::filter(ell_cs, sd_level == 2), pts,
                   "CrimeStat — R-tool vs CrimeStat reference (2×)",
                   file.path(fig_dir, "CrimeStat_R_vs_CSRef_2x.png"))
}

plot_prob(ell_prob, pts, file.path(fig_dir, "Prob_ellipses.png"))

message("All figures regenerated to: ", fig_dir)
