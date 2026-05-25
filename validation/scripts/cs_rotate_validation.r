###########################################################
# CrimeStat synthetic rotation validation
# Compare SDEtool vs CrimeStat using rotated synthetic data
###########################################################

library(sf)
library(SDEtool)

###########################################################
# Generate synthetic rotated data
###########################################################

set.seed(42)

n <- 50

x <- rnorm(n, mean = 0, sd = 3)
y <- rnorm(n, mean = 0, sd = 1)

theta <- pi / 4

x_rot <- x*cos(theta) - y*sin(theta)
y_rot <- x*sin(theta) + y*cos(theta)

df_test <- data.frame(
  longitude = x_rot + 103,
  latitude = y_rot + 18
)

###########################################################
# Generate SDEtool ellipse
###########################################################

sf_pts_proj <- convert_to_sf_utm(df_test)

sde_cs <- generate_sde_ellipses(
  sf_data = sf_pts_proj,
  group_vars = NULL,
  mode = "crimestat"
)

r_rot <- sde_cs[sde_cs$sd_level == 1, ]

###########################################################
# Load CrimeStat output
###########################################################

cs_rot <- st_read(
"C:/Users/.../SDEtool/validation/data/CrimeStatValid/Rotate/SDESDECS_Rotate.shp"
)

st_crs(cs_rot) <- 4326

cs_rot <- st_transform(
  cs_rot,
  st_crs(r_rot)
)

###########################################################
# IoU
###########################################################

inter_area <- st_area(
  st_intersection(
    st_make_valid(cs_rot),
    st_make_valid(r_rot)
  )
)

union_area <- st_area(
  st_union(
    st_make_valid(cs_rot),
    st_make_valid(r_rot)
  )
)

iou <- as.numeric(
  inter_area / union_area
)

print(iou)

###########################################################
# Export figure
###########################################################

png(
"C:/Users/...SDEtool/validation/figures/CrimeStat_RotateValidation.png",
width = 2400,
height = 1800,
res = 300
)

cs_rot_ll <- st_transform(cs_rot, 4326)
r_rot_ll <- st_transform(r_rot, 4326)

plot(
  df_test$longitude,
  df_test$latitude,
  pch = 16,
  col = "grey70",
  asp = 1,
  xlab = "Longitude",
  ylab = "Latitude",
  main = "Synthetic Rotation Validation"
)

grid(col = "grey80")

plot(
  st_geometry(cs_rot_ll),
  add = TRUE,
  border = "red",
  lwd = 2
)

plot(
  st_geometry(r_rot_ll),
  add = TRUE,
  border = "blue",
  lwd = 2
)

legend(
  "bottomright",
  legend = c(
    "Synthetic points",
    "CrimeStat",
    "SDEtool"
  ),
  col = c(
    "grey70",
    "red",
    "blue"
  ),
  pch = c(16, NA, NA),
  lty = c(NA,1,1),
  lwd = c(NA,2,2),
  bty = "n"
)

dev.off()
