library(MASS)

set.seed(123)

results <- data.frame()

for(i in 1:100){

  sim <- MASS::mvrnorm(
    n = 10000,
    mu = c(0,0),
    Sigma = matrix(
      c(4,1,
        1,2),
      nrow = 2
    )
  )

  df_sim <- data.frame(
    longitude = sim[,1],
    latitude  = sim[,2]
  )

  sf_sim <- sf::st_as_sf(
    df_sim,
    coords = c("longitude","latitude"),
    crs = 4326
  )

  prob_sde <- generate_sde_ellipses(
    sf_data = sf_sim,
    group_vars = NULL,
    mode = "prob",
    coverage = c(
      0.6827,
      0.95,
      0.9973
    ),
    compute_in = "input"
  )

  results <- rbind(
    results,
    data.frame(
      iter = i,
      target = prob_sde$target_coverage,
      observed = prob_sde$percent_inside / 100
    )
  )

}

aggregate(
  observed ~ target,
  data = results,
  FUN = function(x)
    c(
      mean = mean(x),
      sd = sd(x)
    )
)
