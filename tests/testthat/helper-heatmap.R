make_heatmap_example <- function(n = 72L, seed = 712L) {
  set.seed(seed)
  coords <- cbind(x = runif(n), y = runif(n))
  x <- cbind(
    income = rnorm(n),
    education = rnorm(n),
    noise = rnorm(n)
  )
  local_income <- -1.1 + 2.2 * coords[, 1L]
  residual_cluster <- 0.45 * (
    coords[, 1L] > 0.65 & coords[, 2L] > 0.65
  )
  y <- 0.5 + local_income * x[, "income"] +
    0.65 * x[, "education"] + residual_cluster +
    rnorm(n, sd = 0.2)
  list(x = x, y = y, coords = coords)
}

make_heatmap_fit <- function(diagnostics = "standard",
                             lambda = 0.25,
                             alpha = 0.85,
                             seed = 712L) {
  example <- make_heatmap_example(seed = seed)
  fit <- gwr_sl_fit(
    example$x,
    example$y,
    coords = example$coords,
    k = 42L,
    lambda = lambda,
    alpha = alpha,
    d = 0.35,
    control = gwrs_control(
      diagnostics = diagnostics,
      keep_data = TRUE,
      n_threads = 1L,
      max_iterations = 3000L
    )
  )
  list(fit = fit, example = example)
}

skip_if_no_heatmap_packages <- function() {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("sf")
}
