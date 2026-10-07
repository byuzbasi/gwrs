# Reproducible core benchmark; run after installing gwrs.
library(gwrs)

set.seed(2026)
n <- 10000L
p <- 12L
k <- 100L
coords <- cbind(runif(n), runif(n))
x <- matrix(rnorm(n * p), nrow = n, ncol = p)
y <- 1 + 1.5 * x[, 1L] - x[, 2L] +
  sin(4 * coords[, 1L]) + rnorm(n, sd = 0.5)

neighbor_time <- system.time({
  neighbors <- gwr_neighbors(coords, k = k)
})

serial_core_time <- system.time({
  serial_core <- gwr_sl_fit(
    x, y,
    neighbors = neighbors,
    lambda = 0.08,
    alpha = 0.7,
    d = 0.5,
    control = gwrs_control(
      n_threads = 1L, keep_data = FALSE, diagnostics = "none"
    )
  )
})

parallel_core_time <- system.time({
  core_fit <- gwr_sl_fit(
    x, y,
    neighbors = neighbors,
    lambda = 0.08,
    alpha = 0.7,
    d = 0.5,
    control = gwrs_control(
      n_threads = -1L, keep_data = FALSE, diagnostics = "none"
    )
  )
})

standard_diagnostic_time <- system.time({
  fit <- gwr_sl_fit(
    x, y,
    neighbors = neighbors,
    lambda = 0.08,
    alpha = 0.7,
    d = 0.5,
    control = gwrs_control(
      n_threads = -1L, keep_data = TRUE, diagnostics = "standard"
    )
  )
})

full_diagnostic_time <- system.time({
  full_fit <- gwr_sl_fit(
    x, y,
    neighbors = neighbors,
    lambda = 0.08,
    alpha = 0.7,
    d = 0.5,
    control = gwrs_control(
      n_threads = -1L, keep_data = FALSE, diagnostics = "full"
    )
  )
})

serial_path_time <- system.time({
  serial_path <- gwr_sl_path(
    x, y,
    neighbors = neighbors,
    n_lambda = 30L,
    alpha = 0.7,
    d = 0.5,
    keep_coefficients = FALSE,
    control = gwrs_control(
      n_threads = 1L, keep_data = FALSE, diagnostics = "none"
    )
  )
})

parallel_path_time <- system.time({
  path <- gwr_sl_path(
    x, y,
    neighbors = neighbors,
    n_lambda = 30L,
    alpha = 0.7,
    d = 0.5,
    keep_coefficients = FALSE,
    control = gwrs_control(
      n_threads = -1L, keep_data = FALSE, diagnostics = "none"
    )
  )
})

diagnostic_path_time <- system.time({
  diagnostic_path <- gwr_sl_path(
    x, y,
    neighbors = neighbors,
    n_lambda = 30L,
    alpha = 0.7,
    d = 0.5,
    keep_coefficients = FALSE,
    return_diagnostics = TRUE,
    control = gwrs_control(
      n_threads = -1L, keep_data = FALSE, diagnostics = "none"
    )
  )
})

collinearity_time <- system.time({
  collinearity <- gwr_local_collinearity(fit)
})

moran_time <- system.time({
  moran <- gwr_moran(fit, permutations = 99L, seed = 2026L)
})

result <- list(
  dimensions = c(n = n, p = p, k = k),
  available_threads = RcppParallel::defaultNumThreads(),
  neighbor_time = neighbor_time,
  serial_core_time = serial_core_time,
  parallel_core_time = parallel_core_time,
  core_elapsed_speedup = unname(
    serial_core_time[["elapsed"]] / parallel_core_time[["elapsed"]]
  ),
  standard_diagnostic_time = standard_diagnostic_time,
  full_diagnostic_time = full_diagnostic_time,
  serial_path_time = serial_path_time,
  parallel_path_time = parallel_path_time,
  path_elapsed_speedup = unname(
    serial_path_time[["elapsed"]] / parallel_path_time[["elapsed"]]
  ),
  diagnostic_path_time = diagnostic_path_time,
  collinearity_time = collinearity_time,
  moran_99_time = moran_time,
  fit_convergence = mean(fit$diagnostics$converged),
  path_convergence = range(path$path_summary$convergence_rate),
  diagnostic_path_convergence = range(
    diagnostic_path$path_summary$convergence_rate
  ),
  object_mib = c(
    core = as.numeric(object.size(core_fit)) / 1024^2,
    standard = as.numeric(object.size(fit)) / 1024^2,
    full = as.numeric(object.size(full_fit)) / 1024^2,
    path = as.numeric(object.size(path)) / 1024^2,
    diagnostic_path = as.numeric(object.size(diagnostic_path)) / 1024^2
  ),
  dense_weight_matrix_gib_avoided = 8 * n^2 / 1024^3,
  diagnostic_sanity = c(
    aicc = gwr_aicc(fit),
    moran_I = moran$I,
    median_condition_number = stats::median(
      collinearity$condition_number$condition_number, na.rm = TRUE
    )
  )
)
print(result)
