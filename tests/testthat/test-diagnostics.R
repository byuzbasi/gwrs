test_that("Gaussian GWR diagnostics equal explicit dense operators", {
  example <- make_spatial_example(n = 32, p = 2, seed = 101)
  neighbors <- gwr_neighbors(example$coords, k = 20)
  fit <- gwr_fit(
    example$x, example$y, neighbors = neighbors, standardize = FALSE,
    control = gwrs_control(
      diagnostics = "full", n_threads = 1, tolerance = 1e-10
    )
  )
  reference <- dense_gwr_operators(example$x, neighbors)
  smoother <- reference$smoother
  diagnostics <- gwr_diagnostics(fit)

  expect_equal(diagnostics$hat_diag, diag(smoother), tolerance = 2e-9)
  expect_equal(
    diagnostics$s2_diag, rowSums(smoother^2), tolerance = 2e-9
  )
  expect_equal(diagnostics$global$trS, sum(diag(smoother)), tolerance = 2e-9)
  expect_equal(
    diagnostics$global$trStS, sum(smoother^2), tolerance = 2e-9
  )
  expect_equal(
    diagnostics$global$enp,
    2 * sum(diag(smoother)) - sum(smoother^2),
    tolerance = 2e-9
  )
  expected_se <- sqrt(
    diagnostics$global$sigma2 * reference$coefficient_map_ss
  )
  expect_equal(
    diagnostics$inference$standard_errors, expected_se,
    tolerance = 2e-9, ignore_attr = TRUE
  )
  expect_equal(gwr_rss(fit), sum(fit$residuals^2))
  expect_equal(gwr_aicc(fit), diagnostics$global$aicc)
  expect_equal(AIC(fit), diagnostics$global$aic)
  expect_equal(as.numeric(logLik(fit)), diagnostics$global$logLik)
  second <- gwr_fit(
    example$x, example$y, neighbors = neighbors, standardize = TRUE,
    control = gwrs_control(diagnostics = "standard", n_threads = 1)
  )
  expect_equal(nrow(AIC(fit, second)), 2L)
})

test_that("conditional penalized leverage agrees with a finite difference", {
  example <- make_spatial_example(n = 30, p = 2, seed = 103)
  neighbors <- gwr_neighbors(example$coords, k = 20)
  control <- gwrs_control(
    diagnostics = "standard", n_threads = 1, tolerance = 1e-11,
    max_iterations = 5000
  )
  fit <- gwr_sl_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.08,
    alpha = 0.65, d = 0.4, standardize = FALSE, control = control
  )
  epsilon <- 1e-6
  for (target in c(2L, 9L, 17L, 26L)) {
    perturbed_y <- example$y
    perturbed_y[target] <- perturbed_y[target] + epsilon
    perturbed <- gwr_sl_fit(
      example$x, perturbed_y, neighbors = neighbors, lambda = 0.08,
      alpha = 0.65, d = 0.4, standardize = FALSE,
      control = gwrs_control(
        diagnostics = "none", n_threads = 1, tolerance = 1e-11,
        max_iterations = 5000
      )
    )
    numerical <- (perturbed$fitted.values[target] -
      fit$fitted.values[target]) / epsilon
    expect_equal(numerical, fit$hat_diag[target], tolerance = 2e-5)
  }
})

test_that("path diagnostics reproduce scalar fits", {
  example <- make_spatial_example(n = 34, p = 3, seed = 107)
  neighbors <- gwr_neighbors(example$coords, k = 22)
  lambda <- c(0.18, 0.09, 0.03)
  path <- gwr_sl_path(
    example$x, example$y, neighbors = neighbors, lambda = lambda,
    alpha = 0.7, d = 0.5, return_diagnostics = TRUE,
    control = serial_control(tolerance = 1e-10)
  )
  for (index in seq_along(lambda)) {
    fit <- gwr_sl_fit(
      example$x, example$y, neighbors = neighbors, lambda = lambda[index],
      alpha = 0.7, d = 0.5,
      control = gwrs_control(
        diagnostics = "standard", n_threads = 1, tolerance = 1e-10,
        max_iterations = 3000
      )
    )
    expect_equal(
      path$hat_diag_path[, index], fit$hat_diag, tolerance = 2e-9
    )
    expect_equal(
      path$s2_diag_path[, index], fit$s2_diag, tolerance = 2e-9
    )
    expect_equal(
      path$path_summary$rss[index], sum(fit$residuals^2), tolerance = 2e-9
    )
  }
  expect_equal(dim(path$post_s2_diag_path), c(34L, 3L))
  expect_true(all(is.finite(path$path_summary$aicc)))
  expect_true(all(is.finite(path$path_summary$ebic)))
})

test_that("serial and parallel statistical diagnostics are deterministic", {
  example <- make_spatial_example(n = 46, p = 3, seed = 109)
  neighbors <- gwr_neighbors(example$coords, k = 24)
  fit_one <- gwr_sl_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.07,
    alpha = 0.6, d = 0.35,
    control = gwrs_control(diagnostics = "full", n_threads = 1)
  )
  fit_two <- gwr_sl_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.07,
    alpha = 0.6, d = 0.35,
    control = gwrs_control(
      diagnostics = "full", n_threads = 2, grain_size = 3
    )
  )
  expect_equal(fit_two$model_diagnostics, fit_one$model_diagnostics,
               tolerance = 1e-12)
})

test_that("diagnostic levels and exploratory inference are explicit", {
  example <- make_spatial_example(n = 30, p = 2, seed = 113)
  neighbors <- gwr_neighbors(example$coords, k = 18)
  none <- gwr_fit(
    example$x, example$y, neighbors = neighbors,
    control = gwrs_control(diagnostics = "none", n_threads = 1)
  )
  expect_null(none$model_diagnostics)
  expect_error(gwr_diagnostics(none), "unavailable")

  standard <- gwr_fit(
    example$x, example$y, neighbors = neighbors,
    control = gwrs_control(diagnostics = "standard", n_threads = 1)
  )
  expect_error(gwr_local_inference(standard), "full")

  full <- gwr_fit(
    example$x, example$y, neighbors = neighbors,
    control = gwrs_control(diagnostics = "full", n_threads = 1)
  )
  inference <- gwr_local_inference(full, adjust = "BH")
  expect_s3_class(inference, "gwrs_local_inference")
  expect_true(inference$settings$exploratory)
  expect_true(all(c("p_value", "p_adjusted", "significant") %in%
                    names(inference$local)))
})
