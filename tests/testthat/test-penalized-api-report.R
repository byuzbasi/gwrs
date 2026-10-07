test_that("unified penalized interfaces map exactly to their core cases", {
  example <- make_spatial_example(n = 40, p = 3, seed = 149)
  neighbors <- gwr_neighbors(example$coords, k = 24)
  control <- gwrs_control(diagnostics = "standard", n_threads = 1)
  specifications <- list(
    lasso = c(alpha = 1, d = 0),
    ridge = c(alpha = 0, d = 0),
    elastic_net = c(alpha = 0.4, d = 0),
    scaled_lasso = c(alpha = 0.4, d = 0.6)
  )
  for (penalty in names(specifications)) {
    parameter <- specifications[[penalty]]
    unified <- gwr_penalized_fit(
      example$x, example$y, neighbors = neighbors, lambda = 0.08,
      penalty = penalty, alpha = parameter[["alpha"]],
      d = parameter[["d"]], control = control
    )
    direct <- gwr_sl_fit(
      example$x, example$y, neighbors = neighbors, lambda = 0.08,
      alpha = parameter[["alpha"]], d = parameter[["d"]], control = control
    )
    expect_equal(unified$coefficients, direct$coefficients, tolerance = 0)
    expect_equal(unified$fitted.values, direct$fitted.values, tolerance = 0)
  }

  data <- data.frame(
    y = example$y, x1 = example$x[, 1], x2 = example$x[, 2],
    x3 = example$x[, 3], cx = example$coords[, 1], cy = example$coords[, 2]
  )
  formula_fit <- gwr_en(
    y ~ x1 + x2 + x3, data, c("cx", "cy"), k = 24,
    lambda = 0.08, alpha = 0.4, control = control
  )
  expect_equal(formula_fit$method, "GWR-Elastic Net")
  expect_equal(colnames(formula_fit$coefficients),
               c("(Intercept)", "x1", "x2", "x3"))
})

test_that("unified paths and spatial CV retain model labels", {
  example <- make_spatial_example(n = 36, p = 2, seed = 151)
  neighbors <- gwr_neighbors(example$coords, k = 20)
  ridge <- gwr_ridge_path(
    example$x, example$y, neighbors = neighbors,
    lambda = c(0.2, 0.1, 0.03), return_diagnostics = TRUE,
    control = gwrs_control(n_threads = 1)
  )
  expect_equal(ridge$method, "GWR-Ridge path")
  expect_equal(ridge$alpha, 0)
  expect_error(
    gwr_ridge_path(example$x, example$y, neighbors = neighbors),
    "explicitly"
  )
  folds <- spatial_folds(example$coords, n_folds = 3, seed = 8)
  cv <- cv_gwr_penalized(
    example$x, example$y, example$coords, k = 14, fold = folds,
    lambda = c(0.12, 0.05), penalty = "lasso", refit = TRUE,
    control = gwrs_control(n_threads = 1, diagnostics = "none")
  )
  expect_equal(cv$penalty, "lasso")
  expect_equal(cv$fit$method, "GWR-Lasso")
  expect_equal(cv$fit$alpha, 1)
  no_refit <- cv_gwr_penalized(
    example$x, example$y, example$coords, k = 14, fold = folds,
    lambda = c(0.12, 0.05), penalty = "lasso", refit = FALSE,
    control = gwrs_control(n_threads = 1, diagnostics = "none")
  )
  expect_error(gwr_rss(no_refit), "refit = FALSE")
})

test_that("non-Gaussian extensions expose only compatible diagnostics", {
  example <- make_spatial_example(n = 24, p = 2, seed = 157)
  neighbors <- gwr_neighbors(example$coords, k = 14)
  lad <- lad_gwr_sl_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.04,
    alpha = 0.7, d = 0.4,
    control = gwrs_control(
      diagnostics = "standard", n_threads = 1, max_iterations = 4000
    )
  )
  expect_true(is.finite(lad$model_diagnostics$global$rss))
  expect_true(is.na(lad$model_diagnostics$global$aicc))
  expect_false(lad$model_diagnostics$validity$information_criteria)
  expect_error(logLik(lad), "not defined")
  expect_error(gwr_local_inference(lad), "Gaussian")

  laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
  ktdd <- gwr_ktdd_fit(
    example$x, example$y, neighbors = neighbors, kernel = "gaussian",
    diffusion_operator = laplacian, source = rep(0, 24),
    coefficient_laplacian = laplacian, lambda_diffusion = 1,
    gamma = 0.05,
    control = gwrs_control(
      diagnostics = "standard", n_threads = 1, max_iterations = 500
    )
  )
  expect_true(is.finite(ktdd$model_diagnostics$global$joint_objective))
  expect_equal(length(ktdd$loss_components), 3)
  expect_false(ktdd$model_diagnostics$validity$information_criteria)
  expect_error(AIC(ktdd), "not defined")
})

test_that("reports and base graphics methods run without optional packages", {
  example <- make_spatial_example(n = 30, p = 2, seed = 163)
  neighbors <- gwr_neighbors(example$coords, k = 18)
  fit <- gwr_lasso_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.05,
    control = gwrs_control(diagnostics = "full", n_threads = 1)
  )
  moran <- gwr_moran(fit, permutations = 9, seed = 4)
  collinearity <- gwr_local_collinearity(fit)
  output <- capture.output(report <- gwrs_report(
    fit, include_inference = TRUE, moran = moran,
    collinearity = collinearity
  ))
  expect_s3_class(report, "gwrs_report")
  expect_true(any(grepl("gwrs model report", output, fixed = TRUE)))

  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot_gwrs_model(fit, type = "coefficient",
                                   coefficient = "x1"))
  expect_invisible(plot_gwrs_inference(fit, type = "p_adjusted",
                                       coefficient = "x1"))
  expect_invisible(plot_gwrs_diagnostics(fit, type = "local_r2"))
})
