guarded_example <- function(p = 6L) {
  set.seed(20260904)
  x <- matrix(rnorm(36 * p), 36, p)
  x[, 2] <- x[, 1] + 0.01 * x[, 2]
  list(x = x, y = drop(x %*% c(2, -1, rep(0, p - 2))) + rnorm(36, sd = .1),
       coords = cbind(seq_len(36), rep(0, 36)))
}

test_that("solver selection is explicit, validated and old controls remain usable", {
  expect_identical(gwrs_control()$nonconvex_solver, "coordinate")
  expect_error(gwrs_control(nonconvex_solver = "unknown"), "arg")
  old <- gwrs_control(); old$nonconvex_solver <- NULL
  expect_identical(gwrs:::validate_control(old)$nonconvex_solver, "coordinate")
  old$nonconvex_solver <- "unknown"
  expect_error(gwrs:::validate_control(old), "arg")
})

test_that("guarded fit, formula, path, prediction and threads share the core", {
  e <- guarded_example()
  control <- gwrs_control(nonconvex_solver = "guarded_block", n_threads = 1,
    max_iterations = 4000, diagnostics = "none")
  for (penalty in c("scad", "mcp")) {
    fit <- gwr_nonconvex_fit(e$x, e$y, coords = e$coords, k = 28,
      lambda = .2, penalty = penalty, control = control)
    path <- gwr_nonconvex_path(e$x, e$y, coords = e$coords, k = 28,
      lambda = .2, penalty = penalty, screening = FALSE, keep_coefficients = TRUE,
      control = control)
    expect_equal(unname(coef(fit)), unname(path$coefficients[, , 1]), tolerance = 1e-12)
    expect_true(all(path$path_state == "converged"))
    expect_identical(path$control$nonconvex_solver, "guarded_block")
    parallel <- control; parallel$n_threads <- 2L
    other <- gwr_nonconvex_path(e$x, e$y, coords = e$coords, k = 28,
      lambda = .2, penalty = penalty, screening = FALSE, keep_coefficients = TRUE,
      control = parallel)
    expect_identical(path$coefficients, other$coefficients)
    expect_identical(path$solver_diagnostics, other$solver_diagnostics)
    colnames(e$x) <- paste0("x", seq_len(ncol(e$x)))
    d <- data.frame(y = e$y, e$x)
    form <- gwr_nonconvex(y ~ ., data = d, coords = e$coords, k = 28,
      lambda = .2, penalty = penalty, control = control)
    expect_equal(unname(coef(form)), unname(coef(fit)), tolerance = 1e-12)
    pred <- predict(form, newdata = d[1:3, -1], newcoords = e$coords[1:3, ], k = 28, type = "both")
    expect_equal(dim(pred$coefficients), c(3L, ncol(e$x) + 1L))
    expect_true(all(is.finite(pred$fit)))
    incomplete <- form; incomplete$control$max_iterations <- 1L
    expect_error(predict(incomplete, newdata = d[1:3, -1],
      newcoords = e$coords[1:3, ], k = 28), class = "gwrs_nonconvex_convergence_error")
  }
})

test_that("HD paths report failure without converting unfinished work into zeros", {
  e <- guarded_example(48L)
  ctl <- gwrs_control(nonconvex_solver = "guarded_block", max_iterations = 1,
    n_threads = 1, diagnostics = "none")
  for (penalty in c("scad", "mcp")) {
    expect_warning(path <- gwr_nonconvex_path(e$x, e$y, coords = e$coords, k = 24,
      lambda = c(100, .01, .001), penalty = penalty, control = ctl,
      keep_coefficients = TRUE, return_diagnostics = TRUE), "failed or unattempted")
    expect_true(all(path$path_state[, 1] == "converged"))
    expect_true(any(path$path_state[, 2] == "failed"))
    failed <- path$path_state[, 2] == "failed"
    expect_true(all(path$path_state[failed, 3] == "unattempted"))
    expect_true(all(is.na(path$fitted.values[failed, 2:3])))
    expect_true(all(is.finite(path$coefficients[failed, , 2])))
    expect_true(all(is.na(path$coefficients[failed, , 3])))
    expect_true(all(is.na(path$path_diagnostics[!path$path_summary$valid, ])))
    error <- tryCatch(gwr_nonconvex_fit(e$x, e$y, coords = e$coords, k = 24,
      lambda = .01, penalty = penalty, control = ctl), error = identity)
    expect_s3_class(error, "gwrs_nonconvex_convergence_error")
    expect_true(any(error$solver_result$path_state == 2L))
  }
})

test_that("guarded CV keeps full denominators and invalidates whole lambda scores", {
  e <- guarded_example(48L)
  ctl <- gwrs_control(nonconvex_solver = "guarded_block", max_iterations = 1,
    n_threads = 1, diagnostics = "none")
  folds <- rep(1:3, each = 12)
  cv <- cv_gwr_mcp(e$x, e$y, e$coords, k = 18, fold = folds,
    lambda = c(100, .01, .001), refit = FALSE, control = ctl)
  expect_equal(cv$cv$observation_count, rep(36L, 3))
  expect_true(cv$cv$valid[1])
  expect_true(any(!cv$cv$valid))
  expect_true(all(is.na(cv$cv[!cv$cv$valid, c("rmse", "mse", "mae")])))
  expect_equal(cv$best$lambda, 100)
  expect_warning(none <- cv_gwr_mcp(e$x, e$y, e$coords, k = 18, fold = folds,
    lambda = c(.01, .001), refit = TRUE, control = ctl), "No lambda converged")
  expect_equal(nrow(none$best), 0L)
  expect_null(none$fit)
  expect_output(print(none), "cross-validation")
  expect_true(all(none$cv$observation_count == 36L))
  ctl$max_iterations <- 1000L
  fitted <- cv_gwr_mcp(e$x, e$y, e$coords, k = 18, fold = folds,
    lambda = c(100, 50), refit = TRUE, control = ctl)
  expect_identical(fitted$fit$control$nonconvex_solver, "guarded_block")
  expect_true(all(fitted$fit$solver_diagnostics$converged))
})

test_that("guarded screening, scaling, kernels and zero-lambda interfaces remain valid", {
  e <- guarded_example()
  ctl <- gwrs_control(nonconvex_solver = "guarded_block", max_iterations = 4000,
    n_threads = 1, diagnostics = "none")
  constant <- e$x; constant[, 6] <- 1
  expect_error(gwr_scad_path(constant, e$y, coords = e$coords, k = 28,
    lambda = c(100, 50), control = ctl), "Constant or numerically constant")
  for (kernel in c("gaussian", "bisquare", "exponential", "tricube", "boxcar")) {
    for (standardize in c(TRUE, FALSE)) {
      path <- gwr_scad_path(e$x, e$y, coords = e$coords, k = 28,
        lambda = c(100, 50), kernel = kernel, standardize = standardize, control = ctl)
      expect_true(all(path$path_summary$valid))
    }
  }
  a <- gwr_mcp_path(e$x, e$y, coords = e$coords, k = 28, lambda = c(10, 5), control = ctl)
  b <- gwr_mcp_path(e$x, e$y, coords = e$coords, k = 28, lambda = c(10, 5), screening = FALSE, control = ctl)
  expect_identical(a$fitted.values, b$fitted.values)
  e$x <- e$x[, 1:4]
  z <- gwr_scad_fit(e$x, e$y, coords = e$coords, k = 28, lambda = 0, control = ctl)
  expect_true(all(z$solver_diagnostics$converged))
})

test_that("guarded statistical diagnostics retain conditional validity labels", {
  e <- make_spatial_example(n = 42, p = 4, seed = 403)
  ctl <- gwrs_control(nonconvex_solver = "guarded_block", max_iterations = 4000,
    n_threads = 1, diagnostics = "full")
  fit <- gwr_mcp_fit(e$x, e$y, coords = e$coords, k = 27, lambda = .1, control = ctl)
  expect_match(fit$model_diagnostics$validity$basis, "nonconvex")
  expect_identical(gwr_local_inference(fit)$settings$estimator, "post_selection")
  tests <- gwr_diagnostic_tests(fit, allow_approximate = TRUE)
  expect_match(tests$validity$basis, "nonconvex")
  path <- gwr_mcp_path(e$x, e$y, coords = e$coords, k = 27, lambda = .1,
    return_diagnostics = TRUE, control = ctl)
  expect_true(all(path$path_summary$valid))
  expect_true(all(is.finite(path$path_diagnostics$rss)))
})
