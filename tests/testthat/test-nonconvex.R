nonconvex_penalty_r <- function(beta, lambda, gamma, penalty) {
  value <- abs(beta)
  if (penalty == "mcp") {
    ifelse(
      value <= gamma * lambda,
      lambda * value - value^2 / (2 * gamma),
      0.5 * gamma * lambda^2
    )
  } else {
    ifelse(
      value <= lambda,
      lambda * value,
      ifelse(
        value <= gamma * lambda,
        (-value^2 + 2 * gamma * lambda * value - lambda^2) /
          (2 * (gamma - 1)),
        0.5 * (gamma + 1) * lambda^2
      )
    )
  }
}

independent_coordinate_minimum <- function(score, curvature, lambda,
                                           gamma, penalty) {
  objective <- function(beta) {
    0.5 * curvature * beta^2 - score * beta +
      nonconvex_penalty_r(beta, lambda, gamma, penalty)
  }
  bound <- abs(score) / curvature + gamma * lambda + 1
  knots <- if (penalty == "mcp") gamma * lambda else c(lambda, gamma * lambda)
  cuts <- sort(unique(c(-bound, -rev(knots), 0, knots, bound)))
  candidates <- cuts
  for (index in seq_len(length(cuts) - 1L)) {
    interval <- cuts[c(index, index + 1L)]
    if (diff(interval) <= 1e-12) next
    candidates <- c(
      candidates,
      stats::optimize(objective, interval = interval, tol = 1e-14)$minimum
    )
  }
  values <- vapply(candidates, objective, numeric(1L))
  candidates[which.min(values)]
}

test_that("exact SCAD/MCP coordinate updates minimize the scalar objective", {
  cases <- expand.grid(
    score = c(-1.7, -0.2, 0.15, 1.4),
    curvature = c(0.08, 0.3, 1.2),
    lambda = c(0.07, 0.4),
    stringsAsFactors = FALSE
  )
  for (penalty in c("mcp", "scad")) {
    code <- if (penalty == "mcp") 0L else 1L
    gamma <- if (penalty == "mcp") 3 else 3.7
    for (row in seq_len(nrow(cases))) {
      current <- cases[row, ]
      observed <- gwrs:::cpp_nonconvex_coordinate_minimum(
        current$score, current$curvature, current$lambda, code, gamma
      )
      reference <- independent_coordinate_minimum(
        current$score, current$curvature, current$lambda, gamma, penalty
      )
      observed_objective <-
        0.5 * current$curvature * observed^2 - current$score * observed +
        nonconvex_penalty_r(observed, current$lambda, gamma, penalty)
      reference_objective <-
        0.5 * current$curvature * reference^2 - current$score * reference +
        nonconvex_penalty_r(reference, current$lambda, gamma, penalty)
      expect_lte(observed_objective, reference_objective + 2e-9)
    }
  }
})

test_that("nonconvex lambda anchor starts at an all-zero local slope surface", {
  example <- make_spatial_example(n = 48, p = 5, seed = 401)
  neighbors <- gwr_neighbors(example$coords, k = 28)
  control <- gwrs_control(n_threads = 1, diagnostics = "none")
  for (penalty in c("mcp", "scad")) {
    gamma <- if (penalty == "mcp") 3 else 3.7
    anchor <- gwr_nonconvex_lambda_max(
      example$x, example$y, neighbors = neighbors,
      penalty = penalty, gamma = gamma, control = control
    )
    path <- gwr_nonconvex_path(
      example$x, example$y, neighbors = neighbors,
      lambda = c(anchor, 0.6 * anchor), penalty = penalty, gamma = gamma,
      keep_coefficients = TRUE, control = control
    )
    expect_true(is.finite(anchor) && anchor >= 0)
    expect_lte(max(abs(path$coefficients[, -1L, 1L])), 1e-9)
    expect_lte(path$path_summary$max_coordinate_gap[1L], 1e-9)
  }
})

test_that("single fits, one-lambda paths, screening, and threads agree", {
  example <- make_spatial_example(n = 55, p = 6, seed = 402)
  neighbors <- gwr_neighbors(example$coords, k = 32)
  serial <- gwrs_control(
    tolerance = 1e-9, max_iterations = 4000, n_threads = 1,
    diagnostics = "none"
  )
  for (penalty in c("mcp", "scad")) {
    gamma <- if (penalty == "mcp") 3 else 3.7
    fit <- gwr_nonconvex_fit(
      example$x, example$y, neighbors = neighbors, lambda = 0.12,
      penalty = penalty, gamma = gamma, control = serial
    )
    path <- gwr_nonconvex_path(
      example$x, example$y, neighbors = neighbors, lambda = 0.12,
      penalty = penalty, gamma = gamma, screening = FALSE,
      keep_coefficients = TRUE, control = serial
    )
    expect_equal(unname(fit$coefficients),
                 unname(path$coefficients[, , 1L]),
                 tolerance = 2e-8)

    lambda <- c(0.35, 0.18, 0.08)
    screened <- gwr_nonconvex_path(
      example$x, example$y, neighbors = neighbors, lambda = lambda,
      penalty = penalty, gamma = gamma, screening = TRUE,
      keep_coefficients = TRUE, control = serial
    )
    full <- gwr_nonconvex_path(
      example$x, example$y, neighbors = neighbors, lambda = lambda,
      penalty = penalty, gamma = gamma, screening = FALSE,
      keep_coefficients = TRUE, control = serial
    )
    expect_equal(screened$fitted.values, full$fitted.values,
                 tolerance = 2e-7)
    expect_equal(screened$coefficients, full$coefficients,
                 tolerance = 2e-7)

    parallel <- gwr_nonconvex_fit(
      example$x, example$y, neighbors = neighbors, lambda = 0.12,
      penalty = penalty, gamma = gamma,
      control = gwrs_control(
        tolerance = 1e-9, max_iterations = 4000, n_threads = 2,
        grain_size = 3, diagnostics = "none"
      )
    )
    expect_equal(parallel$coefficients, fit$coefficients, tolerance = 1e-12)
  }
})

test_that("nonconvex diagnostics and unified APIs retain validity labels", {
  example <- make_spatial_example(n = 42, p = 4, seed = 403)
  neighbors <- gwr_neighbors(example$coords, k = 27)
  control <- gwrs_control(n_threads = 1, diagnostics = "full")
  fit <- gwr_penalized_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.1,
    penalty = "mcp", gamma = 3, control = control
  )
  expect_s3_class(fit, "gwr_mcp_fit")
  expect_identical(fit$penalty, "mcp")
  expect_true(all(c(
    "local_convex", "minimum_active_hessian_eigenvalue",
    "penalty_knot_count"
  ) %in% names(fit$model_diagnostics$local)))
  expect_match(fit$model_diagnostics$validity$basis, "nonconvex")
  expect_identical(gwr_local_inference(fit)$settings$estimator,
                   "post_selection")
  tests <- gwr_diagnostic_tests(fit, allow_approximate = TRUE)
  expect_match(tests$validity$basis, "nonconvex")
})

test_that("spatial CV refits the selected nonconvex model", {
  example <- make_spatial_example(n = 45, p = 4, seed = 404)
  fold <- spatial_folds(example$coords, n_folds = 3, seed = 11)
  cv <- cv_gwr_mcp(
    example$x, example$y, example$coords, k = 18, fold = fold,
    lambda = c(0.25, 0.1), refit = TRUE,
    control = gwrs_control(n_threads = 1, diagnostics = "none")
  )
  expect_s3_class(cv, "gwrs_nonconvex_cv")
  expect_s3_class(cv$fit, "gwr_mcp_fit")
  expect_equal(cv$fit$lambda, cv$best$lambda)
  expect_equal(nrow(cv$cv), 2L)
})

test_that("formula fits and new-location prediction use the nonconvex core", {
  example <- make_spatial_example(n = 46, p = 4, seed = 405)
  predictor_names <- paste0("x", seq_len(ncol(example$x)))
  colnames(example$x) <- predictor_names
  data <- data.frame(y = example$y, example$x, check.names = FALSE)
  control <- gwrs_control(
    n_threads = 1, keep_data = TRUE, diagnostics = "none"
  )
  fit <- gwr_scad(
    y ~ ., data = data, coords = example$coords,
    k = 20, lambda = 0.1, control = control
  )
  expect_s3_class(fit, "gwr_scad_fit")
  expect_identical(fit$penalty, "scad")

  prediction <- predict(
    fit,
    newdata = data[1:5, predictor_names, drop = FALSE],
    newcoords = example$coords[1:5, , drop = FALSE],
    k = 20,
    type = "both"
  )
  expect_length(prediction$fit, 5L)
  expect_equal(dim(prediction$coefficients), c(5L, 5L))
  expect_true(all(is.finite(prediction$fit)))
  expect_true(all(is.finite(prediction$coefficients)))
})
