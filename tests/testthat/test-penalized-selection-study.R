library(gwrs)

study_helper <- function(file) {
  installed <- system.file("simulations", file, package = "gwrs")
  if (nzchar(installed)) return(installed)
  candidates <- c(
    file.path("inst", "simulations", file),
    file.path("..", "..", "inst", "simulations", file)
  )
  candidates[file.exists(candidates)][[1L]]
}

sys.source(study_helper("selection-metrics.R"), envir = environment())
sys.source(study_helper("penalized-study-utils.R"), envir = environment())

test_that("the approved paper design has deterministic dimensions", {
  first <- gwrs_penalized_study_plan("production", 20260826L)
  second <- gwrs_penalized_study_plan("production", 20260826L)

  expect_identical(first, second)
  expect_equal(nrow(first), 1900L)
  expect_equal(length(unique(first$cell_id)), 39L)
  expect_equal(sum(first$section == "main"), 1600L)
  expect_equal(sum(first$section == "reference"), 200L)
  expect_equal(sum(first$section == "stress"), 100L)
  expect_identical(anyDuplicated(first$task_id), 0L)
  expect_identical(anyDuplicated(first$seed), 0L)
  expect_equal(nrow(first) * 8L, 15200L)

  main_structures <- unique(first[
    first$section == "main", c("predictor_rho", "spatial_rho")
  ])
  expect_equal(nrow(main_structures), 4L)
  expect_setequal(
    paste(main_structures$predictor_rho, main_structures$spatial_rho),
    c("0 0", "0.7 0", "0 0.7", "0.7 0.7")
  )

  smoke <- gwrs_penalized_study_plan("smoke", 20260826L)
  expect_equal(nrow(smoke), 2L)

  calibration <- gwrs_penalized_study_plan(
    "max_cell_calibration", 20260826L
  )
  production_match <- first[
    first$section == "stress" &
      first$design_variant == "p100" &
      first$scenario_id == "regional_sparse_discontinuous" &
      first$replication == 1L,
    , drop = FALSE
  ]
  expect_equal(nrow(calibration), 1L)
  expect_identical(calibration$n, 10000L)
  expect_identical(calibration$p, 100L)
  expect_identical(calibration$k_candidates, "250:500:1000")
  expect_identical(calibration$data_seed, production_match$data_seed)
  expect_identical(calibration$fold_seed, production_match$fold_seed)
  expect_identical(calibration$seed, production_match$seed)
  expect_identical(calibration$section, "max_cell_calibration")
})

test_that("the approved k rule respects dimension and density", {
  expect_identical(
    gwrs_penalized_study_k_candidates(2000L, 30L),
    c(62L, 124L, 248L)
  )
  expect_identical(
    gwrs_penalized_study_k_candidates(10000L, 100L),
    c(250L, 500L, 1000L)
  )
  expect_true(all(gwrs_penalized_study_k_candidates(120L, 6L) > 6L))
})

test_that("data generation is deterministic with the declared SNR", {
  task <- gwrs_penalized_study_plan("smoke", 20260826L)[1L, , drop = FALSE]
  config <- list(truth_tolerance = 1e-12)
  first <- gwrs_penalized_study_simulate(task, config)
  second <- gwrs_penalized_study_simulate(task, config)

  expect_identical(first, second)
  expect_equal(first$signal_sd^2 / first$noise_sd^2, task$snr[[1L]])
  expect_lt(max(abs(colMeans(first$x))), 1e-12)
  expect_equal(unname(sqrt(colMeans(first$x^2))), rep(1, task$p[[1L]]),
               tolerance = 1e-12)
  expect_true(any(first$beta == 0))
  expect_true(any(first$beta != 0))
  expect_true(all(first$beta[!first$support] == 0))
  near_edge <- rbind(c(0.30 + 0.40 * (1 - 1e-4), 0.70), c(0.71, 0.70))
  near_edge_bump <- gwrs_penalized_study_plateau_wendland(
    near_edge, center = c(0.30, 0.70), radius = 0.40
  )
  expect_true(near_edge_bump$support[[1L]])
  expect_lt(abs(near_edge_bump$value[[1L]]), 1e-12)
  expect_false(near_edge_bump$support[[2L]])
  expect_equal(nrow(first$x_target), task$target_grid_side[[1L]]^2)
  expect_identical(dim(first$beta_target), dim(first$x_target))
  expect_false(identical(first$y_train, first$y_test))

  null_task <- task
  null_task$scenario_id <- "null_slopes"
  null_data <- gwrs_penalized_study_simulate(null_task, config)
  expect_true(all(null_data$beta == 0))
  expect_false(any(null_data$support))
  expect_equal(null_data$signal_sd^2 / null_data$noise_sd^2,
               null_task$snr[[1L]])
})

test_that("stress designs generate the intended geometry and correlation", {
  plan <- gwrs_penalized_study_plan("production", 20260826L)
  stress_task <- function(variant) {
    task <- plan[
      plan$section == "stress" & plan$design_variant == variant,
      , drop = FALSE
    ][1L, , drop = FALSE]
    task$n <- 400L
    task$p <- min(task$p, 8L)
    task$target_grid_side <- 8L
    task
  }

  clustered <- gwrs_penalized_study_simulate(
    stress_task("clustered_sampling"), list()
  )
  expect_gt(mean(clustered$domain_edge), 0)
  centers <- rbind(c(0.20, 0.22), c(0.30, 0.78), c(0.72, 0.30), c(0.78, 0.76))
  nearest_center <- apply(clustered$coords, 1L, function(point) {
    min(sqrt(rowSums(sweep(centers, 2L, point, FUN = "-")^2)))
  })
  expect_gt(mean(nearest_center < 0.15), 0.75)

  l_shaped <- gwrs_penalized_study_simulate(
    stress_task("l_domain_boundary"), list()
  )
  expect_false(any(
    l_shaped$coords[, 1L] > 0.55 & l_shaped$coords[, 2L] > 0.55
  ))
  expect_lt(nrow(l_shaped$target_coords), 8L^2)

  collinear <- gwrs_penalized_study_simulate(
    stress_task("local_collinearity"), list()
  )
  expect_gt(collinear$realized_local_correlation, 0.8)
  expect_lt(abs(collinear$realized_outside_correlation), 0.25)
  expect_false(any(collinear$support[, 2L]))
})

test_that("oracle-union GWR expands coefficients without selecting noise", {
  task <- gwrs_penalized_study_plan("smoke", 20260826L)[1L, , drop = FALSE]
  config <- list(
    truth_tolerance = 1e-12,
    kernel = "gaussian"
  )
  data <- gwrs_penalized_study_simulate(task, config)
  neighbors <- gwr_neighbors(data$coords, k = 28L)
  control <- gwrs_control(
    n_threads = 1L, keep_data = FALSE, diagnostics = "none"
  )
  fit <- gwrs_penalized_study_oracle_union(
    data, neighbors, config, control
  )

  expect_equal(dim(fit$coefficients), c(task$n[[1L]], task$p[[1L]] + 1L))
  inactive_union <- which(colSums(abs(data$beta) > 1e-12) == 0L)
  expect_true(all(fit$coefficients[, inactive_union + 1L] == 0))
  expect_true(all(fit$diagnostics$converged))

  null_task <- task
  null_task$scenario_id <- "null_slopes"
  null_data <- gwrs_penalized_study_simulate(null_task, config)
  null_fit <- gwrs_penalized_study_oracle_union(
    null_data, neighbors, config, control
  )
  expect_true(all(null_fit$coefficients[, -1L] == 0))
  expect_true(all(null_fit$diagnostics$converged))
})

test_that("oracle-local GWR respects exact location-specific support", {
  task <- gwrs_penalized_study_plan("smoke", 20260826L)[1L, , drop = FALSE]
  config <- list(kernel = "gaussian")
  data <- gwrs_penalized_study_simulate(task, config)
  neighbors <- gwr_neighbors(data$coords, k = 28L)
  target_neighbors <- gwr_neighbors(
    data$coords, query_coords = data$target_coords, k = 28L,
    include_self = FALSE
  )
  control <- gwrs_control(
    n_threads = 1L, keep_data = TRUE, diagnostics = "none"
  )
  fit <- gwrs_penalized_study_oracle_local(
    data, neighbors, target_neighbors, config, control
  )

  expect_true(all(fit$coefficients[, -1L][!data$support] == 0))
  expect_true(all(
    fit$target_prediction$coefficients[, -1L][!data$target_support] == 0
  ))
  expect_equal(length(fit$fitted.values), task$n[[1L]])
  expect_equal(length(fit$target_prediction$fit), nrow(data$x_target))
})

test_that("null-slope summaries retain FPR and undefined MCC", {
  task <- gwrs_penalized_study_plan("smoke", 20260826L)[1L, , drop = FALSE]
  task$scenario_id <- "null_slopes"
  config <- list(
    kernel = "gaussian",
    n_lambda = 4L,
    lambda_min_ratio = 1e-3,
    ridge_lambda = 10^seq(2, -4, length.out = 4L),
    en_alpha = 0.5,
    scad_gamma = 3.7,
    mcp_gamma = 3,
    selection_methods = c("lasso", "elastic_net", "scad", "mcp"),
    selection_tolerance = 1e-8,
    truth_tolerance = 1e-12
  )
  data <- gwrs_penalized_study_simulate(task, config)
  fold <- spatial_folds(data$coords, n_folds = 3L, seed = task$fold_seed)
  neighbors <- gwr_neighbors(data$coords, k = 28L)
  control <- gwrs_control(
    n_threads = 1L, keep_data = TRUE, diagnostics = "none"
  )
  target_neighbors <- gwr_neighbors(
    data$coords, query_coords = data$target_coords, k = 28L,
    include_self = FALSE
  )
  bandwidth_cv <- list(
    k = 28L,
    best = data.frame(k = 28L, rmse = 1, convergence_rate = 1,
                      elapsed_seconds = 0),
    table = data.frame(k = 28L, rmse = 1, convergence_rate = 1,
                       elapsed_seconds = 0)
  )
  methods <- c(
    "gwr", "ridge", "lasso", "elastic_net", "scad", "mcp",
    "oracle_union_gwr", "oracle_local_gwr"
  )
  summaries <- lapply(methods, function(method) {
    fitted <- gwrs_penalized_study_fit_method(
      method, data, neighbors, target_neighbors, 28L, fold, config, control
    )
    gwrs_penalized_study_summary(
      method, fitted, data, task, config, 28L, bandwidth_cv, 0,
      target_neighbors, 0
    )
  })
  summary <- do.call(rbind, summaries)
  selected <- summary$method %in% config$selection_methods

  expect_equal(summary$design_snr, rep(task$snr[[1L]], nrow(summary)))
  expect_true(all(summary$realized_train_noise_variance > 0))
  expect_true(all(summary$realized_test_noise_variance > 0))
  expect_true(all(summary$realized_target_noise_variance > 0))
  expect_true(all(summary$realized_train_snr > 0))
  expect_true(all(summary$realized_test_snr > 0))
  expect_true(all(summary$realized_target_snr > 0))
  expect_equal(
    summary$realized_target_snr,
    rep(
      stats::var(data$target_signal) /
        stats::var(data$y_target - data$target_signal),
      nrow(summary)
    )
  )
  expect_true(all(is.na(summary$active_coefficient_rmse)))
  expect_true(all(is.finite(summary$new_location_signal_rmse)))
  expect_true(all(is.na(summary$mcc[selected])))
  expect_true(all(is.finite(summary$false_positive_rate[selected])))
  expect_true(all(is.na(summary$mcc[!selected])))
  expect_true(all(is.na(summary$false_positive_rate[!selected])))
})
