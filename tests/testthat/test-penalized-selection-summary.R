summary_script_helper <- function() {
  installed <- system.file(
    "simulations", "summarize-penalized-selection-study.R",
    package = "gwrs"
  )
  if (nzchar(installed)) return(installed)
  candidates <- c(
    file.path(
      "inst", "simulations", "summarize-penalized-selection-study.R"
    ),
    file.path(
      "..", "..", "inst", "simulations",
      "summarize-penalized-selection-study.R"
    )
  )
  candidates[file.exists(candidates)][[1L]]
}

summary_script <- summary_script_helper()
sys.source(summary_script, envir = environment())
sys.source(
  file.path(dirname(summary_script), "penalized-study-utils.R"),
  envir = environment()
)

frozen_selection_config <- function() {
  list(
    schema = "gwrs-resumable-simulation-v1",
    study = gwrs_selection_summary_study,
    scope = "paper-level-spatial-cv-selection-and-prediction",
    mode = "production",
    design_version = "penalized-selection-v2",
    package_version = "0.4.0",
    n = 10000L, p = 100L, k = 1000L,
    task_count = 1900L, method_fit_count = 15200L,
    methods = gwrs_selection_summary_methods,
    selection_methods = gwrs_selection_summary_selection_methods,
    master_seed = 20260826L,
    kernel = "gaussian", n_folds = 5L, n_lambda = 30L,
    lambda_min_ratio = 1e-3,
    ridge_lambda = 10^seq(2, -4, length.out = 30L),
    ridge_lambda_definition = "10^seq(2,-4,length.out=n_lambda)",
    en_alpha = 0.5, scad_gamma = 3.7, mcp_gamma = 3,
    k_definition = paste(
      "max(2(p+1),ceil(.025n));",
      "max(4(p+1),ceil(.05n)); max(8(p+1),ceil(.10n))"
    ),
    bandwidth_selection = "common GWR RMSE on spatial folds",
    model_bandwidth_rule =
      "all eight methods use the common GWR-CV selected k",
    predictor_process = paste(
      "AR(1) cross-correlation followed by a sparse row-standardized",
      "8-nearest-neighbor SAR filter"
    ),
    predictor_cross_rho_values = c(0, 0.7),
    predictor_spatial_rho_values = c(0, 0.7),
    spatial_graph_neighbors = 8L,
    support_truth = "explicit location-by-predictor logical masks",
    target_grid_definition =
      "32x32 midpoint grid; restricted to the L-shaped domain when applicable",
    intercept_definition =
      "0.7+0.25*cos(2*pi*east)-0.20*sin(2*pi*north) in every scenario",
    snr_definition =
      "var(signal)/var(noise); sigma=sd(signal)/sqrt(SNR)",
    response_evaluation = paste(
      "independent same-location response plus independent fixed-grid",
      "new-location response and noise-free signal"
    ),
    oracle_scope = paste(
      "global-union and location-specific true supports; simulation only"
    ),
    selection_metric_scope = "pooled location-by-predictor support",
    primary_prediction_metric = "new_location_signal_rmse",
    primary_estimation_metric = "active_coefficient_rmse",
    primary_selection_metric = "mcc",
    selection_tolerance = 1e-8, truth_tolerance = 1e-12,
    solver_tolerance = 1e-7, max_iterations = 2000L,
    n_threads = 56L, grain_size = 16L,
    fit_keep_data_during_prediction = TRUE, keep_task_data = FALSE
  )
}

test_that("resource calibration cannot enter scientific summaries", {
  scientific <- frozen_selection_config()
  expect_silent(gwrs_selection_summary_validate_run_config(scientific))

  calibration <- scientific
  calibration$scope <-
    "resource-calibration-only-never-merge-with-scientific-results"
  calibration$mode <- "max_cell_calibration"
  expect_error(
    gwrs_selection_summary_validate_run_config(calibration),
    "diagnostics only"
  )

  changed_seed <- scientific
  changed_seed$master_seed <- 1L
  expect_error(
    gwrs_selection_summary_validate_run_config(changed_seed),
    "master_seed"
  )
})

test_that("paper analysis accepts only the exact frozen production plan", {
  plan <- gwrs_penalized_study_plan("production", 20260826L)
  expect_silent(gwrs_selection_summary_validate_plan_contract(plan))

  changed <- plan
  changed$predictor_rho[[1L]] <- 0.3
  expect_error(
    gwrs_selection_summary_validate_plan_contract(changed),
    "cell definitions"
  )
})

test_that("zero MCC misses remain in Monte Carlo summaries", {
  summary <- data.frame(
    cell_id = "non-null-cell",
    section = "main",
    scenario_id = "smooth_multiscale_local_sparse",
    design_variant = "baseline",
    sampling_design = "uniform",
    n = 100L,
    p = 6L,
    predictor_rho = 0.7,
    spatial_rho = 0.7,
    local_correlation = 0,
    target_grid_side = 32L,
    target_snr = 3,
    method = "lasso",
    stringsAsFactors = FALSE
  )[rep(1L, 2L), , drop = FALSE]
  for (metric in gwrs_selection_summary_error_metrics) summary[[metric]] <- 1
  for (metric in gwrs_selection_summary_descriptive_metrics) {
    summary[[metric]] <- 0
  }
  for (metric in gwrs_selection_summary_selection_metrics) {
    summary[[metric]] <- 0
  }
  summary$mcc <- c(0, 1)

  metrics <- gwrs_selection_summary_metric_table(summary)
  mcc <- metrics[metrics$metric == "mcc", , drop = FALSE]
  expect_equal(mcc$n_total, 2L)
  expect_equal(mcc$n_finite, 2L)
  expect_equal(mcc$mean, 0.5)
})

test_that("paper summaries retain legitimate NA metrics and paired MCC", {
  task_id <- rep(c("null-rep-00001", "null-rep-00002"), each = 8L)
  method <- rep(gwrs_selection_summary_methods, times = 2L)
  summary <- data.frame(
    task_id = task_id,
    cell_id = "null-cell",
    section = "reference",
    scenario_id = "null_slopes",
    design_variant = "baseline",
    sampling_design = "uniform",
    n = 100L,
    p = 6L,
    predictor_rho = 0.7,
    spatial_rho = 0.7,
    local_correlation = 0,
    target_grid_side = 32L,
    target_snr = 3,
    method = method,
    selection_evaluable = method %in%
      gwrs_selection_summary_selection_methods,
    all_converged = TRUE,
    stringsAsFactors = FALSE
  )
  method_offset <- match(method, gwrs_selection_summary_methods) / 100
  summary$realized_signal_variance <- 3
  summary$noise_sd <- 1
  summary$design_snr <- 3
  summary$realized_train_noise_variance <- 1
  summary$realized_test_noise_variance <- 1
  summary$realized_target_noise_variance <- 1
  summary$realized_train_snr <- 3
  summary$realized_test_snr <- 3
  summary$realized_target_snr <- 3
  for (metric in gwrs_selection_summary_error_metrics) {
    summary[[metric]] <- 1 + method_offset
  }
  active_metrics <- c(
    "active_coefficient_rmse", "core_active_coefficient_rmse",
    "domain_edge_active_coefficient_rmse",
    "domain_interior_active_coefficient_rmse",
    "new_location_active_coefficient_rmse"
  )
  for (metric in active_metrics) summary[[metric]] <- NA_real_
  summary$support_boundary_coefficient_rmse <- NA_real_
  summary$active_attenuation_bias <- NA_real_
  for (metric in gwrs_selection_summary_selection_metrics) {
    summary[[metric]] <- NA_real_
  }
  evaluable <- summary$selection_evaluable
  false_positive <- c(lasso = 10, elastic_net = 8, scad = 5, mcp = 4)
  summary$true_positive[evaluable] <- 0
  summary$false_positive[evaluable] <- false_positive[method[evaluable]]
  summary$false_negative[evaluable] <- 0
  summary$true_negative[evaluable] <- 600 -
    summary$false_positive[evaluable]
  summary$false_positive_rate[evaluable] <-
    summary$false_positive[evaluable] / 600
  summary$specificity[evaluable] <- 1 -
    summary$false_positive_rate[evaluable]
  summary$precision[evaluable] <- 0
  summary$false_discovery_rate[evaluable] <- 1
  summary$f1[evaluable] <- 0
  summary$support_iou[evaluable] <- 0
  summary$selection_accuracy[evaluable] <-
    summary$true_negative[evaluable] / 600
  summary$mean_selected_predictors[evaluable] <-
    summary$false_positive[evaluable] / 100

  expect_silent(gwrs_selection_summary_validate_metrics(summary))
  metrics <- gwrs_selection_summary_metric_table(summary)
  active <- metrics[
    metrics$method == "gwr" &
      metrics$metric == "active_coefficient_rmse", , drop = FALSE
  ]
  mcc <- metrics[
    metrics$method == "lasso" & metrics$metric == "mcc", , drop = FALSE
  ]
  expect_equal(active$n_total, 2L)
  expect_equal(active$n_finite, 0L)
  expect_true(is.na(active$mean))
  expect_equal(mcc$n_total, 2L)
  expect_equal(mcc$n_finite, 0L)

  paired <- gwrs_selection_summary_paired_table(summary)
  scad_fpr <- paired[
    paired$target_method == "scad" &
      paired$reference_method == "lasso" &
      paired$metric == "false_positive_rate", , drop = FALSE
  ]
  paired_mcc <- paired[
    paired$target_method == "scad" & paired$metric == "mcc", , drop = FALSE
  ]
  expect_equal(scad_fpr$n_total, 2L)
  expect_equal(scad_fpr$n_finite, 2L)
  expect_lt(scad_fpr$mean_difference, 0)
  expect_identical(scad_fpr$direction, "lower_is_better")
  expect_identical(scad_fpr$favors, "target")
  expect_equal(paired_mcc$n_total, 2L)
  expect_equal(paired_mcc$n_finite, 0L)
  expect_true(is.na(paired_mcc$mean_difference))

  direct_scad <- paired[
    paired$target_method == "scad" &
      paired$reference_method == "lasso" &
      paired$metric == "test_response_rmse", , drop = FALSE
  ]
  direct_mcp <- paired[
    paired$target_method == "mcp" &
      paired$reference_method == "lasso" &
      paired$metric == "coefficient_rmse", , drop = FALSE
  ]
  expect_equal(nrow(direct_scad), 1L)
  expect_equal(nrow(direct_mcp), 1L)

  sensitivity <- summary
  sensitivity$all_converged[
    sensitivity$method == "scad" &
      sensitivity$task_id == "null-rep-00001"
  ] <- FALSE
  converged_metrics <- gwrs_selection_summary_metric_table(
    sensitivity, "all_converged"
  )
  scad_error <- converged_metrics[
    converged_metrics$method == "scad" &
      converged_metrics$metric == "test_response_rmse", , drop = FALSE
  ]
  expect_identical(scad_error$analysis_set, "all_converged")
  expect_equal(scad_error$n_total, 1L)
  converged_pairs <- gwrs_selection_summary_paired_table(
    sensitivity, "all_converged"
  )
  scad_pair <- converged_pairs[
    converged_pairs$target_method == "scad" &
      converged_pairs$reference_method == "lasso" &
      converged_pairs$metric == "test_response_rmse", , drop = FALSE
  ]
  expect_equal(scad_pair$n_total, 1L)

  fractional <- summary
  fractional$true_positive[which(evaluable)[[1L]]] <- 0.5
  expect_error(
    gwrs_selection_summary_validate_metrics(fractional),
    "nonnegative integers"
  )
  inconsistent <- summary
  inconsistent$selection_accuracy[which(evaluable)[[1L]]] <- 0.5
  expect_error(
    gwrs_selection_summary_validate_metrics(inconsistent),
    "inconsistent"
  )
})
