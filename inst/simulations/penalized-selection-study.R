# Resumable paper-level comparison of penalized geographically weighted models.
# Production is a manual handoff and must not be launched by workflow setup.

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
if (length(file_argument) != 1L) {
  stop("Run this study with Rscript.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", file_argument))
script_directory <- dirname(script_path)
source(file.path(script_directory, "resume-utils.R"))
source(file.path(script_directory, "selection-metrics.R"))
source(file.path(script_directory, "penalized-study-utils.R"))

read_integer <- function(name, default, minimum = 1L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(raw))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop("`", name, "` must be an integer >= ", minimum, ".",
         call. = FALSE)
  }
  value
}

read_number <- function(name, default, lower = -Inf, upper = Inf) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.double(raw))
  if (length(value) != 1L || !is.finite(value) ||
      value < lower || value > upper) {
    stop(
      "`", name, "` must be finite and in [", lower, ", ", upper, "].",
      call. = FALSE
    )
  }
  value
}

read_flag <- function(name, default = FALSE) {
  raw <- Sys.getenv(name, unset = if (default) "YES" else "NO")
  if (!raw %in% c("YES", "NO")) {
    stop("`", name, "` must be YES or NO.", call. = FALSE)
  }
  identical(raw, "YES")
}

output_directory <- Sys.getenv("GWRS_SELECTION_SIM_OUTPUT")
if (!nzchar(output_directory)) {
  stop(
    "Set `GWRS_SELECTION_SIM_OUTPUT` to a new, versioned directory.",
    call. = FALSE
  )
}
action <- match.arg(
  Sys.getenv("GWRS_SELECTION_SIM_ACTION", unset = "run"),
  c("run", "status")
)
if (identical(action, "status")) {
  gwrs_sim_status(output_directory)
  quit(save = "no", status = 0L)
}

library(gwrs)

mode <- match.arg(
  Sys.getenv("GWRS_SELECTION_SIM_MODE", unset = "smoke"),
  c("smoke", "pilot", "max_cell_calibration", "production")
)
mode_defaults <- switch(
  mode,
  smoke = list(n_folds = 3L, n_lambda = 4L),
  pilot = list(n_folds = 5L, n_lambda = 10L),
  max_cell_calibration = list(n_folds = 5L, n_lambda = 30L),
  production = list(n_folds = 5L, n_lambda = 30L)
)
master_seed <- read_integer(
  "GWRS_SELECTION_SIM_SEED", 20260826L, minimum = 0L
)
n_threads <- read_integer("GWRS_SELECTION_SIM_THREADS", 1L)
kernel <- match.arg(
  Sys.getenv("GWRS_SELECTION_SIM_KERNEL", unset = "gaussian"),
  c("gaussian", "bisquare", "exponential", "tricube", "boxcar")
)
solver_tolerance <- read_number(
  "GWRS_SELECTION_SIM_SOLVER_TOLERANCE", 1e-7,
  lower = .Machine$double.eps
)
max_iterations <- read_integer(
  "GWRS_SELECTION_SIM_MAX_ITERATIONS", 2000L
)
grain_size <- read_integer("GWRS_SELECTION_SIM_GRAIN_SIZE", 16L)
selection_tolerance <- read_number(
  "GWRS_SELECTION_SIM_SELECTION_TOLERANCE", 1e-8, lower = 0
)
truth_tolerance <- read_number(
  "GWRS_SELECTION_SIM_TRUTH_TOLERANCE", 1e-12, lower = 0
)
keep_task_data <- read_flag("GWRS_SELECTION_SIM_KEEP_TASK_DATA")
confirm_large <- read_flag("GWRS_SELECTION_SIM_CONFIRM_LARGE")
retry_failed <- read_flag(
  "GWRS_SELECTION_SIM_RETRY_FAILED", default = TRUE
)
recover_lock <- read_flag("GWRS_SELECTION_SIM_RECOVER_LOCK")
max_tasks_raw <- Sys.getenv("GWRS_SELECTION_SIM_MAX_TASKS")
max_tasks <- if (nzchar(max_tasks_raw)) {
  read_integer("GWRS_SELECTION_SIM_MAX_TASKS", 1L)
} else {
  Inf
}
if (mode %in% c("max_cell_calibration", "production") && !confirm_large) {
  stop(
    "This mode requires `GWRS_SELECTION_SIM_CONFIRM_LARGE=YES`.",
    call. = FALSE
  )
}

methods <- c(
  "gwr", "ridge", "lasso", "elastic_net", "scad", "mcp",
  "oracle_union_gwr", "oracle_local_gwr"
)
selection_methods <- c("lasso", "elastic_net", "scad", "mcp")
plan <- gwrs_penalized_study_plan(mode, master_seed)
all_k <- sort(unique(unlist(strsplit(plan$k_candidates, ":", fixed = TRUE))))
all_k <- as.integer(all_k)
ridge_lambda <- 10^seq(2, -4, length.out = mode_defaults$n_lambda)

config <- list(
  schema = gwrs_sim_schema,
  study = "gwr-penalized-selection-paper-v2",
  scope = if (identical(mode, "max_cell_calibration")) {
    "resource-calibration-only-never-merge-with-scientific-results"
  } else {
    "paper-level-spatial-cv-selection-and-prediction"
  },
  mode = mode,
  design_version = "penalized-selection-v2",
  package_version = as.character(utils::packageVersion("gwrs")),
  n = max(plan$n),
  p = max(plan$p),
  k = max(all_k),
  task_count = nrow(plan),
  method_fit_count = nrow(plan) * length(methods),
  methods = methods,
  selection_methods = selection_methods,
  master_seed = master_seed,
  kernel = kernel,
  n_folds = mode_defaults$n_folds,
  n_lambda = mode_defaults$n_lambda,
  lambda_min_ratio = 1e-3,
  ridge_lambda = ridge_lambda,
  ridge_lambda_definition = "10^seq(2,-4,length.out=n_lambda)",
  en_alpha = 0.5,
  scad_gamma = 3.7,
  mcp_gamma = 3,
  k_definition = paste(
    "max(2(p+1),ceil(.025n));",
    "max(4(p+1),ceil(.05n)); max(8(p+1),ceil(.10n))"
  ),
  bandwidth_selection = "common GWR RMSE on spatial folds",
  model_bandwidth_rule = if (identical(mode, "max_cell_calibration")) {
    "bandwidth CV is audited; all eight resource fits force max(k_candidates)"
  } else {
    "all eight methods use the common GWR-CV selected k"
  },
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
  snr_definition = "var(signal)/var(noise); sigma=sd(signal)/sqrt(SNR)",
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
  selection_tolerance = selection_tolerance,
  truth_tolerance = truth_tolerance,
  solver_tolerance = solver_tolerance,
  max_iterations = max_iterations,
  n_threads = n_threads,
  grain_size = grain_size,
  fit_keep_data_during_prediction = TRUE,
  keep_task_data = keep_task_data
)

runtime_signature <- gwrs_sim_runtime_signature("gwrs")
workflow_files <- c(
  resume_utils = file.path(script_directory, "resume-utils.R"),
  selection_metrics = file.path(script_directory, "selection-metrics.R"),
  study_utils = file.path(script_directory, "penalized-study-utils.R"),
  study_script = script_path
)
runtime_signature <- rbind(
  runtime_signature,
  data.frame(
    field = paste0(names(workflow_files), "_sha256"),
    value = gwrs_sim_sha256(unname(workflow_files)),
    stringsAsFactors = FALSE
  )
)

task_function <- function(task, config, heartbeat) {
  heartbeat("simulate-shared-data")
  data <- gwrs_penalized_study_simulate(task, config)
  fold <- spatial_folds(
    data$coords,
    n_folds = config$n_folds,
    seed = task$fold_seed[[1L]]
  )
  control <- gwrs_control(
    tolerance = config$solver_tolerance,
    max_iterations = config$max_iterations,
    n_threads = config$n_threads,
    grain_size = config$grain_size,
    keep_data = TRUE,
    diagnostics = "none"
  )
  k_candidates <- as.integer(strsplit(
    task$k_candidates[[1L]], ":", fixed = TRUE
  )[[1L]])
  bandwidth_cv <- gwrs_penalized_study_select_k(
    data,
    k_candidates,
    fold,
    config,
    control,
    heartbeat
  )
  model_k <- if (identical(config$mode, "max_cell_calibration")) {
    max(k_candidates)
  } else {
    bandwidth_cv$k
  }
  heartbeat("build-selected-neighbors")
  neighbor_timing <- system.time({
    neighbors <- gwr_neighbors(data$coords, k = model_k)
  })
  neighbor_seconds <- unname(neighbor_timing[["elapsed"]])
  heartbeat("build-target-neighbors")
  target_neighbor_timing <- system.time({
    target_neighbors <- gwr_neighbors(
      data$coords,
      query_coords = data$target_coords,
      k = model_k,
      include_self = FALSE
    )
  })
  target_neighbor_seconds <- unname(target_neighbor_timing[["elapsed"]])

  tuning_results <- vector("list", length(config$methods))
  names(tuning_results) <- config$methods
  retained_coefficients <- if (isTRUE(config$keep_task_data)) {
    vector("list", length(config$methods))
  } else NULL
  if (!is.null(retained_coefficients)) {
    names(retained_coefficients) <- config$methods
  }
  summaries <- vector("list", length(config$methods))
  for (index in seq_along(config$methods)) {
    method <- config$methods[[index]]
    fitted <- gwrs_penalized_study_fit_method(
      method,
      data,
      neighbors,
      target_neighbors,
      model_k,
      fold,
      config,
      control,
      heartbeat
    )
    summaries[[index]] <- gwrs_penalized_study_summary(
      method,
      fitted,
      data,
      task,
      config,
      model_k,
      bandwidth_cv,
      neighbor_seconds,
      target_neighbors,
      target_neighbor_seconds
    )
    tuning_results[index] <- list(fitted$tuning)
    if (!is.null(retained_coefficients)) {
      retained_coefficients[[index]] <- fitted$fit$coefficients
    }
    rm(fitted)
  }
  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  result <- list(
    summary = summary,
    bandwidth_cv = bandwidth_cv$table,
    tuning = tuning_results
  )
  if (isTRUE(config$keep_task_data)) {
    result$task_data <- list(
      coords = data$coords,
      target_coords = data$target_coords,
      x = data$x,
      x_target = data$x_target,
      y_train = data$y_train,
      y_test = data$y_test,
      y_target = data$y_target,
      signal = data$signal,
      target_signal = data$target_signal,
      intercept = data$intercept,
      target_intercept = data$target_intercept,
      beta = data$beta,
      beta_target = data$beta_target,
      support = data$support,
      target_support = data$target_support,
      support_boundary = data$support_boundary,
      coefficients = retained_coefficients
    )
  }
  result
}

validate_result <- function(result, task) {
  fail <- function(reason) {
    message("Scientific checkpoint validation failed: ", reason)
    FALSE
  }
  if (!is.list(result) || !is.data.frame(result$summary) ||
      !is.data.frame(result$bandwidth_cv) || !is.list(result$tuning)) {
    return(fail("top-level result structure"))
  }
  summary <- result$summary
  if (nrow(summary) != length(config$methods) ||
      !identical(as.character(summary$method), config$methods) ||
      !identical(names(result$tuning), config$methods)) {
    return(fail("method rows or tuning names"))
  }
  task_matches <- all(summary$cell_id == task$cell_id[[1L]]) &&
    all(summary$section == task$section[[1L]]) &&
    all(summary$scenario_id == task$scenario_id[[1L]]) &&
    all(summary$design_variant == task$design_variant[[1L]]) &&
    all(summary$sampling_design == task$sampling_design[[1L]]) &&
    all(summary$replication == task$replication[[1L]]) &&
    all(summary$data_seed == task$data_seed[[1L]]) &&
    all(summary$fold_seed == task$fold_seed[[1L]]) &&
    all(summary$n == task$n[[1L]]) && all(summary$p == task$p[[1L]]) &&
    all(summary$spatial_rho == task$spatial_rho[[1L]]) &&
    all(summary$local_correlation == task$local_correlation[[1L]]) &&
    all(summary$target_grid_side == task$target_grid_side[[1L]])
  if (!task_matches) return(fail("task metadata"))

  candidates <- as.integer(strsplit(
    task$k_candidates[[1L]], ":", fixed = TRUE
  )[[1L]])
  bandwidth_valid <- identical(
    as.integer(result$bandwidth_cv$k), as.integer(candidates)
  ) && all(is.finite(result$bandwidth_cv$rmse)) &&
    all(is.finite(result$bandwidth_cv$convergence_rate)) &&
    all(is.finite(result$bandwidth_cv$elapsed_seconds)) &&
    length(unique(summary$bandwidth_selected_k)) == 1L &&
    summary$bandwidth_selected_k[[1L]] ==
      result$bandwidth_cv$k[which.min(result$bandwidth_cv$rmse)] &&
    length(unique(summary$chosen_k)) == 1L &&
    summary$chosen_k[[1L]] %in% candidates
  if (identical(config$mode, "max_cell_calibration")) {
    bandwidth_valid <- bandwidth_valid &&
      summary$chosen_k[[1L]] == max(candidates)
  } else {
    bandwidth_valid <- bandwidth_valid &&
      summary$chosen_k[[1L]] == summary$bandwidth_selected_k[[1L]]
  }
  if (!bandwidth_valid) return(fail("bandwidth contract"))

  paired_fields <- c(
    "signal_sum", "signal_sum_squares", "x_sum", "y_train_sum",
    "y_test_sum", "target_signal_sum", "target_signal_sum_squares",
    "x_target_sum", "y_target_sum", "support_count", "target_support_count",
    "noise_sd", "realized_signal_variance", "design_snr",
    "realized_train_noise_variance", "realized_test_noise_variance",
    "realized_target_noise_variance", "realized_train_snr",
    "realized_test_snr", "realized_target_snr",
    "realized_predictor_correlation", "realized_local_correlation",
    "realized_outside_correlation", "realized_spatial_neighbor_product",
    "bandwidth_cv_rmse", "bandwidth_selected_k"
  )
  if (any(vapply(
    summary[paired_fields],
    function(value) length(unique(value)) != 1L || !all(is.finite(value)),
    logical(1L)
  )) || any(abs(summary$design_snr - summary$target_snr) >
              1e-12 * pmax(1, abs(summary$target_snr)))) {
    return(fail("paired data signatures or SNR"))
  }

  finite_fields <- c(
    "predictor_rho", "spatial_rho", "local_correlation", "target_snr",
    "n_target", "bandwidth_selected_k", "chosen_k", "lambda",
    "convergence_rate", "train_response_rmse", "test_response_rmse",
    "signal_rmse", "new_location_response_rmse", "new_location_signal_rmse",
    "intercept_rmse", "coefficient_rmse", "inactive_coefficient_rmse",
    "domain_edge_signal_rmse", "domain_interior_signal_rmse",
    "new_location_coefficient_rmse",
    "new_location_inactive_coefficient_rmse", "bandwidth_elapsed_seconds",
    "neighbor_elapsed_seconds", "target_neighbor_elapsed_seconds",
    "tuning_elapsed_seconds", "fit_elapsed_seconds",
    "prediction_elapsed_seconds", "retained_fit_mib"
  )
  if (!all(is.finite(unlist(summary[finite_fields])))) {
    return(fail("required finite metrics"))
  }
  if (anyNA(summary$all_converged) || anyNA(summary$selection_evaluable)) {
    return(fail("missing convergence or selection flags"))
  }

  active_rmse_fields <- c(
    "active_coefficient_rmse", "core_active_coefficient_rmse",
    "active_attenuation_bias", "domain_edge_active_coefficient_rmse",
    "domain_interior_active_coefficient_rmse",
    "new_location_active_coefficient_rmse"
  )
  active_rmse_valid <- if (identical(task$scenario_id[[1L]], "null_slopes")) {
    all(vapply(summary[active_rmse_fields], function(value) all(is.na(value)),
               logical(1L)))
  } else {
    all(vapply(summary[active_rmse_fields], function(value) all(is.finite(value)),
               logical(1L)))
  }
  if (!active_rmse_valid) return(fail("active-coefficient metrics"))

  boundary_expected <- task$scenario_id[[1L]] %in% c(
    "smooth_multiscale_local_sparse", "regional_sparse_discontinuous",
    "local_collinearity_sparse"
  )
  if (boundary_expected) {
    if (!all(is.finite(summary$support_boundary_coefficient_rmse))) {
      return(fail("support-boundary metric"))
    }
  } else if (!all(is.na(summary$support_boundary_coefficient_rmse))) {
    return(fail("unexpected support-boundary metric"))
  }

  selection_fields <- c(
    "true_positive", "false_positive", "false_negative", "true_negative",
    "true_positive_rate", "false_positive_rate", "specificity", "precision",
    "false_discovery_rate", "f1", "mcc", "support_iou",
    "selection_accuracy", "mean_selected_predictors"
  )
  for (index in seq_len(nrow(summary))) {
    expected <- summary$method[[index]] %in% config$selection_methods
    if (!identical(summary$selection_evaluable[[index]], expected)) {
      return(fail(paste0("selection flag for ", summary$method[[index]])))
    }
    values <- unlist(summary[index, selection_fields, drop = FALSE])
    if (expected) {
      confusion <- values[c(
        "true_positive", "false_positive", "false_negative", "true_negative"
      )]
      required <- values[c("selection_accuracy", "mean_selected_predictors")]
      if (!all(is.finite(c(confusion, required))) || any(confusion < 0) ||
          any(abs(confusion - round(confusion)) > 1e-8) ||
          sum(confusion) != summary$n[[index]] * summary$p[[index]] ||
          !all(is.finite(values) | is.na(values))) {
        return(fail(paste0(
          "selection confusion metrics for ", summary$method[[index]]
        )))
      }
    } else if (!all(is.na(values))) {
      return(fail(paste0(
        "selection metrics on non-selection method ", summary$method[[index]]
      )))
    }
  }

  tuned <- summary$method %in% c(
    "ridge", "lasso", "elastic_net", "scad", "mcp"
  )
  if (!all(is.finite(summary$cv_rmse[tuned])) ||
      !all(is.finite(summary$cv_convergence_rate[tuned])) ||
      !all(is.na(summary$cv_rmse[!tuned])) ||
      !all(is.na(summary$cv_convergence_rate[!tuned]))) {
    return(fail("CV metric contract"))
  }
  nonconvex <- summary$method %in% c("scad", "mcp")
  if (!all(is.finite(summary$maximum_stationarity_violation[nonconvex])) ||
      !all(is.finite(summary$maximum_coordinate_gap[nonconvex]))) {
    return(fail("nonconvex solver audit"))
  }
  TRUE
}

status <- gwrs_sim_run(
  output_directory = output_directory,
  config = config,
  plan = plan,
  task_function = task_function,
  validate_result = validate_result,
  runtime_signature = runtime_signature,
  max_tasks = max_tasks,
  retry_failed = retry_failed,
  recover_lock = recover_lock,
  confirm_large = confirm_large,
  limits = list(n = 20000L, p = 100L, k = 1200L, tasks = 500L)
)
print(status$progress, row.names = FALSE)
cat("Output:", status$output_directory, "\n")
