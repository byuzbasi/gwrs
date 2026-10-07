# Validate and summarize a completed paper-level penalized GWR study.
#
# Usage:
#   Rscript --vanilla summarize-penalized-selection-study.R STUDY_DIR OUTPUT_DIR

gwrs_selection_summary_schema <- "gwrs-penalized-selection-summary-v3"
gwrs_selection_summary_study <- "gwr-penalized-selection-paper-v2"
gwrs_selection_summary_methods <- c(
  "gwr", "ridge", "lasso", "elastic_net", "scad", "mcp",
  "oracle_union_gwr", "oracle_local_gwr"
)
gwrs_selection_summary_selection_methods <- c(
  "lasso", "elastic_net", "scad", "mcp"
)
gwrs_selection_summary_expected_tasks <- 1900L
gwrs_selection_summary_expected_fits <- 15200L
gwrs_selection_summary_error_metrics <- c(
  "train_response_rmse", "test_response_rmse", "signal_rmse",
  "new_location_response_rmse", "new_location_signal_rmse",
  "intercept_rmse", "coefficient_rmse", "active_coefficient_rmse",
  "inactive_coefficient_rmse", "core_active_coefficient_rmse",
  "support_boundary_coefficient_rmse",
  "domain_edge_active_coefficient_rmse",
  "domain_interior_active_coefficient_rmse", "domain_edge_signal_rmse",
  "domain_interior_signal_rmse", "new_location_coefficient_rmse",
  "new_location_active_coefficient_rmse",
  "new_location_inactive_coefficient_rmse"
)
gwrs_selection_summary_descriptive_metrics <- c("active_attenuation_bias")
gwrs_selection_summary_selection_metrics <- c(
  "true_positive", "false_positive", "false_negative", "true_negative",
  "true_positive_rate", "false_positive_rate", "specificity", "precision",
  "false_discovery_rate", "f1", "mcc", "support_iou",
  "selection_accuracy", "mean_selected_predictors"
)
gwrs_selection_summary_cell_fields <- c(
  "cell_id", "section", "scenario_id", "design_variant",
  "sampling_design", "n", "p", "predictor_rho", "spatial_rho",
  "local_correlation", "target_grid_side", "target_snr"
)
gwrs_selection_summary_paired_fields <- c(
  "cell_id", "section", "scenario_id", "design_variant",
  "sampling_design", "replication", "data_seed", "fold_seed", "n", "p",
  "predictor_rho", "spatial_rho", "local_correlation", "target_grid_side",
  "n_target", "target_snr",
  "realized_signal_variance", "noise_sd", "design_snr",
  "realized_train_noise_variance", "realized_test_noise_variance",
  "realized_target_noise_variance", "realized_train_snr",
  "realized_test_snr", "realized_target_snr",
  "realized_predictor_correlation", "realized_local_correlation",
  "realized_outside_correlation", "realized_spatial_neighbor_product",
  "bandwidth_selected_k", "chosen_k",
  "bandwidth_cv_rmse", "bandwidth_elapsed_seconds",
  "neighbor_elapsed_seconds", "target_neighbor_elapsed_seconds",
  "signal_sum", "signal_sum_squares", "x_sum", "y_train_sum", "y_test_sum",
  "target_signal_sum", "target_signal_sum_squares", "x_target_sum",
  "y_target_sum", "support_count", "target_support_count"
)

gwrs_selection_summary_require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      label, " is missing required columns: ", paste(missing, collapse = ", "),
      ".", call. = FALSE
    )
  }
  invisible(data)
}

gwrs_selection_summary_read_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing ", label, ": `", path, "`.", call. = FALSE)
  }
  tryCatch(
    utils::read.csv(
      path, stringsAsFactors = FALSE, check.names = FALSE,
      na.strings = c("NA", "")
    ),
    error = function(error) {
      stop("Could not read ", label, ": ", conditionMessage(error),
           call. = FALSE)
    }
  )
}

gwrs_selection_summary_single_value <- function(data, field, label) {
  value <- data[[field]]
  if (anyNA(value) || length(unique(value)) != 1L) {
    stop(label, " does not have one nonmissing `", field, "` value.",
         call. = FALSE)
  }
  value[[1L]]
}

gwrs_selection_summary_stats <- function(value) {
  value <- as.double(value)
  finite <- is.finite(value)
  observed <- value[finite]
  n_finite <- length(observed)
  standard_deviation <- if (n_finite >= 2L) stats::sd(observed) else NA_real_
  quantiles <- if (n_finite) {
    as.double(stats::quantile(
      observed, probs = c(0.25, 0.75), names = FALSE, type = 7
    ))
  } else {
    c(NA_real_, NA_real_)
  }
  data.frame(
    n_total = length(value),
    n_finite = n_finite,
    mean = if (n_finite) mean(observed) else NA_real_,
    sd = standard_deviation,
    mcse = if (n_finite >= 2L) {
      standard_deviation / sqrt(n_finite)
    } else {
      NA_real_
    },
    median = if (n_finite) stats::median(observed) else NA_real_,
    q25 = quantiles[[1L]],
    q75 = quantiles[[2L]],
    stringsAsFactors = FALSE
  )
}

gwrs_selection_summary_group_indices <- function(summary) {
  split(
    seq_len(nrow(summary)),
    list(
      factor(summary$cell_id, levels = unique(summary$cell_id)),
      factor(
        summary$method, levels = gwrs_selection_summary_methods
      )
    ),
    drop = TRUE, lex.order = TRUE
  )
}

gwrs_selection_summary_group_metadata <- function(data) {
  values <- lapply(gwrs_selection_summary_cell_fields, function(field) {
    gwrs_selection_summary_single_value(data, field, "Cell/method group")
  })
  names(values) <- gwrs_selection_summary_cell_fields
  as.data.frame(values, stringsAsFactors = FALSE)
}

gwrs_selection_summary_metric_table <- function(
    summary,
    analysis_set = c("all_completed", "all_converged")) {
  analysis_set <- match.arg(analysis_set)
  metrics <- c(
    gwrs_selection_summary_error_metrics,
    gwrs_selection_summary_descriptive_metrics,
    gwrs_selection_summary_selection_metrics
  )
  groups <- gwrs_selection_summary_group_indices(summary)
  rows <- vector("list", length(groups) * length(metrics))
  output_index <- 0L
  for (indices in groups) {
    group <- summary[indices, , drop = FALSE]
    metadata <- gwrs_selection_summary_group_metadata(group)
    method <- gwrs_selection_summary_single_value(
      group, "method", "Cell/method group"
    )
    analyzed <- if (identical(analysis_set, "all_converged")) {
      group[as.logical(group$all_converged), , drop = FALSE]
    } else {
      group
    }
    for (metric in metrics) {
      output_index <- output_index + 1L
      rows[[output_index]] <- cbind(
        metadata,
        data.frame(
          analysis_set = analysis_set, method = method, metric = metric,
                   stringsAsFactors = FALSE),
        gwrs_selection_summary_stats(analyzed[[metric]])
      )
    }
  }
  output <- do.call(rbind, rows)
  method_order <- match(output$method, gwrs_selection_summary_methods)
  metric_order <- match(output$metric, metrics)
  output <- output[order(output$cell_id, method_order, metric_order), ,
                   drop = FALSE]
  rownames(output) <- NULL
  output
}

gwrs_selection_summary_pair_stats <- function(value) {
  result <- gwrs_selection_summary_stats(value)
  names(result)[names(result) == "mean"] <- "mean_difference"
  names(result)[names(result) == "sd"] <- "sd_difference"
  if (result$n_finite >= 2L) {
    critical <- stats::qt(0.975, df = result$n_finite - 1L)
    half_width <- critical * result$mcse
    result$ci_lower_95 <- result$mean_difference - half_width
    result$ci_upper_95 <- result$mean_difference + half_width
  } else {
    result$ci_lower_95 <- NA_real_
    result$ci_upper_95 <- NA_real_
  }
  result[c(
    "n_total", "n_finite", "mean_difference", "sd_difference", "mcse",
    "ci_lower_95", "ci_upper_95", "median", "q25", "q75"
  )]
}

gwrs_selection_summary_favors <- function(mean_difference,
                                            ci_lower_95,
                                            ci_upper_95,
                                            direction) {
  if (!is.finite(mean_difference) || !is.finite(ci_lower_95) ||
      !is.finite(ci_upper_95) || identical(direction, "descriptive")) {
    return("not_estimable")
  }
  if (ci_lower_95 <= 0 && ci_upper_95 >= 0) return("inconclusive")
  if (identical(direction, "lower_is_better")) {
    if (ci_upper_95 < 0) "target" else "reference"
  } else {
    if (ci_lower_95 > 0) "target" else "reference"
  }
}

gwrs_selection_summary_paired_table <- function(
    summary,
    analysis_set = c("all_completed", "all_converged")) {
  analysis_set <- match.arg(analysis_set)
  selection_directions <- c(
    true_positive = "higher_is_better",
    false_positive = "lower_is_better",
    false_negative = "lower_is_better",
    true_negative = "higher_is_better",
    true_positive_rate = "higher_is_better",
    false_positive_rate = "lower_is_better",
    specificity = "higher_is_better",
    precision = "higher_is_better",
    false_discovery_rate = "lower_is_better",
    f1 = "higher_is_better",
    mcc = "higher_is_better",
    support_iou = "higher_is_better",
    selection_accuracy = "higher_is_better",
    mean_selected_predictors = "descriptive"
  )
  error_vs_gwr <- expand.grid(
    target_method = setdiff(gwrs_selection_summary_methods, "gwr"),
    reference_method = "gwr",
    metric = gwrs_selection_summary_error_metrics,
    stringsAsFactors = FALSE
  )
  direct_nonconvex <- expand.grid(
    target_method = c("scad", "mcp"),
    reference_method = "lasso",
    metric = gwrs_selection_summary_error_metrics,
    stringsAsFactors = FALSE
  )
  selection_vs_lasso <- expand.grid(
    target_method = setdiff(
      gwrs_selection_summary_selection_methods, "lasso"
    ),
    reference_method = "lasso",
    metric = gwrs_selection_summary_selection_metrics,
    stringsAsFactors = FALSE
  )
  comparisons <- rbind(
    error_vs_gwr,
    direct_nonconvex,
    selection_vs_lasso
  )
  comparisons$direction <- ifelse(
    comparisons$metric %in% gwrs_selection_summary_error_metrics,
    "lower_is_better", unname(selection_directions[comparisons$metric])
  )

  cells <- split(
    seq_len(nrow(summary)),
    factor(summary$cell_id, levels = unique(summary$cell_id)),
    drop = TRUE
  )
  rows <- vector("list", length(cells) * nrow(comparisons))
  output_index <- 0L
  for (indices in cells) {
    cell <- summary[indices, , drop = FALSE]
    metadata <- gwrs_selection_summary_group_metadata(cell)
    for (comparison_index in seq_len(nrow(comparisons))) {
      comparison <- comparisons[comparison_index, , drop = FALSE]
      reference <- cell[
        cell$method == comparison$reference_method, , drop = FALSE
      ]
      target <- cell[
        cell$method == comparison$target_method, , drop = FALSE
      ]
      reference <- reference[order(reference$task_id), , drop = FALSE]
      target <- target[order(target$task_id), , drop = FALSE]
      if (!identical(reference$task_id, target$task_id)) {
        stop("Paired method rows are not aligned by task_id.",
             call. = FALSE)
      }
      include <- rep.int(TRUE, nrow(target))
      if (identical(analysis_set, "all_converged")) {
        include <- as.logical(target$all_converged) &
          as.logical(reference$all_converged)
      }
      difference <- target[[comparison$metric]][include] -
        reference[[comparison$metric]][include]
      statistics <- gwrs_selection_summary_pair_stats(difference)
      output_index <- output_index + 1L
      rows[[output_index]] <- cbind(
        metadata,
        data.frame(
          analysis_set = analysis_set,
          metric = comparison$metric,
          target_method = comparison$target_method,
          reference_method = comparison$reference_method,
          difference_definition = "target - reference",
          direction = comparison$direction,
          favors = gwrs_selection_summary_favors(
            statistics$mean_difference,
            statistics$ci_lower_95,
            statistics$ci_upper_95,
            comparison$direction
          ),
          stringsAsFactors = FALSE
        ),
        statistics
      )
    }
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

gwrs_selection_summary_finite_fields <- function(data, fields) {
  output <- list()
  for (field in fields) {
    value <- as.double(data[[field]])
    finite <- value[is.finite(value)]
    prefix <- paste0(field, "_")
    output[[paste0(prefix, "n_finite")]] <- length(finite)
    output[[paste0(prefix, "mean")]] <- if (length(finite)) {
      mean(finite)
    } else {
      NA_real_
    }
    output[[paste0(prefix, "sd")]] <- if (length(finite) >= 2L) {
      stats::sd(finite)
    } else {
      NA_real_
    }
    output[[paste0(prefix, "maximum")]] <- if (length(finite)) {
      max(finite)
    } else {
      NA_real_
    }
  }
  as.data.frame(output, stringsAsFactors = FALSE)
}

gwrs_selection_summary_convergence_table <- function(summary) {
  fields <- c(
    "convergence_rate", "cv_convergence_rate", "target_convergence_rate",
    "maximum_stationarity_violation", "maximum_coordinate_gap",
    "bandwidth_elapsed_seconds", "neighbor_elapsed_seconds",
    "target_neighbor_elapsed_seconds", "tuning_elapsed_seconds",
    "fit_elapsed_seconds", "prediction_elapsed_seconds", "retained_fit_mib"
  )
  groups <- gwrs_selection_summary_group_indices(summary)
  rows <- lapply(groups, function(indices) {
    group <- summary[indices, , drop = FALSE]
    metadata <- gwrs_selection_summary_group_metadata(group)
    all_converged <- as.logical(group$all_converged)
    cbind(
      metadata,
      data.frame(
        method = gwrs_selection_summary_single_value(
          group, "method", "Cell/method group"
        ),
        n_total = nrow(group),
        n_all_converged = sum(all_converged),
        proportion_all_converged = mean(all_converged),
        stringsAsFactors = FALSE
      ),
      gwrs_selection_summary_finite_fields(group, fields)
    )
  })
  output <- do.call(rbind, rows)
  output <- output[order(
    output$cell_id,
    match(output$method, gwrs_selection_summary_methods)
  ), , drop = FALSE]
  rownames(output) <- NULL
  output
}

gwrs_selection_summary_bandwidth_table <- function(summary) {
  common <- summary[summary$method == "gwr", , drop = FALSE]
  cells <- split(
    seq_len(nrow(common)),
    factor(common$cell_id, levels = unique(common$cell_id)),
    drop = TRUE
  )
  rows <- list()
  output_index <- 0L
  for (indices in cells) {
    cell <- common[indices, , drop = FALSE]
    metadata <- gwrs_selection_summary_group_metadata(cell)
    choices <- sort(unique(cell$chosen_k))
    for (chosen_k in choices) {
      selected <- cell[cell$chosen_k == chosen_k, , drop = FALSE]
      output_index <- output_index + 1L
      cv <- gwrs_selection_summary_stats(selected$bandwidth_cv_rmse)
      rows[[output_index]] <- cbind(
        metadata,
        data.frame(
          chosen_k = chosen_k,
          n_tasks = nrow(selected),
          proportion_of_cell = nrow(selected) / nrow(cell),
          mean_bandwidth_cv_rmse = cv$mean,
          sd_bandwidth_cv_rmse = cv$sd,
          mcse_bandwidth_cv_rmse = cv$mcse,
          mean_bandwidth_elapsed_seconds = mean(
            selected$bandwidth_elapsed_seconds
          ),
          mean_neighbor_elapsed_seconds = mean(
            selected$neighbor_elapsed_seconds
          ),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  output <- do.call(rbind, rows)
  output <- output[order(output$cell_id, output$chosen_k), , drop = FALSE]
  rownames(output) <- NULL
  output
}

gwrs_selection_summary_design_table <- function(plan, summary) {
  cells <- split(
    seq_len(nrow(plan)),
    factor(plan$cell_id, levels = unique(plan$cell_id)),
    drop = TRUE
  )
  rows <- lapply(cells, function(indices) {
    cell <- plan[indices, , drop = FALSE]
    observed <- summary[summary$cell_id == cell$cell_id[[1L]], , drop = FALSE]
    data.frame(
      cell_id = gwrs_selection_summary_single_value(cell, "cell_id", "Plan cell"),
      section = gwrs_selection_summary_single_value(cell, "section", "Plan cell"),
      scenario_id = gwrs_selection_summary_single_value(
        cell, "scenario_id", "Plan cell"
      ),
      design_variant = gwrs_selection_summary_single_value(
        cell, "design_variant", "Plan cell"
      ),
      sampling_design = gwrs_selection_summary_single_value(
        cell, "sampling_design", "Plan cell"
      ),
      n = gwrs_selection_summary_single_value(cell, "n", "Plan cell"),
      p = gwrs_selection_summary_single_value(cell, "p", "Plan cell"),
      predictor_rho = gwrs_selection_summary_single_value(
        cell, "predictor_rho", "Plan cell"
      ),
      spatial_rho = gwrs_selection_summary_single_value(
        cell, "spatial_rho", "Plan cell"
      ),
      local_correlation = gwrs_selection_summary_single_value(
        cell, "local_correlation", "Plan cell"
      ),
      target_grid_side = gwrs_selection_summary_single_value(
        cell, "target_grid_side", "Plan cell"
      ),
      target_snr = gwrs_selection_summary_single_value(
        cell, "snr", "Plan cell"
      ),
      k_candidates = gwrs_selection_summary_single_value(
        cell, "k_candidates", "Plan cell"
      ),
      replication_min = min(cell$replication),
      replication_max = max(cell$replication),
      planned_tasks = nrow(cell),
      completed_tasks = length(unique(observed$task_id)),
      expected_method_count = length(gwrs_selection_summary_methods),
      observed_fit_rows = nrow(observed),
      expected_fit_rows = nrow(cell) * length(gwrs_selection_summary_methods),
      methods = paste(gwrs_selection_summary_methods, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

gwrs_selection_summary_validate_metrics <- function(summary) {
  all_converged <- as.logical(summary$all_converged)
  if (anyNA(all_converged)) {
    stop("Convergence flags contain missing values.", call. = FALSE)
  }
  snr_fields <- c(
    "realized_signal_variance", "noise_sd", "design_snr",
    "realized_train_noise_variance", "realized_test_noise_variance",
    "realized_target_noise_variance", "realized_train_snr",
    "realized_test_snr", "realized_target_snr"
  )
  snr_values <- as.matrix(summary[snr_fields])
  storage.mode(snr_values) <- "double"
  if (!all(is.finite(snr_values)) || any(snr_values <= 0) ||
      any(abs(summary$design_snr - summary$target_snr) >
            1e-12 * pmax(1, abs(summary$target_snr)))) {
    stop("Design and empirical SNR diagnostics are invalid.",
         call. = FALSE)
  }
  error_metrics <- gwrs_selection_summary_error_metrics
  active_metrics <- c(
    "active_coefficient_rmse", "core_active_coefficient_rmse",
    "domain_edge_active_coefficient_rmse",
    "domain_interior_active_coefficient_rmse",
    "new_location_active_coefficient_rmse"
  )
  conditional_metrics <- c(active_metrics, "support_boundary_coefficient_rmse")
  for (metric in setdiff(error_metrics, conditional_metrics)) {
    if (!all(is.finite(summary[[metric]]))) {
      stop("Unexpected nonfinite values in `", metric, "`.", call. = FALSE)
    }
  }
  null <- summary$scenario_id == "null_slopes"
  active_valid <- all(vapply(active_metrics, function(metric) {
    all(is.na(summary[[metric]][null])) &&
      all(is.finite(summary[[metric]][!null]))
  }, logical(1L)))
  bias_valid <- all(is.na(summary$active_attenuation_bias[null])) &&
    all(is.finite(summary$active_attenuation_bias[!null]))
  if (!active_valid || !bias_valid) {
    stop(
      "Active-coefficient metrics must be NA only for null-slope scenarios.",
      call. = FALSE
    )
  }
  boundary_scenario <- summary$scenario_id %in% c(
    "smooth_multiscale_local_sparse", "regional_sparse_discontinuous",
    "local_collinearity_sparse"
  )
  if (!all(is.finite(
    summary$support_boundary_coefficient_rmse[boundary_scenario]
  )) || !all(is.na(
    summary$support_boundary_coefficient_rmse[!boundary_scenario]
  ))) {
    stop("Support-boundary RMSE does not match the frozen scenarios.",
         call. = FALSE)
  }

  evaluable <- summary$method %in%
    gwrs_selection_summary_selection_methods
  if (!identical(as.logical(summary$selection_evaluable), evaluable)) {
    stop("Selection-evaluable flags do not match the frozen methods.",
         call. = FALSE)
  }
  selection <- summary[gwrs_selection_summary_selection_metrics]
  if (!all(vapply(selection[!evaluable, , drop = FALSE], function(value) {
    all(is.na(value))
  }, logical(1L)))) {
    stop("Non-selection methods contain selection metrics.", call. = FALSE)
  }
  count_fields <- c(
    "true_positive", "false_positive", "false_negative", "true_negative"
  )
  counts <- as.matrix(summary[evaluable, count_fields, drop = FALSE])
  storage.mode(counts) <- "double"
  dimensions <- as.double(summary$n[evaluable]) *
    as.double(summary$p[evaluable])
  if (!all(is.finite(counts)) || any(counts < 0) ||
      any(abs(counts - round(counts)) > 1e-8) ||
      any(!is.finite(dimensions)) || any(dimensions <= 0) ||
      any(abs(dimensions - round(dimensions)) > 1e-8) ||
      any(rowSums(counts) != dimensions)) {
    stop(
      "Selection confusion counts must be nonnegative integers summing to n*p.",
      call. = FALSE
    )
  }

  tp <- counts[, "true_positive"]
  fp <- counts[, "false_positive"]
  fn <- counts[, "false_negative"]
  tn <- counts[, "true_negative"]
  safe_ratio <- function(numerator, denominator) {
    output <- rep.int(NA_real_, length(denominator))
    defined <- is.finite(denominator) & denominator != 0
    output[defined] <- numerator[defined] / denominator[defined]
    output
  }
  predicted_positive <- tp + fp
  predicted_negative <- tn + fn
  truth_positive <- tp + fn
  truth_negative <- tn + fp
  mcc_denominator <- sqrt(
    predicted_positive * predicted_negative * truth_positive * truth_negative
  )
  expected_mcc <- rep.int(NA_real_, length(tp))
  constant_prediction <- truth_positive > 0 & truth_negative > 0 &
    (predicted_positive == 0 | predicted_negative == 0)
  expected_mcc[constant_prediction] <- 0
  regular_mcc <- is.finite(mcc_denominator) & mcc_denominator > 0
  expected_mcc[regular_mcc] <-
    (tp[regular_mcc] * tn[regular_mcc] -
       fp[regular_mcc] * fn[regular_mcc]) / mcc_denominator[regular_mcc]
  expected <- list(
    true_positive_rate = safe_ratio(tp, truth_positive),
    false_positive_rate = safe_ratio(fp, truth_negative),
    specificity = safe_ratio(tn, truth_negative),
    precision = safe_ratio(tp, predicted_positive),
    false_discovery_rate = safe_ratio(fp, predicted_positive),
    f1 = safe_ratio(2 * tp, 2 * tp + fp + fn),
    mcc = expected_mcc,
    support_iou = safe_ratio(tp, tp + fp + fn),
    selection_accuracy = (tp + tn) / dimensions,
    mean_selected_predictors = predicted_positive /
      as.double(summary$n[evaluable])
  )
  equal_metric <- function(recorded, reconstructed) {
    recorded <- as.double(recorded)
    same_missing <- identical(
      unname(is.na(recorded)), unname(is.na(reconstructed))
    )
    finite <- is.finite(reconstructed)
    same_missing &&
      all(is.finite(recorded[finite])) &&
      all(abs(recorded[finite] - reconstructed[finite]) <=
            1e-10 * pmax(1, abs(reconstructed[finite])))
  }
  for (metric in names(expected)) {
    if (!equal_metric(summary[[metric]][evaluable], expected[[metric]])) {
      stop(
        "Selection metric `", metric,
        "` is inconsistent with its confusion counts.", call. = FALSE
      )
    }
  }
  invisible(summary)
}

gwrs_selection_summary_validate_pairing <- function(summary, plan) {
  expected <- gwrs_selection_summary_methods
  if (anyNA(summary$task_id) || anyNA(summary$method) ||
      !setequal(unique(summary$task_id), plan$task_id) ||
      !setequal(unique(summary$method), expected)) {
    stop("Task IDs or frozen method names do not match the plan.",
         call. = FALSE)
  }
  method_counts <- table(
    factor(summary$task_id, levels = plan$task_id),
    factor(summary$method, levels = expected)
  )
  if (!all(method_counts == 1L)) {
    stop("Every task_id must contain each of the eight frozen methods once.",
         call. = FALSE)
  }

  plan_index <- match(summary$task_id, plan$task_id)
  mappings <- c(
    cell_id = "cell_id", section = "section", scenario_id = "scenario_id",
    design_variant = "design_variant", sampling_design = "sampling_design",
    replication = "replication", n = "n", p = "p",
    predictor_rho = "predictor_rho", spatial_rho = "spatial_rho",
    local_correlation = "local_correlation",
    target_grid_side = "target_grid_side", target_snr = "snr",
    data_seed = "data_seed", fold_seed = "fold_seed", seed = "seed"
  )
  for (summary_field in names(mappings)) {
    plan_field <- mappings[[summary_field]]
    left <- summary[[summary_field]]
    right <- plan[[plan_field]][plan_index]
    if (anyNA(left) || anyNA(right) || !all(left == right)) {
      stop("Task summary does not match the plan for `", summary_field,
           "`.", call. = FALSE)
    }
  }

  tasks <- split(
    seq_len(nrow(summary)),
    factor(summary$task_id, levels = plan$task_id),
    drop = TRUE
  )
  for (indices in tasks) {
    task <- summary[indices, , drop = FALSE]
    for (field in gwrs_selection_summary_paired_fields) {
      value <- task[[field]]
      if (anyNA(value) || length(unique(value)) != 1L) {
        stop(
          "Paired data signature differs across methods for task `",
          task$task_id[[1L]], "` in field `", field, "`.", call. = FALSE
        )
      }
    }
  }
  invisible(summary)
}

gwrs_selection_summary_validate_run_config <- function(config) {
  if (!is.list(config) ||
      !identical(as.character(config$study), gwrs_selection_summary_study)) {
    stop("The source is not `", gwrs_selection_summary_study, "`.",
         call. = FALSE)
  }
  if (!identical(
    as.character(config$scope),
    "paper-level-spatial-cv-selection-and-prediction"
  ) || !identical(as.character(config$mode), "production")) {
    stop(
      "Paper summaries require the frozen production scope/mode; smoke, ",
      "pilot, and resource calibration outputs are diagnostics only.",
      call. = FALSE
    )
  }
  recorded_methods <- as.character(unlist(config$methods, use.names = FALSE))
  if (!identical(recorded_methods, gwrs_selection_summary_methods)) {
    stop("The run configuration does not contain the eight frozen methods.",
         call. = FALSE)
  }
  recorded_selection <- as.character(unlist(
    config$selection_methods, use.names = FALSE
  ))
  if (!identical(
    recorded_selection, gwrs_selection_summary_selection_methods
  )) {
    stop("The run configuration has different selection methods.",
         call. = FALSE)
  }

  expected <- list(
    schema = "gwrs-resumable-simulation-v1",
    design_version = "penalized-selection-v2",
    package_version = "0.4.0",
    n = 10000L,
    p = 100L,
    k = 1000L,
    task_count = gwrs_selection_summary_expected_tasks,
    method_fit_count = gwrs_selection_summary_expected_fits,
    master_seed = 20260826L,
    kernel = "gaussian",
    n_folds = 5L,
    n_lambda = 30L,
    lambda_min_ratio = 1e-3,
    ridge_lambda = 10^seq(2, -4, length.out = 30L),
    ridge_lambda_definition = "10^seq(2,-4,length.out=n_lambda)",
    en_alpha = 0.5,
    scad_gamma = 3.7,
    mcp_gamma = 3,
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
    selection_tolerance = 1e-8,
    truth_tolerance = 1e-12,
    solver_tolerance = 1e-7,
    max_iterations = 2000L,
    n_threads = 56L,
    grain_size = 16L,
    fit_keep_data_during_prediction = TRUE,
    keep_task_data = FALSE
  )
  missing <- setdiff(names(expected), names(config))
  if (length(missing)) {
    stop(
      "The production configuration is missing frozen fields: ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  same_value <- function(recorded, frozen) {
    if (is.character(frozen)) {
      identical(as.character(recorded), frozen)
    } else if (is.logical(frozen)) {
      identical(as.logical(recorded), frozen)
    } else {
      isTRUE(all.equal(
        as.double(recorded), as.double(frozen), tolerance = 0,
        check.attributes = FALSE
      ))
    }
  }
  mismatched <- names(expected)[!vapply(
    names(expected),
    function(field) same_value(config[[field]], expected[[field]]),
    logical(1L)
  )]
  if (length(mismatched)) {
    stop(
      "The source differs from the frozen production protocol in: ",
      paste(mismatched, collapse = ", "), ".", call. = FALSE
    )
  }
  invisible(config)
}

gwrs_selection_summary_validate_plan_contract <- function(plan) {
  if (nrow(plan) != gwrs_selection_summary_expected_tasks ||
      length(unique(plan$cell_id)) != 39L) {
    stop("The production plan must contain 1,900 tasks in 39 cells.",
         call. = FALSE)
  }
  core <- expand.grid(
    scenario_id = c(
      "smooth_multiscale_local_sparse",
      "regional_sparse_discontinuous"
    ),
    n = c(2000L, 10000L),
    predictor_rho = c(0, 0.7),
    spatial_rho = c(0, 0.7),
    snr = c(1L, 3L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  core$section <- "main"
  core$design_variant <- "baseline"
  core$sampling_design <- "uniform"
  core$p <- 30L
  core$local_correlation <- 0
  core$target_grid_side <- 32L
  reference <- data.frame(
    scenario_id = c("stationary_sparse", "null_slopes"),
    design_variant = "baseline", sampling_design = "uniform",
    n = 5000L, predictor_rho = 0.7, spatial_rho = 0.7,
    local_correlation = 0, snr = 3L, target_grid_side = 32L,
    section = "reference", p = 30L,
    stringsAsFactors = FALSE
  )
  stress <- data.frame(
    scenario_id = c(
      "smooth_multiscale_local_sparse",
      "regional_sparse_discontinuous",
      "local_collinearity_sparse",
      "smooth_multiscale_local_sparse",
      "regional_sparse_discontinuous"
    ),
    design_variant = c(
      "p100", "p100", "local_collinearity", "clustered_sampling",
      "l_domain_boundary"
    ),
    sampling_design = c(
      "uniform", "uniform", "uniform", "clustered", "l_shaped"
    ),
    n = 10000L, p = c(100L, 100L, 30L, 30L, 30L),
    predictor_rho = c(0.7, 0.7, 0, 0.7, 0.7), spatial_rho = 0.7,
    local_correlation = c(0, 0, 0.9, 0, 0), snr = 3L,
    target_grid_side = 32L, section = "stress",
    stringsAsFactors = FALSE
  )
  expected_cells <- rbind(
    core[c(
      "scenario_id", "design_variant", "sampling_design", "n", "p",
      "predictor_rho", "spatial_rho", "local_correlation", "snr",
      "target_grid_side", "section"
    )],
    reference,
    stress
  )
  cell_fields <- c(
    "section", "scenario_id", "design_variant", "sampling_design", "n", "p",
    "predictor_rho", "spatial_rho", "local_correlation", "snr",
    "target_grid_side"
  )
  expected_cells <- expected_cells[cell_fields]
  observed_cells <- unique(plan[cell_fields])
  row_key <- function(value) do.call(paste, c(value, sep = "|"))
  observed_cells <- observed_cells[order(row_key(observed_cells)), , drop = FALSE]
  expected_cells <- expected_cells[order(row_key(expected_cells)), , drop = FALSE]
  rownames(observed_cells) <- NULL
  rownames(expected_cells) <- NULL
  if (!isTRUE(all.equal(
    observed_cells, expected_cells, tolerance = 0, check.attributes = FALSE
  ))) {
    stop("The 39 production-cell definitions differ from the protocol.",
         call. = FALSE)
  }
  expected_replications <- ifelse(
    plan$section == "main", 50L,
    ifelse(plan$section == "reference", 100L, 20L)
  )
  cells <- split(seq_len(nrow(plan)), plan$cell_id)
  valid_replications <- vapply(cells, function(indices) {
    expected <- unique(expected_replications[indices])
    length(expected) == 1L &&
      identical(as.integer(plan$replication[indices]), seq_len(expected))
  }, logical(1L))
  allowed_sections <- c(main = 1600L, reference = 200L, stress = 100L)
  section_counts <- table(factor(plan$section, levels = names(allowed_sections)))
  if (!all(valid_replications) || anyNA(plan$replication) ||
      !identical(as.integer(section_counts), unname(allowed_sections))) {
    stop("The production sections or replication ranges are not frozen.",
         call. = FALSE)
  }

  cross_code <- sprintf("%03d", round(100 * plan$predictor_rho))
  spatial_code <- sprintf("%03d", round(100 * plan$spatial_rho))
  snr_code <- sprintf("%02d", round(10 * plan$snr))
  expected_cell <- sprintf(
    "%s-%s-%s-n%05d-p%03d-rx%s-rs%s-s%s",
    plan$section, plan$scenario_id, plan$design_variant, plan$n, plan$p,
    cross_code, spatial_code, snr_code
  )
  expected_task <- sprintf(
    "%s-rep-%05d", expected_cell, plan$replication
  )
  task_index <- seq_len(nrow(plan))
  expected_seed <- as.integer(
    ((20260826L + task_index * 104729) %% 2147483646) + 1
  )
  expected_data_seed <- as.integer(
    ((20260826L + task_index * 130363) %% 2147483646) + 1
  )
  expected_fold_seed <- as.integer(
    ((20260826L + task_index * 169087) %% 2147483646) + 1
  )
  expected_k <- vapply(seq_len(nrow(plan)), function(index) {
    n <- as.integer(plan$n[[index]])
    p <- as.integer(plan$p[[index]])
    paste(sort(unique(as.integer(c(
      max(2L * (p + 1L), ceiling(0.025 * n)),
      max(4L * (p + 1L), ceiling(0.050 * n)),
      max(8L * (p + 1L), ceiling(0.100 * n))
    )))), collapse = ":")
  }, character(1L))
  exact <- identical(as.character(plan$cell_id), expected_cell) &&
    identical(as.character(plan$task_id), expected_task) &&
    identical(as.integer(plan$seed), expected_seed) &&
    identical(as.integer(plan$data_seed), expected_data_seed) &&
    identical(as.integer(plan$fold_seed), expected_fold_seed) &&
    identical(as.character(plan$k_candidates), expected_k)
  if (!exact) {
    stop("Task IDs, frozen seeds, or k candidates differ from the protocol.",
         call. = FALSE)
  }
  invisible(plan)
}

gwrs_selection_summary_validate_source <- function(study_directory) {
  required_files <- c(
    "run-config.rds", "run-config.csv", "task-plan.csv",
    "runtime-signature.csv", "scientific-signature.csv",
    "task-manifest.csv", "task-summary.csv", "manifest-sha256.csv",
    "COMPLETED"
  )
  missing <- required_files[
    !file.exists(file.path(study_directory, required_files))
  ]
  if (length(missing)) {
    stop("Completed study is missing: ", paste(missing, collapse = ", "),
         ".", call. = FALSE)
  }
  status <- gwrs_sim_status(study_directory, quiet = TRUE)
  if (!isTRUE(status$finalized) || !isTRUE(status$manifest_valid) ||
      !isTRUE(status$completed) || !isTRUE(status$completion_valid) ||
      !all(status$tasks$status == "completed")) {
    stop(
      "The study must have valid completion, task, scientific, and SHA-256 ",
      "manifests before analysis.", call. = FALSE
    )
  }
  if (!gwrs_sim_verify_scientific_signature(study_directory) ||
      !gwrs_sim_verify_manifest(study_directory) ||
      !gwrs_sim_verify_completion(study_directory)) {
    stop("Independent source-manifest verification failed.", call. = FALSE)
  }

  config <- readRDS(file.path(study_directory, "run-config.rds"))
  gwrs_selection_summary_validate_run_config(config)

  plan <- gwrs_selection_summary_read_csv(
    file.path(study_directory, "task-plan.csv"), "task plan"
  )
  manifest <- gwrs_selection_summary_read_csv(
    file.path(study_directory, "task-manifest.csv"), "task manifest"
  )
  summary <- gwrs_selection_summary_read_csv(
    file.path(study_directory, "task-summary.csv"), "task summary"
  )
  gwrs_selection_summary_require_columns(
    plan,
    c(
      "task_id", "cell_id", "section", "scenario_id", "n", "p",
      "design_variant", "sampling_design", "predictor_rho", "spatial_rho",
      "local_correlation", "snr", "target_grid_side", "replication",
      "k_candidates", "data_seed", "fold_seed", "seed"
    ),
    "Task plan"
  )
  gwrs_selection_summary_require_columns(
    manifest, c("task_id", "status", "bytes", "sha256"), "Task manifest"
  )
  gwrs_selection_summary_require_columns(
    summary,
    unique(c(
      "task_id", "seed", "method", "selection_evaluable", "all_converged",
      gwrs_selection_summary_cell_fields,
      gwrs_selection_summary_paired_fields,
      gwrs_selection_summary_error_metrics,
      gwrs_selection_summary_descriptive_metrics,
      gwrs_selection_summary_selection_metrics,
      "cv_convergence_rate", "convergence_rate",
      "target_convergence_rate", "maximum_stationarity_violation",
      "maximum_coordinate_gap", "tuning_elapsed_seconds",
      "fit_elapsed_seconds", "prediction_elapsed_seconds", "retained_fit_mib"
    )),
    "Task summary"
  )
  gwrs_selection_summary_validate_plan_contract(plan)
  if (anyDuplicated(plan$task_id) || anyDuplicated(manifest$task_id) ||
      !setequal(plan$task_id, manifest$task_id) ||
      !all(manifest$status == "completed")) {
    stop("Task plan and completed task manifest do not agree.",
         call. = FALSE)
  }
  scanned <- status$tasks[match(plan$task_id, status$tasks$task_id), ,
                          drop = FALSE]
  recorded <- manifest[match(plan$task_id, manifest$task_id), , drop = FALSE]
  if (anyNA(scanned$task_id) ||
      !identical(as.character(recorded$sha256), as.character(scanned$sha256)) ||
      !identical(as.double(recorded$bytes), as.double(scanned$bytes))) {
    stop("The task manifest does not match the checkpoint shards.",
         call. = FALSE)
  }
  if (!identical(as.integer(config$task_count), as.integer(nrow(plan))) ||
      !identical(
        as.integer(config$method_fit_count), as.integer(nrow(summary))
      ) || nrow(plan) != gwrs_selection_summary_expected_tasks ||
      nrow(summary) != gwrs_selection_summary_expected_fits) {
    stop("Configured task or method-fit counts do not match the run.",
         call. = FALSE)
  }
  gwrs_selection_summary_validate_pairing(summary, plan)
  gwrs_selection_summary_validate_metrics(summary)
  list(
    status = status,
    config = config,
    plan = plan,
    task_manifest = manifest,
    summary = summary
  )
}

gwrs_selection_summary_output_manifest <- function(directory) {
  relative <- list.files(
    directory, recursive = TRUE, all.files = TRUE, no.. = TRUE,
    include.dirs = FALSE
  )
  relative <- relative[!relative %in% c("manifest-sha256.csv", "COMPLETED")]
  paths <- file.path(directory, relative)
  manifest <- data.frame(
    file = relative,
    bytes = as.double(unname(file.info(paths)$size)),
    sha256 = gwrs_sim_sha256(paths),
    stringsAsFactors = FALSE
  )
  manifest[order(manifest$file), , drop = FALSE]
}

gwrs_selection_summary_verify_output <- function(directory) {
  manifest_path <- file.path(directory, "manifest-sha256.csv")
  if (!file.exists(manifest_path)) return(FALSE)
  recorded <- tryCatch(
    utils::read.csv(manifest_path, stringsAsFactors = FALSE),
    error = function(error) NULL
  )
  if (is.null(recorded) || !identical(
    names(recorded), c("file", "bytes", "sha256")
  ) || anyDuplicated(recorded$file)) return(FALSE)
  files <- file.path(directory, recorded$file)
  all(file.exists(files)) && identical(
    as.double(unname(file.info(files)$size)), as.double(recorded$bytes)
  ) && identical(
    unname(gwrs_sim_sha256(files)), as.character(recorded$sha256)
  )
}

gwrs_selection_summary_analysis_config <- function(study_directory,
                                                     source,
                                                     script_path) {
  source_files <- file.path(study_directory, c(
    "COMPLETED", "manifest-sha256.csv", "scientific-signature.csv",
    "task-manifest.csv", "task-summary.csv"
  ))
  data.frame(
    field = c(
      "schema", "created_utc", "source_study", "source_mode",
      "source_design_version", "source_directory",
      "source_task_count", "source_method_fit_count", "frozen_methods",
      "selection_methods", "metric_summary_metrics",
      "analysis_sets", "error_comparison", "selection_comparison",
      "difference_definition", "mc_ci", "favors_rule", "na_handling",
      "convergence_handling", "paired_signature_fields", "script_sha256",
      paste0("source_", basename(source_files), "_sha256")
    ),
    value = c(
      gwrs_selection_summary_schema,
      gwrs_sim_utc(),
      gwrs_selection_summary_study,
      source$config$mode,
      source$config$design_version,
      study_directory,
      nrow(source$plan),
      nrow(source$summary),
      paste(gwrs_selection_summary_methods, collapse = ","),
      paste(gwrs_selection_summary_selection_methods, collapse = ","),
      paste(c(
        gwrs_selection_summary_error_metrics,
        gwrs_selection_summary_descriptive_metrics,
        gwrs_selection_summary_selection_metrics
      ), collapse = ","),
      "all_completed,all_converged",
      paste(
        "each non-GWR method minus GWR;",
        "direct SCAD minus Lasso and MCP minus Lasso"
      ),
      "elastic_net,scad,mcp minus lasso",
      "target - reference",
      "two-sided 95% paired Monte Carlo t interval",
      "target/reference only when the 95% interval excludes zero",
      paste(
        "n_total retains every replication; n_finite records estimable values;",
        "constant prediction with two-class truth has MCC=0;",
        "one-class-truth MCC and null active-coefficient RMSE remain NA"
      ),
      paste(
        "all_completed is primary intention-to-fit; all_converged is a",
        "pre-specified sensitivity analysis and paired rows require both fits"
      ),
      paste(gwrs_selection_summary_paired_fields, collapse = ","),
      gwrs_sim_sha256(script_path),
      gwrs_sim_sha256(source_files)
    ),
    stringsAsFactors = FALSE
  )
}

gwrs_selection_summary_resolve_output <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("OUTPUT_DIR must be one nonempty path.", call. = FALSE)
  }
  if (file.exists(path) || dir.exists(path)) {
    stop("Refusing existing OUTPUT_DIR: `", path, "`.", call. = FALSE)
  }
  parent <- normalizePath(dirname(path), mustWork = TRUE)
  name <- basename(path)
  if (!nzchar(name) || name %in% c(".", "..")) {
    stop("OUTPUT_DIR must name a new child directory.", call. = FALSE)
  }
  file.path(parent, name)
}

gwrs_selection_summary_script_path <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", arguments, value = TRUE)
  if (length(file_argument) != 1L) {
    stop("Run this analysis with Rscript.", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
}

gwrs_selection_summary_main <- function(arguments = commandArgs(
                                           trailingOnly = TRUE
                                         ),
                                         script_path = NULL) {
  if (length(arguments) != 2L || any(!nzchar(arguments))) {
    stop(
      "Usage: Rscript --vanilla summarize-penalized-selection-study.R ",
      "STUDY_DIR OUTPUT_DIR", call. = FALSE
    )
  }
  if (is.null(script_path)) script_path <- gwrs_selection_summary_script_path()
  script_path <- normalizePath(script_path, mustWork = TRUE)
  source(file.path(dirname(script_path), "resume-utils.R"))

  study_directory <- normalizePath(arguments[[1L]], mustWork = TRUE)
  output_directory <- gwrs_selection_summary_resolve_output(arguments[[2L]])
  if (identical(study_directory, output_directory)) {
    stop("STUDY_DIR and OUTPUT_DIR must differ.", call. = FALSE)
  }
  source <- gwrs_selection_summary_validate_source(study_directory)

  metric_summary <- rbind(
    gwrs_selection_summary_metric_table(source$summary, "all_completed"),
    gwrs_selection_summary_metric_table(source$summary, "all_converged")
  )
  rownames(metric_summary) <- NULL
  paired_differences <- rbind(
    gwrs_selection_summary_paired_table(source$summary, "all_completed"),
    gwrs_selection_summary_paired_table(source$summary, "all_converged")
  )
  rownames(paired_differences) <- NULL
  convergence_summary <- gwrs_selection_summary_convergence_table(
    source$summary
  )
  bandwidth_summary <- gwrs_selection_summary_bandwidth_table(source$summary)
  design_summary <- gwrs_selection_summary_design_table(
    source$plan, source$summary
  )
  analysis_config <- gwrs_selection_summary_analysis_config(
    study_directory, source, script_path
  )

  staging <- tempfile(
    paste0(".", basename(output_directory), ".partial-"),
    tmpdir = dirname(output_directory)
  )
  if (!dir.create(staging, showWarnings = FALSE)) {
    stop("Could not create the analysis staging directory.", call. = FALSE)
  }
  published <- FALSE
  on.exit({
    if (!published && dir.exists(staging)) {
      unlink(staging, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  outputs <- list(
    "metric-summary.csv" = metric_summary,
    "paired-differences.csv" = paired_differences,
    "convergence-summary.csv" = convergence_summary,
    "bandwidth-summary.csv" = bandwidth_summary,
    "design-summary.csv" = design_summary,
    "analysis-config.csv" = analysis_config
  )
  for (name in names(outputs)) {
    gwrs_sim_atomic_csv(outputs[[name]], file.path(staging, name))
  }
  gwrs_sim_atomic_text(
    gwrs_sim_session_info(), file.path(staging, "session-info.txt")
  )
  manifest <- gwrs_selection_summary_output_manifest(staging)
  gwrs_sim_atomic_csv(manifest, file.path(staging, "manifest-sha256.csv"))
  if (!gwrs_selection_summary_verify_output(staging)) {
    stop("The analysis size/SHA-256 manifest failed verification.",
         call. = FALSE)
  }

  completion <- c(
    paste0("Schema: ", gwrs_selection_summary_schema),
    paste0("Completed-UTC: ", gwrs_sim_utc()),
    paste0("Source-Study: ", gwrs_selection_summary_study),
    paste0("Source-Tasks: ", nrow(source$plan)),
    paste0(
      "Source-Completion-SHA256: ",
      gwrs_sim_sha256(file.path(study_directory, "COMPLETED"))
    ),
    paste0(
      "Manifest-SHA256: ",
      gwrs_sim_sha256(file.path(staging, "manifest-sha256.csv"))
    )
  )
  gwrs_sim_atomic_text(completion, file.path(staging, "COMPLETED"))
  if (!gwrs_selection_summary_verify_output(staging)) {
    stop("Analysis verification failed after completion marking.",
         call. = FALSE)
  }
  completion_record <- tryCatch(
    as.list(as.data.frame(
      read.dcf(file.path(staging, "COMPLETED")), stringsAsFactors = FALSE
    )[1L, ]),
    error = function(error) NULL
  )
  if (is.null(completion_record) ||
      !identical(completion_record$Schema, gwrs_selection_summary_schema) ||
      !identical(
        completion_record$`Manifest-SHA256`,
        gwrs_sim_sha256(file.path(staging, "manifest-sha256.csv"))
      )) {
    stop("The analysis completion marker failed verification.",
         call. = FALSE)
  }
  if (!file.rename(staging, output_directory)) {
    stop("Could not publish the completed analysis directory.",
         call. = FALSE)
  }
  published <- TRUE
  cat("Validated study:", study_directory, "\n")
  cat("Analysis output:", output_directory, "\n")
  invisible(output_directory)
}

if (identical(sys.nframe(), 0L)) {
  gwrs_selection_summary_main()
}
