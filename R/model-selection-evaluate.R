selection_fit <- function(candidate, prep, settings, neighbors) {
  components <- selection_components(selection_decode(candidate$subset_id), prep)
  if (candidate$family == "OLS") {
    design <- cbind("(Intercept)" = 1, components$x)
    fit <- stats::lm.fit(design, components$y)
    fit$terms <- components$terms
    fit$xlevels <- components$xlevels
    fit$contrasts <- components$contrasts
    fit$model <- components$model
    fit$call <- call("lm", formula = stats::formula(components$terms))
    class(fit) <- "lm"
  } else {
    nb <- neighbors
    nb$k <- as.integer(candidate$k)
    nb$index <- neighbors$index[seq_len(nb$k), , drop = FALSE]
    nb$distance <- neighbors$distance[seq_len(nb$k), , drop = FALSE]
    fit <- gwr_fit(components$x, components$y, neighbors = nb,
      kernel = candidate$kernel,
      bandwidth = if (candidate$bandwidth_type == "fixed") candidate$bw else NULL,
      standardize = settings$standardize, control = settings$control)
    fit <- decorate_formula_fit(fit, components,
      call("gwr_fit", formula = stats::formula(components$terms)))
  }
  fit
}

selection_evaluate <- function(candidate, prep, settings, neighbors, weights) {
  row <- candidate
  row$n <- length(prep$full$y)
  for (key in c(selection_criteria(), "rss", "trS", "trStS", "edf", "r_squared",
                "min_rank", "min_neff", "min_positive_weights", "max_leverage",
                "min_abs_loo_denominator", "max_condition_number", "max_vif",
                "f1_p", "f2_p", "moran_I", "moran_p")) row[[key]] <- NA_real_
  row$eligible <- FALSE
  row$status <- "failed"
  row$reason <- ""
  row$flags <- ""
  row$diagnostic_notes <- ""
  row$f_status <- if (candidate$family == "OLS") "not_applicable" else
    if (settings$f_tests) "not_calculated" else "disabled"
  row$moran_status <- if (is.null(weights)) "not_calculated: weights_not_supplied" else "not_calculated"
  f <- NULL
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  outcome <- tryCatch(withCallingHandlers({
    fit <- selection_fit(candidate, prep, settings, neighbors)
    p <- candidate$n_columns
    if (candidate$family == "OLS") {
      h <- rowSums(qr.Q(fit$qr)[, seq_len(fit$rank), drop = FALSE]^2)
      g <- gaussian_smoother_diagnostics(prep$full$y, fit$fitted.values, fit$residuals, h, h)
      valid <- fit$rank == p + 1L && all(is.finite(fit$coefficients))
      row$min_rank <- fit$rank - 1L
      row$min_neff <- row$n
      row$min_positive_weights <- row$n
      if (!valid) row$reason <- "rank_deficient_global_design"
      if (p > 0L) {
        x <- prep$full$x[, prep$assignment %in% selection_decode(candidate$subset_id), drop = FALSE]
        col <- tryCatch(gwr_global_collinearity(x), error = identity)
        if (!inherits(col, "error")) {
          row$max_condition_number <- col$condition_number$condition_number
          row$max_vif <- max(col$vif$vif)
        }
      }
    } else {
      dg <- gwr_diagnostics(fit)
      g <- dg$global
      h <- dg$hat_diag
      local <- dg$local
      row$min_rank <- min(local$local_rank)
      row$min_neff <- min(local$effective_sample_size)
      row$min_positive_weights <- min(local$positive_weights)
      reasons <- character()
      if (!all(fit$solver_diagnostics$converged)) reasons <- c(reasons, "nonconverged")
      if (!all(fit$solver_diagnostics$solver == "direct_solve") ||
          !all(local$operator_solver == "direct_solve")) reasons <- c(reasons, "non_direct_solver")
      if (!all(local$local_rank == p)) reasons <- c(reasons, "rank_deficient_local_design")
      if (!all(is.finite(fit$coefficients)) || !all(is.finite(fit$fitted.values)))
        reasons <- c(reasons, "nonfinite_fit")
      valid <- !length(reasons)
      row$reason <- paste(reasons, collapse = ";")
      col <- tryCatch(gwr_local_collinearity(fit), error = identity)
      if (inherits(col, "error")) {
        row$diagnostic_notes <- paste("collinearity:", conditionMessage(col))
      } else {
        row$max_condition_number <- max(col$condition_number$condition_number)
        row$max_vif <- max(col$vif)
      }
    }
    for (key in c(selection_criteria(), "rss", "trS", "trStS", "edf", "r_squared"))
      row[[key]] <- g[[key]]
    row$max_leverage <- max(h)
    row$min_abs_loo_denominator <- min(abs(1 - h))
    if (!is.finite(g$edf) || g$edf <= 0) {
      valid <- FALSE
      row$reason <- paste(c(row$reason[nzchar(row$reason)], "nonpositive_residual_edf"), collapse = ";")
    }
    row$eligible <- isTRUE(valid)
    row$status <- if (row$eligible) "eligible" else "invalid"
    if (row$eligible && !is.finite(row[[settings$criterion]]))
      row$status <- "criterion_unavailable"
    flag <- c(
      if (isTRUE(row$max_condition_number > 30)) "condition_number_above_30",
      if (isTRUE(row$max_vif > 10)) "vif_above_10",
      if (isTRUE(row$min_neff <= p + 1L)) "effective_sample_at_most_parameter_count",
      if (isTRUE(row$max_leverage > 0.99)) "leverage_above_0.99",
      if (isTRUE(row$edf < 5)) "residual_edf_below_5",
      if (isTRUE(row$boundary) && candidate$family == "GWR") "grid_boundary")
    row$flags <- paste(flag, collapse = ";")
    if (row$eligible && candidate$family == "GWR" && settings$f_tests) {
      f <- tryCatch(gwr_diagnostic_tests(fit, memory_limit_mb = settings$memory_limit_mb), error = identity)
      if (inherits(f, "error")) {
        row$f_status <- paste("unavailable:", conditionMessage(f))
        f <- NULL
      } else {
        row$f1_p <- f$global_tests$p_value[match("F1", f$global_tests$test)]
        row$f2_p <- f$global_tests$p_value[match("F2", f$global_tests$test)]
        row$f_status <- if (anyNA(c(f$global_tests$p_value, f$coefficient_tests$p_value)))
          "partially_unavailable" else "conditional_exploratory"
        f <- list(global_tests = f$global_tests, coefficient_tests = f$coefficient_tests,
          validity = f$validity)
      }
    }
    if (row$eligible && !is.null(weights)) {
      m <- tryCatch(compute_gwrs_weighted_moran(fit$residuals, weights,
        settings$moran_permutations, settings$moran_row_standardize, settings$moran_seed,
        settings$moran_zero_policy, settings$control$n_threads,
        settings$control$grain_size, settings$memory_limit_mb), error = identity)
      if (inherits(m, "error")) row$moran_status <- paste("unavailable:", conditionMessage(m)) else {
        row$moran_I <- m$I
        row$moran_p <- m$p_value
        row$moran_status <- if (settings$moran_permutations == 0) "statistic_only" else "exploratory_permutation"
      }
    }
    NULL
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  if (inherits(outcome, "error")) {
    row$eligible <- FALSE
    row$status <- "failed"
    row$reason <- conditionMessage(outcome)
  }
  row$diagnostic_notes <- paste(c(row$diagnostic_notes[nzchar(row$diagnostic_notes)],
    unique(warnings)), collapse = ";")
  row$elapsed_seconds <- proc.time()[["elapsed"]] - started
  list(row = row, f_tests = f)
}

selection_record_table <- function(records) {
  if (!length(records)) return(data.frame())
  out <- do.call(rbind, lapply(records, function(x) x$row))
  rownames(out) <- NULL
  out
}

selection_order <- function(table, criterion, family = NULL) {
  if (!nrow(table)) return(integer())
  good <- table$eligible & is.finite(table[[criterion]])
  if (!is.null(family)) good <- good & table$family == family
  at <- which(good)
  b <- ifelse(table$bandwidth_type == "adaptive", table$k, table$bw)
  b[is.na(b)] <- 0
  at[order(table[[criterion]][at], table$n_columns[at],
    table$family[at] != "OLS", -b[at], table$kernel[at],
    table$subset_id[at], table$candidate_id[at], na.last = TRUE, method = "radix")]
}

selection_result <- function(records, prep, settings, neighbors, provenance, call,
    complete, reason, evaluated, expected, max_models) {
  table <- selection_record_table(records)
  winners <- list()
  for (metric in selection_criteria()) for (scope in c("GWR", "overall")) {
    at <- selection_order(table, metric, if (scope == "GWR") "GWR" else NULL)
    if (length(at)) {
      win <- table[at[1L], , drop = FALSE]
      win$criterion <- metric
      win$comparison <- scope
      win$criterion_value <- win[[metric]]
      winners[[length(winners) + 1L]] <- win
    }
  }
  best <- selection_order(table, settings$criterion, "GWR")
  overall <- selection_order(table, settings$criterion)
  selected <- if (length(best)) table[best[1L], , drop = FALSE] else NULL
  selected_overall <- if (length(overall)) table[overall[1L], , drop = FALSE] else NULL
  refit <- function(row, full = FALSE) {
    if (is.null(row)) return(NULL)
    final_settings <- settings
    if (full) final_settings$control$diagnostic_level <- "full"
    fit <- selection_fit(row, prep, final_settings, neighbors)
    if (row$family == "GWR") {
      score <- gwr_diagnostics(fit)$global[[settings$criterion]]
      if (!isTRUE(all.equal(score, row[[settings$criterion]], tolerance = 1e-10)))
        stop("Selected-model refit did not reproduce its criterion.", call. = FALSE)
    }
    fit
  }
  best_model <- refit(selected, full = settings$local_inference)
  local_inference <- NULL
  local_status <- if (!settings$local_inference) "disabled" else
    if (is.null(best_model)) "not_calculated: no_valid_gwr" else "not_calculated"
  if (!is.null(best_model)) {
    best_model$selection_context <- list(candidate_id = selected$candidate_id,
      criterion = settings$criterion, terms = selected$terms,
      kernel = selected$kernel, bandwidth_type = selected$bandwidth_type,
      k = selected$k, bw = selected$bw, search_complete = complete,
      sample_rows = provenance$sample_rows)
    if (settings$local_inference) {
      local_inference <- tryCatch(gwr_local_inference(best_model), error = identity)
      if (inherits(local_inference, "error")) {
        local_status <- paste("unavailable:", conditionMessage(local_inference))
        local_inference <- NULL
      } else local_status <- local_inference$status
    }
  }
  best_overall <- if (!is.null(selected) && !is.null(selected_overall) &&
    selected$candidate_id == selected_overall$candidate_id) best_model else refit(selected_overall)
  if (nrow(table)) {
    table$rank <- NA_integer_
    table$rank[overall] <- seq_along(overall)
    table$delta <- NA_real_
    if (length(overall)) table$delta[overall] <- table[[settings$criterion]][overall] -
      table[[settings$criterion]][overall[1L]]
    table <- table[c(overall, setdiff(seq_len(nrow(table)), overall)), , drop = FALSE]
    rownames(table) <- NULL
  }
  combine_tests <- function(field) {
    rows <- lapply(records, function(r) {
      x <- r$f_tests[[field]]
      if (is.null(x)) return(NULL)
      cbind(candidate_id = r$row$candidate_id, x, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }
  result <- list(call = call, best_model = best_model, best_overall = best_overall,
    selected = selected, selected_overall = selected_overall,
    ranking = table, best_by_criterion = if (length(winners)) do.call(rbind, winners) else data.frame(),
    global_tests = combine_tests("global_tests"), coefficient_tests = combine_tests("coefficient_tests"),
    local_inference = local_inference, local_inference_status = local_status,
    candidate_diagnostics = table, settings = settings, provenance = provenance,
    search_status = list(complete = complete,
      status = if (!complete) "incomplete" else if (is.null(best_model)) "no_valid_gwr" else "complete",
      stop_reason = reason, method = settings$search,
      exhaustive_optimum = complete && settings$search == "exhaustive" && !is.null(best_model),
      evaluated = length(records), new_evaluations = evaluated,
      expected = expected, max_models = max_models,
      final_refits = as.integer(!is.null(best_model)) +
        as.integer(!is.null(selected_overall) && (is.null(selected) ||
          selected_overall$candidate_id != selected$candidate_id))),
    notes = c("Best means best among eligible evaluated candidates under the selected criterion.",
      "Adaptive Gaussian/exponential weights are truncated at the retained k neighbors.",
      "F and Moran p-values are exploratory, unadjusted and not ranking objectives.",
      "F1/F2 transform the same RSS ratio; separate reference p-values are not independent evidence or calibrated tests.",
      "Local p-values and per-coefficient BH adjustment are exploratory; they do not account for subset/kernel/bandwidth selection.",
      "Conditional LOO is not external spatial validation.",
      if (!complete) "Provisional result: search is incomplete.",
      if (settings$search == "stepwise") "Stepwise search has no global optimum guarantee."))
  class(result) <- "gwrs_selection"
  result
}
