# Install output-validation wrappers around the frozen native model paths.
# The statistical objectives, paths, anchors, folds and retry budgets are not
# changed. Full coefficient cubes used for inline certificates are discarded
# after the compact diagnostics have been computed, except for selected outer
# fits whose coefficients are required scientific outputs.

usa_study_v1_install_models <- function(e) {
  original_path <- e$usa_path
  original_score <- e$usa_score
  original_evaluate <- e$usa_evaluate

  e$usa_path <- function(ctx, train, candidate, until = NULL, keep = FALSE,
                         full_diagnostics = FALSE) {
    capture <- ctx$native_capture
    capture$nonconvex_calls <- list()
    capture$convex_calls <- list()
    capture$working_metrics <- list()
    result <- original_path(
      ctx, train, candidate, until, keep, full_diagnostics
    )
    result$native_states <- capture$nonconvex_calls
    result$convex_native_states <- capture$convex_calls
    result$working_diagnostics <- if (length(capture$working_metrics)) {
      value <- do.call(rbind, capture$working_metrics)
      rownames(value) <- NULL
      value
    } else NULL

    if (candidate$method %in% c("scad", "mcp")) {
      stopifnot(length(result$native_states) == length(result$attempts))
      for (j in seq_along(result$native_states)) {
        state <- result$native_states[[j]]
        stopifnot(
          all(state$path_state %in% 0:2),
          identical(state$path_state == 1L, state$converged == 1L),
          all(state$iterations <= result$attempts[[j]]$iterations_cap)
        )
      }
      result$valid <- result$valid & apply(
        tail(result$native_states, 1)[[1]]$path_state == 1L, 2, all
      )
    }
    if (candidate$method %in% c("lasso", "en", "ridge")) {
      stopifnot(length(result$convex_native_states) == length(result$attempts))
      for (j in seq_along(result$convex_native_states)) {
        state <- result$convex_native_states[[j]]
        stopifnot(
          ncol(state$converged) == length(result$parameter),
          all(state$converged %in% 0:1),
          all(state$iterations <= result$attempts[[j]]$iterations_cap)
        )
      }
      result$valid <- result$valid & apply(
        tail(result$convex_native_states, 1)[[1]]$converged == 1L, 2, all
      )
    }
    if (candidate$method %in% c("lasso", "en")) {
      stopifnot(length(capture$working_metrics) == length(result$attempts))
    }
    capture$last_nonconvex <- result$native_states
    capture$last_convex <- result$convex_native_states
    capture$last_working <- result$working_diagnostics
    result
  }

  e$usa_score <- function(fit, y) {
    value <- original_score(fit, y)
    value$sse[!value$valid] <- NA_real_
    value$sae[!value$valid] <- NA_real_
    value
  }

  e$usa_evaluate <- function(ctx, train, data, split, candidate,
                             until = NULL, keep = FALSE,
                             full_diagnostics = FALSE) {
    capture <- ctx$native_capture
    capture$last_nonconvex <- list()
    capture$last_convex <- list()
    capture$last_working <- NULL
    value <- original_evaluate(
      ctx, train, data, split, candidate, until, keep, full_diagnostics
    )
    value$native_states <- capture$last_nonconvex
    value$convex_native_states <- capture$last_convex
    value$working_diagnostics <- capture$last_working
    if (!is.null(value$score)) {
      value$score$sse[!value$score$valid] <- NA_real_
      value$score$sae[!value$score$valid] <- NA_real_
    }
    value
  }
}
