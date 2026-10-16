#' @export
print.gwrs_selection <- function(x, ...) {
  cat("Joint unpenalized GWR selection\n")
  cat("  search:", x$search_status$method, " / ", x$search_status$status, "\n", sep = "")
  cat("  evaluated:", x$search_status$evaluated, "  criterion:", x$settings$criterion, "\n")
  if (is.null(x$selected)) cat("  no eligible GWR selected\n") else {
    z <- x$selected
    cat("  best GWR:", z$terms, "\n")
    cat("  kernel:", z$kernel, "  bandwidth:", z$bandwidth_type,
      if (z$bandwidth_type == "adaptive") z$k else z$bw, "\n")
    cat("  criterion:", format(z[[x$settings$criterion]]), "\n")
    cat("  Moran:", z$moran_status, "\n")
    if (!is.null(x$local_inference_status))
      cat("  local coefficient inference:", x$local_inference_status, "\n")
    if (nzchar(z$flags)) cat("  flags:", z$flags, "\n")
  }
  if (!is.null(x$selected_overall))
    cat("  best including OLS:", x$selected_overall$family, "-", x$selected_overall$terms, "\n")
  if (!x$search_status$complete) cat("  PROVISIONAL:", x$search_status$stop_reason, "\n")
  invisible(x)
}

#' Summarize joint GWR selection
#' @param object A result from [gwr_select()].
#' @param n Number of leading candidate rows to show.
#' @param ... Additional arguments.
#' @return A summary.gwrs_selection object retaining the selected tests,
#'   local coefficient inference/status, criterion winners, top candidates
#'   and search status.
#' @export
summary.gwrs_selection <- function(object, n = 10L, ...) {
  selection_integer(n, "n", 1)
  select_tests <- function(x) {
    if (is.null(x) || is.null(object$selected)) return(NULL)
    x[x$candidate_id == object$selected$candidate_id, , drop = FALSE]
  }
  out <- list(selection = object, top = utils::head(object$ranking, n),
    global_tests = select_tests(object$global_tests),
    coefficient_tests = select_tests(object$coefficient_tests),
    local_inference = object$local_inference,
    local_inference_status = object$local_inference_status)
  class(out) <- "summary.gwrs_selection"
  out
}

#' @export
print.summary.gwrs_selection <- function(x, ...) {
  print(x$selection)
  cols <- c("terms", "family", "kernel", "bandwidth_type", "k", "bw",
    selection_criteria(), "status", "flags")
  if (nrow(x$top)) {
    cat("\nLeading candidates\n")
    print(x$top[, cols, drop = FALSE], row.names = FALSE)
  }
  if (!is.null(x$global_tests)) {
    cat("\nSelected GWR: conditional, exploratory F1/F2 (uncalibrated references)\n")
    cat("F1/F2 transform the same RSS ratio; their p-values are not independent evidence.\n")
    print(x$global_tests[, setdiff(names(x$global_tests), "candidate_id"), drop = FALSE], row.names = FALSE)
  }
  if (!is.null(x$coefficient_tests)) {
    cat("\nSelected GWR: coefficient spatial variability (F3, raw p-values)\n")
    print(x$coefficient_tests[, setdiff(names(x$coefficient_tests), "candidate_id"), drop = FALSE], row.names = FALSE)
  }
  if (!is.null(x$local_inference)) {
    cat("\nSelected GWR: exploratory local coefficient inference\n")
    print(x$local_inference)
  }
  invisible(x)
}

#' @export
print.gwrs_selection_plan <- function(x, ...) {
  cat("GWR selection design (no fits)\n")
  cat("  n:", x$n, "  scope terms:", length(x$scope), "\n")
  cat("  search:", x$search, "  candidates:",
    if (is.na(x$expected_candidates)) "adaptive search; count unknown" else x$expected_candidates, "\n")
  cat("  neighbor workspace estimate (MiB):", x$neighbor_workspace_mb, "\n")
  cat("  Moran:", x$moran, "\n")
  invisible(x)
}

#' @export
coef.gwrs_selection <- function(object, ...) {
  if (is.null(object$best_model)) stop("No eligible GWR was selected.", call. = FALSE)
  stats::coef(object$best_model, ...)
}

#' @export
predict.gwrs_selection <- function(object, ...) {
  if (is.null(object$best_model)) stop("No eligible GWR was selected.", call. = FALSE)
  stats::predict(object$best_model, ...)
}
