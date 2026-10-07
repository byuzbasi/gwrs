#' Structured gwrs model report
#'
#' Combines calibration settings, a global OLS baseline when available,
#' statistical and numerical diagnostics, coefficient-surface summaries, and
#' optional local inference, F-tests, Moran, and collinearity objects. Expensive
#' spatial diagnostics are not computed implicitly.
#'
#' @param object A fitted `gwrs_fit`.
#' @param include_global Include a global OLS baseline when training data were
#'   retained.
#' @param include_inference Attempt local inference when full diagnostics exist.
#' @param include_f_tests Include GWR-style F-tests.
#' @param inference Optional precomputed [gwr_local_inference()] result.
#' @param f_tests Optional precomputed [gwr_diagnostic_tests()] result.
#' @param moran Optional precomputed [gwr_moran()] result.
#' @param collinearity Optional precomputed [gwr_local_collinearity()] result.
#' @param global_collinearity Optional [gwr_global_collinearity()] result.
#' @param local_moran Optional [gwr_local_moran()] result. A compact cluster
#'   summary is printed; all local results are retained in the report.
#' @param allow_approximate_tests Permit labelled conditional F-tests for a
#'   penalized model.
#' @param alpha,adjust Settings used if local inference must be computed.
#' @param digits Significant digits used by the print method.
#' @param ... Unused.
#'
#' @return A printed `gwrs_report` object, invisibly.
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_lasso_fit(
#'   x, quake$mag, coords = coords, k = 20, lambda = 0.05,
#'   control = gwrs_control(n_threads = 1)
#' )
#' report <- gwrs_report(fit, include_inference = FALSE)
#' report$calibration
#'
#' @export
gwrs_report <- function(object,
                        include_global = TRUE,
                        include_inference = TRUE,
                        include_f_tests = FALSE,
                        inference = NULL,
                        f_tests = NULL,
                        moran = NULL,
                        collinearity = NULL,
                        allow_approximate_tests = FALSE,
                        alpha = 0.05,
                        adjust = "BH",
                        digits = 4,
                        ...,
                        global_collinearity = NULL,
                        local_moran = NULL) {
  if (!inherits(object, "gwrs_fit")) {
    stop("`object` must inherit from 'gwrs_fit'.", call. = FALSE)
  }
  diagnostics <- tryCatch(gwr_diagnostics(object), error = identity)
  notes <- character()
  if (inherits(diagnostics, "error")) {
    notes <- c(notes, paste("Statistical diagnostics unavailable:",
                            conditionMessage(diagnostics)))
    diagnostics <- NULL
  }
  global_baseline <- NULL
  if (isTRUE(include_global)) {
    if (!is.null(object$x) && !is.null(object$y)) {
      design <- cbind(`(Intercept)` = 1, object$x)
      fit <- stats::lm.fit(design, object$y)
      rss <- sum(fit$residuals^2)
      tss <- sum((object$y - mean(object$y))^2)
      global_baseline <- data.frame(
        metric = c("RSS", "Residual df", "R-squared"),
        value = c(rss, length(object$y) - fit$rank,
                  if (tss > 0) 1 - rss / tss else NA_real_)
      )
    } else {
      notes <- c(notes, "Global OLS unavailable because training data were not retained.")
    }
  }
  if (isTRUE(include_inference) && is.null(inference)) {
    inference <- tryCatch(
      gwr_local_inference(object, alpha = alpha, adjust = adjust),
      error = identity
    )
    if (inherits(inference, "error")) {
      notes <- c(notes, paste("Local inference not shown:",
                              conditionMessage(inference)))
      inference <- NULL
    }
  }
  if (isTRUE(include_f_tests) && is.null(f_tests)) {
    f_tests <- tryCatch(
      gwr_diagnostic_tests(
        object, allow_approximate = allow_approximate_tests
      ),
      error = identity
    )
    if (inherits(f_tests, "error")) {
      notes <- c(notes, paste("F-tests not shown:", conditionMessage(f_tests)))
      f_tests <- NULL
    }
  }
  coefficient_summary <- t(vapply(
    seq_len(ncol(object$coefficients)),
    function(column) {
      value <- object$coefficients[, column]
      c(
        mean = mean(value), sd = stats::sd(value), min = min(value),
        q25 = unname(stats::quantile(value, 0.25)),
        median = stats::median(value),
        q75 = unname(stats::quantile(value, 0.75)), max = max(value),
        zero_fraction = mean(abs(value) <= 1e-10)
      )
    },
    numeric(8L)
  ))
  rownames(coefficient_summary) <- colnames(object$coefficients)
  calibration <- data.frame(
    setting = c(
      "Method", "Observations", "Predictors", "Kernel", "Bandwidth type",
      "Bandwidth/k", "Standardized", "Diagnostic level"
    ),
    value = c(
      object$method, nrow(object$coefficients), ncol(object$coefficients) - 1L,
      object$kernel, if (isTRUE(object$adaptive)) "adaptive" else "fixed",
      if (isTRUE(object$adaptive)) object$neighbors$k else object$bandwidth,
      object$standardize,
      object$control$diagnostic_level %||% "legacy"
    ),
    stringsAsFactors = FALSE
  )
  result <- list(
    call = object$call,
    calibration = calibration,
    global_baseline = global_baseline,
    diagnostics = diagnostics,
    solver_diagnostics = object$solver_diagnostics %||% object$diagnostics,
    coefficients = coefficient_summary,
    inference = inference,
    f_tests = f_tests,
    moran = moran,
    collinearity = collinearity,
    global_collinearity = global_collinearity,
    local_moran = local_moran,
    notes = notes,
    digits = digits
  )
  class(result) <- "gwrs_report"
  print(result, digits = digits)
  invisible(result)
}

#' @export
print.gwrs_report <- function(x, digits = x$digits %||% 4, ...) {
  rule <- paste(rep("-", 72), collapse = "")
  cat(rule, "\n", sep = "")
  cat("gwrs model report\n")
  cat(rule, "\n", sep = "")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  cat("\nCalibration\n")
  print(x$calibration, row.names = FALSE, right = FALSE)
  if (!is.null(x$global_baseline)) {
    cat("\nGlobal OLS baseline\n")
    baseline <- x$global_baseline
    baseline$value <- signif(baseline$value, digits)
    print(baseline, row.names = FALSE)
  }
  if (!is.null(x$diagnostics)) {
    cat("\nModel diagnostics\n")
    global <- x$diagnostics$global
    keys <- intersect(
      c("rss", "rmse", "mae", "r_squared", "adjusted_r_squared", "trS",
        "trStS", "enp", "edf", "aic", "aicc", "bic", "gcv", "loo_cv",
        "median_absolute_error", "joint_objective"),
      names(global)
    )
    table <- data.frame(
      metric = keys,
      value = vapply(global[keys], as.numeric, numeric(1L))
    )
    table$value <- signif(table$value, digits)
    print(table, row.names = FALSE)
    cat("Basis:", x$diagnostics$validity$basis, "\n")
  }
  cat("\nCoefficient-surface summaries\n")
  print(round(x$coefficients, digits))
  if (!is.null(x$inference)) {
    cat("\n")
    print(x$inference)
  }
  if (!is.null(x$f_tests)) {
    cat("\n")
    print(x$f_tests)
  }
  if (!is.null(x$moran)) {
    cat("\n")
    print(x$moran)
  }
  if (!is.null(x$collinearity)) {
    cat("\n")
    print(x$collinearity)
  }
  if (!is.null(x$global_collinearity)) {
    cat("\n")
    print(x$global_collinearity)
  }
  if (!is.null(x$local_moran)) {
    cat("\nResidual LISA (total randomization; unadjusted p-values)\n")
    assessed <- is.finite(x$local_moran$p_value)
    cat("  assessed/unassessed:", sum(assessed), "/", sum(!assessed), "\n")
    print(table(x$local_moran$cluster[assessed]))
  }
  if (length(x$notes)) {
    cat("\nNotes\n")
    for (note in x$notes) cat(" -", note, "\n")
  }
  cat(rule, "\n", sep = "")
  invisible(x)
}
