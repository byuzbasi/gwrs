#' Exploratory local coefficient inference
#'
#' Computes local standard errors, test statistics, raw p-values, adjusted
#' p-values, and significance flags from a fit created with
#' `gwrs_control(diagnostics = "full")`. Sparse fits default to an unpenalized
#' refit on each selected local active set; all reported p-values remain
#' exploratory because selection and overlapping neighborhoods induce
#' dependence.
#'
#' @param object A Gaussian least-squares `gwrs_fit` with full diagnostics, or
#'   a `gwrs_selection` whose selected GWR retains full diagnostics.
#' @param alpha Significance level.
#' @param df Reference degrees-of-freedom rule.
#' @param adjust Multiplicity adjustment passed to [stats::p.adjust()].
#' @param adjust_scope Adjust separately by coefficient or over all local tests.
#' @param estimator Use the penalized active-set derivative, the local
#'   post-selection refit, or choose automatically.
#' @param ... Unused.
#'
#' @details The two-sided reference calculation concerns whether each local
#'   coefficient is zero. It is distinct from the F3 spatial-constancy
#'   diagnostic. By default, BH adjustment is applied separately to the
#'   locations of each design-matrix coefficient, including the intercept.
#'   Factor contrasts are separate coefficients. `adjust_scope = "global"`
#'   instead adjusts over all local coefficient tests in the fitted model.
#'   Neither scope includes tests from other candidate models.
#'
#'   The default `df = "edf"` uses effective residual degrees of freedom;
#'   `"residual"` uses n minus the smoother trace. If the requested degrees
#'   of freedom are unavailable or nonpositive, p-values are NA with a warning.
#'   There is no automatic normal fallback; `df = "normal"` requests it
#'   explicitly. Unavailable tests also have NA significance flags and NA
#'   summary values, rather than apparent nonsignificance or infinite minima.
#'
#'   All results are exploratory. BH does not establish calibration of the
#'   underlying reference p-values or account for response-based predictor,
#'   kernel, bandwidth or penalty selection, smoothing bias, or dependence
#'   between overlapping neighborhoods. Selection context is retained when
#'   available; no calibrated post-selection or FDR guarantee is claimed.
#' @return A `gwrs_local_inference` object with local estimates, standard
#'   errors, statistics, raw `p_value`, adjusted `p_adjusted`, `significant`,
#'   coefficient summaries, settings, notes and a calculation status.
#'
#' @examples
#' quake <- datasets::quakes[1:36, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_lasso_fit(
#'   x, quake$mag, coords = coords, k = 24, lambda = 0.05,
#'   control = gwrs_control(n_threads = 1, diagnostics = "full")
#' )
#' inference <- gwr_local_inference(fit, adjust = "BH")
#' inference$summary
#' head(inference$local)
#'
#' @export
gwr_local_inference <- function(object,
                                alpha = 0.05,
                                df = c("edf", "residual", "normal"),
                                adjust = c("BH", "none", "holm", "bonferroni",
                                           "BY", "fdr"),
                                adjust_scope = c("coefficient", "global"),
                                estimator = c("auto", "post_selection",
                                              "penalized"),
                                ...) {
  if (inherits(object, "gwrs_selection")) {
    if (is.null(object$best_model))
      stop("No eligible GWR was selected.", call. = FALSE)
    object <- object$best_model
  }
  if (!inherits(object, "gwr_sl_fit") &&
      !inherits(object, "gwr_nonconvex_fit")) {
    stop(
      "Local coefficient inference is currently defined for Gaussian least-squares GWR fits only.",
      call. = FALSE
    )
  }
  df <- match.arg(df)
  adjust <- match.arg(adjust)
  adjust_scope <- match.arg(adjust_scope)
  estimator <- match.arg(estimator)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  diagnostics <- gwr_diagnostics(object)
  if (!identical(diagnostics$level, "full") ||
      is.null(diagnostics$inference)) {
    stop(
      paste0("Local inference requires `gwrs_control(diagnostics = \"full\")`; ",
        "for gwr_select(), use local_inference = TRUE."),
      call. = FALSE
    )
  }
  if (estimator == "auto") {
    estimator <- if (gwrs_is_sparse_selection_fit(object)) {
      "post_selection"
    } else {
      "penalized"
    }
  }
  source <- if (estimator == "post_selection") {
    diagnostics$post_selection
  } else {
    diagnostics$inference
  }
  if (is.null(source)) {
    stop("The requested inference estimator is unavailable.", call. = FALSE)
  }
  reference_df <- switch(
    df,
    edf = if (estimator == "post_selection") {
      source$global$edf
    } else {
      diagnostics$global$edf
    },
    residual = if (estimator == "post_selection") {
      length(object$residuals) - source$global$trS
    } else {
      length(object$residuals) - diagnostics$global$trS
    },
    normal = Inf
  )
  reference_available <- identical(df, "normal") ||
    (length(reference_df) == 1L && is.finite(reference_df) && reference_df > 0)
  if (!reference_available) {
    warning(
      "Positive reference degrees of freedom are unavailable; p-values are NA. Use df = \"normal\" only for an explicit normal approximation.",
      call. = FALSE
    )
    reference_df <- NA_real_
  }
  coefficient_names <- colnames(source$coefficients)
  rows <- lapply(seq_along(coefficient_names), function(column) {
    statistic <- source$statistics[, column]
    p_value <- if (!reference_available) {
      rep(NA_real_, length(statistic))
    } else if (is.infinite(reference_df)) {
      2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
    } else {
      2 * stats::pt(abs(statistic), df = reference_df, lower.tail = FALSE)
    }
    data.frame(
      location = seq_len(nrow(source$coefficients)),
      coefficient = coefficient_names[column],
      estimate = source$coefficients[, column],
      standard_error = source$standard_errors[, column],
      statistic = statistic,
      p_value = p_value,
      selected = column == 1L | abs(source$coefficients[, column]) > 1e-10,
      stringsAsFactors = FALSE
    )
  })
  local <- do.call(rbind, rows)
  if (adjust_scope == "global") {
    local$p_adjusted <- stats::p.adjust(local$p_value, method = adjust)
  } else {
    local$p_adjusted <- stats::ave(
      local$p_value,
      local$coefficient,
      FUN = function(value) stats::p.adjust(value, method = adjust)
    )
  }
  local$significant <- ifelse(is.finite(local$p_adjusted),
    local$p_adjusted < alpha, NA)
  coefficient_summary <- do.call(rbind, lapply(
    split(local, local$coefficient),
    function(value) {
      usable <- is.finite(value$p_adjusted)
      statistics <- abs(value$statistic[is.finite(value$statistic)])
      data.frame(
        coefficient = value$coefficient[1L],
        selected_percent = 100 * mean(value$selected),
        significant_percent = if (any(usable))
          100 * mean(value$significant[usable]) else NA_real_,
        median_abs_statistic = if (length(statistics))
          stats::median(statistics) else NA_real_,
        min_adjusted_p = if (any(usable)) min(value$p_adjusted[usable]) else NA_real_,
        available_tests = sum(usable),
        unavailable_tests = sum(!usable),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(coefficient_summary) <- NULL
  selection <- object$selection_context
  basis <- if (identical(object$method, "GWR"))
    "unpenalized GWR linear smoother" else source$basis
  result <- list(
    local = local,
    summary = coefficient_summary,
    status = if (!any(is.finite(local$p_value))) "unavailable" else
      if (any(!is.finite(local$p_value))) "partially_available" else "calculated",
    settings = list(
      alpha = alpha,
      df_rule = df,
      reference_df = reference_df,
      adjust = adjust,
      adjust_scope = adjust_scope,
      estimator = estimator,
      basis = basis,
      exploratory = TRUE,
      reference_available = reference_available,
      reference_calibrated = FALSE,
      selection_conditioned = !is.null(selection) || gwrs_is_sparse_selection_fit(object),
      selection = selection
    ),
    notes = c(
      "Two-sided local coefficient-zero reference calculations, distinct from F3 spatial constancy.",
      "Multiplicity adjustment does not calibrate the reference p-values or account for response-based model/tuning selection.",
      "Overlapping neighborhoods and smoothing bias limit interpretation; no calibrated post-selection or FDR guarantee is claimed.",
      if (!is.null(selection) && !isTRUE(selection$search_complete))
        "Provisional selected model: the search is incomplete.")
  )
  class(result) <- "gwrs_local_inference"
  result
}

#' @export
print.gwrs_local_inference <- function(x, ...) {
  cat("gwrs exploratory local coefficient inference\n")
  cat("  estimator:", x$settings$estimator, "\n")
  cat("  basis:", x$settings$basis, "\n")
  cat("  reference:", if (is.na(x$settings$reference_df)) {
    "unavailable"
  } else if (is.infinite(x$settings$reference_df)) {
    "normal"
  } else {
    paste0("t(df = ", signif(x$settings$reference_df, 6), ")")
  }, "\n")
  cat("  adjustment:", x$settings$adjust,
      "(", x$settings$adjust_scope, ")\n")
  print(x$summary, row.names = FALSE)
  cat("Note: p-values are exploratory and conditional on overlapping local neighborhoods.\n")
  cat("Multiplicity adjustment does not establish calibration or account for model/tuning selection.\n")
  invisible(x)
}
