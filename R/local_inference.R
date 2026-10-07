#' Exploratory local coefficient inference
#'
#' Computes local standard errors, test statistics, raw p-values, adjusted
#' p-values, and significance flags from a fit created with
#' `gwrs_control(diagnostics = "full")`. Sparse fits default to an unpenalized
#' refit on each selected local active set; all reported p-values remain
#' exploratory because selection and overlapping neighborhoods induce
#' dependence.
#'
#' @param object A least-squares `gwrs_fit`.
#' @param alpha Significance level.
#' @param df Reference degrees-of-freedom rule.
#' @param adjust Multiplicity adjustment passed to [stats::p.adjust()].
#' @param adjust_scope Adjust separately by coefficient or over all local tests.
#' @param estimator Use the penalized active-set derivative, the local
#'   post-selection refit, or choose automatically.
#' @param ... Unused.
#'
#' @return A `gwrs_local_inference` object.
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
      "Local inference requires `gwrs_control(diagnostics = \"full\")`.",
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
  if (!is.infinite(reference_df) &&
      (!is.finite(reference_df) || reference_df <= 0)) {
    warning(
      "A positive reference degree of freedom was unavailable; using the normal approximation.",
      call. = FALSE
    )
    reference_df <- Inf
  }
  coefficient_names <- colnames(source$coefficients)
  rows <- lapply(seq_along(coefficient_names), function(column) {
    statistic <- source$statistics[, column]
    p_value <- if (is.infinite(reference_df)) {
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
  local$significant <- is.finite(local$p_adjusted) &
    local$p_adjusted < alpha
  coefficient_summary <- do.call(rbind, lapply(
    split(local, local$coefficient),
    function(value) data.frame(
      coefficient = value$coefficient[1L],
      selected_percent = 100 * mean(value$selected),
      significant_percent = 100 * mean(value$significant, na.rm = TRUE),
      median_abs_statistic = stats::median(abs(value$statistic), na.rm = TRUE),
      min_adjusted_p = suppressWarnings(
        min(value$p_adjusted, na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    )
  ))
  rownames(coefficient_summary) <- NULL
  result <- list(
    local = local,
    summary = coefficient_summary,
    settings = list(
      alpha = alpha,
      df_rule = df,
      reference_df = reference_df,
      adjust = adjust,
      adjust_scope = adjust_scope,
      estimator = estimator,
      basis = source$basis,
      exploratory = TRUE
    )
  )
  class(result) <- "gwrs_local_inference"
  result
}

#' @export
print.gwrs_local_inference <- function(x, ...) {
  cat("gwrs exploratory local coefficient inference\n")
  cat("  estimator:", x$settings$estimator, "\n")
  cat("  basis:", x$settings$basis, "\n")
  cat("  reference:", if (is.infinite(x$settings$reference_df)) {
    "normal"
  } else {
    paste0("t(df = ", signif(x$settings$reference_df, 6), ")")
  }, "\n")
  cat("  adjustment:", x$settings$adjust,
      "(", x$settings$adjust_scope, ")\n")
  print(x$summary, row.names = FALSE)
  cat("Note: p-values are exploratory and conditional on overlapping local neighborhoods.\n")
  invisible(x)
}
