#' Global predictor collinearity diagnostics
#'
#' Computes one global condition number, VIFs and the predictor correlation
#' matrix, using all training observations with equal weights. No model fit is
#' needed. The RcppArmadillo core stores an n-by-p centered design and p-by-p
#' matrices, never an n-by-n spatial matrix.
#'
#' @param object A finite numeric predictor matrix/data frame without an
#'   intercept, or a `gwrs_fit` retaining its design matrix.
#' @param threshold_cn,threshold_vif Positive flag thresholds, matching
#'   [gwr_local_collinearity()].
#' @param ... Unused.
#' @details Predictors are centered and scaled by their root mean square
#'   deviations. The condition number is the square root of the largest/smallest
#'   eigenvalue ratio of the correlation matrix; the intercept is excluded.
#'   This equals the singular-value condition number of the standardized design.
#'   VIFs are the diagonal of the inverse correlation matrix. Eigenvalues at or
#'   below `max(1, largest eigenvalue) * 1e-10` are numerically rank deficient.
#'   Rank-deficient designs return infinite CN and VIFs. Constant columns make
#'   the correlation matrix undefined: CN/VIF are infinite, correlation and
#'   eigenvalues are missing, and rank is unavailable. No observations or
#'   predictors are silently removed. Threshold flags are diagnostics, not tests.
#' @return A `gwrs_global_collinearity` object including `condition_number`,
#'   `vif`, `correlation`, `eigenvalues`, `rank` and `constant_columns`.
#' @references Wheeler, D. and Tiefelsdorf, M. (2005). Multicollinearity and
#'   correlation among local regression coefficients in geographically weighted
#'   regression. Journal of Geographical Systems, 7, 161-187.
#'   \doi{10.1007/s10109-005-0155-6}.
#' @examples
#' global <- gwr_global_collinearity(datasets::mtcars[, c("wt", "hp", "disp")])
#' global
#' global$correlation
#' @export
gwr_global_collinearity <- function(object, threshold_cn = 30,
                                    threshold_vif = 10, ...) {
  if (inherits(object, "gwrs_fit")) {
    validate_gwrs_diagnostic_object(object, require_x = TRUE)
    object <- object$x
  }
  if (!is.numeric(as.matrix(object))) {
    stop("Predictors must be numeric; encode factors explicitly.", call. = FALSE)
  }
  x <- as_design_matrix(object)
  if (nrow(x) < 2L) stop("At least two rows are required.", call. = FALSE)
  for (value in list(threshold_cn, threshold_vif)) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop("Collinearity thresholds must be positive finite scalars.", call. = FALSE)
    }
  }
  raw <- cpp_gwr_global_collinearity(x)
  predictors <- colnames(x) %||% paste0("x", seq_len(ncol(x)))
  dimnames(raw$correlation) <- list(predictors, predictors)
  result <- list(
    condition_number = data.frame(condition_number = raw$condition_number,
                                  flagged = raw$condition_number > threshold_cn),
    vif = data.frame(predictor = predictors, vif = as.double(raw$vif),
                     flagged = as.double(raw$vif) > threshold_vif),
    correlation = raw$correlation, eigenvalues = as.double(raw$eigenvalues),
    rank = raw$rank, full_rank = raw$full_rank,
    constant_columns = predictors[as.integer(raw$constant_columns)],
    max_abs_correlation = raw$max_abs_correlation, n = nrow(x), p = ncol(x),
    thresholds = list(condition_number = threshold_cn, vif = threshold_vif),
    settings = list(center = TRUE, scale = TRUE, intercept = FALSE,
                    weighting = "equal", relative_rank_tolerance = 1e-10)
  )
  class(result) <- "gwrs_global_collinearity"
  result
}

#' @export
print.gwrs_global_collinearity <- function(x, ...) {
  cat("gwrs global predictor collinearity\n")
  cat("  observations/predictors:", x$n, "/", x$p, "\n")
  cat("  condition number:", format(x$condition_number$condition_number), "\n")
  cat("  predictor correlation rank:", x$rank, "\n")
  print(x$vif, row.names = FALSE)
  if (length(x$constant_columns)) {
    cat("  constant columns:", paste(x$constant_columns, collapse = ", "), "\n")
  }
  invisible(x)
}
