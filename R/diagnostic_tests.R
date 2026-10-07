#' GWR-style diagnostic F-tests
#'
#' Computes F1 (local model versus global OLS), F2 (joint coefficient
#' stationarity), and coefficient-specific F3 spatial-variability diagnostics.
#' Exact first and second operator moments are computed in neighbor blocks.
#' The F reference distributions are moment-matched approximations for a fixed,
#' identified Gaussian GWR smoother, not exact finite-sample distributions.
#' Penalized fits require `allow_approximate = TRUE` and remain exploratory.
#'
#' @param object A least-squares `gwrs_fit` retaining `x` and `y`.
#' @param allow_approximate Permit conditional active-set tests for a penalized
#'   fit.
#' @param ... Unused.
#' @param block_size Positive integer number of operator columns per block.
#'   Smaller blocks reduce workspace without changing the target moments.
#' @param memory_limit_mb Positive finite workspace budget in MiB. A conservative
#'   estimate is checked before operator allocation. This is not a hard process
#'   RSS limit; the R session, allocator and numerical libraries use other memory.
#'   No stochastic or diagonal-only fallback is used when the budget is exceeded.
#'
#' @return A `gwrs_diagnostic_tests` object. `diagnostics` includes full operator
#'   moments, coefficient moments and workspace settings. Unidentified global
#'   OLS references or degenerate variability produce unavailable (`NA`) tests.
#'   A singular/non-positive local Hessian or a nonconvex penalty knot raises
#'   an informative error instead of silently substituting a pseudoinverse.
#'
#' @details F3 tests spatial variation of an entire coefficient surface; it is
#'   not a pointwise test of a zero coefficient. Exact moments do not correct
#'   selection bias, response-dependent tuning or residual spatial dependence.
#'   Selection maps should label gray as "not selected", with inferential
#'   decisions displayed separately. Exact all-coefficient second moments can
#'   be expensive even though no dense spatial smoother is stored.
#' @references Leung, Y., Mei, C.-L. and Zhang, W.-X. (2000). Statistical Tests
#'   for Spatial Nonstationarity Based on the Geographically Weighted Regression
#'   Model. Environment and Planning A, 32, 9-32. \doi{10.1068/a3162}.
#'
#' @examples
#' quake <- datasets::quakes[1:36, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_fit(
#'   x, quake$mag, coords = coords, k = 24,
#'   control = gwrs_control(n_threads = 1)
#' )
#' tests <- gwr_diagnostic_tests(fit)
#' tests$global_tests
#' tests$coefficient_tests
#'
#' @export
gwr_diagnostic_tests <- function(object,
                                 allow_approximate = FALSE,
                                 ...,
                                 block_size = 64L,
                                 memory_limit_mb = 512) {
  if (!inherits(object, "gwr_sl_fit") &&
      !inherits(object, "gwr_nonconvex_fit")) {
    stop("F1/F2/F3 diagnostics are defined for Gaussian least-squares GWR fits.",
         call. = FALSE)
  }
  validate_gwrs_diagnostic_object(object, require_x = TRUE)
  if (is.null(object$y)) {
    stop("Refit with `gwrs_control(keep_data = TRUE)` for diagnostic F-tests.",
         call. = FALSE)
  }
  if (object$lambda > 0 && !isTRUE(allow_approximate)) {
    stop(
      "Classical F1/F2/F3 reference distributions do not account for penalization or local selection; set `allow_approximate = TRUE` for labelled conditional diagnostics.",
      call. = FALSE
    )
  }
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      !is.finite(block_size) || block_size < 1 ||
      block_size != floor(block_size) || block_size > .Machine$integer.max) {
    stop("`block_size` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(memory_limit_mb) || length(memory_limit_mb) != 1L ||
      !is.finite(memory_limit_mb) || memory_limit_mb <= 0) {
    stop("`memory_limit_mb` must be positive and finite.", call. = FALSE)
  }
  n <- length(object$y)
  p <- ncol(object$x)
  input_mb <- 8 * (8 * as.double(n) * (p + 1) +
                    4 * length(object$neighbors$index)) / 1024^2
  if (input_mb > memory_limit_mb) {
    stop("F-test input workspace exceeds `memory_limit_mb`; no approximation was used.",
         call. = FALSE)
  }
  gwr_diagnostics(object) # Preserve the requirement for retained diagnostics.
  global_design <- cbind(1, object$x)
  global_fit <- stats::lm.fit(global_design, object$y)
  rss_ols <- sum(global_fit$residuals^2)
  df_ols <- n - global_fit$rank
  global_basis <- matrix(0, n, if (df_ols > 0) global_fit$rank else 0L)
  if (ncol(global_basis)) {
    at <- seq_len(ncol(global_basis))
    global_basis[cbind(at, at)] <- 1
    global_basis <- base::qr.qy(global_fit$qr, global_basis)
  }
  threads <- object$control$n_threads
  if (threads < 0L) threads <- RcppParallel::defaultNumThreads()
  standardized_x <- apply_standardization(
    object$x, object$x_center, object$x_scale
  )
  kernel <- kernel_code(object$kernel)
  bandwidth <- diagnostic_bandwidth(object)
  penalty_mode <- switch(
    object$penalty %||% "convex",
    mcp = 1L,
    scad = 2L,
    0L
  )
  components <- cpp_gwr_f_test_components(
    standardized_x,
    t(object$coefficients),
    object$neighbors$index,
    object$neighbors$distance,
    object$x_center,
    object$x_scale,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    object$lambda,
    object[["alpha"]] %||% 1,
    object[["d"]] %||% 0,
    penalty_mode,
    object$gamma %||% 0,
    threads,
    object$control$grain_size,
    global_basis,
    as.integer(block_size),
    as.double(memory_limit_mb)
  )
  rss_local <- sum(object$residuals^2)
  l_delta1 <- components$delta1
  l_delta2 <- components$delta2
  delta1 <- l_delta1

  f1 <- (rss_local / l_delta1) / (rss_ols / df_ols)
  f1_df1 <- l_delta1^2 / l_delta2
  f1_df2 <- df_ols
  if (l_delta1 <= 64 * .Machine$double.eps * n) {
    f1 <- f1_df1 <- NA_real_
  }
  f1_p <- gwr_f_reference_p(f1, f1_df1, f1_df2, lower.tail = TRUE)
  f2 <- ((rss_ols - rss_local) / (df_ols - l_delta1)) /
    (rss_ols / df_ols)
  f2_df1 <- (df_ols - l_delta1)^2 /
    (df_ols - 2 * (l_delta1 - components$trQH) + l_delta2)
  f2_df2 <- df_ols
  if (df_ols - l_delta1 <= 64 * .Machine$double.eps * max(n, l_delta1)) {
    f2 <- f2_df1 <- NA_real_
  }
  f2_p <- gwr_f_reference_p(f2, f2_df1, f2_df2, lower.tail = FALSE)

  coefficient_names <- colnames(object$coefficients)
  sigma2_delta1 <- if (is.finite(f1_df1)) rss_local / delta1 else NA_real_
  f3 <- lapply(seq_along(coefficient_names), function(column) {
    beta <- object$coefficients[, column]
    variability <- sum((beta - mean(beta))^2) / length(beta)
    gamma1 <- components$gamma1[column]
    gamma2 <- components$gamma2[column]
    df1 <- gamma1^2 / gamma2
    statistic <- (variability / gamma1) / sigma2_delta1
    data.frame(
      coefficient = coefficient_names[column],
      statistic = statistic,
      df1 = df1,
      df2 = f1_df1,
      p_value = gwr_f_reference_p(statistic, df1, f1_df1, lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  })
  global_tests <- data.frame(
    test = c("F1", "F2"),
    null_hypothesis = c(
      "The local model does not improve on global OLS",
      "All coefficient surfaces are spatially stationary"
    ),
    statistic = c(f1, f2),
    df1 = c(f1_df1, f2_df1),
    df2 = c(f1_df2, f2_df2),
    p_value = c(f1_p, f2_p),
    stringsAsFactors = FALSE
  )
  global_tests <- sanitize_gwr_f_tests(global_tests)
  coefficient_tests <- sanitize_gwr_f_tests(do.call(rbind, f3))
  result <- list(
    global_tests = global_tests,
    coefficient_tests = coefficient_tests,
    diagnostics = list(
      rss_ols = rss_ols,
      rss_local = rss_local,
      df_ols = df_ols,
      trS = components$trS,
      trStS = components$trStS,
      delta1 = delta1,
      L_delta1 = l_delta1,
      L_delta2 = l_delta2,
      trQH = components$trQH,
      coefficient_moments = data.frame(
        coefficient = coefficient_names,
        gamma1 = as.double(components$gamma1),
        gamma2 = as.double(components$gamma2)
      ),
      moment_method = components$moment_method,
      block_size = components$block_size,
      threads = components$threads,
      workspace_estimate_mb = components$workspace_estimate_mb,
      memory_limit_mb = memory_limit_mb
    ),
    validity = list(
      classical = object$lambda <= 0,
      basis = if (object$lambda <= 0) {
        "classical Gaussian GWR smoother"
      } else if (inherits(object, "gwr_nonconvex_fit")) {
        "conditional local-stationary approximation for a nonconvex penalized smoother"
      } else {
        "conditional active-set approximation for a penalized smoother"
      },
      exploratory = object$lambda > 0,
      reference_distribution = "moment-matched F approximation",
      exact_operator_moments = TRUE,
      post_selection_validated = FALSE
    )
  )
  class(result) <- "gwrs_diagnostic_tests"
  result
}

gwr_f_reference_p <- function(statistic, df1, df2, lower.tail) {
  if (!is.finite(statistic) || statistic < 0 ||
      !is.finite(df1) || df1 <= 0 || !is.finite(df2) || df2 <= 0) {
    return(NA_real_)
  }
  stats::pf(statistic, df1, df2, lower.tail = lower.tail)
}

sanitize_gwr_f_tests <- function(table) {
  bad <- !is.finite(table$statistic) | table$statistic < 0 |
    !is.finite(table$df1) | table$df1 <= 0 |
    !is.finite(table$df2) | table$df2 <= 0
  table$statistic[bad] <- NA_real_
  table$df1[bad] <- NA_real_
  table$df2[bad] <- NA_real_
  table$p_value[bad] <- NA_real_
  table
}

#' @export
print.gwrs_diagnostic_tests <- function(x, ...) {
  cat("gwrs GWR-style diagnostic F-tests\n")
  cat("  basis:", x$validity$basis, "\n\n")
  print(x$global_tests, row.names = FALSE)
  cat("\nCoefficient-specific spatial variability tests\n")
  print(x$coefficient_tests, row.names = FALSE)
  if (isTRUE(x$validity$exploratory)) {
    cat("\nNote: penalized-fit tests are conditional exploratory approximations.\n")
  }
  invisible(x)
}
