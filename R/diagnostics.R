#' Extract statistical diagnostics from a gwrs fit
#'
#' Statistical diagnostics are computed without constructing a dense spatial
#' weight or smoother matrix. For sparse penalties, effective-degrees-of-
#' freedom quantities are conditional on the fitted active sets.
#'
#' @param object A fitted `gwrs_fit` object.
#' @param ... Unused.
#'
#' @return A `gwrs_diagnostics` object.
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_fit(
#'   x, quake$mag, coords = coords, k = 20,
#'   control = gwrs_control(n_threads = 1)
#' )
#' diagnostics <- gwr_diagnostics(fit)
#' diagnostics$global[c("rss", "r_squared", "aicc", "gcv")]
#' head(diagnostics$local)
#'
#' @export
gwr_diagnostics <- function(object, ...) {
  if (!inherits(object, "gwrs_fit")) {
    stop("`object` must inherit from 'gwrs_fit'.", call. = FALSE)
  }
  diagnostics <- object$model_diagnostics
  if (is.null(diagnostics)) {
    stop(
      "Statistical diagnostics are unavailable; refit with diagnostics enabled in `gwrs_control()`.",
      call. = FALSE
    )
  }
  diagnostics
}

#' Compact diagnostic extractors
#'
#' @inheritParams gwr_diagnostics
#' @return A single numeric diagnostic value.
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_fit(
#'   x, quake$mag, coords = coords, k = 20,
#'   control = gwrs_control(n_threads = 1)
#' )
#' c(rss = gwr_rss(fit), aicc = gwr_aicc(fit), cv = gwr_cv_score(fit))
#'
#' @name gwr_diagnostic_extractors
NULL

#' @rdname gwr_diagnostic_extractors
#' @export
gwr_rss <- function(object, ...) {
  if (inherits(object, "gwrs_cv")) {
    if (is.null(object$fit)) {
      stop("This cross-validation object was created with `refit = FALSE`.",
           call. = FALSE)
    }
    return(sum(object$fit$residuals^2))
  }
  gwr_diagnostics(object)$global$rss
}

#' @rdname gwr_diagnostic_extractors
#' @export
gwr_aicc <- function(object, ...) {
  gwr_diagnostics(object)$global$aicc
}

#' @rdname gwr_diagnostic_extractors
#' @export
gwr_cv_score <- function(object, ...) {
  if (inherits(object, "gwrs_cv")) {
    return(as.double(object$best[[object$metric]]))
  }
  gwr_diagnostics(object)$global$loo_cv
}

#' @export
print.gwrs_diagnostics <- function(x, ...) {
  cat("gwrs statistical diagnostics\n")
  cat("  basis:", x$validity$basis, "\n")
  values <- x$global[c(
    "rss", "rmse", "mae", "r_squared", "adjusted_r_squared",
    "trS", "trStS", "enp", "edf", "sigma2", "aic", "aicc", "bic",
    "gcv", "loo_cv"
  )]
  print(round(unlist(values), 6L))
  if (!is.null(x$validity$note)) cat("  note:", x$validity$note, "\n")
  invisible(x)
}

#' @export
logLik.gwrs_fit <- function(object, ...) {
  diagnostics <- gwr_diagnostics(object)
  if (!isTRUE(diagnostics$validity$information_criteria) ||
      !is.finite(diagnostics$global$deviance)) {
    stop(
      paste0(
        "A Gaussian likelihood/information criterion is not defined for ",
        object$method, "."
      ),
      call. = FALSE
    )
  }
  value <- -0.5 * diagnostics$global$deviance
  structure(
    value,
    df = diagnostics$global$model_df,
    nobs = length(object$residuals),
    class = "logLik"
  )
}

compute_gwr_model_diagnostics <- function(prepared_x,
                                          y,
                                          fitted,
                                          residuals,
                                          raw_coefficients,
                                          coefficients,
                                          neighbors,
                                          kernel,
                                          bandwidth,
                                          lambda,
                                          alpha,
                                          d,
                                          x_center,
                                          x_scale,
                                          control) {
  level <- control$diagnostic_level %||%
    if (isTRUE(control$diagnostics)) "standard" else "none"
  if (identical(level, "none")) return(NULL)
  full <- identical(level, "full")
  core <- cpp_gwr_model_diagnostics(
    prepared_x,
    y,
    fitted,
    raw_coefficients,
    neighbors$index,
    neighbors$distance,
    x_center,
    x_scale,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    as.double(lambda),
    as.double(alpha),
    as.double(d),
    0L,
    0,
    full,
    control$n_threads,
    control$grain_size
  )
  assemble_gwr_model_diagnostics(
    core = core,
    y = y,
    fitted = fitted,
    residuals = residuals,
    coefficients = coefficients,
    x_center = x_center,
    x_scale = x_scale,
    lambda = lambda,
    alpha = alpha,
    d = d,
    level = level
  )
}

assemble_gwr_model_diagnostics <- function(core,
                                           y,
                                           fitted,
                                           residuals,
                                           coefficients,
                                           x_center,
                                           x_scale,
                                           lambda,
                                           alpha,
                                           d,
                                           level) {
  global <- gaussian_smoother_diagnostics(
    y, fitted, residuals, core$hat_diag, core$s2_diag
  )
  basis <- if (lambda <= 0) {
    "classical GWR linear smoother"
  } else if (alpha <= 0) {
    "exact linear penalized smoother"
  } else {
    "conditional active-set smoother"
  }
  note <- if (lambda > 0 && alpha > 0) {
    paste(
      "Degrees of freedom and information criteria condition on the selected",
      "local active sets; they are intended for model comparison, not exact",
      "post-selection inference."
    )
  } else if (lambda > 0) {
    "Trace quantities are exact for the fitted linear smoother; penalized coefficient inference remains bias-aware and exploratory."
  } else {
    "Local coefficient tests remain exploratory because neighboring local fits are dependent."
  }
  local <- data.frame(
    location = seq_along(y),
    leverage = as.double(core$hat_diag),
    smoother_row_norm2 = as.double(core$s2_diag),
    local_r_squared = as.double(core$local_r2),
    effective_sample_size = as.double(core$local_neff),
    positive_weights = as.integer(core$local_n_positive),
    condition_number = as.double(core$local_kappa),
    local_rank = as.integer(core$local_rank),
    operator_solver = status_labels(core$operator_status)
  )
  result <- list(
    level = level,
    global = global,
    local = local,
    hat_diag = as.double(core$hat_diag),
    s2_diag = as.double(core$s2_diag),
    validity = list(
      basis = basis,
      exact_linear = lambda <= 0 || alpha <= 0,
      conditional = lambda > 0 && alpha > 0,
      information_criteria = TRUE,
      note = note
    ),
    penalty = list(lambda = lambda, alpha = alpha, d = d),
    inference = NULL,
    post_selection = NULL
  )
  if (!identical(level, "full")) {
    class(result) <- "gwrs_diagnostics"
    return(result)
  }

  standard_variances <- global$sigma2 * as.matrix(core$coefficient_map_ss)
  standard_variances[is.finite(standard_variances) & standard_variances < 0] <- 0
  standard_errors <- sqrt(standard_variances)
  colnames(standard_errors) <- colnames(coefficients)
  statistics <- coefficients / standard_errors
  statistics[!is.finite(statistics)] <- NA_real_
  result$inference <- list(
    basis = "penalized active-set derivative",
    coefficients = coefficients,
    standard_errors = standard_errors,
    statistics = statistics
  )

  post_coefficients <- backtransform_coefficients(
    core$post_coefficients, x_center, x_scale
  )
  colnames(post_coefficients) <- colnames(coefficients)
  post_residuals <- y - as.double(core$post_fitted)
  post_global <- gaussian_smoother_diagnostics(
    y, as.double(core$post_fitted), post_residuals,
    as.double(core$post_hat_diag), as.double(core$post_s2_diag)
  )
  post_variances <- post_global$sigma2 *
    as.matrix(core$post_coefficient_map_ss)
  post_variances[is.finite(post_variances) & post_variances < 0] <- 0
  post_standard_errors <- sqrt(post_variances)
  colnames(post_standard_errors) <- colnames(coefficients)
  post_statistics <- post_coefficients / post_standard_errors
  post_statistics[!is.finite(post_statistics)] <- NA_real_
  result$post_selection <- list(
    basis = "unpenalized refit on each selected local active set",
    coefficients = post_coefficients,
    fitted.values = as.double(core$post_fitted),
    residuals = post_residuals,
    global = post_global,
    hat_diag = as.double(core$post_hat_diag),
    s2_diag = as.double(core$post_s2_diag),
    standard_errors = post_standard_errors,
    statistics = post_statistics
  )
  class(result) <- "gwrs_diagnostics"
  result
}

gaussian_smoother_diagnostics <- function(y,
                                          fitted,
                                          residuals,
                                          hat_diag,
                                          s2_diag) {
  n <- length(y)
  rss <- sum(residuals^2)
  tss <- sum((y - mean(y))^2)
  trS <- sum(hat_diag)
  trStS <- sum(s2_diag)
  enp <- 2 * trS - trStS
  edf <- n - enp
  sigma2 <- if (is.finite(edf) && edf > .Machine$double.eps) {
    rss / edf
  } else {
    NA_real_
  }
  safe_rss <- max(rss, .Machine$double.xmin)
  deviance <- n * (log(safe_rss / n) + log(2 * pi) + 1)
  model_df <- trS + 1
  aic <- deviance + 2 * model_df
  aicc_denominator <- n - trS - 2
  aicc <- if (aicc_denominator > 0) {
    deviance + 2 * n * model_df / aicc_denominator
  } else {
    Inf
  }
  bic <- deviance + log(n) * model_df
  gcv_denominator <- 1 - trS / n
  gcv <- if (gcv_denominator > 0) {
    (rss / n) / gcv_denominator^2
  } else {
    Inf
  }
  loo_denominator <- 1 - hat_diag
  valid_loo <- is.finite(loo_denominator) & abs(loo_denominator) > 1e-10
  loo_cv <- if (all(valid_loo)) {
    sum((residuals / loo_denominator)^2)
  } else {
    Inf
  }
  r_squared <- if (tss > 0) 1 - rss / tss else NA_real_
  adjusted_r_squared <- if (is.finite(r_squared) && edf > 0) {
    1 - (1 - r_squared) * (n - 1) / edf
  } else {
    NA_real_
  }
  list(
    n = n,
    rss = rss,
    tss = tss,
    rmse = sqrt(mean(residuals^2)),
    mae = mean(abs(residuals)),
    bias = mean(residuals),
    r_squared = r_squared,
    adjusted_r_squared = adjusted_r_squared,
    trS = trS,
    trStS = trStS,
    enp = enp,
    edf = edf,
    sigma2 = sigma2,
    deviance = deviance,
    logLik = -0.5 * deviance,
    model_df = model_df,
    aic = aic,
    aicc = aicc,
    bic = bic,
    gcv = gcv,
    loo_cv = loo_cv
  )
}

compute_descriptive_model_diagnostics <- function(y,
                                                   fitted,
                                                   residuals,
                                                   neighbors,
                                                   kernel,
                                                   bandwidth,
                                                   control,
                                                   method,
                                                   extra = list()) {
  level <- control$diagnostic_level %||%
    if (isTRUE(control$diagnostics)) "standard" else "none"
  if (identical(level, "none")) return(NULL)
  local <- cpp_gwr_local_residual_diagnostics(
    y,
    fitted,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    control$n_threads,
    control$grain_size
  )
  rss <- sum(residuals^2)
  tss <- sum((y - mean(y))^2)
  global <- c(
    list(
      n = length(y),
      rss = rss,
      tss = tss,
      rmse = sqrt(mean(residuals^2)),
      mae = mean(abs(residuals)),
      median_absolute_error = stats::median(abs(residuals)),
      bias = mean(residuals),
      r_squared = if (tss > 0) 1 - rss / tss else NA_real_,
      adjusted_r_squared = NA_real_,
      trS = NA_real_,
      trStS = NA_real_,
      enp = NA_real_,
      edf = NA_real_,
      sigma2 = NA_real_,
      deviance = NA_real_,
      logLik = NA_real_,
      model_df = NA_real_,
      aic = NA_real_,
      aicc = NA_real_,
      bic = NA_real_,
      gcv = NA_real_,
      loo_cv = NA_real_
    ),
    extra
  )
  result <- list(
    level = level,
    global = global,
    local = data.frame(
      location = seq_along(y),
      local_r_squared = as.double(local$local_r2),
      local_rmse = as.double(local$local_rmse),
      local_mae = as.double(local$local_mae),
      effective_sample_size = as.double(local$effective_sample_size),
      positive_weights = as.integer(local$positive_weights)
    ),
    validity = list(
      basis = paste("descriptive residual diagnostics for", method),
      exact_linear = FALSE,
      conditional = FALSE,
      information_criteria = FALSE,
      note = paste(
        "Classical Gaussian smoother traces, AICc, local t-tests, and F-tests",
        "are not reported because they do not match this estimator's loss or",
        "joint operator."
      )
    ),
    inference = NULL,
    post_selection = NULL
  )
  class(result) <- "gwrs_diagnostics"
  result
}
