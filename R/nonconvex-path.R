#' Lambda anchor for GWR-SCAD and GWR-MCP
#'
#' Computes an all-zero coordinatewise-minimum anchor after the requested
#' global predictor standardization. Unlike the Lasso stationarity anchor,
#' this function accounts for location-specific curvature and directly checks
#' the exact SCAD/MCP one-coordinate objective. This distinction matters when
#' a local curvature falls inside a concave penalty region.
#'
#' @inheritParams gwr_nonconvex_fit
#'
#' @return One nonnegative number.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' gwr_nonconvex_lambda_max(
#'   x, quake$mag, coords = coords, k = 18,
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#'
#' @export
gwr_nonconvex_lambda_max <- function(x,
                                     y,
                                     coords = NULL,
                                     neighbors = NULL,
                                     k = NULL,
                                     kernel = c("bisquare", "gaussian",
                                                "exponential", "tricube",
                                                "boxcar"),
                                     bandwidth = NULL,
                                     penalty = c("scad", "mcp"),
                                     gamma = NULL,
                                     standardize = TRUE,
                                     control = gwrs_control()) {
  x <- as_design_matrix(x)
  y <- as_response(y, nrow(x))
  control <- validate_control(control)
  neighbors <- resolve_neighbors(x, coords, neighbors, k)
  kernel <- kernel_code(kernel)
  bandwidth <- resolve_bandwidth(bandwidth)
  specification <- resolve_nonconvex_penalty(penalty, gamma)
  prepared <- standardize_design(x, standardize)
  cpp_gwr_nonconvex_lambda_max(
    prepared$x,
    y,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    specification$code,
    specification$gamma,
    control$n_threads,
    control$grain_size
  )
}

#' Warm-started GWR-SCAD and GWR-MCP paths
#'
#' Fits a decreasing nonconvex regularization path. Local neighborhoods and
#' centered design columns are reused. Warm starts, persistent active sets,
#' sequential candidate screening, and a complete inactive-coordinate
#' objective-improvement scan reduce work while retaining a coordinatewise
#' minimum check at every accepted solution.
#'
#' In guarded-block mode, a failed target stops its remaining lambda path.
#' `path_state` labels converged, failed and unattempted target/lambda entries;
#' `solver_diagnostics` retains objectives, stopping diagnostics and compact
#' pair/block counters. Failed coefficients, if requested, are diagnostic
#' iterates only; invalid fitted values and incomplete-path statistical
#' summaries are unavailable. A warning identifies partial paths. The
#' coordinate mode retains the original behavior.
#'
#' @inheritParams gwr_nonconvex_fit
#' @inheritParams gwr_sl_path
#' @param ... Additional arguments passed by [gwr_scad_path()] or
#'   [gwr_mcp_path()] to `gwr_nonconvex_path()`.
#'
#' @return An object inheriting from `gwrs_nonconvex_path` and `gwrs_path`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' path <- gwr_mcp_path(
#'   x, quake$mag, coords = coords, k = 18,
#'   lambda = c(0.2, 0.08, 0.02),
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#' path
#' path$path_summary
#'
#' @export
gwr_nonconvex_path <- function(x,
                               y,
                               coords = NULL,
                               neighbors = NULL,
                               k = NULL,
                               kernel = c("bisquare", "gaussian", "exponential",
                                          "tricube", "boxcar"),
                               bandwidth = NULL,
                               lambda = NULL,
                               n_lambda = 50L,
                               lambda_min_ratio = 1e-3,
                               penalty = c("scad", "mcp"),
                               gamma = NULL,
                               standardize = TRUE,
                               screening = TRUE,
                               keep_coefficients = FALSE,
                               return_diagnostics = FALSE,
                               ebic_gamma = 0.5,
                               control = gwrs_control()) {
  call <- match.call()
  x <- as_design_matrix(x)
  y <- as_response(y, nrow(x))
  control <- validate_control(control)
  neighbors <- resolve_neighbors(x, coords, neighbors, k)
  kernel <- kernel_code(kernel)
  bandwidth <- resolve_bandwidth(bandwidth)
  specification <- resolve_nonconvex_penalty(penalty, gamma)
  if (length(ebic_gamma) != 1L || !is.finite(ebic_gamma) ||
      ebic_gamma < 0 || ebic_gamma > 1) {
    stop("`ebic_gamma` must be one number in [0, 1].", call. = FALSE)
  }
  prepared <- standardize_design(x, standardize)

  lambda_max <- NULL
  if (is.null(lambda)) {
    lambda_max <- cpp_gwr_nonconvex_lambda_max(
      prepared$x,
      y,
      neighbors$index,
      neighbors$distance,
      kernel$code,
      bandwidth$adaptive,
      bandwidth$fixed,
      specification$code,
      specification$gamma,
      control$n_threads,
      control$grain_size
    )
    lambda <- make_lambda_sequence(lambda_max, n_lambda, lambda_min_ratio)
  } else {
    lambda <- validate_lambda_path(lambda)
  }

  raw <- cpp_gwr_nonconvex_path_predict(
    prepared$x,
    y,
    prepared$x,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    lambda,
    specification$code,
    specification$gamma,
    control$tolerance,
    control$max_iterations,
    isTRUE(screening),
    isTRUE(keep_coefficients),
    isTRUE(return_diagnostics),
    control$n_threads,
    control$grain_size,
    nonconvex_guarded(control)
  )
  guarded <- nonconvex_guarded(control)
  if (any(!is.finite(raw$predictions)) && !guarded) {
    stop("At least one target has no valid local weighted fit.", call. = FALSE)
  }
  if (guarded && !all(nonconvex_valid(raw))) {
    warning("Guarded path contains failed or unattempted solutions; see path_state and solver_diagnostics.", call. = FALSE)
  }

  coefficient_path <- NULL
  predictor_names <- colnames(x) %||% paste0("x", seq_len(ncol(x)))
  if (isTRUE(keep_coefficients)) {
    coefficient_path <- backtransform_coefficient_path(
      raw$coefficients, prepared$center, prepared$scale
    )
    dimnames(coefficient_path) <- list(
      location = seq_len(nrow(x)),
      coefficient = c("(Intercept)", predictor_names),
      lambda = format(lambda, digits = 7L)
    )
  }

  path_summary <- data.frame(
    lambda = lambda,
    mean_nonzero = colMeans(raw$nonzero),
    mean_iterations = colMeans(raw$iterations),
    convergence_rate = colMeans(raw$converged != 0L),
    max_stationarity_violation = if (guarded) nonconvex_column_max(raw$stationarity) else apply(
      raw$stationarity, 2L, max, na.rm = TRUE
    ),
    max_coordinate_gap = if (guarded) nonconvex_column_max(raw$coordinate_gap) else apply(
      raw$coordinate_gap, 2L, max, na.rm = TRUE
    ),
    mean_objective = colMeans(raw$objective),
    mean_screened = colMeans(raw$screened_count),
    mean_active_coordinates = colMeans(raw$active_count),
    mean_full_scans = colMeans(raw$full_scans)
  )
  if (guarded) path_summary <- cbind(path_summary, nonconvex_solver_summary(raw))
  path_diagnostics <- NULL
  if (isTRUE(return_diagnostics)) {
    path_diagnostics <- assemble_gwr_path_diagnostics(
      y = y,
      fitted_path = raw$predictions,
      hat_diag_path = raw$hat_diag_path,
      s2_diag_path = raw$s2_diag_path,
      post_fitted_path = raw$post_fitted_path,
      post_hat_diag_path = raw$post_hat_diag_path,
      post_s2_diag_path = raw$post_s2_diag_path,
      surface_active_count = raw$surface_active_count,
      ebic_gamma = ebic_gamma
    )
    path_diagnostics$locally_convex_fraction <-
      colMeans(raw$local_convex != 0L)
    path_diagnostics$minimum_active_hessian_eigenvalue <- if (guarded) -nonconvex_column_max(-raw$minimum_hessian) else apply(
      raw$minimum_hessian, 2L, min, na.rm = TRUE
    )
    path_diagnostics$mean_penalty_knot_count <-
      colMeans(raw$penalty_knot_count)
    if (guarded) path_diagnostics[!path_summary$valid, ] <- NA_real_
    path_summary <- cbind(path_summary, path_diagnostics)
  }

  object <- list(
    call = call,
    method = paste(specification$method, "path"),
    penalty = specification$name,
    gamma = specification$gamma,
    lambda = lambda,
    lambda_max = lambda_max,
    fitted.values = raw$predictions,
    coefficients = coefficient_path,
    path_summary = path_summary,
    path_diagnostics = path_diagnostics,
    hat_diag_path = raw$hat_diag_path,
    s2_diag_path = raw$s2_diag_path,
    post_fitted_path = raw$post_fitted_path,
    post_hat_diag_path = raw$post_hat_diag_path,
    post_s2_diag_path = raw$post_s2_diag_path,
    local_convex_path = raw$local_convex,
    minimum_hessian_path = raw$minimum_hessian,
    surface_active_count = raw$surface_active_count,
    ebic_gamma = ebic_gamma,
    status = nonconvex_status_labels(raw$status),
    weight_sum = raw$weight_sum,
    local_bandwidth = raw$bandwidth,
    kernel = kernel$name,
    adaptive = bandwidth$adaptive,
    bandwidth = if (bandwidth$adaptive) NULL else bandwidth$fixed,
    neighbors = neighbors,
    standardize = isTRUE(standardize),
    x_center = prepared$center,
    x_scale = prepared$scale,
    predictor_names = predictor_names,
    x = if (control$keep_data) x else NULL,
    y = if (control$keep_data) y else NULL,
    control = control
  )
  if (guarded) {
    object$path_state <- nonconvex_state_labels(raw$path_state)
    object$solver_diagnostics <- raw[c("path_state", "converged", "iterations",
      "objective", "stationarity", "coordinate_gap", "pair_attempts", "pair_accepted",
      "pair_gain", "block_checks", "block_accepted", "block_gain", "active_block_policy")]
    object$fitted.values[!nonconvex_valid(raw)] <- NA_real_
    # Diagnostic coefficients retain failed iterates; path_state is mandatory.
    # Public fitted values never present these iterates as accepted predictions.
  }
  class(object) <- c("gwrs_nonconvex_path", "gwrs_path")
  object
}

#' @rdname gwr_nonconvex_path
#' @export
gwr_scad_path <- function(x, y, ..., gamma = 3.7) {
  gwr_nonconvex_path(
    x, y, ..., penalty = "scad", gamma = gamma
  )
}

#' @rdname gwr_nonconvex_path
#' @export
gwr_mcp_path <- function(x, y, ..., gamma = 3) {
  gwr_nonconvex_path(
    x, y, ..., penalty = "mcp", gamma = gamma
  )
}
