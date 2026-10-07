#' Fast nonconvex geographically weighted regression
#'
#' Fits Gaussian geographically weighted regression with either the smoothly
#' clipped absolute deviation (SCAD) or minimax concave penalty (MCP). At
#' target location `i`, the profiled local objective is
#' \deqn{\frac{1}{2c_i}(y-X\beta_i)^T W_i(y-X\beta_i)+
#' \sum_{m=1}^p p_{\lambda,\gamma}(|\beta_{im}|),}
#' where `c_i` is the sum of local kernel weights. The intercept is locally
#' profiled and is never penalized. Predictor coefficients are optimized in a
#' parallel RcppArmadillo core using exact one-coordinate SCAD/MCP minimizers.
#'
#' With `control = gwrs_control(nonconvex_solver = "guarded_block")`, the
#' same objective also admits guarded pair and small flat-tail QR proposals.
#' This is an optimizer choice, not another penalty. Solutions can differ
#' because the objective is nonconvex. A failed guarded fit raises a
#' `gwrs_nonconvex_convergence_error` condition containing `solver_result`;
#' it is not used to construct statistical diagnostics or predictions.
#'
#' @inheritParams gwr_sl_fit
#' @param penalty Nonconvex penalty, either `"scad"` or `"mcp"`.
#' @param gamma Concavity parameter. The default is `3.7` for SCAD and `3` for
#'   MCP. SCAD requires `gamma > 2`; MCP requires `gamma > 1`.
#' @param ... Additional arguments passed by [gwr_scad_fit()] or
#'   [gwr_mcp_fit()] to `gwr_nonconvex_fit()`.
#'
#' @return An object inheriting from `gwr_nonconvex_fit` and `gwrs_fit`.
#'
#' @references
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. *Journal of the American Statistical
#' Association*, 96, 1348--1360. \doi{10.1198/016214501753382273}.
#'
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. *The Annals of Statistics*, 38, 894--942.
#' \doi{10.1214/09-AOS729}.
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_mcp_fit(
#'   x, quake$mag, coords = coords, k = 20, lambda = 0.05,
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#' fit
#' head(coef(fit))
#'
#' @export
gwr_nonconvex_fit <- function(x,
                              y,
                              coords = NULL,
                              neighbors = NULL,
                              k = NULL,
                              kernel = c("bisquare", "gaussian", "exponential",
                                         "tricube", "boxcar"),
                              bandwidth = NULL,
                              lambda,
                              penalty = c("scad", "mcp"),
                              gamma = NULL,
                              standardize = TRUE,
                              control = gwrs_control()) {
  call <- match.call()
  x <- as_design_matrix(x)
  y <- as_response(y, nrow(x))
  control <- validate_control(control)
  neighbors <- resolve_neighbors(x, coords, neighbors, k)
  kernel <- kernel_code(kernel)
  bandwidth <- resolve_bandwidth(bandwidth)
  specification <- resolve_nonconvex_penalty(penalty, gamma)
  if (missing(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda < 0) {
    stop("`lambda` must be one nonnegative finite number.", call. = FALSE)
  }

  prepared <- standardize_design(x, standardize)
  raw <- cpp_gwr_nonconvex_path_predict(
    prepared$x,
    y,
    prepared$x,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    as.double(lambda),
    specification$code,
    specification$gamma,
    control$tolerance,
    control$max_iterations,
    FALSE,
    TRUE,
    FALSE,
    control$n_threads,
    control$grain_size,
    nonconvex_guarded(control)
  )
  nonconvex_require_fit(raw, control)
  fitted <- as.double(raw$predictions[, 1L])
  if (any(!is.finite(fitted))) {
    stop("At least one target has no valid local weighted fit.", call. = FALSE)
  }
  raw_coefficients <- raw$coefficients[, , 1L, drop = FALSE][, , 1L]
  coefficients <- backtransform_coefficients(
    raw_coefficients, prepared$center, prepared$scale
  )
  predictor_names <- colnames(x) %||% paste0("x", seq_len(ncol(x)))
  colnames(coefficients) <- c("(Intercept)", predictor_names)
  residuals <- y - fitted

  solver_diagnostics <- data.frame(
    location = seq_len(nrow(x)),
    converged = as.logical(raw$converged[, 1L]),
    iterations = as.integer(raw$iterations[, 1L]),
    solver = nonconvex_status_labels(raw$status),
    weight_sum = as.double(raw$weight_sum),
    bandwidth = as.double(raw$bandwidth),
    objective = as.double(raw$objective[, 1L]),
    max_stationarity_violation = as.double(raw$stationarity[, 1L]),
    max_coordinate_gap = as.double(raw$coordinate_gap[, 1L]),
    active_coordinates = as.integer(raw$active_count[, 1L]),
    stringsAsFactors = FALSE
  )
  if (nonconvex_guarded(control)) {
    if (lambda > 0) solver_diagnostics$solver[raw$status == 0L] <- "guarded_block"
    solver_diagnostics$solver_mode <- "guarded_block"
    solver_diagnostics$state <- drop(nonconvex_state_labels(raw$path_state))
    solver_diagnostics$pair_accepted <- raw$pair_accepted[, 1L]
    solver_diagnostics$block_accepted <- raw$block_accepted[, 1L]
    solver_diagnostics$block_checks <- raw$block_checks[, 1L]
  }

  model_diagnostics <- compute_gwr_nonconvex_model_diagnostics(
    prepared_x = prepared$x,
    y = y,
    fitted = fitted,
    residuals = residuals,
    raw_coefficients = raw_coefficients,
    coefficients = coefficients,
    neighbors = neighbors,
    kernel = kernel,
    bandwidth = bandwidth,
    lambda = as.double(lambda),
    penalty = specification$name,
    gamma = specification$gamma,
    x_center = prepared$center,
    x_scale = prepared$scale,
    control = control
  )

  object <- list(
    call = call,
    method = if (lambda == 0) "GWR" else specification$method,
    penalty = specification$name,
    gamma = specification$gamma,
    coefficients = coefficients,
    fitted.values = fitted,
    residuals = residuals,
    lambda = as.double(lambda),
    kernel = kernel$name,
    adaptive = bandwidth$adaptive,
    bandwidth = if (bandwidth$adaptive) NULL else bandwidth$fixed,
    neighbors = neighbors,
    standardize = isTRUE(standardize),
    x_center = prepared$center,
    x_scale = prepared$scale,
    diagnostics = solver_diagnostics,
    solver_diagnostics = solver_diagnostics,
    model_diagnostics = model_diagnostics,
    hat_diag = model_diagnostics$hat_diag %||% NULL,
    s2_diag = model_diagnostics$s2_diag %||% NULL,
    local_r2 = model_diagnostics$local$local_r_squared %||% NULL,
    terms = NULL,
    predictor_names = predictor_names,
    x = if (control$keep_data) x else NULL,
    y = if (control$keep_data) y else NULL,
    control = control
  )
  class(object) <- c(
    paste0("gwr_", specification$name, "_fit"),
    "gwr_nonconvex_fit", "gwrs_fit"
  )
  object
}

#' @rdname gwr_nonconvex_fit
#' @export
gwr_scad_fit <- function(x, y, ..., gamma = 3.7) {
  call <- match.call()
  fit <- gwr_nonconvex_fit(
    x, y, ..., penalty = "scad", gamma = gamma
  )
  fit$call <- call
  fit
}

#' @rdname gwr_nonconvex_fit
#' @export
gwr_mcp_fit <- function(x, y, ..., gamma = 3) {
  call <- match.call()
  fit <- gwr_nonconvex_fit(
    x, y, ..., penalty = "mcp", gamma = gamma
  )
  fit$call <- call
  fit
}

#' Formula interfaces for nonconvex GWR
#'
#' @param formula Model formula with an intercept.
#' @param data Data frame containing model variables.
#' @param coords Coordinate matrix or two column names in `data`.
#' @param ... Additional arguments passed to [gwr_nonconvex_fit()].
#' @inheritParams gwr_nonconvex_fit
#'
#' @return An object inheriting from `gwr_nonconvex_fit` and `gwrs_fit`.
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' fit <- gwr_scad(
#'   mag ~ depth + stations, quake, coords = c("long", "lat"),
#'   k = 20, lambda = 0.05,
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#' fit
#'
#' @export
gwr_nonconvex <- function(formula,
                          data,
                          coords,
                          ...,
                          penalty = c("scad", "mcp"),
                          gamma = NULL) {
  call <- match.call()
  penalty <- match.arg(penalty)
  components <- formula_model_components(formula, data, coords)
  fit <- gwr_nonconvex_fit(
    components$x,
    components$y,
    coords = components$coords,
    ...,
    penalty = penalty,
    gamma = gamma
  )
  decorate_formula_fit(fit, components, call)
}

#' @rdname gwr_nonconvex
#' @export
gwr_scad <- function(formula, data, coords, ..., gamma = 3.7) {
  gwr_nonconvex(
    formula, data, coords, ..., penalty = "scad", gamma = gamma
  )
}

#' @rdname gwr_nonconvex
#' @export
gwr_mcp <- function(formula, data, coords, ..., gamma = 3) {
  gwr_nonconvex(
    formula, data, coords, ..., penalty = "mcp", gamma = gamma
  )
}

resolve_nonconvex_penalty <- function(penalty, gamma = NULL) {
  penalty <- match.arg(penalty, c("scad", "mcp"))
  if (is.null(gamma)) gamma <- if (penalty == "scad") 3.7 else 3
  gamma <- as.double(gamma)
  if (length(gamma) != 1L || !is.finite(gamma) ||
      (penalty == "scad" && gamma <= 2) ||
      (penalty == "mcp" && gamma <= 1)) {
    stop(
      if (penalty == "scad") {
        "SCAD requires one finite `gamma > 2`."
      } else {
        "MCP requires one finite `gamma > 1`."
      },
      call. = FALSE
    )
  }
  list(
    name = penalty,
    code = if (penalty == "mcp") 0L else 1L,
    diagnostic_code = if (penalty == "mcp") 1L else 2L,
    gamma = gamma,
    method = if (penalty == "mcp") "GWR-MCP" else "GWR-SCAD"
  )
}

nonconvex_status_labels <- function(status) {
  status <- as.integer(status)
  labels <- rep("coordinate_descent", length(status))
  labels[status == 1L] <- "general_solve"
  labels[status == 2L] <- "pseudoinverse"
  labels[status == 3L] <- "zero_weights_or_failed_operator"
  labels[status == 4L] <- "invalid_neighbor"
  labels
}

compute_gwr_nonconvex_model_diagnostics <- function(prepared_x,
                                                    y,
                                                    fitted,
                                                    residuals,
                                                    raw_coefficients,
                                                    coefficients,
                                                    neighbors,
                                                    kernel,
                                                    bandwidth,
                                                    lambda,
                                                    penalty,
                                                    gamma,
                                                    x_center,
                                                    x_scale,
                                                    control) {
  level <- control$diagnostic_level %||%
    if (isTRUE(control$diagnostics)) "standard" else "none"
  if (identical(level, "none")) return(NULL)
  full <- identical(level, "full")
  specification <- resolve_nonconvex_penalty(penalty, gamma)
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
    1,
    0,
    specification$diagnostic_code,
    specification$gamma,
    full,
    control$n_threads,
    control$grain_size
  )
  result <- assemble_gwr_model_diagnostics(
    core = core,
    y = y,
    fitted = fitted,
    residuals = residuals,
    coefficients = coefficients,
    x_center = x_center,
    x_scale = x_scale,
    lambda = lambda,
    alpha = 1,
    d = 0,
    level = level
  )
  result$local$local_convex <- as.logical(core$local_convex)
  result$local$minimum_active_hessian_eigenvalue <-
    as.double(core$local_min_hessian)
  result$local$penalty_knot_count <-
    as.integer(core$local_penalty_knot_count)
  result$validity <- list(
    basis = "conditional local-stationary nonconvex smoother",
    exact_linear = FALSE,
    conditional = TRUE,
    locally_convex_fraction = mean(core$local_convex != 0L),
    information_criteria = TRUE,
    note = paste(
      "Trace, information-criterion, and coefficient-derivative quantities",
      "condition on the selected local stationary solution. Non-locally",
      "convex fits and coefficients at penalty knots require additional",
      "caution; these are not exact post-selection inference."
    )
  )
  result$penalty <- list(
    name = specification$name,
    lambda = lambda,
    gamma = specification$gamma
  )
  if (!is.null(result$inference)) {
    result$inference$basis <-
      "nonconvex penalized active-set derivative at the local stationary solution"
  }
  result
}
