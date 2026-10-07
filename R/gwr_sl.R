#' Fit a fast geographically weighted scaled-lasso model
#'
#' The local objective at target location `i` is
#' \deqn{\frac{1}{2c_i}(y-X\beta_i)^T W_i(y-X\beta_i)
#' + \lambda\alpha\|\beta_i\|_1
#' + \frac{\lambda(1-\alpha)}{2}
#' \|\beta_i-d\widehat\beta_i^{GWR}\|_2^2,}
#' where `c_i` is the sum of the local kernel weights. The intercept is
#' unpenalized and profiled using weighted local means.
#'
#' @param x Numeric design matrix without an intercept.
#' @param y Numeric response vector.
#' @param coords Optional coordinate matrix used to construct neighbors.
#' @param neighbors Optional reusable object from [gwr_neighbors()].
#' @param k Number of neighbors when `neighbors` is not supplied.
#' @param kernel Spatial kernel.
#' @param bandwidth Fixed distance bandwidth. The default `NULL` uses the
#'   farthest retained neighbor as an adaptive bandwidth at each target.
#' @param lambda Nonnegative regularization strength.
#' @param alpha Mixing parameter in `[0, 1]`.
#' @param d Center-location parameter in `[0, 1]`.
#' @param standardize Whether predictors are standardized globally before the
#'   local fits. Coefficients are returned on the original scale.
#' @param control Numerical controls from [gwrs_control()].
#'
#' @return An object of class `gwrs_fit`.
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_sl_fit(
#'   x, quake$mag, coords = coords, k = 20,
#'   lambda = 0.05, alpha = 0.7, d = 0.5,
#'   control = gwrs_control(n_threads = 1)
#' )
#' fit
#' head(coef(fit))
#'
#' @export
gwr_sl_fit <- function(x,
                       y,
                       coords = NULL,
                       neighbors = NULL,
                       k = NULL,
                       kernel = c("bisquare", "gaussian", "exponential",
                                  "tricube", "boxcar"),
                       bandwidth = NULL,
                       lambda,
                       alpha = 0.5,
                       d = 0.5,
                       standardize = TRUE,
                       control = gwrs_control()) {
  call <- match.call()
  x <- as_design_matrix(x)
  y <- as_response(y, nrow(x))
  control <- validate_control(control)
  neighbors <- resolve_neighbors(x, coords, neighbors, k)
  kernel <- kernel_code(kernel)

  if (missing(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda < 0) {
    stop("`lambda` must be one nonnegative finite number.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be one number in [0, 1].", call. = FALSE)
  }
  if (length(d) != 1L || !is.finite(d) || d < 0 || d > 1) {
    stop("`d` must be one number in [0, 1].", call. = FALSE)
  }
  adaptive <- is.null(bandwidth)
  fixed_bandwidth <- if (adaptive) 1 else as.double(bandwidth)
  if (!adaptive && (length(fixed_bandwidth) != 1L ||
                    !is.finite(fixed_bandwidth) || fixed_bandwidth <= 0)) {
    stop("`bandwidth` must be one positive finite distance.", call. = FALSE)
  }

  prepared <- standardize_design(x, standardize)
  raw <- cpp_gwr_sl_fit(
    prepared$x,
    y,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    adaptive,
    fixed_bandwidth,
    as.double(lambda),
    as.double(alpha),
    as.double(d),
    control$tolerance,
    control$max_iterations,
    control$n_threads,
    control$grain_size
  )

  coefficients <- backtransform_coefficients(
    raw$coefficients,
    prepared$center,
    prepared$scale
  )
  predictor_names <- colnames(x) %||% paste0("x", seq_len(ncol(x)))
  colnames(coefficients) <- c("(Intercept)", predictor_names)
  gwr_center <- sweep(t(raw$gwr_center), 2L, prepared$scale, FUN = "/")
  colnames(gwr_center) <- predictor_names

  diagnostics <- data.frame(
    location = seq_len(nrow(x)),
    converged = as.logical(raw$converged),
    iterations = as.integer(raw$iterations),
    solver = status_labels(raw$status),
    weight_sum = as.double(raw$weight_sum),
    bandwidth = as.double(raw$bandwidth),
    objective = as.double(raw$objective),
    max_kkt_violation = as.double(raw$kkt)
  )

  bandwidth_specification <- list(
    adaptive = adaptive,
    fixed = fixed_bandwidth
  )
  model_diagnostics <- compute_gwr_model_diagnostics(
    prepared_x = prepared$x,
    y = y,
    fitted = as.double(raw$fitted),
    residuals = as.double(raw$residuals),
    raw_coefficients = raw$coefficients,
    coefficients = coefficients,
    neighbors = neighbors,
    kernel = kernel,
    bandwidth = bandwidth_specification,
    lambda = as.double(lambda),
    alpha = as.double(alpha),
    d = as.double(d),
    x_center = prepared$center,
    x_scale = prepared$scale,
    control = control
  )

  object <- list(
    call = call,
    method = if (lambda == 0) "GWR" else "GWR-SL",
    coefficients = coefficients,
    fitted.values = as.double(raw$fitted),
    residuals = as.double(raw$residuals),
    lambda = as.double(lambda),
    alpha = as.double(alpha),
    d = as.double(d),
    kernel = kernel$name,
    adaptive = adaptive,
    bandwidth = if (adaptive) NULL else fixed_bandwidth,
    neighbors = neighbors,
    standardize = isTRUE(standardize),
    x_center = prepared$center,
    x_scale = prepared$scale,
    diagnostics = diagnostics,
    solver_diagnostics = diagnostics,
    model_diagnostics = model_diagnostics,
    hat_diag = model_diagnostics$hat_diag %||% NULL,
    s2_diag = model_diagnostics$s2_diag %||% NULL,
    local_r2 = model_diagnostics$local$local_r_squared %||% NULL,
    gwr_center = gwr_center,
    terms = NULL,
    predictor_names = predictor_names,
    x = if (control$keep_data) x else NULL,
    y = if (control$keep_data) y else NULL,
    control = control
  )
  class(object) <- c("gwr_sl_fit", "gwrs_fit")
  object
}

#' Fit an unpenalized geographically weighted regression
#'
#' @inheritParams gwr_sl_fit
#' @return An object of class `gwrs_fit`.
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_fit(
#'   x, quake$mag, coords = coords, k = 20,
#'   control = gwrs_control(n_threads = 1)
#' )
#' summary(fit)
#'
#' @export
gwr_fit <- function(x,
                    y,
                    coords = NULL,
                    neighbors = NULL,
                    k = NULL,
                    kernel = c("bisquare", "gaussian", "exponential",
                               "tricube", "boxcar"),
                    bandwidth = NULL,
                    standardize = TRUE,
                    control = gwrs_control()) {
  gwr_sl_fit(
    x = x,
    y = y,
    coords = coords,
    neighbors = neighbors,
    k = k,
    kernel = kernel,
    bandwidth = bandwidth,
    lambda = 0,
    alpha = 1,
    d = 0,
    standardize = standardize,
    control = control
  )
}

#' Formula interface for geographically weighted scaled lasso
#'
#' @param formula Model formula.
#' @param data Data frame containing model variables.
#' @param coords Coordinate matrix or two column names in `data`.
#' @param ... Additional arguments passed to [gwr_sl_fit()].
#'
#' @return An object of class `gwrs_fit`.
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' fit <- gwr_sl(
#'   mag ~ depth + stations,
#'   data = quake,
#'   coords = c("long", "lat"),
#'   k = 20,
#'   lambda = 0.05,
#'   alpha = 0.7,
#'   d = 0.5,
#'   control = gwrs_control(n_threads = 1)
#' )
#' fit
#'
#' @export
gwr_sl <- function(formula, data, coords, ...) {
  call <- match.call()
  components <- formula_model_components(formula, data, coords)
  fit <- gwr_sl_fit(
    components$x, components$y, coords = components$coords, ...
  )
  decorate_formula_fit(fit, components, call)
}
