#' Fit robust LAD geographically weighted scaled lasso
#'
#' At location `i`, the numerical core minimizes
#' \deqn{\sum_j \bar w_{ij}|y_j-\beta_{i0}-x_j^T\beta_i|
#' +\lambda\alpha\|\beta_i\|_1
#' +\frac{\lambda(1-\alpha)}{2}
#' \|\beta_i-d\widehat\beta_i^{GWR}\|_2^2.}
#' The intercept is unpenalized. A two-split ADMM solver uses a factorization
#' cached once per target location; target locations are processed in parallel.
#'
#' @inheritParams gwr_sl_fit
#' @param rho_residual Positive ADMM parameter for the residual split.
#' @param rho_penalty Positive ADMM parameter for the lasso split.
#'
#' @return An object of class `lad_gwr_sl_fit` and `gwrs_fit`.
#'
#' @examples
#' n <- 24
#' coords <- cbind(east = seq_len(n), north = rep(0, n))
#' x <- cbind(trend = seq(-1, 1, length.out = n), wave = cos(seq_len(n)))
#' y <- 0.5 + 1.2 * x[, "trend"] - 0.3 * x[, "wave"]
#' y[12] <- y[12] + 2
#' fit <- lad_gwr_sl_fit(
#'   x, y, coords = coords, k = 18,
#'   lambda = 0.05, alpha = 0.7, d = 0.4,
#'   control = gwrs_control(
#'     tolerance = 1e-5, max_iterations = 1500, n_threads = 1
#'   )
#' )
#' fit
#' head(coef(fit))
#'
#' @export
lad_gwr_sl_fit <- function(x,
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
                           rho_residual = 0.1,
                           rho_penalty = 0.1,
                           control = gwrs_control(max_iterations = 2000L)) {
  call <- match.call()
  x <- as_design_matrix(x)
  y <- as_response(y, nrow(x))
  control <- validate_control(control)
  neighbors <- resolve_neighbors(x, coords, neighbors, k)
  kernel <- kernel_code(kernel)
  bandwidth <- resolve_bandwidth(bandwidth)
  validate_penalty_parameters(alpha, d)
  if (length(d) != 1L) {
    stop("`d` must be one number in [0, 1].", call. = FALSE)
  }
  if (missing(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda < 0) {
    stop("`lambda` must be one nonnegative finite number.", call. = FALSE)
  }
  if (length(rho_residual) != 1L || !is.finite(rho_residual) ||
      rho_residual <= 0 || length(rho_penalty) != 1L ||
      !is.finite(rho_penalty) || rho_penalty <= 0) {
    stop("ADMM rho parameters must be positive finite numbers.",
         call. = FALSE)
  }

  prepared <- standardize_design(x, standardize)
  raw <- cpp_lad_gwr_sl_fit(
    prepared$x,
    y,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    as.double(lambda),
    as.double(alpha),
    as.double(d),
    as.double(rho_residual),
    as.double(rho_penalty),
    control$tolerance,
    control$max_iterations,
    control$n_threads,
    control$grain_size
  )

  coefficients <- backtransform_coefficients(
    raw$coefficients, prepared$center, prepared$scale
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
    primal_residual = as.double(raw$primal_residual),
    dual_residual = as.double(raw$dual_residual)
  )
  model_diagnostics <- compute_descriptive_model_diagnostics(
    y = y,
    fitted = as.double(raw$fitted),
    residuals = as.double(raw$residuals),
    neighbors = neighbors,
    kernel = kernel,
    bandwidth = bandwidth,
    control = control,
    method = "LAD-GWR-SL",
    extra = list(
      objective_sum = sum(raw$objective),
      absolute_residual_sum = sum(abs(raw$residuals))
    )
  )

  object <- list(
    call = call,
    method = "LAD-GWR-SL",
    loss = "weighted absolute deviation",
    coefficients = coefficients,
    fitted.values = as.double(raw$fitted),
    residuals = as.double(raw$residuals),
    lambda = as.double(lambda),
    alpha = as.double(alpha),
    d = as.double(d),
    rho_residual = as.double(rho_residual),
    rho_penalty = as.double(rho_penalty),
    kernel = kernel$name,
    adaptive = bandwidth$adaptive,
    bandwidth = if (bandwidth$adaptive) NULL else bandwidth$fixed,
    neighbors = neighbors,
    standardize = isTRUE(standardize),
    x_center = prepared$center,
    x_scale = prepared$scale,
    diagnostics = diagnostics,
    solver_diagnostics = diagnostics,
    model_diagnostics = model_diagnostics,
    local_r2 = model_diagnostics$local$local_r_squared %||% NULL,
    gwr_center = gwr_center,
    terms = NULL,
    predictor_names = predictor_names,
    x = if (control$keep_data) x else NULL,
    y = if (control$keep_data) y else NULL,
    control = control
  )
  class(object) <- c("lad_gwr_sl_fit", "gwrs_fit")
  object
}

#' Formula interface for robust LAD-GWR-SL
#'
#' @param formula Model formula.
#' @param data Data frame containing model variables.
#' @param coords Coordinate matrix or two column names in `data`.
#' @param ... Additional arguments passed to [lad_gwr_sl_fit()].
#'
#' @return An object of class `lad_gwr_sl_fit` and `gwrs_fit`.
#'
#' @examples
#' n <- 24
#' dat <- data.frame(
#'   east = seq_len(n),
#'   north = rep(0, n),
#'   trend = seq(-1, 1, length.out = n),
#'   wave = cos(seq_len(n))
#' )
#' dat$y <- 0.5 + 1.2 * dat$trend - 0.3 * dat$wave
#' dat$y[12] <- dat$y[12] + 2
#' fit <- lad_gwr_sl(
#'   y ~ trend + wave, data = dat, coords = c("east", "north"),
#'   k = 18, lambda = 0.05, alpha = 0.7, d = 0.4,
#'   control = gwrs_control(
#'     tolerance = 1e-5, max_iterations = 1500, n_threads = 1
#'   )
#' )
#' fit
#'
#' @export
lad_gwr_sl <- function(formula, data, coords, ...) {
  call <- match.call()
  components <- formula_model_components(formula, data, coords)
  fit <- lad_gwr_sl_fit(
    components$x, components$y, coords = components$coords, ...
  )
  decorate_formula_fit(fit, components, call)
}
