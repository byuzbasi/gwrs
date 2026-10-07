#' Fit geographically weighted lasso, ridge, or elastic net
#'
#' The three wrappers use the same normalized local least-squares loss,
#' neighbor structure and preprocessing convention. Ridge applies a quadratic
#' penalty, lasso an absolute-value penalty, and elastic net a mixture controlled
#' by alpha. Positive ridge uses a direct Cholesky solve; an unsuccessful solve
#' is reported as a failed fit. Lasso and elastic net use coordinate updates.
#' Coefficients are returned on the original predictor scale.
#'
#' @inheritParams gwr_sl_fit
#'
#' @return An object of class `gwr_sl_fit` and `gwrs_fit`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' neighbors <- gwr_neighbors(coords, k = 18)
#' control <- gwrs_control(n_threads = 1)
#'
#' lasso <- gwr_lasso_fit(
#'   x, quake$mag, neighbors = neighbors, lambda = 0.05, control = control
#' )
#' ridge <- gwr_ridge_fit(
#'   x, quake$mag, neighbors = neighbors, lambda = 0.05, control = control
#' )
#' elastic_net <- gwr_en_fit(
#'   x, quake$mag, neighbors = neighbors, lambda = 0.05,
#'   alpha = 0.6, control = control
#' )
#' c(lasso = lasso$alpha, ridge = ridge$alpha, elastic_net = elastic_net$alpha)
#'
#' @name gwr_reference_fits
NULL

#' @rdname gwr_reference_fits
#' @export
gwr_lasso_fit <- function(x,
                          y,
                          coords = NULL,
                          neighbors = NULL,
                          k = NULL,
                          kernel = c("bisquare", "gaussian", "exponential",
                                     "tricube", "boxcar"),
                          bandwidth = NULL,
                          lambda,
                          standardize = TRUE,
                          control = gwrs_control()) {
  call <- match.call()
  fit <- gwr_sl_fit(
    x, y, coords = coords, neighbors = neighbors, k = k, kernel = kernel,
    bandwidth = bandwidth, lambda = lambda, alpha = 1, d = 0,
    standardize = standardize, control = control
  )
  fit$call <- call
  fit$method <- "GWR-Lasso"
  fit
}

#' @rdname gwr_reference_fits
#' @export
gwr_ridge_fit <- function(x,
                          y,
                          coords = NULL,
                          neighbors = NULL,
                          k = NULL,
                          kernel = c("bisquare", "gaussian", "exponential",
                                     "tricube", "boxcar"),
                          bandwidth = NULL,
                          lambda,
                          standardize = TRUE,
                          control = gwrs_control()) {
  call <- match.call()
  fit <- gwr_sl_fit(
    x, y, coords = coords, neighbors = neighbors, k = k, kernel = kernel,
    bandwidth = bandwidth, lambda = lambda, alpha = 0, d = 0,
    standardize = standardize, control = control
  )
  fit$call <- call
  fit$method <- "GWR-Ridge"
  fit
}

#' @rdname gwr_reference_fits
#' @export
gwr_en_fit <- function(x,
                       y,
                       coords = NULL,
                       neighbors = NULL,
                       k = NULL,
                       kernel = c("bisquare", "gaussian", "exponential",
                                  "tricube", "boxcar"),
                       bandwidth = NULL,
                       lambda,
                       alpha = 0.5,
                       standardize = TRUE,
                       control = gwrs_control()) {
  call <- match.call()
  fit <- gwr_sl_fit(
    x, y, coords = coords, neighbors = neighbors, k = k, kernel = kernel,
    bandwidth = bandwidth, lambda = lambda, alpha = alpha, d = 0,
    standardize = standardize, control = control
  )
  fit$call <- call
  fit$method <- "GWR-Elastic Net"
  fit
}

#' Unified penalized geographically weighted regression
#'
#' `gwr_penalized_fit()` is the matrix interface and `gwr_penalized()` is the
#' formula interface for the least-squares models that share the fast
#' RcppArmadillo numerical cores. Specify the intended penalty explicitly.
#' Lasso, ridge and elastic net use the convex core. SCAD and MCP use the
#' dedicated nonconvex pathwise coordinate-descent core.
#'
#' @inheritParams gwr_sl_fit
#' @param penalty Penalized GWR specification.
#' @param gamma Concavity parameter used only for SCAD or MCP. A `NULL`
#'   default selects `3.7` for SCAD and `3` for MCP.
#' @param ... Additional arguments passed from the formula interface to
#'   `gwr_penalized_fit()`.
#'
#' @return An object inheriting from `gwrs_fit`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' control <- gwrs_control(n_threads = 1)
#'
#' matrix_fit <- gwr_penalized_fit(
#'   x, quake$mag, coords = coords, k = 18, lambda = 0.05,
#'   penalty = "elastic_net", alpha = 0.6, control = control
#' )
#' formula_fit <- gwr_penalized(
#'   mag ~ depth + stations, quake, coords = c("long", "lat"),
#'   k = 18, lambda = 0.05, penalty = "lasso", control = control
#' )
#' c(matrix_fit$method, formula_fit$method)
#'
#' @export
gwr_penalized_fit <- function(x,
                              y,
                              coords = NULL,
                              neighbors = NULL,
                              k = NULL,
                              kernel = c("bisquare", "gaussian", "exponential",
                                         "tricube", "boxcar"),
                              bandwidth = NULL,
                              lambda,
                              penalty = c("scaled_lasso", "lasso", "ridge",
                                          "elastic_net", "scad", "mcp"),
                              alpha = 0.5,
                              d = 0.5,
                              gamma = NULL,
                              standardize = TRUE,
                              control = gwrs_control()) {
  call <- match.call()
  penalty <- match.arg(penalty)
  if (penalty %in% c("scad", "mcp")) {
    fit <- gwr_nonconvex_fit(
      x = x,
      y = y,
      coords = coords,
      neighbors = neighbors,
      k = k,
      kernel = kernel,
      bandwidth = bandwidth,
      lambda = lambda,
      penalty = penalty,
      gamma = gamma,
      standardize = standardize,
      control = control
    )
    fit$call <- call
    return(fit)
  }
  parameters <- penalized_gwr_parameters(penalty, alpha, d)
  fit <- gwr_sl_fit(
    x = x,
    y = y,
    coords = coords,
    neighbors = neighbors,
    k = k,
    kernel = kernel,
    bandwidth = bandwidth,
    lambda = lambda,
    alpha = parameters$alpha,
    d = parameters$d,
    standardize = standardize,
    control = control
  )
  fit$call <- call
  fit$penalty <- penalty
  fit$method <- parameters$method
  class(fit) <- unique(c(paste0("gwr_", penalty, "_fit"), class(fit)))
  fit
}

#' @rdname gwr_penalized_fit
#' @param formula Model formula with an intercept.
#' @param data Data frame containing model variables.
#' @export
gwr_penalized <- function(formula,
                          data,
                          coords,
                          ...,
                          penalty = c("scaled_lasso", "lasso", "ridge",
                                      "elastic_net", "scad", "mcp")) {
  call <- match.call()
  penalty <- match.arg(penalty)
  components <- formula_model_components(formula, data, coords)
  fit <- gwr_penalized_fit(
    components$x,
    components$y,
    coords = components$coords,
    penalty = penalty,
    ...
  )
  decorate_formula_fit(fit, components, call)
}

#' Formula interfaces for reference penalized GWR models
#'
#' These are concise formula wrappers around [gwr_penalized()].
#'
#' @inheritParams gwr_penalized
#' @param ... Additional arguments passed to [gwr_penalized()].
#' @return An object inheriting from `gwrs_fit`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' control <- gwrs_control(n_threads = 1)
#' lasso <- gwr_lasso(
#'   mag ~ depth + stations, quake, c("long", "lat"),
#'   k = 18, lambda = 0.05, control = control
#' )
#' ridge <- gwr_ridge(
#'   mag ~ depth + stations, quake, c("long", "lat"),
#'   k = 18, lambda = 0.05, control = control
#' )
#' elastic_net <- gwr_en(
#'   mag ~ depth + stations, quake, c("long", "lat"),
#'   k = 18, lambda = 0.05, alpha = 0.6, control = control
#' )
#' c(lasso$method, ridge$method, elastic_net$method)
#'
#' @name gwr_reference_formulas
NULL

#' @rdname gwr_reference_formulas
#' @export
gwr_lasso <- function(formula, data, coords, ...) {
  gwr_penalized(formula, data, coords, ..., penalty = "lasso")
}

#' @rdname gwr_reference_formulas
#' @export
gwr_ridge <- function(formula, data, coords, ...) {
  gwr_penalized(formula, data, coords, ..., penalty = "ridge")
}

#' @rdname gwr_reference_formulas
#' @param alpha Elastic-net mixing parameter in `[0, 1]`.
#' @export
gwr_en <- function(formula, data, coords, ..., alpha = 0.5) {
  gwr_penalized(
    formula, data, coords, ..., penalty = "elastic_net", alpha = alpha
  )
}

penalized_gwr_parameters <- function(penalty, alpha, d) {
  switch(
    penalty,
    lasso = list(alpha = 1, d = 0, method = "GWR-Lasso"),
    ridge = list(alpha = 0, d = 0, method = "GWR-Ridge"),
    elastic_net = {
      validate_penalty_parameters(alpha)
      list(alpha = as.double(alpha), d = 0, method = "GWR-Elastic Net")
    },
    scaled_lasso = {
      validate_penalty_parameters(alpha, d)
      if (length(d) != 1L) {
        stop("`d` must be one number in [0, 1].", call. = FALSE)
      }
      list(
        alpha = as.double(alpha), d = as.double(d), method = "GWR-SL"
      )
    }
  )
}
