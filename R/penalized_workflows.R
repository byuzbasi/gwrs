#' Unified regularization paths for penalized GWR
#'
#' GWR-Lasso, GWR-Ridge, GWR-Elastic Net, GWR-SCAD and GWR-MCP use the
#' same normalized local loss convention. Convex and nonconvex models
#' dispatch to their corresponding RcppArmadillo path cores.
#' Always specify the intended penalty when using the unified interface.
#'
#' @inheritParams gwr_sl_path
#' @param penalty Penalized GWR specification.
#' @param gamma Concavity parameter used only for SCAD or MCP.
#' @param ... Additional arguments passed to [gwr_sl_path()].
#'
#' @return A `gwrs_path` object.
#'
#' @examples
#' quake <- datasets::quakes[1:28, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' neighbors <- gwr_neighbors(coords, k = 17)
#' lambda <- c(0.2, 0.08, 0.02)
#' control <- gwrs_control(n_threads = 1, diagnostics = "none")
#'
#' lasso <- gwr_lasso_path(
#'   x, quake$mag, neighbors = neighbors, lambda = lambda,
#'   control = control
#' )
#' ridge <- gwr_ridge_path(
#'   x, quake$mag, neighbors = neighbors, lambda = lambda,
#'   control = control
#' )
#' elastic_net <- gwr_en_path(
#'   x, quake$mag, neighbors = neighbors, lambda = lambda,
#'   alpha = 0.6, control = control
#' )
#' vapply(
#'   list(lasso = lasso, ridge = ridge, elastic_net = elastic_net),
#'   function(path) length(path$lambda), integer(1)
#' )
#'
#' @export
gwr_penalized_path <- function(x,
                               y,
                               ...,
                               penalty = c("scaled_lasso", "lasso", "ridge",
                                           "elastic_net", "scad", "mcp"),
                               alpha = 0.5,
                               d = 0.5,
                               gamma = NULL) {
  call <- match.call()
  penalty <- match.arg(penalty)
  if (penalty %in% c("scad", "mcp")) {
    result <- gwr_nonconvex_path(
      x = x,
      y = y,
      ...,
      penalty = penalty,
      gamma = gamma
    )
    result$call <- call
    return(result)
  }
  parameters <- penalized_gwr_parameters(penalty, alpha, d)
  result <- gwr_sl_path(
    x = x,
    y = y,
    ...,
    alpha = parameters$alpha,
    d = parameters$d
  )
  result$call <- call
  result$penalty <- penalty
  result$method <- paste(parameters$method, "path")
  result
}

#' @rdname gwr_penalized_path
#' @export
gwr_lasso_path <- function(x, y, ...) {
  gwr_penalized_path(x, y, ..., penalty = "lasso")
}

#' @rdname gwr_penalized_path
#' @export
gwr_ridge_path <- function(x, y, ...) {
  gwr_penalized_path(x, y, ..., penalty = "ridge")
}

#' @rdname gwr_penalized_path
#' @export
gwr_en_path <- function(x, y, ..., alpha = 0.5) {
  gwr_penalized_path(
    x, y, ..., penalty = "elastic_net", alpha = alpha
  )
}

#' Unified spatial cross-validation for penalized GWR
#'
#' @inheritParams cv_gwr_sl
#' @param penalty Penalized GWR specification.
#' @param gamma Concavity parameter used only for SCAD or MCP.
#' @param ... Additional arguments passed to [cv_gwr_sl()].
#'
#' @return A `gwrs_cv` object.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fold <- spatial_folds(coords, n_folds = 3, seed = 42)
#' cv <- cv_gwr_penalized(
#'   x, quake$mag, coords, k = 8, fold = fold,
#'   lambda = c(0.12, 0.04), penalty = "lasso", refit = FALSE,
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#' cv$best
#'
#' @export
cv_gwr_penalized <- function(x,
                             y,
                             coords,
                             k,
                             fold,
                             ...,
                             penalty = c("scaled_lasso", "lasso", "ridge",
                                         "elastic_net", "scad", "mcp"),
                             alpha = 0.5,
                             d = c(0, 0.5, 1),
                             gamma = NULL) {
  call <- match.call()
  penalty <- match.arg(penalty)
  if (penalty %in% c("scad", "mcp")) {
    result <- cv_gwr_nonconvex(
      x = x,
      y = y,
      coords = coords,
      k = k,
      fold = fold,
      ...,
      penalty = penalty,
      gamma = gamma
    )
    result$call <- call
    return(result)
  }
  parameters <- penalized_gwr_parameters(
    penalty,
    alpha,
    if (penalty == "scaled_lasso") d[[1L]] else 0
  )
  selected_d <- if (penalty == "scaled_lasso") d else 0
  result <- cv_gwr_sl(
    x = x,
    y = y,
    coords = coords,
    k = k,
    fold = fold,
    ...,
    alpha = parameters$alpha,
    d = selected_d
  )
  result$call <- call
  result$penalty <- penalty
  if (!is.null(result$fit)) {
    result$fit$penalty <- penalty
    result$fit$method <- parameters$method
  }
  result
}
