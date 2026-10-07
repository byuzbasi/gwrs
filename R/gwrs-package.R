#' Fast Geographically Weighted Regression and Selection
#'
#' The package provides geographically weighted ridge, lasso and
#' elastic-net regression, together with GWR-SCAD and GWR-MCP for local
#' variable selection. Matrix and formula interfaces support fitting,
#' regularization paths, spatial cross-validation, prediction and coefficient
#' summaries. RcppArmadillo and RcppParallel provide the compiled kernels.
#' Use vignette("gwrs-introduction", package = "gwrs") for the five-method guide.
#'
#' @author Bahadir Yuzbasi
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_lasso_fit(
#'   x, quake$mag, coords = coords, k = 20, lambda = 0.05,
#'   control = gwrs_control(n_threads = 1)
#' )
#' head(coef(fit))
#'
#' @keywords internal
#' @useDynLib gwrs, .registration = TRUE
#' @importFrom Matrix Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
"_PACKAGE"

utils::globalVariables(".data")
