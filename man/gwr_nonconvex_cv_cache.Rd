\name{gwr_nonconvex_cv_cache}
\alias{gwr_nonconvex_cv_cache}
\title{Prepare reusable training-fold inputs for SCAD/MCP cross-validation}
\usage{gwr_nonconvex_cv_cache(x, coords, k, fold, standardize = TRUE)}
\arguments{
\item{x}{Numeric predictor matrix with observations in rows.}
\item{coords}{Numeric coordinate matrix or an sf point object.}
\item{k}{Retained neighbor count, no larger than any training-fold sample.}
\item{fold}{One non-missing spatial-fold label per observation.}
\item{standardize}{Whether to center and scale predictors within each training fold.}
}
\value{A read-only \code{gwrs_nonconvex_cv_cache} environment with locked bindings.}
\description{Prepares full-data and per-training-fold scaling and neighbor structures for reuse by SCAD and MCP. Pass it as \code{cv_cache} to \code{\link{cv_gwr_nonconvex}}, \code{\link{cv_gwr_scad}} or \code{\link{cv_gwr_mcp}}.}
\details{Predictors, coordinates, neighbor count, fold labels and standardization must match exactly, including row order. Each fold uses its own training-only transformation. Responses, penalties, coefficients, warm-start state and lambda paths are never cached. Automatic lambda anchors retain their existing method-specific full-data definition. Kernel weights are not cached, allowing preparation to be reused across kernels and bandwidths. Holding all prepared folds costs additional memory; \code{cv_cache = NULL} retains the one-fold-at-a-time route. No global cache is created.}
