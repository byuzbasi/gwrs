#' Prepare reusable training-fold inputs for SCAD/MCP cross-validation
#'
#' Computes full-data and per-training-fold scaling and neighbor structures.
#' Pass the returned object as `cv_cache` to [cv_gwr_nonconvex()],
#' [cv_gwr_scad()] or [cv_gwr_mcp()] with the same inputs. No response,
#' penalty, lambda sequence, coefficients, or warm-start state is cached.
#' Each training fold retains its own scaling; held-out responses are never
#' used by this preparation. Automatic lambda anchors remain method-specific
#' and retain their existing full-data definition.
#'
#' The returned environment has locked bindings and is read-only. A cache
#' for different `x`, coordinates, `k`, fold labels or standardization is
#' rejected, including changes to row order. It is never stored globally.
#' Holding all prepared folds increases memory use; pass `NULL` to retain
#' the uncached, one-fold-at-a-time route. Kernel weights are not cached,
#' so the same preparation can be used with different kernels/bandwidths.
#'
#' @inheritParams cv_gwr_nonconvex
#' @return A read-only `gwrs_nonconvex_cv_cache` environment.
#' @export
gwr_nonconvex_cv_cache <- function(x, coords, k, fold, standardize = TRUE) {
  x <- as_design_matrix(x)
  coords <- as_coordinate_matrix(coords)
  if (nrow(coords) != nrow(x)) stop("`coords` and `x` must have the same number of rows.", call. = FALSE)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L) stop("`k` must be a positive integer.", call. = FALSE)
  if (missing(fold) || length(fold) != nrow(x) || anyNA(fold)) {
    stop("Supply one non-missing spatial `fold` label per observation.", call. = FALSE)
  }
  factor_fold <- factor(fold)
  if (nlevels(factor_fold) < 2L) stop("At least two spatial folds are required.", call. = FALSE)
  fold_id <- as.integer(factor_fold)
  if (any(vapply(seq_len(nlevels(factor_fold)), function(i) sum(fold_id != i), integer(1L)) < k)) {
    stop("`k` cannot exceed the training size of any fold.", call. = FALSE)
  }
  cache <- new.env(parent = emptyenv())
  cache$schema <- "gwrs-nonconvex-cv-cache-v1"
  cache$inputs <- list(x = x, coords = coords, k = k, fold = fold, standardize = isTRUE(standardize))
  cache$full_neighbors <- gwr_neighbors(coords, k = k)
  cache$full_prepared <- standardize_design(x, standardize)
  cache$folds <- lapply(seq_len(nlevels(factor_fold)), function(i) {
    test <- fold_id == i
    train <- !test
    prepared <- standardize_design(x[train, , drop = FALSE], standardize)
    list(prepared = prepared,
      x_test = apply_standardization(x[test, , drop = FALSE], prepared$center, prepared$scale),
      neighbors = gwr_neighbors(coords[train, , drop = FALSE],
        query_coords = coords[test, , drop = FALSE], k = k, include_self = FALSE))
  })
  class(cache) <- "gwrs_nonconvex_cv_cache"
  lockEnvironment(cache, bindings = TRUE)
  cache
}

validate_nonconvex_cv_cache <- function(cache, x, coords, k, fold, standardize) {
  fields <- c("schema", "inputs", "full_neighbors", "full_prepared", "folds")
  if (!inherits(cache, "gwrs_nonconvex_cv_cache") || !is.environment(cache) ||
      !environmentIsLocked(cache) || !setequal(ls(cache), fields) ||
      !all(vapply(fields, bindingIsLocked, logical(1L), env = cache)) ||
      any(vapply(fields, bindingIsActive, logical(1L), env = cache)) ||
      !identical(cache$schema, "gwrs-nonconvex-cv-cache-v1")) {
    stop("`cv_cache` must be an unmodified object from `gwr_nonconvex_cv_cache()`.", call. = FALSE)
  }
  inputs <- list(x = x, coords = coords, k = k, fold = fold, standardize = isTRUE(standardize))
  if (!identical(cache$inputs, inputs)) {
    stop("`cv_cache` does not match x, coords, k, fold or standardize; rebuild it for these inputs.", call. = FALSE)
  }
  invisible(cache)
}
