#' Sparse residual Moran diagnostics
#'
#' Computes global Moran's I and local Moran (LISA) statistics directly from a
#' fitted model's nearest-neighbor graph. Self-neighbor weights are removed and
#' no dense spatial matrix is formed. Permutations use deterministic independent
#' C++ random streams and can be evaluated in parallel.
#'
#' @param object A fitted `gwrs_fit`.
#' @param permutations Number of random permutations; zero skips inference.
#' @param row_standardize Whether non-self spatial weights are row-standardized.
#' @param seed Nonnegative deterministic C++ permutation seed.
#' @param alpha Significance level used for LISA cluster labels.
#' @param weights Optional nonnegative numeric or sparse spatial weight matrix,
#'   in observation order, with a zero diagonal. When supplied it replaces the
#'   fitted kernel graph. Directed weights are allowed; no symmetrization occurs.
#' @param zero_policy For explicit weights, `"error"` rejects islands; `"keep"`
#'   retains them in global centering and permutations and reports local NA.
#' @param memory_limit_mb Workspace guard for the explicit-weight calculation.
#' @param ... Unused.
#'
#' @return `gwr_moran()` returns a `gwrs_moran` object. `gwr_local_moran()`
#'   returns a `gwrs_local_moran` data frame. `gwr_moran_diagnostics()` returns
#'   both as `global` and `local`, reusing one C++ permutation calculation.
#' @details The existing permutation null randomly relabels the entire centered
#'   residual vector (total randomization); it is not conditional LISA with the
#'   focal residual fixed. Two-sided global exceedances are centered at
#'   `-1/(n-1)`; local exceedances are centered at each location's simulated mean.
#'   Both use the plus-one permutation correction. LISA cluster labels use raw,
#'   unadjusted p-values. Regression residuals need not be exchangeable, so these
#'   permutation diagnostics are exploratory, not exact residual-model tests.
#'   No fitted model or tuning parameters are re-estimated during permutations.
#' @references Anselin, L. (1995). Local Indicators of Spatial Association--LISA.
#'   Geographical Analysis, 27, 93-115. \doi{10.1111/j.1538-4632.1995.tb00338.x}.
#'
#' @examples
#' quake <- datasets::quakes[1:36, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_fit(
#'   x, quake$mag, coords = coords, k = 24,
#'   control = gwrs_control(n_threads = 1)
#' )
#' global <- gwr_moran(fit, permutations = 19, seed = 42)
#' local <- gwr_local_moran(fit, permutations = 19, seed = 42)
#' global
#' head(local)
#' plot(local, significant_only = TRUE)
#'
#' @export
gwr_moran <- function(object,
                      permutations = 999L,
                      row_standardize = TRUE,
                      seed = 1L,
                      ..., weights = NULL,
                      zero_policy = c("error", "keep"),
                      memory_limit_mb = 512) {
  raw <- compute_gwr_moran(
    object, permutations = permutations,
    row_standardize = row_standardize, seed = seed, weights = weights,
    zero_policy = zero_policy, memory_limit_mb = memory_limit_mb
  )
  assemble_gwr_moran(raw, object, permutations, row_standardize, seed)
}

assemble_gwr_moran <- function(raw, object, permutations, row_standardize, seed) {
  result <- list(
    I = raw$I,
    expected_I = raw$expected_I,
    p_value = raw$p_value,
    permutations = as.integer(permutations),
    n = length(object$residuals),
    kernel = object$kernel,
    adaptive = object$adaptive,
    bandwidth = object$bandwidth,
    row_standardize = isTRUE(row_standardize),
    seed = as.integer(seed),
    permuted_I = raw$permuted_I
  )
  if (!is.null(raw$weight_specification)) {
    result$weights <- raw$weight_specification
    result$kernel <- "user_supplied"
    result$adaptive <- NA
    result$bandwidth <- NULL
  }
  class(result) <- "gwrs_moran"
  result
}

#' @rdname gwr_moran
#' @export
gwr_local_moran <- function(object,
                            permutations = 999L,
                            row_standardize = TRUE,
                            seed = 1L,
                            alpha = 0.05,
                            ..., weights = NULL,
                            zero_policy = c("error", "keep"),
                            memory_limit_mb = 512) {
  validate_lisa_alpha(alpha)
  raw <- compute_gwr_moran(
    object, permutations = permutations,
    row_standardize = row_standardize, seed = seed, weights = weights,
    zero_policy = zero_policy, memory_limit_mb = memory_limit_mb
  )
  assemble_gwr_local_moran(raw, object, permutations, row_standardize, seed, alpha)
}

assemble_gwr_local_moran <- function(raw, object, permutations,
                                      row_standardize, seed, alpha) {
  significant <- is.finite(raw$local_p) & raw$local_p <= alpha
  cluster <- rep("Not Significant", length(raw$local_I))
  cluster[significant & raw$centered_residual > 0 & raw$spatial_lag > 0] <-
    "High-High"
  cluster[significant & raw$centered_residual < 0 & raw$spatial_lag < 0] <-
    "Low-Low"
  cluster[significant & raw$centered_residual > 0 & raw$spatial_lag < 0] <-
    "High-Low"
  cluster[significant & raw$centered_residual < 0 & raw$spatial_lag > 0] <-
    "Low-High"
  result <- data.frame(
    location = seq_along(raw$local_I),
    local_I = as.double(raw$local_I),
    z_score = as.double(raw$local_z),
    p_value = as.double(raw$local_p),
    cluster = factor(
      cluster,
      levels = c("Not Significant", "High-High", "Low-Low", "High-Low",
                 "Low-High")
    ),
    residual = as.double(object$residuals),
    residual_centered = as.double(raw$centered_residual),
    spatial_lag = as.double(raw$spatial_lag)
  )
  attr(result, "settings") <- list(
    permutations = permutations,
    row_standardize = row_standardize,
    seed = seed,
    alpha = alpha,
    kernel = object$kernel
  )
  if (!is.null(raw$weight_specification)) {
    attr(result, "settings")$weights <- raw$weight_specification
    attr(result, "settings")$kernel <- "user_supplied"
  }
  class(result) <- c("gwrs_local_moran", "data.frame")
  result
}

#' @rdname gwr_moran
#' @export
gwr_moran_diagnostics <- function(object, permutations = 999L,
                                   row_standardize = TRUE, seed = 1L,
                                   alpha = 0.05, ..., weights = NULL,
                                   zero_policy = c("error", "keep"),
                                   memory_limit_mb = 512) {
  validate_lisa_alpha(alpha)
  raw <- compute_gwr_moran(object, permutations, row_standardize, seed,
    weights, zero_policy, memory_limit_mb)
  list(global = assemble_gwr_moran(raw, object, permutations, row_standardize, seed),
       local = assemble_gwr_local_moran(raw, object, permutations,
                                       row_standardize, seed, alpha))
}

validate_lisa_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one number strictly between zero and one.", call. = FALSE)
  }
}

compute_gwr_moran <- function(object, permutations, row_standardize, seed,
                              weights = NULL, zero_policy = c("error", "keep"),
                              memory_limit_mb = 512) {
  if (!is.null(weights)) {
    if (!inherits(object, "gwrs_fit")) stop("`object` must be a gwrs_fit.", call. = FALSE)
    return(compute_gwrs_weighted_moran(object$residuals, weights, permutations,
      row_standardize, seed, zero_policy, object$control$n_threads,
      object$control$grain_size, memory_limit_mb))
  }
  validate_gwrs_diagnostic_object(object)
  permutations <- as.integer(permutations)
  seed <- as.integer(seed)
  if (length(permutations) != 1L || is.na(permutations) || permutations < 0L ||
      length(seed) != 1L || is.na(seed) || seed < 0L) {
    stop("`permutations` and `seed` must be nonnegative integers.",
         call. = FALSE)
  }
  kernel <- kernel_code(object$kernel)
  bandwidth <- diagnostic_bandwidth(object)
  cpp_gwr_moran(
    as.double(object$residuals),
    object$neighbors$index,
    object$neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    isTRUE(row_standardize),
    permutations,
    as.double(seed),
    object$control$n_threads,
    object$control$grain_size
  )
}

#' @export
print.gwrs_moran <- function(x, ...) {
  cat("gwrs residual Moran's I\n")
  cat("  I:", format(signif(x$I, 7)), "\n")
  cat("  expected I:", format(signif(x$expected_I, 7)), "\n")
  cat("  permutation p-value:", format.pval(x$p_value, digits = 5), "\n")
  cat("  permutations:", x$permutations, "\n")
  invisible(x)
}

#' Plot residual LISA statistics or a cluster map
#' @param x A result from [gwr_local_moran()].
#' @param significant_only Gray out nonsignificant locations at `alpha`.
#' @param alpha Display significance level; defaults to the calculation setting.
#' @param coords Optional finite two-column coordinates in result row order.
#'   When supplied, draw a cluster map. Otherwise preserve the location-ID plot.
#' @param ... Additional base graphics arguments.
#' @return The input object invisibly. Missing permutation p-values are shown
#'   as unassessed on coordinate maps. Labels use raw p-values.
#' @export
plot.gwrs_local_moran <- function(x,
                                  significant_only = FALSE,
                                  alpha = attr(x, "settings")$alpha %||% 0.05,
                                  ...,
                                  coords = NULL) {
  validate_lisa_alpha(alpha)
  palette <- c(
    "Not Significant" = "grey75", "High-High" = "#D7191C",
    "Low-Low" = "#2C00A0", "High-Low" = "#FDA0A0",
    "Low-High" = "#ABD9E9"
  )
  colors <- unname(palette[as.character(x$cluster)])
  if (isTRUE(significant_only)) {
    colors[!is.finite(x$p_value) | x$p_value > alpha] <- "grey80"
  }
  if (is.null(coords)) {
    graphics::plot(
      x$location, x$local_I, pch = 19, col = colors,
      xlab = "Location", ylab = "Local Moran's I",
      main = "Residual LISA", ...
    )
    graphics::abline(h = 0, col = "grey50", lty = 2)
  } else {
    coords <- validate_diagnostic_plot_coords(coords, nrow(x))
    colors[!is.finite(x$p_value)] <- "black"
    graphics::plot(coords, pch = ifelse(is.finite(x$p_value), 19, 4), col = colors,
      xlab = "Coordinate 1", ylab = "Coordinate 2", asp = 1,
      main = "Residual LISA (unadjusted p-values)", ...)
    if (any(!is.finite(x$p_value))) palette <- c(palette, "Unassessed" = "black")
  }
  graphics::legend(
    "topright", legend = names(palette), col = palette,
    pch = ifelse(names(palette) == "Unassessed", 4, 19),
    bty = if (is.null(coords)) "n" else "o", bg = "white", cex = 0.8
  )
  invisible(x)
}

# Shared by the explicit-weight public diagnostic and GWR/OLS selection records.
validate_gwrs_moran_weights <- function(weights, n, zero_policy = "error") {
  zero_policy <- match.arg(zero_policy, c("error", "keep"))
  if (!(is.matrix(weights) && is.numeric(weights)) && !inherits(weights, "Matrix"))
    stop("`weights` must be a numeric or Matrix spatial weight matrix.", call. = FALSE)
  if (!identical(as.integer(dim(weights)), c(as.integer(n), as.integer(n))))
    stop("`weights` must have one row and column per observation.", call. = FALSE)
  weights <- methods::as(Matrix::Matrix(weights, sparse = TRUE), "generalMatrix")
  weights <- methods::as(weights, "dMatrix")
  weights <- methods::as(weights, "CsparseMatrix")
  values <- methods::slot(weights, "x")
  if (any(!is.finite(values)) || any(values < 0) || any(Matrix::diag(weights) != 0))
    stop("`weights` must be finite, nonnegative and have a zero diagonal.", call. = FALSE)
  mass <- Matrix::rowSums(weights)
  if (any(!is.finite(mass)) || sum(mass) <= 0)
    stop("`weights` must contain finite positive total weight.", call. = FALSE)
  if (zero_policy == "error" && any(mass == 0))
    stop("`weights` contain islands; explicitly use zero_policy = 'keep'.", call. = FALSE)
  weights
}

compute_gwrs_weighted_moran <- function(residuals, weights, permutations = 999L,
    row_standardize = TRUE, seed = 1L, zero_policy = c("error", "keep"),
    n_threads = 1L, grain_size = 16L, memory_limit_mb = 512) {
  zero_policy <- match.arg(zero_policy)
  weights <- validate_gwrs_moran_weights(weights, length(residuals), zero_policy)
  for (value in list(permutations, seed)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value != floor(value) || value >= .Machine$integer.max)
      stop("`permutations` and `seed` must be nonnegative integers.", call. = FALSE)
  }
  if (!is.logical(row_standardize) || length(row_standardize) != 1L || is.na(row_standardize))
    stop("`row_standardize` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(memory_limit_mb) || length(memory_limit_mb) != 1L ||
      !is.finite(memory_limit_mb) || memory_limit_mb <= 0)
    stop("`memory_limit_mb` must be positive and finite.", call. = FALSE)
  raw <- cpp_gwr_moran_weights(as.double(residuals), weights, row_standardize,
    zero_policy == "keep", as.integer(permutations), as.double(seed),
    n_threads, grain_size, memory_limit_mb)
  raw$weight_specification <- list(source = "user_supplied", n = length(residuals),
    row_standardize = row_standardize, zero_policy = zero_policy,
    islands = as.integer(raw$islands), directed = !Matrix::isSymmetric(weights),
    inference = "exploratory total residual randomization; no refitting or selection adjustment")
  raw
}
