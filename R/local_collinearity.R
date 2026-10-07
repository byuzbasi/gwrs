#' Fast local multicollinearity diagnostics
#'
#' Computes local condition numbers, variance-inflation factors, maximum
#' absolute local correlations, effective sample sizes, and positive kernel
#' counts directly from the reusable nearest-neighbor structure.
#'
#' @param object A fitted `gwrs_fit` retaining its design matrix.
#' @param threshold_cn Condition-number flag threshold.
#' @param threshold_vif VIF flag threshold.
#' @param center Whether predictors are locally centered.
#' @param scale Whether predictors are locally scaled.
#' @param return_pairwise Whether all pairwise local correlations are returned.
#' @param ... Unused.
#'
#' @return A `gwrs_local_collinearity` object.
#' @details With the default centering, CN is the square root of the ratio of
#'   largest to smallest eigenvalues of the weighted predictor correlation
#'   matrix. The intercept is excluded. For comparable global diagnostics use
#'   [gwr_global_collinearity()]. The rank tolerance is the largest eigenvalue
#'   (or one, if larger) times `1e-10`. Singular designs return infinite CN/VIF.
#'   Constant local predictors have undefined correlations and infinite CN/VIF.
#'   Even with `scale = FALSE`, the covariance is converted to a correlation
#'   matrix before computing CN/VIF. These predictor diagnostics differ from
#'   the local design condition number in [gwr_diagnostics()], which includes
#'   the intercept and uses the model's working scale.
#' @seealso [summary.gwrs_local_collinearity()], [plot.gwrs_local_collinearity()]
#'
#' @examples
#' quake <- datasets::quakes[1:36, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_fit(
#'   x, quake$mag, coords = coords, k = 24,
#'   control = gwrs_control(n_threads = 1)
#' )
#' collinearity <- gwr_local_collinearity(fit, return_pairwise = TRUE)
#' head(collinearity$condition_number)
#' head(collinearity$vif_long)
#'
#' @export
gwr_local_collinearity <- function(object,
                                   threshold_cn = 30,
                                   threshold_vif = 10,
                                   center = TRUE,
                                   scale = TRUE,
                                   return_pairwise = FALSE,
                                   ...) {
  validate_gwrs_diagnostic_object(object, require_x = TRUE)
  if (!is.numeric(threshold_cn) || length(threshold_cn) != 1L ||
      !is.finite(threshold_cn) || threshold_cn <= 0 ||
      !is.numeric(threshold_vif) || length(threshold_vif) != 1L ||
      !is.finite(threshold_vif) || threshold_vif <= 0) {
    stop("Collinearity thresholds must be positive finite scalars.",
         call. = FALSE)
  }
  kernel <- kernel_code(object$kernel)
  bandwidth <- diagnostic_bandwidth(object)
  raw <- cpp_gwr_local_collinearity(
    object$x,
    object$neighbors$index,
    object$neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    isTRUE(center),
    isTRUE(scale),
    isTRUE(return_pairwise),
    object$control$n_threads,
    object$control$grain_size
  )
  predictor_names <- object$predictor_names %||%
    colnames(object$x) %||% paste0("x", seq_len(ncol(object$x)))
  vif <- as.matrix(raw$vif)
  colnames(vif) <- predictor_names
  condition_number <- data.frame(
    location = seq_len(nrow(object$x)),
    condition_number = as.double(raw$condition_number),
    flagged = as.double(raw$condition_number) > threshold_cn
  )
  vif_long <- data.frame(
    location = rep(seq_len(nrow(vif)), times = ncol(vif)),
    predictor = rep(colnames(vif), each = nrow(vif)),
    vif = as.double(vif),
    stringsAsFactors = FALSE
  )
  vif_long$flagged <- vif_long$vif > threshold_vif
  pairwise <- NULL
  if (isTRUE(return_pairwise) && length(predictor_names) > 1L) {
    correlations <- raw$correlations
    pairs <- utils::combn(seq_along(predictor_names), 2L)
    pairwise <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(index) {
      first <- pairs[1L, index]
      second <- pairs[2L, index]
      data.frame(
        location = seq_len(dim(correlations)[3L]),
        predictor_1 = predictor_names[first],
        predictor_2 = predictor_names[second],
        correlation = correlations[first, second, ],
        stringsAsFactors = FALSE
      )
    }))
  }
  result <- list(
    condition_number = condition_number,
    vif = vif,
    vif_long = vif_long,
    max_abs_correlation = as.double(raw$max_abs_correlation),
    effective_sample_size = as.double(raw$effective_sample_size),
    positive_weights = as.integer(raw$positive_weights),
    pairwise_correlations = pairwise,
    thresholds = list(condition_number = threshold_cn, vif = threshold_vif),
    settings = list(center = center, scale = scale, kernel = object$kernel)
  )
  class(result) <- "gwrs_local_collinearity"
  result
}

#' @export
print.gwrs_local_collinearity <- function(x, ...) {
  cat("gwrs local multicollinearity diagnostics\n")
  cat("  condition number median/max:",
      format(stats::median(x$condition_number$condition_number, na.rm = TRUE)),
      "/", format(max(x$condition_number$condition_number, na.rm = TRUE)),
      "\n")
  cat("  locations above threshold:",
      sum(x$condition_number$flagged, na.rm = TRUE), "\n")
  cat("  predictor-location VIF flags:", sum(x$vif_long$flagged, na.rm = TRUE),
      "\n")
  invisible(x)
}

#' Summarize local predictor condition numbers
#' @param object A result from [gwr_local_collinearity()].
#' @param ... Unused.
#' @return One row with location counts, median, 90th percentile, maximum,
#'   and counts above the object's CN threshold. Infinite CN values are retained;
#'   unavailable values are counted separately. Quantiles use empirical order
#'   statistics (`type = 1`) to handle infinite values without interpolation.
#' @export
summary.gwrs_local_collinearity <- function(object, ...) {
  cn <- object$condition_number$condition_number
  valid <- cn[!is.na(cn)]
  q <- if (length(valid)) stats::quantile(valid, c(.5, .9, 1), type = 1) else
    rep(NA_real_, 3L)
  data.frame(locations = length(cn), unavailable = sum(is.na(cn)),
    infinite = sum(is.infinite(cn)), median = unname(q[1]),
    p90 = unname(q[2]), maximum = unname(q[3]),
    above_threshold = sum(cn > object$thresholds$condition_number, na.rm = TRUE))
}

#' Plot local predictor collinearity
#' @param x A result from [gwr_local_collinearity()].
#' @param coords Optional finite two-column coordinates in result row order.
#'   Supply coordinates for a spatial map; otherwise plot against location ID.
#' @param type Diagnostic to display.
#' @param ... Additional base graphics arguments.
#' @return The input object invisibly. Infinite/undefined map values have black
#'   crosses; finite map colors are explained by a value legend.
#' @export
plot.gwrs_local_collinearity <- function(x, coords = NULL,
                                        type = c("condition_number", "max_vif",
                                                 "max_abs_correlation"), ...) {
  type <- match.arg(type)
  value <- switch(type, condition_number = x$condition_number$condition_number,
                  max_vif = apply(x$vif, 1, max),
                  max_abs_correlation = x$max_abs_correlation)
  title <- paste("Local predictor", gsub("_", " ", type))
  if (is.null(coords)) {
    good <- is.finite(value)
    ylim <- if (any(good)) range(value[good]) else c(0, 1)
    graphics::plot(seq_along(value), replace(value, !good, NA_real_),
      xlab = "Location", ylab = title, ylim = ylim, pch = 19, ...)
    if (any(!good)) graphics::mtext("Infinite/undefined values omitted", side = 3)
  } else {
    coords <- validate_diagnostic_plot_coords(coords, length(value))
    colors <- gwrs_value_colors(replace(value, !is.finite(value), NA_real_))
    good <- is.finite(value)
    colors[!good] <- "black"
    graphics::plot(coords, col = colors, pch = ifelse(good, 19, 4),
                   xlab = "Coordinate 1", ylab = "Coordinate 2", asp = 1,
                   main = title, ...)
    if (any(good)) {
      limits <- range(value[good])
      labels <- if (diff(limits) == 0) limits[1] else seq(limits[1], limits[2], length.out = 5)
      legend_colors <- if (length(labels) == 1) colors[which(good)[1]] else
        gwrs_value_colors(labels)
      graphics::legend("topright", legend = format(signif(labels, 3)),
                       col = legend_colors, pch = 19, bty = "o", bg = "white",
                       title = "Value")
    }
    if (any(!good)) graphics::legend("bottomleft", "Infinite/undefined",
                                     col = "black", pch = 4, bty = "o", bg = "white")
  }
  invisible(x)
}

validate_diagnostic_plot_coords <- function(coords, n) {
  coords <- as_coordinate_matrix(coords)
  if (nrow(coords) != n || ncol(coords) != 2L) {
    stop("`coords` must have two columns and one row per location.", call. = FALSE)
  }
  coords
}

validate_gwrs_diagnostic_object <- function(object, require_x = FALSE) {
  if (!inherits(object, "gwrs_fit")) {
    stop("`object` must inherit from 'gwrs_fit'.", call. = FALSE)
  }
  if (is.null(object$neighbors)) {
    stop("The fitted object does not retain a neighbor structure.",
         call. = FALSE)
  }
  if (require_x && is.null(object$x)) {
    stop("Refit with `gwrs_control(keep_data = TRUE)` for this diagnostic.",
         call. = FALSE)
  }
  invisible(TRUE)
}

diagnostic_bandwidth <- function(object) {
  adaptive <- isTRUE(object$adaptive) || is.null(object$bandwidth)
  list(
    adaptive = adaptive,
    fixed = if (adaptive) 1 else as.double(object$bandwidth)
  )
}
