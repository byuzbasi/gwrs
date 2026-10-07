#' Plot fitted gwrs surfaces
#'
#' Draws coefficient, fitted-value, residual, local-fit, inference, and local
#' conditioning surfaces on the first two stored coordinate dimensions.
#'
#' @param x A fitted `gwrs_fit`.
#' @param type Surface to draw.
#' @param coefficient Coefficient name for coefficient-level surfaces.
#' @param coords Optional plotting coordinates, overriding stored coordinates.
#' @param estimator Inference estimator passed to [gwr_local_inference()].
#' @param alpha,p_adjust Local-inference display settings.
#' @param palette Base R HCL palette.
#' @param point_size Point expansion factor.
#' @param main Optional plot title.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_lasso_fit(
#'   x, quake$mag, coords = coords, k = 18, lambda = 0.05,
#'   control = gwrs_control(n_threads = 1)
#' )
#' plot(fit, type = "coefficient", coefficient = "depth")
#' plot(fit, type = "residual")
#'
#' @export
plot.gwrs_fit <- function(x,
                          type = c(
                            "coefficient", "standard_error", "statistic",
                            "p_value", "p_adjusted", "significance", "fitted",
                            "residual", "local_r2", "condition_number",
                            "effective_sample_size", "process"
                          ),
                          coefficient = NULL,
                          coords = NULL,
                          estimator = c("auto", "post_selection", "penalized"),
                          alpha = 0.05,
                          p_adjust = "BH",
                          palette = "Viridis",
                          point_size = 1,
                          main = NULL,
                          ...) {
  type <- match.arg(type)
  estimator <- match.arg(estimator)
  coordinates <- resolve_gwrs_plot_coords(x, coords)
  value_info <- gwrs_plot_values(
    x, type, coefficient, estimator, alpha, p_adjust
  )
  values <- value_info$values
  colors <- gwrs_value_colors(values, palette)
  graphics::plot(
    coordinates[, 1L], coordinates[, 2L], pch = 19, col = colors,
    cex = point_size, asp = 1,
    xlab = colnames(coordinates)[1L] %||% "Coordinate 1",
    ylab = colnames(coordinates)[2L] %||% "Coordinate 2",
    main = main %||% value_info$title,
    ...
  )
  invisible(x)
}

#' Convenience gwrs surface plots
#'
#' @inheritParams plot.gwrs_fit
#' @return Invisibly returns `object`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_lasso_fit(
#'   x, quake$mag, coords = coords, k = 18, lambda = 0.05,
#'   control = gwrs_control(n_threads = 1, diagnostics = "full")
#' )
#' plot_gwrs_model(fit, type = "coefficient", coefficient = "depth")
#' plot_gwrs_inference(fit, type = "significance", coefficient = "depth")
#' plot_gwrs_diagnostics(fit, type = "local_r2")
#'
#' @name gwrs_plot_helpers
NULL

#' @rdname gwrs_plot_helpers
#' @param object A fitted `gwrs_fit`.
#' @export
plot_gwrs_model <- function(object,
                            type = c("coefficient", "fitted", "residual",
                                     "local_r2", "process"),
                            ...) {
  plot(object, type = match.arg(type), ...)
  invisible(object)
}

#' @rdname gwrs_plot_helpers
#' @export
plot_gwrs_inference <- function(object,
                                type = c("standard_error", "statistic",
                                         "p_value", "p_adjusted",
                                         "significance"),
                                ...) {
  plot(object, type = match.arg(type), ...)
  invisible(object)
}

#' @rdname gwrs_plot_helpers
#' @export
plot_gwrs_diagnostics <- function(object,
                                  type = c("residual", "local_r2",
                                           "condition_number",
                                           "effective_sample_size"),
                                  ...) {
  plot(object, type = match.arg(type), ...)
  invisible(object)
}

#' Plot a gwrs regularization path
#'
#' @param x A `gwrs_path`.
#' @param type Plot an information criterion, mean active count, or convergence.
#' @param criterion Criterion column in `x$path_summary`.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return Invisibly returns `x`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' path <- gwr_lasso_path(
#'   x, quake$mag, coords = coords, k = 18,
#'   lambda = c(0.2, 0.08, 0.02),
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#' plot(path, type = "active")
#'
#' @export
plot.gwrs_path <- function(x,
                           type = c("criterion", "active", "convergence"),
                           criterion = c("aicc", "bic", "gcv", "ebic",
                                         "ebic_refit"),
                           ...) {
  type <- match.arg(type)
  criterion <- match.arg(criterion)
  if (type == "criterion") {
    if (!criterion %in% names(x$path_summary)) {
      stop(
        "Criterion plots require `gwr_sl_path(..., return_diagnostics = TRUE)`.",
        call. = FALSE
      )
    }
    values <- x$path_summary[[criterion]]
    ylab <- toupper(criterion)
  } else if (type == "active") {
    values <- x$path_summary$mean_nonzero
    ylab <- "Mean active slopes"
  } else {
    values <- x$path_summary$convergence_rate
    ylab <- "Convergence rate"
  }
  graphics::plot(
    log(x$lambda), values, type = "b", pch = 19,
    xlab = "log(lambda)", ylab = ylab, ...
  )
  invisible(x)
}

resolve_gwrs_plot_coords <- function(object, coords = NULL) {
  if (is.null(coords)) {
    coords <- object$neighbors$query_coords %||% object$neighbors$coords
  }
  coords <- as_coordinate_matrix(coords)
  if (nrow(coords) != nrow(object$coefficients) || ncol(coords) < 2L) {
    stop("Plot coordinates must have at least two columns and one row per location.",
         call. = FALSE)
  }
  coords[, 1:2, drop = FALSE]
}

gwrs_plot_values <- function(object,
                             type,
                             coefficient,
                             estimator,
                             alpha,
                             p_adjust) {
  coefficient_types <- c(
    "coefficient", "standard_error", "statistic", "p_value",
    "p_adjusted", "significance"
  )
  if (type %in% coefficient_types) {
    coefficient <- coefficient %||% colnames(object$coefficients)[1L]
    column <- match(coefficient, colnames(object$coefficients))
    if (is.na(column)) stop("Unknown coefficient: ", coefficient, call. = FALSE)
    if (type == "coefficient") {
      return(list(
        values = object$coefficients[, column],
        title = paste("Local coefficient:", coefficient)
      ))
    }
    inference <- gwr_local_inference(
      object, alpha = alpha, adjust = p_adjust, estimator = estimator
    )
    local <- inference$local[inference$local$coefficient == coefficient, ]
    values <- switch(
      type,
      standard_error = local$standard_error,
      statistic = local$statistic,
      p_value = local$p_value,
      p_adjusted = local$p_adjusted,
      significance = as.numeric(local$significant)
    )
    return(list(
      values = values,
      title = paste(gsub("_", " ", type), coefficient)
    ))
  }
  diagnostics <- gwr_diagnostics(object)
  values <- switch(
    type,
    fitted = object$fitted.values,
    residual = object$residuals,
    local_r2 = diagnostics$local$local_r_squared,
    condition_number = {
      if ("condition_number" %in% names(diagnostics$local)) {
        diagnostics$local$condition_number
      } else {
        gwr_local_collinearity(object)$condition_number$condition_number
      }
    },
    effective_sample_size = diagnostics$local$effective_sample_size,
    process = {
      if (is.null(object$process)) {
        stop("A process surface is available only for GWR-KTDD fits.",
             call. = FALSE)
      }
      object$process
    }
  )
  list(values = values, title = gsub("_", " ", type))
}

gwrs_value_colors <- function(values, palette = "Viridis", n = 100L) {
  colors <- grDevices::hcl.colors(n, palette)
  finite <- is.finite(values)
  if (!any(finite)) return(rep("grey75", length(values)))
  range <- base::range(values[finite])
  if (diff(range) <= .Machine$double.eps) {
    result <- rep(colors[ceiling(n / 2)], length(values))
  } else {
    index <- cut(values, breaks = n, include.lowest = TRUE, labels = FALSE)
    result <- colors[pmax(1L, pmin(n, index))]
  }
  result[!finite] <- "grey75"
  result
}
