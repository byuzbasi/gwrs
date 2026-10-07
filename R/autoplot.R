#' ggplot2 autoplot method for fitted gwrs models
#'
#' Routes continuous surface and bivariate requests to [plot_heatmap()] and
#' [plot_bivariate()]. The existing base-R [plot.gwrs_fit()] method is left
#' unchanged. A polygon choropleth is not fabricated because fitted `gwrs`
#' objects store point coordinates rather than observation-level polygons.
#'
#' @param object A fitted `gwrs_fit` object.
#' @param type Heat-map, contour, combined, or bivariate output.
#' @param variable Primary spatial variable.
#' @param variable2 Secondary variable for bivariate output.
#' @param ... Arguments passed to [plot_heatmap()] or [plot_bivariate()].
#'
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("sf", quietly = TRUE)) {
#'   set.seed(11)
#'   coords <- cbind(runif(30), runif(30))
#'   x <- cbind(income = rnorm(30), education = rnorm(30))
#'   y <- 1 + x[, "income"] - 0.5 * x[, "education"] + rnorm(30, sd = 0.2)
#'   fit <- gwr_fit(
#'     x, y, coords = coords, k = 20,
#'     control = gwrs_control(n_threads = 1)
#'   )
#'   ggplot2::autoplot(
#'     fit, type = "heatmap", variable = "income", resolution = 30
#'   )
#' }
#' @exportS3Method ggplot2::autoplot
autoplot.gwrs_fit <- function(object,
                              type = c(
                                "heatmap", "contour", "heatmap_contour",
                                "bivariate", "choropleth"
                              ),
                              variable = NULL,
                              variable2 = NULL,
                              ...) {
  type <- match.arg(type)
  if (!inherits(object, "gwrs_fit")) {
    stop("`object` must inherit from 'gwrs_fit'.", call. = FALSE)
  }
  if (type == "choropleth") {
    stop(
      paste(
        "A fitted `gwrs_fit` stores point coordinates, not one polygon per",
        "observation, so a choropleth cannot be constructed safely.",
        "Use `type = \"heatmap\"` or `plot_heatmap()` with a boundary."
      ),
      call. = FALSE
    )
  }
  if (is.null(variable)) {
    coefficient_names <- colnames(object$coefficients)
    variable <- if (length(coefficient_names) > 1L) {
      coefficient_names[2L]
    } else {
      coefficient_names[1L]
    }
  }
  if (type == "bivariate") {
    if (is.null(variable2)) {
      stop("`variable2` is required for a bivariate autoplot.",
           call. = FALSE)
    }
    return(plot_bivariate(
      object, variable1 = variable, variable2 = variable2, ...
    ))
  }
  plot_heatmap(object, variable = variable, type = type, ...)
}
