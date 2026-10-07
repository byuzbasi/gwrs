#' Publication maps from fitted gwrs models
#'
#' Compose existing local results into a six-panel dashboard or a coefficient
#' atlas. A dashboard shows observed coefficients, a bivariate coefficient and
#' local-R-squared surface, coefficient heat map, heat map with contours, local
#' R-squared and residuals. All continuous panels reuse one visual bandwidth.
#' No regression is fitted by these functions.
#'
#' @param object A fitted `gwrs_fit` with at least standard diagnostics for
#'   the dashboard. The atlas only requires stored coefficients and coordinates.
#' @param variable One coefficient name or alias, such as `"beta_1"`.
#' @param variables Coefficient names for the atlas. Defaults to all slopes,
#'   excluding the intercept. Select a subset for high-dimensional fits.
#' @param boundary Plotting boundary passed to [plot_heatmap()].
#' @param overlays Optional `sf`/`sfc` boundaries or a list of such layers.
#'   Each must have a known CRS, as must the model/plot coordinates. Layers
#'   are transformed to the plot CRS and drawn as outlines, without aggregation.
#' @param resolution Grid resolution passed to [plot_heatmap()].
#' @param bandwidth Visual interpolation bandwidth, distinct from model bandwidth.
#'   With kernel interpolation and `NULL`, the first map supplies the shared value.
#' @param selection Display locally unselected support as `"gray"`, `"hide"`,
#'   or `"show"`. Observed support comes from the stored coefficient extractor;
#'   grid support uses the existing `selection_threshold` interpolation rule.
#' @param unselected_color Colour for locally unselected coefficients.
#' @param coefficient_limits Optional increasing finite limits, containing the
#'   midpoint, shared by observed coefficients and both continuous coefficient
#'   panels. By default they are symmetric about zero (or `midpoint` in `...`).
#' @param scales Atlas colour scales: `"free"` per coefficient or `"fixed"`
#'   across all requested coefficients. Fixed scales are appropriate only when
#'   coefficient units and ranges are comparable.
#' @param labels Optional named character vector mapping coefficient names to
#'   display labels. Unspecified labels retain the coefficient name.
#' @param ncol Number of panel columns.
#' @param title Overall title; `NULL` chooses a model-labelled title.
#' @param return_data Return panels, underlying surfaces and settings alongside
#'   the composite when `TRUE`. Defaults to returning a `patchwork` plot.
#' @param ... Named [plot_heatmap()] options, such as `crs`, `max_neighbors`,
#'   `chunk_size`, `mask`, `selection_threshold`, `palette`, or `significance`.
#'   Significance overlays apply to continuous coefficient panels only.
#'
#' @return A `patchwork` object accepted by [ggplot2::ggsave()], or a list with
#'   `plot`, `panels`, `surfaces`, and `settings` when `return_data = TRUE`.
#'   Plot attributes include `gwrs_publication_settings`.
#'
#' @details Grey means locally unselected, not statistically nonsignificant.
#'   Optional significance overlays remain conditional exploratory displays.
#'   Interpolated surfaces are visual summaries, not new-location predictions.
#'   The observed panel retains individual observations, including colocated
#'   observations. Polygon layers describe boundaries, not observation areas.
#'   The default atlas is deliberately not given a common coefficient scale:
#'   a one-unit slope has different meaning for differently scaled predictors.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("sf", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   set.seed(42)
#'   coords <- cbind(east = runif(40), north = runif(40))
#'   x <- cbind(income = rnorm(40), access = rnorm(40))
#'   y <- 1 + (2 * coords[, 1] - 1) * x[, 1] + rnorm(40, sd = 0.2)
#'   fit <- gwr_lasso_fit(x, y, coords = coords, k = 25, lambda = 0.1,
#'                       control = gwrs_control(n_threads = 1))
#'   plot_gwrs_dashboard(fit, "income", resolution = 25)
#'   plot_gwrs_atlas(fit, resolution = 25)
#' }
#' @name gwrs_publication_plots
NULL

#' @rdname gwrs_publication_plots
#' @export
plot_gwrs_dashboard <- function(object, variable, boundary = NULL,
                                overlays = NULL, resolution = 150L,
                                bandwidth = NULL, selection = "gray",
                                unselected_color = "grey80",
                                coefficient_limits = NULL, ncol = 3L,
                                title = NULL, return_data = FALSE, ...) {
  options <- gwrs_publication_options(ncol, return_data, list(...))
  selection <- match.arg(selection, c("gray", "show", "hide"))
  extracted <- extract_gwrs_spatial_values(object, variable)
  if (extracted$value_type != "coefficient") {
    stop("`variable` must identify a coefficient.", call. = FALSE)
  }
  common <- c(list(object = object, boundary = boundary,
                   resolution = resolution, bandwidth = bandwidth,
                   selection = selection, unselected_color = unselected_color,
                   return_data = TRUE), options)
  beta <- do.call(plot_heatmap, c(common, list(variable = variable)))
  common$bandwidth <- beta$bandwidth
  contours <- do.call(plot_heatmap, c(common, list(
    variable = variable, type = "heatmap_contour"
  )))
  diagnostic_options <- common
  diagnostic_options$selection <- "show"
  diagnostic_options$significance <- FALSE
  diagnostic_options$palette <- NULL
  diagnostic_options$midpoint <- NULL
  local_r2 <- do.call(plot_heatmap, c(diagnostic_options,
                                    list(variable = "local_r2")))
  residual <- do.call(plot_heatmap, c(diagnostic_options,
                                    list(variable = "residuals")))
  bivariate_options <- common[intersect(names(common), names(formals(plot_bivariate)))]
  bivariate <- do.call(plot_bivariate, c(bivariate_options, list(
    variable1 = variable, variable2 = "local_r2"
  )))
  # Missing/outside cells stay transparent, rather than becoming a legend class.
  bivariate_scale <- bivariate$plot$scales$get_scales("fill")
  bivariate_scale$na.value <- NA
  bivariate_scale$na.translate <- FALSE
  # Compact bin pairs retain the existing low-to-high class order. Keep this
  # categorical legend below its panel so it cannot displace adjacent maps.
  bivariate_scale$labels <- function(labels) {
    gsub("V1:|V2:| ", "", labels)
  }
  bivariate_scale$name <- "Coefficient / local R-squared bins\n1 = low, 3 = high"
  bivariate$plot$layers[[1L]]$show.legend <- TRUE
  bivariate$plot <- bivariate$plot +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3, byrow = TRUE))
  limits <- gwrs_publication_limits(beta$points$value_1,
                                    coefficient_limits, options$midpoint)
  beta$plot <- gwrs_publication_scale(beta$plot, limits, options)
  contours$plot <- gwrs_publication_scale(contours$plot, limits, options)
  valid_contours <- contours$grid$value[contours$grid$selected &
                                       contours$grid$interpolation_valid]
  if (length(unique(valid_contours[is.finite(valid_contours)])) < 2L) {
    contours$plot$layers <- Filter(function(layer) {
      !inherits(layer$geom, "GeomContour")
    }, contours$plot$layers)
  }
  points <- beta$points
  points$value <- points$value_1
  points$selected <- points$selected_1 >= 0.5
  displayed <- if (selection == "show") rep(TRUE, nrow(points)) else points$selected
  observed <- ggplot2::ggplot()
  if (selection == "gray") {
    observed <- observed + ggplot2::geom_point(
      data = points[!points$selected, , drop = FALSE],
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = unselected_color, size = 1.2
    )
  }
  observed <- observed + ggplot2::geom_point(
    data = points[displayed, , drop = FALSE],
    ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["value"]]),
    size = 1.2
  ) + ggplot2::geom_sf(data = gwrs_boundary_sf(beta$boundary), fill = NA,
                      colour = "grey45", linewidth = 0.25, inherit.aes = FALSE) +
    ggplot2::coord_sf(datum = NA, expand = FALSE)
  observed <- gwrs_publication_scale(observed, limits, options, "colour")
  surfaces <- list(beta = beta, bivariate = bivariate, contour = contours,
                   local_r2 = local_r2, residuals = residual)
  panels <- c(list(observed = observed), lapply(surfaces[c(
    "bivariate", "beta", "contour", "local_r2", "residuals"
  )], `[[`, "plot"))
  panel_titles <- c("1. Observed coefficients", "2. Coefficient and local R-squared",
                    "3. Coefficient heat map", "4. Coefficients with contours",
                    "5. Local R-squared", "6. Residuals")
  overlay_layers <- gwrs_publication_overlays(overlays, beta$settings$crs)
  panels <- Map(function(plot, label) {
    gwrs_publication_style(plot, label, overlay_layers)
  }, panels, panel_titles)
  panels$bivariate <- panels$bivariate + ggplot2::theme(
    legend.position = "bottom", legend.key.size = grid::unit(3, "mm")
  )
  settings <- list(variable = extracted$variable, bandwidth = beta$bandwidth,
                   coefficient_limits = limits, selection = selection,
                   observed_rows = nrow(points),
                   surface_settings = lapply(surfaces, `[[`, "settings"),
                   visual_interpolation_only = TRUE)
  gwrs_publication_result(panels, surfaces, settings, ncol,
                          title %||% paste(object$method, "spatial results:",
                                           extracted$variable), return_data,
                          isTRUE(options$significance))
}

#' @rdname gwrs_publication_plots
#' @export
plot_gwrs_atlas <- function(object, variables = NULL, boundary = NULL,
                            overlays = NULL, resolution = 150L,
                            bandwidth = NULL, selection = "gray",
                            unselected_color = "grey80", scales = c("free", "fixed"),
                            labels = NULL, ncol = 2L, title = NULL,
                            return_data = FALSE, ...) {
  options <- gwrs_publication_options(ncol, return_data, list(...))
  selection <- match.arg(selection, c("gray", "show", "hide"))
  scales <- match.arg(scales)
  if (!inherits(object, "gwrs_fit")) stop("`object` must be a gwrs_fit.", call. = FALSE)
  available <- colnames(object$coefficients)
  if (is.null(variables)) variables <- setdiff(available, "(Intercept)")
  if (!is.character(variables) || !length(variables) || anyNA(variables) ||
      anyDuplicated(variables) || any(!variables %in% available)) {
    stop("`variables` must be distinct coefficient names present in the fit.", call. = FALSE)
  }
  if (!is.null(labels) && (!is.character(labels) || anyNA(labels) ||
      is.null(names(labels)) || anyNA(names(labels)) || any(!nzchar(names(labels))) ||
      anyDuplicated(names(labels)) || any(!names(labels) %in% variables))) {
    stop("`labels` must be a named character vector for requested variables.", call. = FALSE)
  }
  common <- c(list(object = object, boundary = boundary, resolution = resolution,
                   bandwidth = bandwidth, selection = selection,
                   unselected_color = unselected_color, return_data = TRUE), options)
  fixed_limits <- if (scales == "fixed") {
    gwrs_publication_limits(object$coefficients[, variables, drop = FALSE],
                            NULL, options$midpoint)
  } else NULL
  surfaces <- panels <- stats::setNames(vector("list", length(variables)), variables)
  for (variable in variables) {
    surface <- do.call(plot_heatmap, c(common, list(variable = variable)))
    if (is.null(common$bandwidth)) common$bandwidth <- surface$bandwidth
    limits <- fixed_limits %||% gwrs_publication_limits(
      surface$points$value_1, NULL, options$midpoint
    )
    surface$plot <- gwrs_publication_scale(surface$plot, limits, options)
    label <- if (variable %in% names(labels)) labels[[variable]] else variable
    panels[[variable]] <- gwrs_publication_style(surface$plot, label,
      gwrs_publication_overlays(overlays, surface$settings$crs))
    surfaces[[variable]] <- surface
  }
  settings <- list(variables = variables, bandwidth = common$bandwidth,
                   scales = scales, selection = selection,
                   surface_settings = lapply(surfaces, `[[`, "settings"),
                   visual_interpolation_only = TRUE)
  gwrs_publication_result(panels, surfaces, settings, ncol,
                          title %||% paste(object$method, "coefficient atlas"),
                          return_data, isTRUE(options$significance))
}

gwrs_publication_options <- function(ncol, return_data, options) {
  require_gwrs_visual_packages()
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install the optional 'patchwork' package to compose publication maps.", call. = FALSE)
  }
  if (!is.numeric(ncol) || length(ncol) != 1L || !is.finite(ncol) ||
      ncol < 1 || ncol != floor(ncol)) stop("`ncol` must be a positive integer.", call. = FALSE)
  if (!is.logical(return_data) || length(return_data) != 1L || is.na(return_data)) {
    stop("`return_data` must be TRUE or FALSE.", call. = FALSE)
  }
  reserved <- c("object", "variable", "type", "boundary", "resolution", "bandwidth",
                "selection", "unselected_color", "return_data", "...")
  allowed <- setdiff(names(formals(plot_heatmap)), reserved)
  if (length(options) && (is.null(names(options)) || anyNA(names(options)) ||
      any(!names(options) %in% allowed) || anyDuplicated(names(options)))) {
    stop("Supply distinct, named plot_heatmap options in `...`; type is set by the panel.",
         call. = FALSE)
  }
  options
}

gwrs_publication_limits <- function(values, limits, midpoint = NULL) {
  midpoint <- midpoint %||% 0
  if (is.null(limits)) {
    values <- values[is.finite(values)]
    radius <- if (length(values)) max(abs(values - midpoint)) else 0
    if (radius == 0) radius <- 1
    limits <- midpoint + c(-radius, radius)
  }
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) ||
      limits[1L] >= limits[2L] || midpoint < limits[1L] || midpoint > limits[2L]) {
    stop("`coefficient_limits` must increase and contain the midpoint.", call. = FALSE)
  }
  limits
}

gwrs_publication_scale <- function(plot, limits, options, aesthetic = "fill") {
  colors <- gwrs_heatmap_palette(options$palette, TRUE)
  plot$scales$scales <- Filter(function(scale) {
    !aesthetic %in% scale$aesthetics
  }, plot$scales$scales)
  scale <- if (aesthetic == "fill") ggplot2::scale_fill_gradient2 else
    ggplot2::scale_colour_gradient2
  plot + scale(low = colors[1L], mid = colors[ceiling(length(colors) / 2)],
                high = colors[length(colors)], midpoint = options$midpoint %||% 0,
                limits = limits, na.value = NA, name = "Coefficient",
                oob = function(x, range, ...) pmin(pmax(x, range[1L]), range[2L]))
}

gwrs_publication_overlays <- function(overlays, crs) {
  if (is.null(overlays)) return(list())
  if (inherits(overlays, c("sf", "sfc"))) overlays <- list(overlays)
  if (!is.list(overlays) || is.na(sf::st_crs(crs))) {
    stop("`overlays` require sf/sfc layers and a known plot CRS.", call. = FALSE)
  }
  lapply(overlays, function(layer) {
    if (!inherits(layer, c("sf", "sfc")) || is.na(sf::st_crs(layer))) {
      stop("Every overlay must be an sf/sfc object with a known CRS.", call. = FALSE)
    }
    ggplot2::geom_sf(data = sf::st_transform(layer, crs), fill = NA,
                     colour = "grey45", linewidth = 0.2, inherit.aes = FALSE)
  })
}

gwrs_publication_style <- function(plot, label, overlays) {
  plot + overlays + ggplot2::labs(title = label, subtitle = NULL, caption = NULL) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                   legend.title = ggplot2::element_text(size = 8),
                   legend.text = ggplot2::element_text(size = 7),
                   plot.margin = ggplot2::margin(6, 6, 6, 6))
}

gwrs_publication_result <- function(panels, surfaces, settings, ncol, title,
                                    return_data, significance) {
  caption <- paste(
    "Observed values are point estimates; surfaces are visual interpolations.",
    if (settings$selection == "gray") "Grey: locally unselected support." else "",
    if (significance) "Pale overlay: selected but not significant (exploratory)." else ""
  )
  plot <- patchwork::wrap_plots(panels, ncol = ncol) +
    patchwork::plot_annotation(title = title, caption = trimws(caption),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                            plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                            plot.caption = ggplot2::element_text(size = 8)))
  attr(plot, "gwrs_publication_settings") <- settings
  if (!return_data) return(plot)
  list(plot = plot, panels = panels, surfaces = surfaces, settings = settings)
}
