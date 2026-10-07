gwrs_heatmap_palette <- function(palette, diverging) {
  if (is.null(palette)) {
    return(if (diverging) {
      c("#3B4CC0", "#F7F7F7", "#B40426")
    } else {
      grDevices::hcl.colors(256L, "Viridis")
    })
  }
  if (!is.character(palette) || !length(palette) || anyNA(palette)) {
    stop("`palette` must be NULL, an HCL palette name, or color values.",
         call. = FALSE)
  }
  if (length(palette) == 1L) {
    colors <- tryCatch(
      grDevices::hcl.colors(256L, palette),
      error = function(error) NULL
    )
    if (is.null(colors)) {
      stop("Unknown HCL palette name: ", palette, call. = FALSE)
    }
    return(colors)
  }
  invalid <- vapply(
    palette,
    function(value) inherits(try(grDevices::col2rgb(value), silent = TRUE),
                             "try-error"),
    logical(1L)
  )
  if (any(invalid)) stop("`palette` contains invalid colors.", call. = FALSE)
  palette
}

gwrs_add_fill_scale <- function(plot,
                                extracted,
                                values,
                                palette,
                                midpoint) {
  diverging <- extracted$value_type %in%
    c("coefficient", "residuals", "process")
  colors <- gwrs_heatmap_palette(palette, diverging)
  if (diverging) {
    midpoint <- midpoint %||% 0
    if (!is.numeric(midpoint) || length(midpoint) != 1L ||
        !is.finite(midpoint)) {
      stop("`midpoint` must be one finite number.", call. = FALSE)
    }
    low <- colors[1L]
    middle <- colors[ceiling(length(colors) / 2)]
    high <- colors[length(colors)]
    return(plot + ggplot2::scale_fill_gradient2(
      low = low, mid = middle, high = high, midpoint = midpoint,
      na.value = NA, name = extracted$label
    ))
  }
  limits <- if (extracted$value_type == "local_r2" &&
                all(values[is.finite(values)] >= 0 &
                    values[is.finite(values)] <= 1)) {
    c(0, 1)
  } else {
    NULL
  }
  plot + ggplot2::scale_fill_gradientn(
    colors = colors, limits = limits, na.value = NA,
    name = extracted$label
  )
}

gwrs_validate_heatmap_flags <- function(mask,
                                        significance,
                                        show_points,
                                        show_boundaries,
                                        return_data,
                                        alpha,
                                        contour_bins,
                                        nonsignificant_alpha) {
  flags <- list(
    mask = mask, significance = significance, show_points = show_points,
    show_boundaries = show_boundaries, return_data = return_data
  )
  if (any(!vapply(flags, function(value) {
    is.logical(value) && length(value) == 1L && !is.na(value)
  }, logical(1L)))) {
    stop("Plot switches must each be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.numeric(contour_bins) || length(contour_bins) != 1L ||
      !is.finite(contour_bins) || contour_bins != floor(contour_bins) ||
      contour_bins < 2L || contour_bins > 100L) {
    stop("`contour_bins` must be an integer between 2 and 100.",
         call. = FALSE)
  }
  contour_bins <- as.integer(contour_bins)
  if (!is.numeric(nonsignificant_alpha) ||
      length(nonsignificant_alpha) != 1L ||
      !is.finite(nonsignificant_alpha) || nonsignificant_alpha < 0 ||
      nonsignificant_alpha > 1) {
    stop("`nonsignificant_alpha` must be one number in [0, 1].",
         call. = FALSE)
  }
  contour_bins
}

gwrs_heatmap_caption <- function(extracted, selection, significance) {
  notes <- character()
  if (selection == "gray" && any(!extracted$points$selected)) {
    notes <- c(notes, "Grey: locally unselected coefficient support")
  }
  if (isTRUE(significance)) {
    notes <- c(
      notes,
      "Pale grey: selected but not significant at the requested level"
    )
  }
  if (!length(notes)) return(NULL)
  paste(notes, collapse = "; ")
}

prepare_gwrs_heatmap_result <- function(object,
                                        variable,
                                        interpolation,
                                        distance,
                                        kernel,
                                        bandwidth,
                                        resolution,
                                        boundary,
                                        mask,
                                        significance,
                                        alpha,
                                        selection_threshold,
                                        crs,
                                        idw_power,
                                        max_neighbors,
                                        chunk_size,
                                        extrapolation,
                                        p_adjust,
                                        estimator,
                                        na.rm) {
  extracted <- extract_gwrs_spatial_values(
    object, variable, na.rm = na.rm
  )
  if (isTRUE(significance)) {
    extracted <- add_gwrs_spatial_significance(
      extracted, object, alpha, p_adjust, estimator
    )
  }
  surface <- prepare_gwrs_spatial_surface(
    object = object,
    extractions = list(extracted),
    interpolation = interpolation,
    distance = distance,
    kernel = kernel,
    bandwidth = bandwidth,
    resolution = resolution,
    boundary = boundary,
    mask = mask,
    crs = crs,
    max_neighbors = max_neighbors,
    chunk_size = chunk_size,
    idw_power = idw_power,
    selection_threshold = selection_threshold,
    extrapolation = extrapolation
  )
  surface$grid$value <- surface$grid$value_1
  surface$grid$selection_support <- surface$grid$selected_1
  surface$grid$selected <- is.finite(surface$grid$selected_1) &
    surface$grid$selected_1 >= selection_threshold
  surface$plot_grid$value <- surface$plot_grid$value_1
  surface$plot_grid$selection_support <- surface$plot_grid$selected_1
  surface$plot_grid$selected <- is.finite(surface$plot_grid$selected_1) &
    surface$plot_grid$selected_1 >= selection_threshold
  if (isTRUE(significance)) {
    surface$grid$significance_support <- surface$grid$significant_1
    surface$grid$significant <- is.finite(surface$grid$significant_1) &
      surface$grid$significant_1 >= selection_threshold
    surface$plot_grid$significance_support <-
      surface$plot_grid$significant_1
    surface$plot_grid$significant <-
      is.finite(surface$plot_grid$significant_1) &
      surface$plot_grid$significant_1 >= selection_threshold
  } else {
    surface$grid$significant <- TRUE
    surface$plot_grid$significant <- TRUE
  }
  surface$extracted <- extracted
  surface
}

#' Continuous spatial heat maps for fitted gwrs results
#'
#' Interpolates observed local coefficients, diagnostics, fitted values, or
#' residuals onto a regular geographic grid and returns a `ggplot` surface.
#' A heat map is a visually smoothed representation of fitted local results;
#' it is not a model prediction at new locations. Use [predict.gwrs_fit()] for
#' formal new-location predictions.
#'
#' A continuous heat map is a visually smoothed representation of observed
#' local `gwrs` model results. It must not be interpreted as an official
#' `gwrs` prediction at new locations. Predictions at new locations should be
#' obtained through the model's [predict.gwrs_fit()] method.
#'
#' Sparse coefficient surfaces retain a separate selection-support layer.
#' Unselected regions are grey by default. Optional significance masking is a
#' distinct exploratory layer derived from [gwr_local_inference()] and requires
#' a Gaussian fit created with `gwrs_control(diagnostics = "full")`.
#'
#' @param object A fitted object inheriting from `gwrs_fit`.
#' @param variable Coefficient name, `beta_0`, `beta_1`, ..., or one of
#'   `"local_r2"`, `"residuals"`, `"abs_residuals"`, `"fitted"`,
#'   `"response"`, and `"process"` when available.
#' @param type Draw a heat map, contour lines, or both.
#' @param interpolation Kernel smoothing or inverse-distance weighting.
#' @param distance Distance used only for visual interpolation. Euclidean is
#'   appropriate for projected/Cartesian coordinates; geographic uses
#'   earth-surface metres. Hyperbolic interpolation errors unless a valid
#'   new-location mapping exists, which current `gwrs_fit` objects do not
#'   provide.
#' @param kernel Kernel for kernel interpolation.
#' @param bandwidth Positive kernel bandwidth. `NULL` uses the median distance
#'   to a deterministic k-nearest distinct location.
#' @param resolution Number of grid points on the longer spatial axis. The
#'   shorter axis is scaled to the bounding-box aspect ratio, so 300 produces
#'   at most approximately 90,000 cells.
#' @param boundary Optional `sf` polygon/multipolygon or numeric x/y polygon.
#'   A boundary with a different known CRS is transformed to the coordinate
#'   CRS before grid construction. Multiple supplied polygon outlines are
#'   retained for the boundary overlay.
#' @param mask Whether to retain only cells inside the boundary or default
#'   buffered convex hull. `FALSE` shows the complete bounding box. With
#'   masking, cells farther than the extrapolation threshold remain blank.
#' @param significance Whether to fade locally nonsignificant coefficient
#'   support.
#' @param alpha Significance level.
#' @param selection Treatment of sparse-model regions not selected locally.
#' @param selection_threshold Interpolated binary-support threshold.
#' @param unselected_color Color for unselected coefficient regions.
#' @param nonsignificant_color,nonsignificant_alpha Overlay for selected but
#'   nonsignificant regions.
#' @param show_points Whether to overlay observed locations.
#' @param show_boundaries Whether to draw the plotting boundary.
#' @param palette `NULL`, an HCL palette name, or color vector.
#' @param midpoint Diverging-scale midpoint; coefficients, residuals, and KTDD
#'   process fields default to zero.
#' @param na.rm Whether observations with non-finite requested values are
#'   removed.
#' @param return_data Return the plot and interpolation data in a list.
#' @param crs Optional CRS for the model coordinate matrix, for example 4326.
#' @param idw_power Positive inverse-distance exponent.
#' @param max_neighbors Maximum source locations used per grid cell. This
#'   nearest-neighbor truncation bounds work and memory on large data. Use
#'   `Inf` for the exact all-observation smoother.
#' @param chunk_size Number of grid cells processed in each distance block.
#' @param extrapolation Maximum nearest-observation distance as a multiple of
#'   the median nearest-neighbor distance. `NULL` uses 3 with masking and `Inf`
#'   without masking.
#' @param contour_bins Number of contour levels.
#' @param p_adjust Multiplicity adjustment used for significance masking.
#' @param estimator Inference estimator passed to [gwr_local_inference()].
#' @param ... Reserved for future graphical controls.
#'
#' @return A `ggplot` object, or when `return_data = TRUE`, a list containing
#'   `plot`, `grid`, `points`, `boundary`, `bandwidth`, and `settings`.
#'
#' @details
#' With `bandwidth = NULL`, the kernel bandwidth is the median distance to a
#' deterministic nearest-neighbor order: the smaller of 20 and the square
#' root of the number of distinct locations, subject to the available sample.
#' Geographic distance transforms known source coordinates to longitude and
#' latitude and reports distances in metres. Euclidean distance uses the
#' supplied coordinate units. No artificial hyperbolic coordinates are
#' constructed for grid cells.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("sf", quietly = TRUE)) {
#'   set.seed(42)
#'   coords <- cbind(runif(36), runif(36))
#'   x <- cbind(
#'     income = rnorm(36), education = rnorm(36), noise = rnorm(36)
#'   )
#'   y <- 0.5 + (2 * coords[, 1] - 1) * x[, "income"] +
#'     0.6 * x[, "education"] + rnorm(36, sd = 0.2)
#'   fit <- gwr_en_fit(
#'     x, y, coords = coords, k = 24, lambda = 0.35,
#'     alpha = 0.85,
#'     control = gwrs_control(n_threads = 1)
#'   )
#'   mapped <- plot_heatmap(
#'     fit, "income", resolution = 30, selection = "gray",
#'     unselected_color = "grey75", return_data = TRUE
#'   )
#'   mapped$plot
#'   table(selected = mapped$grid$selected, useNA = "ifany")
#'
#'   # Coefficient aliases and names resolve to the same observed values.
#'   alias <- plot_heatmap(
#'     fit, "beta_1", resolution = 30, return_data = TRUE
#'   )
#'   all.equal(mapped$points$value, alias$points$value)
#'   plot_heatmap(
#'     fit, "local_r2", type = "heatmap_contour", resolution = 30
#'   )
#' }
#' @export
plot_heatmap <- function(object,
                         variable,
                         type = c("heatmap", "contour", "heatmap_contour"),
                         interpolation = c("kernel", "idw"),
                         distance = c("euclidean", "geographic", "hyperbolic"),
                         kernel = c("gaussian", "bisquare", "exponential"),
                         bandwidth = NULL,
                         resolution = 300L,
                         boundary = NULL,
                         mask = TRUE,
                         significance = FALSE,
                         alpha = 0.05,
                         selection = c("gray", "show", "hide"),
                         selection_threshold = 0.5,
                         unselected_color = "grey80",
                         nonsignificant_color = "grey94",
                         nonsignificant_alpha = 0.65,
                         show_points = FALSE,
                         show_boundaries = TRUE,
                         palette = NULL,
                         midpoint = NULL,
                         na.rm = TRUE,
                         return_data = FALSE,
                         crs = NULL,
                         idw_power = 2,
                         max_neighbors = 200L,
                         chunk_size = 2000L,
                         extrapolation = NULL,
                         contour_bins = 10L,
                         p_adjust = "BH",
                         estimator = c("auto", "post_selection", "penalized"),
                         ...) {
  type <- match.arg(type)
  interpolation <- match.arg(interpolation)
  distance <- match.arg(distance)
  kernel <- match.arg(kernel)
  selection <- match.arg(selection)
  estimator <- match.arg(estimator)
  contour_bins <- gwrs_validate_heatmap_flags(
    mask, significance, show_points, show_boundaries, return_data,
    alpha, contour_bins, nonsignificant_alpha
  )
  surface <- prepare_gwrs_heatmap_result(
    object, variable, interpolation, distance, kernel, bandwidth,
    resolution, boundary, mask, significance, alpha,
    selection_threshold, crs, idw_power, max_neighbors, chunk_size,
    extrapolation, p_adjust, estimator, na.rm
  )
  grid <- surface$plot_grid
  grid$display_value <- grid$value
  if (selection != "show") {
    grid$display_value[!grid$selected] <- NA_real_
  }
  grid$contour_value <- grid$display_value
  plot <- ggplot2::ggplot()
  if (type %in% c("heatmap", "heatmap_contour")) {
    plot <- plot + ggplot2::geom_raster(
      data = grid,
      mapping = ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], fill = .data[["display_value"]]
      ),
      na.rm = TRUE
    )
    plot <- gwrs_add_fill_scale(
      plot, surface$extracted, grid$display_value, palette, midpoint
    )
  }
  if (type %in% c("contour", "heatmap_contour")) {
    plot <- plot + ggplot2::geom_contour(
      data = grid,
      mapping = ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]], z = .data[["contour_value"]]
      ),
      bins = contour_bins, color = "grey20", linewidth = 0.35,
      na.rm = TRUE
    )
  }
  if (selection == "gray") {
    unselected <- grid[
      grid$inside_boundary & grid$interpolation_valid & !grid$selected,
      , drop = FALSE
    ]
    if (nrow(unselected)) {
      plot <- plot + ggplot2::geom_raster(
        data = unselected,
        mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
        fill = unselected_color, inherit.aes = FALSE
      )
    }
  }
  if (isTRUE(significance)) {
    nonsignificant <- grid[
      grid$inside_boundary & grid$interpolation_valid & grid$selected &
        !grid$significant,
      , drop = FALSE
    ]
    if (nrow(nonsignificant)) {
      plot <- plot + ggplot2::geom_raster(
        data = nonsignificant,
        mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
        fill = nonsignificant_color, alpha = nonsignificant_alpha,
        inherit.aes = FALSE
      )
    }
  }
  boundary_color <- if (isTRUE(surface$settings$boundary_provided)) {
    "white"
  } else {
    "grey35"
  }
  if (isTRUE(show_boundaries)) {
    plot <- plot + ggplot2::geom_sf(
      data = gwrs_boundary_sf(surface$boundary), fill = NA,
      color = boundary_color,
      linewidth = 0.45, inherit.aes = FALSE
    )
  }
  if (isTRUE(show_points)) {
    plot <- plot + ggplot2::geom_point(
      data = surface$points,
      mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      size = 0.55, color = "black", alpha = 0.65, inherit.aes = FALSE
    )
  }
  map_title <- switch(
    type,
    heatmap = "Heat map",
    contour = "Contour map",
    heatmap_contour = "Heat map + contours"
  )
  plot <- plot +
    ggplot2::coord_sf(datum = NA, expand = FALSE) +
    ggplot2::labs(
      title = paste(surface$extracted$label, map_title),
      subtitle = paste(
        paste0(toupper(substr(interpolation, 1L, 1L)),
               substring(interpolation, 2L)), "visual interpolation;",
        "not a new-location model prediction"
      ),
      caption = gwrs_heatmap_caption(
        surface$extracted, selection, significance
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  surface$grid$value <- surface$grid$value_1
  surface$settings$variable <- surface$extracted$variable
  surface$settings$value_type <- surface$extracted$value_type
  surface$settings$selection <- selection
  surface$settings$significance <- isTRUE(significance)
  surface$settings$alpha <- alpha
  surface$settings$p_adjust <- if (significance) p_adjust else NULL
  surface$settings$visual_interpolation_only <- TRUE
  attr(plot, "gwrs_bandwidth") <- surface$bandwidth
  attr(plot, "gwrs_heatmap_settings") <- surface$settings
  if (!isTRUE(return_data)) return(plot)
  list(
    plot = plot,
    grid = surface$grid,
    points = surface$points,
    boundary = surface$boundary,
    bandwidth = surface$bandwidth,
    settings = surface$settings,
    interpolation_points = surface$interpolation_points
  )
}
