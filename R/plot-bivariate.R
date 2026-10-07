gwrs_quantile_bin <- function(values, bins) {
  finite <- is.finite(values)
  result <- rep(NA_integer_, length(values))
  if (!any(finite)) return(result)
  ranks <- rank(values[finite], ties.method = "average")
  result[finite] <- pmin(
    bins,
    pmax(1L, as.integer(ceiling(ranks / sum(finite) * bins)))
  )
  result
}

gwrs_bivariate_palette <- function(bins) {
  corner <- grDevices::col2rgb(c(
    low = "#E8E8E8", variable1 = "#BE64AC",
    variable2 = "#5AC8C8", both = "#3B4994"
  )) / 255
  colors <- character(bins * bins)
  index <- 0L
  for (second in seq_len(bins)) {
    y <- if (bins == 1L) 0 else (second - 1) / (bins - 1)
    for (first in seq_len(bins)) {
      x <- if (bins == 1L) 0 else (first - 1) / (bins - 1)
      rgb <- (1 - x) * (1 - y) * corner[, "low"] +
        x * (1 - y) * corner[, "variable1"] +
        (1 - x) * y * corner[, "variable2"] +
        x * y * corner[, "both"]
      index <- index + 1L
      colors[index] <- grDevices::rgb(rgb[1L], rgb[2L], rgb[3L])
    }
  }
  colors
}

gwrs_significance_extraction <- function(object,
                                         extracted,
                                         alpha,
                                         p_adjust,
                                         estimator) {
  enriched <- add_gwrs_spatial_significance(
    extracted, object, alpha, p_adjust, estimator
  )
  result <- enriched
  result$points$value <- as.double(enriched$points$significant)
  result$points$selected <- TRUE
  result$variable <- "significance"
  result$label <- paste("Significance of", extracted$label)
  result$value_type <- "significance"
  result$coefficient <- FALSE
  result
}

#' Bivariate continuous spatial maps for gwrs results
#'
#' Interpolates two observed local-result surfaces onto a common grid and
#' combines their quantile classes. This is a visual summary, not a formal
#' bivariate prediction at new locations.
#'
#' @inheritParams plot_heatmap
#' @param variable1,variable2 Spatial variables accepted by [plot_heatmap()].
#'   `variable2 = "significance"` visualizes the adjusted significance support
#'   of a coefficient supplied as `variable1`.
#' @param bins Number of quantile classes per variable. The default produces a
#'   3-by-3 legend.
#'
#' @return A `ggplot` object, or a list with plot and interpolation data when
#'   `return_data = TRUE`.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("sf", quietly = TRUE)) {
#'   set.seed(7)
#'   coords <- cbind(runif(36), runif(36))
#'   x <- cbind(income = rnorm(36), education = rnorm(36))
#'   y <- 0.4 + (1.5 * coords[, 1] - 0.5) * x[, "income"] +
#'     0.5 * x[, "education"] + rnorm(36, sd = 0.2)
#'   fit <- gwr_fit(
#'     x, y, coords = coords, k = 24,
#'     control = gwrs_control(n_threads = 1)
#'   )
#'   plot_bivariate(
#'     fit, "income", "local_r2", bins = 3, resolution = 30
#'   )
#' }
#' @export
plot_bivariate <- function(object,
                           variable1,
                           variable2,
                           bins = 3L,
                           boundary = NULL,
                           interpolation = c("kernel", "idw"),
                           distance = c("euclidean", "geographic", "hyperbolic"),
                           kernel = c("gaussian", "bisquare", "exponential"),
                           bandwidth = NULL,
                           resolution = 300L,
                           mask = TRUE,
                           alpha = 0.05,
                           selection = c("gray", "show", "hide"),
                           selection_threshold = 0.5,
                           unselected_color = "grey80",
                           show_points = FALSE,
                           show_boundaries = TRUE,
                           return_data = FALSE,
                           crs = NULL,
                           idw_power = 2,
                           max_neighbors = 200L,
                           chunk_size = 2000L,
                           extrapolation = NULL,
                           p_adjust = "BH",
                           estimator = c("auto", "post_selection", "penalized"),
                           na.rm = TRUE,
                           ...) {
  require_gwrs_visual_packages()
  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) ||
      bins != floor(bins) || bins < 2L || bins > 5L) {
    stop("`bins` must be an integer between 2 and 5.", call. = FALSE)
  }
  bins <- as.integer(bins)
  interpolation <- match.arg(interpolation)
  distance <- match.arg(distance)
  kernel <- match.arg(kernel)
  selection <- match.arg(selection)
  estimator <- match.arg(estimator)
  gwrs_validate_heatmap_flags(
    mask, FALSE, show_points, show_boundaries, return_data,
    alpha, 10L, 0.65
  )
  first <- extract_gwrs_spatial_values(object, variable1, na.rm = na.rm)
  significance_second <- is.character(variable2) && length(variable2) == 1L &&
    tolower(trimws(variable2)) %in% c("significance", "significant")
  second <- if (significance_second) {
    gwrs_significance_extraction(
      object, first, alpha, p_adjust, estimator
    )
  } else {
    extract_gwrs_spatial_values(object, variable2, na.rm = na.rm)
  }
  surface <- prepare_gwrs_spatial_surface(
    object = object,
    extractions = list(first, second),
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
  grid <- surface$plot_grid
  grid$selected <- is.finite(grid$selected_1) &
    grid$selected_1 >= selection_threshold &
    is.finite(grid$selected_2) & grid$selected_2 >= selection_threshold
  grid$bin_1 <- gwrs_quantile_bin(grid$value_1, bins)
  grid$bin_2 <- if (significance_second) {
    ifelse(is.finite(grid$value_2) & grid$value_2 >= selection_threshold,
           bins, 1L)
  } else {
    gwrs_quantile_bin(grid$value_2, bins)
  }
  code <- (grid$bin_2 - 1L) * bins + grid$bin_1
  if (selection != "show") code[!grid$selected] <- NA_integer_
  levels <- seq_len(bins * bins)
  labels <- unlist(lapply(seq_len(bins), function(second_bin) {
    paste0("V1:", seq_len(bins), " / V2:", second_bin)
  }))
  grid$bivariate_class <- factor(code, levels = levels, labels = labels)
  colors <- gwrs_bivariate_palette(bins)
  names(colors) <- labels

  plot <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = grid,
      mapping = ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        fill = .data[["bivariate_class"]]
      ),
      na.rm = TRUE, show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = colors, drop = FALSE, na.value = NA, na.translate = FALSE,
      name = paste0(
        second$label, " (rows)\n", first$label, " (columns)"
      ),
      guide = ggplot2::guide_legend(
        nrow = bins, byrow = TRUE, title.position = "top"
      )
    )
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
  if (isTRUE(show_boundaries)) {
    boundary_color <- if (isTRUE(surface$settings$boundary_provided)) {
      "white"
    } else {
      "grey35"
    }
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
  plot <- plot +
    ggplot2::coord_sf(datum = NA, expand = FALSE) +
    ggplot2::labs(
      title = paste(first$label, "and", second$label),
      subtitle = "Bivariate visual interpolation; not a new-location model prediction",
      caption = if (selection == "gray" &&
                    (any(!first$points$selected) ||
                     any(!second$points$selected))) {
        "Grey: locally unselected in at least one coefficient surface"
      } else {
        NULL
      },
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  active_match <- match(surface$grid$grid_id, grid$grid_id)
  surface$grid$bin_1 <- grid$bin_1[active_match]
  surface$grid$bin_2 <- grid$bin_2[active_match]
  surface$grid$bivariate_class <- grid$bivariate_class[active_match]
  surface$grid$selected <- grid$selected[active_match]
  surface$settings$variables <- c(first$variable, second$variable)
  surface$settings$bins <- bins
  surface$settings$selection <- selection
  surface$settings$significance_second <- significance_second
  surface$settings$visual_interpolation_only <- TRUE
  attr(plot, "gwrs_bandwidth") <- surface$bandwidth
  attr(plot, "gwrs_bivariate_settings") <- surface$settings
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
