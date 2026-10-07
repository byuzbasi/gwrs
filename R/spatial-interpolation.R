gwrs_source_location_key <- function(x, y) {
  paste(sprintf("%.17g", x), sprintf("%.17g", y), sep = "\r")
}

combine_gwrs_spatial_extractions <- function(extractions) {
  if (!length(extractions)) stop("No spatial variables were supplied.",
                                 call. = FALSE)
  common_id <- Reduce(
    intersect,
    lapply(extractions, function(value) value$points$observation_id)
  )
  if (length(common_id) < 3L) {
    stop("Fewer than three observations have finite values for all variables.",
         call. = FALSE)
  }
  common_id <- sort(common_id)
  first_match <- match(common_id, extractions[[1L]]$points$observation_id)
  combined <- data.frame(
    x = extractions[[1L]]$points$x[first_match],
    y = extractions[[1L]]$points$y[first_match],
    observation_id = common_id,
    stringsAsFactors = FALSE
  )
  for (index in seq_along(extractions)) {
    point <- extractions[[index]]$points
    point_match <- match(common_id, point$observation_id)
    if (anyNA(point_match)) {
      stop("Spatial variable observations could not be aligned.",
           call. = FALSE)
    }
    if (any(point$x[point_match] != combined$x |
            point$y[point_match] != combined$y)) {
      stop("Spatial variables use inconsistent coordinates.", call. = FALSE)
    }
    combined[[paste0("value_", index)]] <- point$value[point_match]
    combined[[paste0("selected_", index)]] <-
      as.double(point$selected[point_match])
    if ("significant" %in% names(point)) {
      combined[[paste0("significant_", index)]] <-
        as.double(point$significant[point_match])
    }
  }
  combined
}

aggregate_duplicate_gwrs_locations <- function(points) {
  key <- gwrs_source_location_key(points$x, points$y)
  unique_key <- unique(key)
  group <- match(key, unique_key)
  if (length(unique_key) == nrow(points)) {
    points$duplicate_count <- 1L
    return(points)
  }
  value_columns <- grep(
    "^(value|selected|significant)_[0-9]+$", names(points), value = TRUE
  )
  result <- data.frame(
    x = vapply(split(points$x, group), function(value) value[1L], numeric(1L)),
    y = vapply(split(points$y, group), function(value) value[1L], numeric(1L)),
    observation_id = vapply(
      split(points$observation_id, group),
      function(value) paste(value, collapse = ","),
      character(1L)
    ),
    duplicate_count = as.integer(tabulate(group)),
    stringsAsFactors = FALSE
  )
  for (column in value_columns) {
    result[[column]] <- vapply(
      split(points[[column]], group),
      mean, numeric(1L), na.rm = TRUE
    )
  }
  if (nrow(result) < 3L) {
    stop("At least three unique locations are required after aggregating duplicates.",
         call. = FALSE)
  }
  result
}

validate_gwrs_interpolation_controls <- function(max_neighbors,
                                                 chunk_size,
                                                 idw_power,
                                                 selection_threshold,
                                                 extrapolation,
                                                 mask) {
  if (!is.numeric(max_neighbors) || length(max_neighbors) != 1L ||
      is.na(max_neighbors) || max_neighbors <= 0) {
    stop("`max_neighbors` must be a positive integer or Inf.",
         call. = FALSE)
  }
  if (is.infinite(max_neighbors)) {
    max_neighbors <- Inf
  } else {
    if (max_neighbors != floor(max_neighbors)) {
      stop("`max_neighbors` must be a positive integer or Inf.",
           call. = FALSE)
    }
    max_neighbors <- as.integer(max_neighbors)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 10L ||
      chunk_size != floor(chunk_size)) {
    stop("`chunk_size` must be an integer of at least 10.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)
  if (!is.numeric(idw_power) || length(idw_power) != 1L ||
      !is.finite(idw_power) || idw_power <= 0) {
    stop("`idw_power` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(selection_threshold) ||
      length(selection_threshold) != 1L ||
      !is.finite(selection_threshold) || selection_threshold < 0 ||
      selection_threshold > 1) {
    stop("`selection_threshold` must be one number in [0, 1].",
         call. = FALSE)
  }
  if (is.null(extrapolation)) extrapolation <- if (isTRUE(mask)) 3 else Inf
  if (!is.numeric(extrapolation) || length(extrapolation) != 1L ||
      is.na(extrapolation) || extrapolation <= 0) {
    stop("`extrapolation` must be positive, Inf, or NULL.", call. = FALSE)
  }
  list(
    max_neighbors = max_neighbors,
    chunk_size = chunk_size,
    idw_power = as.double(idw_power),
    selection_threshold = as.double(selection_threshold),
    extrapolation = as.double(extrapolation)
  )
}

gwrs_kernel_interpolation_weights <- function(distance, kernel, bandwidth) {
  ratio <- distance / bandwidth
  switch(
    kernel,
    gaussian = exp(-0.5 * ratio^2),
    bisquare = ifelse(ratio < 1, (1 - ratio^2)^2, 0),
    exponential = exp(-ratio)
  )
}

gwrs_interpolate_value_matrix <- function(source,
                                          query,
                                          values,
                                          interpolation,
                                          distance,
                                          kernel,
                                          bandwidth,
                                          idw_power,
                                          max_neighbors,
                                          chunk_size) {
  n_source <- nrow(source)
  n_query <- nrow(query)
  k <- if (is.infinite(max_neighbors)) {
    n_source
  } else {
    min(n_source, as.integer(max_neighbors))
  }
  result <- matrix(NA_real_, nrow = n_query, ncol = ncol(values))
  nearest_distance <- rep(NA_real_, n_query)
  zero_weight <- logical(n_query)
  chunks <- split(seq_len(n_query), ceiling(seq_len(n_query) / chunk_size))
  for (rows in chunks) {
    nearest <- gwrs_knn_query(
      source, query[rows, , drop = FALSE], k = k, distance = distance
    )
    index <- as.matrix(nearest$index)
    local_distance <- as.matrix(nearest$distance)
    nearest_distance[rows] <- local_distance[, 1L]
    exact <- local_distance <= 1e-10
    has_exact <- rowSums(exact) > 0L
    weights <- if (interpolation == "idw") {
      1 / pmax(local_distance, .Machine$double.eps)^idw_power
    } else {
      gwrs_kernel_interpolation_weights(local_distance, kernel, bandwidth)
    }
    weights[exact] <- 0
    weight_sum <- rowSums(weights)
    zero_weight[rows] <- !has_exact &
      (!is.finite(weight_sum) | weight_sum <= .Machine$double.eps)
    for (column in seq_len(ncol(values))) {
      local_value <- matrix(
        values[index, column], nrow = nrow(index), ncol = ncol(index)
      )
      estimate <- rep(NA_real_, length(rows))
      if (any(has_exact)) {
        estimate[has_exact] <- rowSums(
          local_value[has_exact, , drop = FALSE] *
            exact[has_exact, , drop = FALSE]
        ) / rowSums(exact[has_exact, , drop = FALSE])
      }
      regular <- !has_exact & is.finite(weight_sum) &
        weight_sum > .Machine$double.eps
      if (any(regular)) {
        estimate[regular] <- rowSums(
          local_value[regular, , drop = FALSE] *
            weights[regular, , drop = FALSE]
        ) / weight_sum[regular]
      }
      result[rows, column] <- estimate
    }
  }
  colnames(result) <- colnames(values)
  list(
    values = result,
    nearest_distance = nearest_distance,
    zero_weight = zero_weight,
    neighbors_used = k,
    truncated = k < n_source
  )
}

prepare_gwrs_spatial_surface <- function(object,
                                         extractions,
                                         interpolation,
                                         distance,
                                         kernel,
                                         bandwidth,
                                         resolution,
                                         boundary,
                                         mask,
                                         crs,
                                         max_neighbors,
                                         chunk_size,
                                         idw_power,
                                         selection_threshold,
                                         extrapolation) {
  require_gwrs_visual_packages()
  interpolation <- match.arg(interpolation, c("kernel", "idw"))
  distance <- match.arg(distance, c("euclidean", "geographic", "hyperbolic"))
  kernel <- match.arg(kernel, c("gaussian", "bisquare", "exponential"))
  controls <- validate_gwrs_interpolation_controls(
    max_neighbors, chunk_size, idw_power, selection_threshold,
    extrapolation, mask
  )
  original_points <- combine_gwrs_spatial_extractions(extractions)
  source_points <- aggregate_duplicate_gwrs_locations(original_points)
  domain <- gwrs_prepare_spatial_domain(
    source_points, object, boundary, mask, resolution, crs
  )
  full_grid <- domain$grid$grid
  active <- which(full_grid$inside_boundary)
  if (length(active) < 3L) {
    stop("The boundary and resolution leave fewer than three grid cells.",
         call. = FALSE)
  }
  active_grid <- full_grid[active, , drop = FALSE]
  distance_coordinates <- gwrs_distance_coordinates(
    source_points, active_grid, domain$crs, distance
  )

  bandwidth_info <- NULL
  if (interpolation == "kernel") {
    bandwidth_info <- validate_gwrs_bandwidth(
      bandwidth, distance_coordinates$source, distance
    )
    selected_bandwidth <- bandwidth_info$bandwidth
  } else {
    if (!is.null(bandwidth)) {
      warning("`bandwidth` is ignored when `interpolation = \"idw\"`.",
              call. = FALSE)
    }
    selected_bandwidth <- NULL
  }

  value_columns <- grep(
    "^(value|selected|significant)_[0-9]+$",
    names(source_points), value = TRUE
  )
  source_values <- as.matrix(source_points[, value_columns, drop = FALSE])
  interpolated <- gwrs_interpolate_value_matrix(
    source = distance_coordinates$source,
    query = distance_coordinates$query,
    values = source_values,
    interpolation = interpolation,
    distance = distance,
    kernel = kernel,
    bandwidth = selected_bandwidth,
    idw_power = controls$idw_power,
    max_neighbors = controls$max_neighbors,
    chunk_size = controls$chunk_size
  )
  for (column in colnames(interpolated$values)) {
    active_grid[[column]] <- interpolated$values[, column]
  }
  active_grid$nearest_observation_distance <- interpolated$nearest_distance
  typical_distance <- gwrs_typical_neighbor_distance(
    distance_coordinates$source, distance
  )
  maximum_distance <- if (is.infinite(controls$extrapolation)) {
    Inf
  } else {
    controls$extrapolation * typical_distance
  }
  active_grid$within_extrapolation <- is.finite(
    active_grid$nearest_observation_distance
  ) & active_grid$nearest_observation_distance <= maximum_distance
  active_grid$interpolation_valid <- active_grid$within_extrapolation &
    !interpolated$zero_weight
  invalid <- !active_grid$interpolation_valid
  if (any(invalid)) {
    for (column in value_columns) active_grid[[column]][invalid] <- NA_real_
  }
  if (any(interpolated$zero_weight)) {
    warning(
      sprintf(
        "%s grid cells received zero interpolation weight and were left blank.",
        sum(interpolated$zero_weight)
      ),
      call. = FALSE
    )
  }

  for (column in c(value_columns, "nearest_observation_distance",
                   "within_extrapolation", "interpolation_valid")) {
    full_grid[[column]] <- NA
    full_grid[[column]][active] <- active_grid[[column]]
  }
  list(
    grid = active_grid,
    plot_grid = full_grid,
    points = original_points,
    interpolation_points = source_points,
    boundary = domain$boundary,
    bandwidth = selected_bandwidth,
    bandwidth_info = bandwidth_info,
    settings = list(
      interpolation = interpolation,
      distance = distance,
      distance_units = distance_coordinates$units,
      kernel = if (interpolation == "kernel") kernel else NULL,
      idw_power = if (interpolation == "idw") controls$idw_power else NULL,
      resolution = domain$grid$resolution,
      grid_dimensions = c(nx = domain$grid$nx, ny = domain$grid$ny),
      grid_cells = domain$grid$n_total,
      inside_cells = domain$grid$n_inside,
      mask = isTRUE(mask),
      extrapolation_factor = controls$extrapolation,
      extrapolation_distance = maximum_distance,
      max_neighbors = controls$max_neighbors,
      neighbors_used = interpolated$neighbors_used,
      truncated_neighbors = interpolated$truncated,
      chunk_size = controls$chunk_size,
      selection_threshold = controls$selection_threshold,
      crs = domain$crs,
      boundary_provided = domain$boundary_provided,
      duplicate_locations = sum(source_points$duplicate_count > 1L)
    )
  )
}
