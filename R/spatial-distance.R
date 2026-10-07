require_gwrs_visual_packages <- function() {
  missing <- c(
    if (!requireNamespace("ggplot2", quietly = TRUE)) "ggplot2",
    if (!requireNamespace("sf", quietly = TRUE)) "sf"
  )
  if (length(missing)) {
    stop(
      sprintf(
        "Spatial heat maps require suggested package%s: %s.",
        if (length(missing) > 1L) "s" else "",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

gwrs_resolve_crs <- function(object, crs = NULL) {
  stored <- object$neighbors$query_crs %||%
    object$neighbors$crs %||%
    attr(object$neighbors$query_coords, "crs") %||%
    attr(object$neighbors$coords, "crs")
  candidate <- crs %||% stored
  if (is.null(candidate)) return(sf::st_crs(NA))
  resolved <- tryCatch(sf::st_crs(candidate), error = identity)
  if (inherits(resolved, "error") || is.na(resolved)) {
    stop("`crs` is not a valid coordinate reference system.", call. = FALSE)
  }
  resolved
}

gwrs_crs_missing <- function(crs) {
  is.null(crs) || isTRUE(is.na(crs))
}

gwrs_distance_coordinates <- function(points, grid, crs, distance) {
  if (distance == "hyperbolic") {
    stop(
      paste(
        "Hyperbolic interpolation is unavailable because `gwrs_fit` does",
        "not provide a valid mapping from new geographic grid locations to",
        "a fitted hyperbolic/Lorentz representation."
      ),
      call. = FALSE
    )
  }
  source <- as.matrix(points[, c("x", "y"), drop = FALSE])
  query <- as.matrix(grid[, c("x", "y"), drop = FALSE])
  storage.mode(source) <- "double"
  storage.mode(query) <- "double"
  if (distance == "euclidean") {
    if (!gwrs_crs_missing(crs) && isTRUE(sf::st_is_longlat(crs))) {
      warning(
        paste(
          "Euclidean interpolation was requested for longitude/latitude",
          "coordinates; use `distance = \"geographic\"` for earth-surface",
          "distances."
        ),
        call. = FALSE
      )
    }
    return(list(source = source, query = query, units = "coordinate units"))
  }
  if (gwrs_crs_missing(crs)) {
    stop(
      paste(
        "Geographic distance requires a known CRS. Supply `crs = 4326`",
        "for longitude/latitude coordinates or the correct source CRS."
      ),
      call. = FALSE
    )
  }
  source_sf <- sf::st_as_sf(
    data.frame(x = source[, 1L], y = source[, 2L]),
    coords = c("x", "y"), crs = crs
  )
  query_sf <- sf::st_as_sf(
    data.frame(x = query[, 1L], y = query[, 2L]),
    coords = c("x", "y"), crs = crs
  )
  source_ll <- sf::st_coordinates(sf::st_transform(source_sf, 4326))[, 1:2,
                                                                          drop = FALSE]
  query_ll <- sf::st_coordinates(sf::st_transform(query_sf, 4326))[, 1:2,
                                                                        drop = FALSE]
  list(source = source_ll, query = query_ll, units = "metres")
}

gwrs_lonlat_unit_sphere <- function(coordinates) {
  longitude <- coordinates[, 1L] * pi / 180
  latitude <- coordinates[, 2L] * pi / 180
  cbind(
    cos(latitude) * cos(longitude),
    cos(latitude) * sin(longitude),
    sin(latitude)
  )
}

gwrs_knn_query <- function(source,
                           query,
                           k,
                           distance = c("euclidean", "geographic")) {
  distance <- match.arg(distance)
  k <- min(as.integer(k), nrow(source))
  if (k < 1L) stop("At least one interpolation neighbor is required.",
                   call. = FALSE)
  if (distance == "euclidean") {
    nearest <- FNN::get.knnx(source, query, k = k, algorithm = "kd_tree")
    return(list(index = nearest$nn.index, distance = nearest$nn.dist))
  }
  source_sphere <- gwrs_lonlat_unit_sphere(source)
  query_sphere <- gwrs_lonlat_unit_sphere(query)
  nearest <- FNN::get.knnx(
    source_sphere, query_sphere, k = k, algorithm = "kd_tree"
  )
  chord <- nearest$nn.dist
  chord[] <- pmin(2, pmax(0, chord))
  earth_radius <- 6371008.8
  list(
    index = nearest$nn.index,
    distance = 2 * earth_radius * asin(chord / 2)
  )
}

gwrs_observation_neighbor_distances <- function(source,
                                                k,
                                                distance) {
  n <- nrow(source)
  k <- min(as.integer(k), n - 1L)
  if (k < 1L) return(matrix(numeric(), nrow = n, ncol = 0L))
  if (distance == "euclidean") {
    return(FNN::get.knn(source, k = k, algorithm = "kd_tree")$nn.dist)
  }
  sphere <- gwrs_lonlat_unit_sphere(source)
  chord <- FNN::get.knn(sphere, k = k, algorithm = "kd_tree")$nn.dist
  chord[] <- pmin(2, pmax(0, chord))
  2 * 6371008.8 * asin(chord / 2)
}

gwrs_auto_bandwidth <- function(source, distance) {
  n <- nrow(source)
  bandwidth_neighbors <- min(n - 1L, max(2L, min(20L, ceiling(sqrt(n)))))
  distances <- gwrs_observation_neighbor_distances(
    source, bandwidth_neighbors, distance
  )
  bandwidth <- stats::median(
    distances[, ncol(distances)], na.rm = TRUE
  )
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    stop("Automatic bandwidth selection failed for the supplied locations.",
         call. = FALSE)
  }
  list(
    bandwidth = as.double(bandwidth),
    neighbors = bandwidth_neighbors,
    method = sprintf(
      "median distance to the %s nearest distinct location",
      bandwidth_neighbors
    )
  )
}

gwrs_typical_neighbor_distance <- function(source, distance) {
  distances <- gwrs_observation_neighbor_distances(source, 1L, distance)
  stats::median(distances[, 1L], na.rm = TRUE)
}

gwrs_spatial_diameter <- function(source, distance) {
  n <- nrow(source)
  sample_size <- min(n, 1000L)
  index <- unique(as.integer(round(seq(1, n, length.out = sample_size))))
  sampled <- source[index, , drop = FALSE]
  if (distance == "euclidean") {
    return(max(stats::dist(sampled)))
  }
  sphere <- gwrs_lonlat_unit_sphere(sampled)
  cosine <- tcrossprod(sphere)
  cosine[] <- pmin(1, pmax(-1, cosine))
  6371008.8 * max(acos(cosine))
}

validate_gwrs_bandwidth <- function(bandwidth,
                                    source,
                                    distance) {
  automatic <- is.null(bandwidth)
  selection <- if (automatic) {
    gwrs_auto_bandwidth(source, distance)
  } else {
    if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
        !is.finite(bandwidth) || bandwidth <= 0) {
      stop("`bandwidth` must be one positive finite number.", call. = FALSE)
    }
    list(
      bandwidth = as.double(bandwidth), neighbors = NA_integer_,
      method = "user supplied"
    )
  }
  typical <- gwrs_typical_neighbor_distance(source, distance)
  diameter <- gwrs_spatial_diameter(source, distance)
  if (is.finite(typical) && selection$bandwidth < 0.5 * typical) {
    warning(
      "The selected bandwidth is less than half the typical nearest-neighbor distance; gaps are likely.",
      call. = FALSE
    )
  }
  if (is.finite(diameter) && selection$bandwidth > diameter) {
    warning(
      "The selected bandwidth exceeds the approximate spatial diameter; the surface may be over-smoothed.",
      call. = FALSE
    )
  }
  selection$typical_neighbor_distance <- typical
  selection$approximate_diameter <- diameter
  selection$automatic <- automatic
  selection
}
