gwrs_polygon_from_coordinates <- function(coordinates, crs) {
  coordinates <- as.matrix(coordinates)
  storage.mode(coordinates) <- "double"
  if (ncol(coordinates) < 2L || nrow(coordinates) < 3L ||
      any(!is.finite(coordinates[, 1:2, drop = FALSE]))) {
    stop("A numeric boundary must contain at least three finite x/y rows.",
         call. = FALSE)
  }
  coordinates <- coordinates[, 1:2, drop = FALSE]
  if (!all(coordinates[1L, ] == coordinates[nrow(coordinates), ])) {
    coordinates <- rbind(coordinates, coordinates[1L, ])
  }
  sf::st_sfc(sf::st_polygon(list(coordinates)), crs = crs)
}

gwrs_as_boundary <- function(boundary, source_crs) {
  if (inherits(boundary, "sf")) {
    geometry <- sf::st_geometry(boundary)
  } else if (inherits(boundary, "sfc")) {
    geometry <- boundary
  } else if (is.matrix(boundary) || is.data.frame(boundary)) {
    if (ncol(boundary) < 2L ||
        !all(vapply(boundary[, 1:2, drop = FALSE], is.numeric, logical(1L)))) {
      converted <- tryCatch(sf::st_as_sf(boundary), error = identity)
      if (inherits(converted, "error")) {
        stop("`boundary` could not be converted to an sf polygon.",
             call. = FALSE)
      }
      geometry <- sf::st_geometry(converted)
    } else {
      geometry <- gwrs_polygon_from_coordinates(boundary, source_crs)
    }
  } else {
    converted <- tryCatch(sf::st_as_sf(boundary), error = identity)
    if (inherits(converted, "error")) {
      stop("`boundary` must be an sf polygon or a convertible spatial object.",
           call. = FALSE)
    }
    geometry <- sf::st_geometry(converted)
  }

  geometry_type <- as.character(sf::st_geometry_type(geometry))
  if (any(!geometry_type %in% c("POLYGON", "MULTIPOLYGON"))) {
    geometry <- suppressWarnings(sf::st_collection_extract(
      sf::st_make_valid(geometry), "POLYGON"
    ))
  }
  if (!length(geometry)) {
    stop("`boundary` does not contain polygon geometry.", call. = FALSE)
  }
  boundary_crs <- sf::st_crs(geometry)
  if (gwrs_crs_missing(source_crs) && !gwrs_crs_missing(boundary_crs)) {
    stop(
      paste(
        "The boundary has a CRS but the model coordinates do not.",
        "Supply the model-coordinate CRS through `crs`."
      ),
      call. = FALSE
    )
  }
  if (!gwrs_crs_missing(source_crs) && gwrs_crs_missing(boundary_crs)) {
    warning(
      "The boundary has no CRS; it is being assigned the model-coordinate CRS.",
      call. = FALSE
    )
    geometry <- sf::st_set_crs(geometry, source_crs)
  } else if (!gwrs_crs_missing(source_crs) &&
             !gwrs_crs_missing(boundary_crs) &&
             !isTRUE(source_crs == boundary_crs)) {
    geometry <- sf::st_transform(geometry, source_crs)
  } else if (gwrs_crs_missing(source_crs) && gwrs_crs_missing(boundary_crs)) {
    warning(
      "Model coordinates and boundary have no CRS; masking uses their shared numeric coordinate units.",
      call. = FALSE
    )
  }
  geometry <- sf::st_make_valid(geometry)
  geometry_type <- as.character(sf::st_geometry_type(geometry))
  if (any(!geometry_type %in% c("POLYGON", "MULTIPOLYGON"))) {
    geometry <- suppressWarnings(sf::st_collection_extract(
      geometry, "POLYGON"
    ))
  }
  geometry
}

gwrs_bbox_polygon <- function(coordinates, crs, expansion = 0.02) {
  x_range <- range(coordinates[, 1L])
  y_range <- range(coordinates[, 2L])
  x_width <- diff(x_range)
  y_width <- diff(y_range)
  scale <- max(x_width, y_width)
  if (!is.finite(scale) || scale <= 0) {
    stop("Spatial coordinates have no usable two-dimensional extent.",
         call. = FALSE)
  }
  x_padding <- if (x_width > 0) expansion * x_width else expansion * scale
  y_padding <- if (y_width > 0) expansion * y_width else expansion * scale
  x_range <- x_range + c(-x_padding, x_padding)
  y_range <- y_range + c(-y_padding, y_padding)
  ring <- rbind(
    c(x_range[1L], y_range[1L]), c(x_range[2L], y_range[1L]),
    c(x_range[2L], y_range[2L]), c(x_range[1L], y_range[2L]),
    c(x_range[1L], y_range[1L])
  )
  sf::st_sfc(sf::st_polygon(list(ring)), crs = crs)
}

gwrs_default_boundary <- function(points, crs, expansion = 0.02) {
  coordinates <- as.matrix(points[, c("x", "y"), drop = FALSE])
  hull_index <- unique(grDevices::chull(coordinates))
  if (length(hull_index) < 3L) {
    return(gwrs_bbox_polygon(coordinates, crs, expansion))
  }
  hull <- coordinates[hull_index, , drop = FALSE]
  center <- colMeans(hull)
  hull <- sweep(hull, 2L, center, FUN = "-") * (1 + expansion)
  hull <- sweep(hull, 2L, center, FUN = "+")
  hull <- rbind(hull, hull[1L, ])
  polygon <- sf::st_sfc(sf::st_polygon(list(hull)), crs = crs)
  if (isTRUE(any(sf::st_is_empty(polygon))) ||
      !isTRUE(all(sf::st_is_valid(polygon)))) {
    return(gwrs_bbox_polygon(coordinates, crs, expansion))
  }
  polygon
}

validate_gwrs_resolution <- function(resolution) {
  if (!is.numeric(resolution) || length(resolution) != 1L ||
      !is.finite(resolution) || resolution != floor(resolution) ||
      resolution < 25L || resolution > 1000L) {
    stop("`resolution` must be an integer between 25 and 1000.",
         call. = FALSE)
  }
  resolution <- as.integer(resolution)
  if (resolution > 500L) {
    warning(
      "A resolution above 500 can create more than 250,000 grid cells.",
      call. = FALSE
    )
  }
  resolution
}

gwrs_regular_grid <- function(boundary,
                              resolution,
                              mask = TRUE) {
  resolution <- validate_gwrs_resolution(resolution)
  bounds <- sf::st_bbox(boundary)
  width <- as.double(bounds[["xmax"]] - bounds[["xmin"]])
  height <- as.double(bounds[["ymax"]] - bounds[["ymin"]])
  if (!is.finite(width) || !is.finite(height) || width <= 0 || height <= 0) {
    stop("The plotting boundary has an invalid spatial extent.",
         call. = FALSE)
  }
  if (width >= height) {
    nx <- resolution
    ny <- max(3L, as.integer(round(resolution * height / width)))
  } else {
    ny <- resolution
    nx <- max(3L, as.integer(round(resolution * width / height)))
  }
  dx <- width / nx
  dy <- height / ny
  x <- as.double(bounds[["xmin"]]) + (seq_len(nx) - 0.5) * dx
  y <- as.double(bounds[["ymin"]]) + (seq_len(ny) - 0.5) * dy
  grid <- expand.grid(x = x, y = y, KEEP.OUT.ATTRS = FALSE)
  grid$grid_id <- seq_len(nrow(grid))
  grid$inside_boundary <- TRUE
  if (isTRUE(mask)) {
    grid_sf <- sf::st_as_sf(
      grid, coords = c("x", "y"), crs = sf::st_crs(boundary),
      remove = FALSE
    )
    grid$inside_boundary <- lengths(
      sf::st_intersects(grid_sf, boundary)
    ) > 0L
  }
  list(
    grid = grid,
    nx = nx,
    ny = ny,
    dx = dx,
    dy = dy,
    resolution = resolution,
    n_total = nrow(grid),
    n_inside = sum(grid$inside_boundary)
  )
}

gwrs_prepare_spatial_domain <- function(points,
                                        object,
                                        boundary,
                                        mask,
                                        resolution,
                                        crs) {
  source_crs <- gwrs_resolve_crs(object, crs)
  boundary_provided <- !is.null(boundary)
  plotting_boundary <- if (boundary_provided) {
    gwrs_as_boundary(boundary, source_crs)
  } else {
    gwrs_default_boundary(points, source_crs)
  }
  grid <- gwrs_regular_grid(plotting_boundary, resolution, mask = mask)
  list(
    boundary = plotting_boundary,
    boundary_provided = boundary_provided,
    crs = source_crs,
    grid = grid
  )
}

gwrs_boundary_sf <- function(boundary) {
  sf::st_sf(geometry = boundary)
}
