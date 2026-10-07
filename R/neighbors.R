#' Build a reusable nearest-neighbor structure
#'
#' `gwr_neighbors()` computes nearest neighbors once and stores them in the
#' target-contiguous layout used by the C++ fitting kernels. Reusing this object
#' avoids repeating spatial searches across models and tuning parameters.
#'
#' @param coords Numeric coordinate matrix with observations in rows, or an
#'   `sf` point object. Optional matrix `"crs"` metadata is retained.
#' @param query_coords Optional matrix or `sf` point coordinates at new target
#'   locations. When
#'   supplied, neighbors are selected from `coords` for every row of
#'   `query_coords`.
#' @param k Number of retained neighbors, including the target itself when
#'   `include_self = TRUE`.
#' @param include_self Whether each observation is inserted as its own first
#'   neighbor.
#' @param algorithm Nearest-neighbor algorithm passed to [FNN::get.knn()].
#'
#' @return An object of class `gwrs_neighbors`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#'
#' neighbors <- gwr_neighbors(coords, k = 18)
#' neighbors
#' head(neighbors$index[, 1:3])
#'
#' query <- coords[1:3, , drop = FALSE] + 0.01
#' cross_neighbors <- gwr_neighbors(
#'   coords, query_coords = query, k = 12, include_self = FALSE
#' )
#' dim(cross_neighbors$index)
#'
#' @export
gwr_neighbors <- function(coords,
                          query_coords = NULL,
                          k,
                          include_self = is.null(query_coords),
                          algorithm = c("kd_tree", "cover_tree", "CR", "brute")) {
  if (missing(include_self)) {
    include_self <- is.null(query_coords)
  } else {
    include_self <- isTRUE(include_self)
  }
  coords_crs <- if (inherits(coords, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package `sf` is required for sf coordinate inputs.", call. = FALSE)
    }
    sf::st_crs(coords)
  } else {
    attr(coords, "crs")
  }
  query_crs <- if (inherits(query_coords, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package `sf` is required for sf coordinate inputs.", call. = FALSE)
    }
    sf::st_crs(query_coords)
  } else {
    attr(query_coords, "crs")
  }
  coords <- as_coordinate_matrix(coords)
  n_train <- nrow(coords)
  same_locations <- is.null(query_coords)
  if (same_locations) {
    query_coords <- coords
    query_crs <- coords_crs
  } else {
    query_coords <- as_coordinate_matrix(query_coords)
    if (ncol(query_coords) != ncol(coords)) {
      stop("Training and target coordinates must have the same dimension.",
           call. = FALSE)
    }
  }
  n_target <- nrow(query_coords)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L || k > n_train) {
    stop("`k` must be between 1 and the number of training rows.",
         call. = FALSE)
  }
  algorithm <- match.arg(algorithm)
  if (!same_locations && include_self) {
    stop("`include_self` is only available when targets equal training data.",
         call. = FALSE)
  }

  external_k <- if (include_self) k - 1L else k
  if (same_locations && external_k > n_train - 1L) {
    stop("At most n - 1 external neighbors can be requested.", call. = FALSE)
  }

  if (external_k == 0L) {
    index <- matrix(seq_len(n_train), nrow = 1L)
    distance <- matrix(0, nrow = 1L, ncol = n_train)
  } else if (same_locations) {
    nearest <- FNN::get.knn(coords, k = external_k, algorithm = algorithm)
    if (include_self) {
      index <- rbind(seq_len(n_train), t(nearest$nn.index))
      distance <- rbind(rep.int(0, n_train), t(nearest$nn.dist))
    } else {
      index <- t(nearest$nn.index)
      distance <- t(nearest$nn.dist)
    }
  } else {
    nearest <- FNN::get.knnx(coords, query_coords, k = external_k,
                            algorithm = algorithm)
    index <- t(nearest$nn.index)
    distance <- t(nearest$nn.dist)
  }

  storage.mode(index) <- "integer"
  storage.mode(distance) <- "double"
  structure(
    list(
      index = index,
      distance = distance,
      coords = coords,
      query_coords = query_coords,
      crs = coords_crs,
      query_crs = query_crs,
      n = n_train,
      n_train = n_train,
      n_target = n_target,
      k = nrow(index),
      include_self = include_self,
      algorithm = algorithm,
      metric = "euclidean"
    ),
    class = "gwrs_neighbors"
  )
}

#' @export
print.gwrs_neighbors <- function(x, ...) {
  cat("Reusable gwrs neighbor structure\n")
  cat("  training observations:", x$n_train %||% x$n, "\n")
  cat("  target locations:", x$n_target %||% x$n, "\n")
  cat("  retained neighbors:", x$k, "\n")
  cat("  algorithm:", x$algorithm, "\n")
  invisible(x)
}

validate_neighbors <- function(neighbors, n) {
  if (!inherits(neighbors, "gwrs_neighbors")) {
    stop("`neighbors` must be created with `gwr_neighbors()`.", call. = FALSE)
  }
  if (!identical(as.integer(neighbors$n_train %||% neighbors$n), as.integer(n)) ||
      !identical(as.integer(neighbors$n_target %||% neighbors$n), as.integer(n))) {
    stop("The neighbor structure and model data have different row counts.",
         call. = FALSE)
  }
  if (!is.matrix(neighbors$index) || !is.matrix(neighbors$distance) ||
      !identical(dim(neighbors$index), dim(neighbors$distance))) {
    stop("The neighbor structure is malformed.", call. = FALSE)
  }
  neighbors
}

validate_cross_neighbors <- function(neighbors, n_train, n_target) {
  if (!inherits(neighbors, "gwrs_neighbors")) {
    stop("`neighbors` must be created with `gwr_neighbors()`.", call. = FALSE)
  }
  observed_train <- as.integer(neighbors$n_train %||% neighbors$n)
  observed_target <- as.integer(neighbors$n_target %||% neighbors$n)
  if (!identical(observed_train, as.integer(n_train)) ||
      !identical(observed_target, as.integer(n_target))) {
    stop("The neighbor structure does not match the training and target data.",
         call. = FALSE)
  }
  if (!is.matrix(neighbors$index) || !is.matrix(neighbors$distance) ||
      !identical(dim(neighbors$index), dim(neighbors$distance))) {
    stop("The neighbor structure is malformed.", call. = FALSE)
  }
  neighbors
}
