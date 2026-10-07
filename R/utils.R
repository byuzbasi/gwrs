as_coordinate_matrix <- function(coords) {
  if (inherits(coords, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package `sf` is required for sf coordinate inputs.", call. = FALSE)
    }
    coords <- sf::st_coordinates(coords)
  }
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  if (nrow(coords) < 1L || ncol(coords) < 1L || any(!is.finite(coords))) {
    stop("`coords` must be a finite numeric matrix with observations in rows.",
         call. = FALSE)
  }
  coords
}

as_design_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) < 1L || ncol(x) < 1L || any(!is.finite(x))) {
    stop("`x` must be a finite numeric matrix with at least one column.",
         call. = FALSE)
  }
  x
}

as_response <- function(y, n) {
  y <- as.double(y)
  if (length(y) != n || any(!is.finite(y))) {
    stop("`y` must contain one finite value for every row of `x`.",
         call. = FALSE)
  }
  y
}

standardize_design <- function(x, standardize) {
  p <- ncol(x)
  if (!isTRUE(standardize)) {
    return(list(x = x, center = rep.int(0, p), scale = rep.int(1, p)))
  }
  center <- colMeans(x)
  centered <- sweep(x, 2L, center, FUN = "-")
  scale <- sqrt(colMeans(centered^2))
  bad <- !is.finite(scale) | scale <= 1e-12
  if (any(bad)) {
    stop(
      sprintf(
        "Constant or numerically constant predictors are not supported: %s.",
        paste(colnames(x)[bad] %||% which(bad), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  list(
    x = sweep(centered, 2L, scale, FUN = "/"),
    center = center,
    scale = scale
  )
}

apply_standardization <- function(x, center, scale) {
  sweep(sweep(x, 2L, center, FUN = "-"), 2L, scale, FUN = "/")
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

backtransform_coefficients <- function(coefficients, center, scale) {
  coefficients <- t(coefficients)
  slopes <- sweep(coefficients[, -1L, drop = FALSE], 2L, scale, FUN = "/")
  intercept <- coefficients[, 1L] - drop(slopes %*% center)
  out <- cbind(`(Intercept)` = intercept, slopes)
  out
}

backtransform_coefficient_path <- function(coefficients, center, scale) {
  dimensions <- dim(coefficients)
  if (length(dimensions) != 3L || dimensions[1L] != length(center) + 1L) {
    stop("Malformed coefficient path returned by the numerical core.",
         call. = FALSE)
  }
  result <- array(
    NA_real_,
    dim = c(dimensions[2L], dimensions[1L], dimensions[3L])
  )
  for (path_index in seq_len(dimensions[3L])) {
    result[, , path_index] <- backtransform_coefficients(
      coefficients[, , path_index, drop = FALSE][, , 1L],
      center,
      scale
    )
  }
  result
}

validate_penalty_parameters <- function(alpha, d = NULL) {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be one number in [0, 1].", call. = FALSE)
  }
  if (!is.null(d) &&
      (length(d) < 1L || any(!is.finite(d)) || any(d < 0 | d > 1))) {
    stop("Every value of `d` must be in [0, 1].", call. = FALSE)
  }
  invisible(TRUE)
}

resolve_bandwidth <- function(bandwidth) {
  adaptive <- is.null(bandwidth)
  fixed_bandwidth <- if (adaptive) 1 else as.double(bandwidth)
  if (!adaptive && (length(fixed_bandwidth) != 1L ||
                    !is.finite(fixed_bandwidth) || fixed_bandwidth <= 0)) {
    stop("`bandwidth` must be one positive finite distance.", call. = FALSE)
  }
  list(adaptive = adaptive, fixed = fixed_bandwidth)
}

validate_lambda_path <- function(lambda) {
  lambda <- as.double(lambda)
  if (!length(lambda) || any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("`lambda` must contain finite nonnegative values.", call. = FALSE)
  }
  if (is.unsorted(-lambda, strictly = FALSE)) {
    stop("`lambda` must be ordered from largest to smallest.", call. = FALSE)
  }
  lambda
}

make_lambda_sequence <- function(lambda_max, n_lambda, lambda_min_ratio) {
  n_lambda <- as.integer(n_lambda)
  if (length(n_lambda) != 1L || is.na(n_lambda) || n_lambda < 1L) {
    stop("`n_lambda` must be a positive integer.", call. = FALSE)
  }
  if (length(lambda_min_ratio) != 1L || !is.finite(lambda_min_ratio) ||
      lambda_min_ratio <= 0 || lambda_min_ratio > 1) {
    stop("`lambda_min_ratio` must be in (0, 1].", call. = FALSE)
  }
  if (!is.finite(lambda_max) || lambda_max < 0) {
    stop("The computed lambda maximum is invalid.", call. = FALSE)
  }
  if (lambda_max <= .Machine$double.eps || n_lambda == 1L) {
    return(as.double(lambda_max))
  }
  exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
          length.out = n_lambda))
}

resolve_neighbors <- function(x, coords, neighbors, k) {
  if (is.null(neighbors)) {
    if (is.null(coords)) {
      stop("Supply either `neighbors` or both `coords` and `k`.", call. = FALSE)
    }
    if (is.null(k)) stop("`k` is required when neighbors are constructed.",
                         call. = FALSE)
    neighbors <- gwr_neighbors(coords, k = k)
  }
  validate_neighbors(neighbors, nrow(x))
}

status_labels <- function(status) {
  labels <- c(
    "direct_solve",
    "pseudoinverse",
    "failed_solve",
    "zero_weights",
    "invalid_neighbor"
  )
  labels[pmin(pmax(as.integer(status), 0L), 4L) + 1L]
}

formula_model_components <- function(formula, data, coords) {
  model_frame <- stats::model.frame(
    formula, data = data, na.action = stats::na.fail
  )
  terms <- stats::terms(model_frame)
  if (attr(terms, "intercept") != 1L) {
    stop("The formula must include an intercept; it is locally unpenalized.",
         call. = FALSE)
  }
  y <- stats::model.response(model_frame)
  design <- stats::model.matrix(terms, model_frame)
  intercept <- match("(Intercept)", colnames(design), nomatch = 0L)
  if (intercept > 0L) design <- design[, -intercept, drop = FALSE]
  if (is.character(coords)) {
    if (length(coords) != 2L || !all(coords %in% names(data))) {
      stop("Character `coords` must name two columns in `data`.",
           call. = FALSE)
    }
    coords <- as.matrix(data[, coords, drop = FALSE])
  }
  list(
    x = design,
    y = y,
    coords = coords,
    terms = terms,
    contrasts = attr(design, "contrasts"),
    xlevels = stats::.getXlevels(terms, model_frame),
    model = model_frame
  )
}

decorate_formula_fit <- function(fit, components, call) {
  fit$call <- call
  fit$terms <- components$terms
  fit$contrasts <- components$contrasts
  fit$xlevels <- components$xlevels
  fit$model <- components$model
  fit
}
