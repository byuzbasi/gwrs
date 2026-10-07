#' Construct reproducible spatially contiguous folds
#'
#' Uses k-means clustering of standardized coordinates to form spatial blocks.
#' The spatial CV routines never substitute ordinary random folds when a
#' spatial fold assignment is absent.
#'
#' @param coords Numeric coordinate matrix.
#' @param n_folds Number of spatial folds.
#' @param seed Integer seed used by k-means initialization. The caller's random
#'   number state is restored before returning.
#' @param nstart Number of k-means starts.
#'
#' @return An integer fold vector with one value per coordinate row.
#'
#' @examples
#' coords <- as.matrix(datasets::quakes[1:30, c("long", "lat")])
#' fold <- spatial_folds(coords, n_folds = 3, seed = 42)
#' table(fold)
#'
#' @export
spatial_folds <- function(coords, n_folds = 5L, seed = 1L, nstart = 20L) {
  coords <- as_coordinate_matrix(coords)
  n_folds <- as.integer(n_folds)
  nstart <- as.integer(nstart)
  seed <- as.integer(seed)
  if (length(n_folds) != 1L || is.na(n_folds) || n_folds < 2L ||
      n_folds >= nrow(coords)) {
    stop("`n_folds` must be between 2 and n - 1.", call. = FALSE)
  }
  if (length(nstart) != 1L || is.na(nstart) || nstart < 1L ||
      length(seed) != 1L || is.na(seed)) {
    stop("`seed` and `nstart` must be valid scalar integers.", call. = FALSE)
  }
  scaled <- scale(coords)
  if (any(!is.finite(scaled))) {
    stop("Coordinate dimensions must have positive finite variation.",
         call. = FALSE)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  as.integer(stats::kmeans(scaled, centers = n_folds,
                           nstart = nstart)$cluster)
}

run_cross_gwr_sl_path <- function(x_train,
                                  y_train,
                                  x_target,
                                  neighbors,
                                  kernel,
                                  bandwidth,
                                  lambda,
                                  alpha,
                                  d,
                                  screening,
                                  keep_coefficients,
                                  control) {
  validate_cross_neighbors(neighbors, nrow(x_train), nrow(x_target))
  cpp_gwr_sl_path_predict(
    x_train,
    y_train,
    x_target,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    lambda,
    as.double(alpha),
    as.double(d),
    control$tolerance,
    control$max_iterations,
    isTRUE(screening),
    isTRUE(keep_coefficients),
    FALSE,
    control$n_threads,
    control$grain_size
  )
}

# Internal serial grid for low-dimensional CV. Columns are lambda within paired
# alpha/d blocks; no state or preprocessing is shared across CV training folds.
run_cross_gwr_sl_grid <- function(x_train, y_train, x_target, neighbors,
                                  kernel, bandwidth, lambda, alpha, d,
                                  screening, keep_coefficients, control,
                                  active_solve = FALSE) {
  validate_cross_neighbors(neighbors, nrow(x_train), nrow(x_target))
  cpp_gwr_sl_grid_predict(x_train, y_train, x_target,
    neighbors$index, neighbors$distance, kernel$code, bandwidth$adaptive,
    bandwidth$fixed, lambda, as.double(alpha), as.double(d),
    control$tolerance, control$max_iterations, isTRUE(screening),
    isTRUE(keep_coefficients), isTRUE(active_solve))
}

#' Spatial cross-validation for GWR-SL
#'
#' Evaluates a common decreasing lambda path using held-out spatial blocks.
#' Predictor centering and scaling are estimated separately inside every
#' training fold, preventing transformation leakage. Nearest-neighbor searches
#' are also restricted to each training fold.
#'
#' @inheritParams gwr_sl_path
#' @param fold User-supplied spatial fold labels, one per observation. Create a
#'   reproducible assignment with [spatial_folds()] or supply application-
#'   specific blocks/buffers.
#' @param d One or more center-location values in `[0, 1]`.
#' @param metric Criterion used for selecting the final pair: root mean squared
#'   error (`"rmse"`), mean squared error (`"mse"`), or mean absolute error
#'   (`"mae"`). All three are returned.
#' @param refit Whether to fit the selected model on all observations.
#'
#' @return An object of class `gwrs_cv` containing the error surface and,
#'   when `refit = TRUE`, the selected `gwrs_fit`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fold <- spatial_folds(coords, n_folds = 3, seed = 42)
#' cv <- cv_gwr_sl(
#'   x, quake$mag, coords, k = 8, fold = fold,
#'   lambda = c(0.12, 0.04), alpha = 0.7, d = c(0, 0.5),
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#' cv
#' cv$best
#'
#' @export
cv_gwr_sl <- function(x,
                      y,
                      coords,
                      k,
                      fold,
                      kernel = c("bisquare", "gaussian", "exponential",
                                 "tricube", "boxcar"),
                      bandwidth = NULL,
                      lambda = NULL,
                      n_lambda = 50L,
                      lambda_min_ratio = 1e-3,
                      alpha = 0.5,
                      d = c(0, 0.5, 1),
                      standardize = TRUE,
                      screening = TRUE,
                      metric = c("rmse", "mse", "mae"),
                      refit = TRUE,
                      control = gwrs_control()) {
  call <- match.call()
  x <- as_design_matrix(x)
  y <- as_response(y, nrow(x))
  coords <- as_coordinate_matrix(coords)
  if (nrow(coords) != nrow(x)) {
    stop("`coords` and `x` must have the same number of rows.", call. = FALSE)
  }
  control <- validate_control(control)
  kernel <- kernel_code(kernel)
  bandwidth <- resolve_bandwidth(bandwidth)
  metric <- match.arg(metric)
  validate_penalty_parameters(alpha, d)
  d <- as.double(d)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }
  if (missing(fold) || length(fold) != nrow(x) || anyNA(fold)) {
    stop("Supply one non-missing spatial `fold` label per observation.",
         call. = FALSE)
  }
  fold_factor <- factor(fold)
  if (nlevels(fold_factor) < 2L) {
    stop("At least two spatial folds are required.", call. = FALSE)
  }
  fold_id <- as.integer(fold_factor)
  fold_labels <- levels(fold_factor)
  training_sizes <- vapply(
    seq_len(nlevels(fold_factor)),
    function(index) sum(fold_id != index),
    integer(1L)
  )
  if (any(training_sizes < k)) {
    stop("`k` cannot exceed the training size of any fold.", call. = FALSE)
  }

  full_neighbors <- gwr_neighbors(coords, k = k)
  full_prepared <- standardize_design(x, standardize)
  lambda_max <- NULL
  if (is.null(lambda)) {
    if (alpha <= 0) {
      stop("Supply `lambda` explicitly when `alpha = 0`.", call. = FALSE)
    }
    lambda_max <- cpp_gwr_sl_lambda_max(
      full_prepared$x, y, full_neighbors$index, full_neighbors$distance,
      kernel$code, bandwidth$adaptive, bandwidth$fixed, as.double(alpha),
      control$n_threads, control$grain_size
    )
    lambda <- make_lambda_sequence(lambda_max, n_lambda, lambda_min_ratio)
  } else {
    lambda <- validate_lambda_path(lambda)
  }

  n_lambda <- length(lambda)
  n_d <- length(d)
  squared_error <- matrix(0, nrow = n_lambda, ncol = n_d)
  absolute_error <- matrix(0, nrow = n_lambda, ncol = n_d)
  observation_count <- matrix(0L, nrow = n_lambda, ncol = n_d)
  convergence_count <- matrix(0, nrow = n_lambda, ncol = n_d)
  fold_results <- vector("list", nlevels(fold_factor) * n_d)
  result_index <- 0L

  for (fold_index in seq_len(nlevels(fold_factor))) {
    test <- fold_id == fold_index
    train <- !test
    prepared <- standardize_design(x[train, , drop = FALSE], standardize)
    x_test <- apply_standardization(
      x[test, , drop = FALSE], prepared$center, prepared$scale
    )
    cross_neighbors <- gwr_neighbors(
      coords[train, , drop = FALSE],
      query_coords = coords[test, , drop = FALSE],
      k = k,
      include_self = FALSE
    )

    batch <- if (control$n_threads == 1L && ncol(x) <= 64L && alpha > 0)
      run_cross_gwr_sl_grid(prepared$x, y[train], x_test, cross_neighbors,
        kernel, bandwidth, lambda, rep(alpha, n_d), d, screening, FALSE, control)
      else NULL
    for (d_index in seq_along(d)) {
      columns <- (d_index - 1L) * n_lambda + seq_len(n_lambda)
      raw <- if (!is.null(batch)) list(
        predictions=batch$predictions[, columns, drop=FALSE],
        converged=batch$converged[, columns, drop=FALSE]) else run_cross_gwr_sl_path(
        prepared$x,
        y[train],
        x_test,
        cross_neighbors,
        kernel,
        bandwidth,
        lambda,
        alpha,
        d[d_index],
        screening,
        FALSE,
        control
      )
      if (any(!is.finite(raw$predictions))) {
        stop(
          sprintf(
            "Fold %s contains a target with zero or invalid local weights; increase `k` or `bandwidth`.",
            fold_labels[fold_index]
          ),
          call. = FALSE
        )
      }
      errors <- sweep(raw$predictions, 1L, y[test], FUN = "-")
      squared_error[, d_index] <- squared_error[, d_index] +
        colSums(errors^2)
      absolute_error[, d_index] <- absolute_error[, d_index] +
        colSums(abs(errors))
      observation_count[, d_index] <- observation_count[, d_index] +
        nrow(errors)
      convergence_count[, d_index] <- convergence_count[, d_index] +
        colSums(raw$converged != 0L)
      result_index <- result_index + 1L
      fold_results[[result_index]] <- data.frame(
        fold = fold_labels[fold_index],
        d = d[d_index],
        lambda = lambda,
        rmse = sqrt(colMeans(errors^2)),
        mae = colMeans(abs(errors)),
        convergence_rate = colMeans(raw$converged != 0L)
      )
    }
  }

  mse <- squared_error / observation_count
  mae <- absolute_error / observation_count
  cv <- expand.grid(
    lambda_index = seq_len(n_lambda),
    d_index = seq_len(n_d),
    KEEP.OUT.ATTRS = FALSE
  )
  cv$lambda <- lambda[cv$lambda_index]
  cv$d <- d[cv$d_index]
  cv$mse <- mse[cbind(cv$lambda_index, cv$d_index)]
  cv$rmse <- sqrt(cv$mse)
  cv$mae <- mae[cbind(cv$lambda_index, cv$d_index)]
  cv$convergence_rate <- convergence_count[
    cbind(cv$lambda_index, cv$d_index)
  ] / observation_count[cbind(cv$lambda_index, cv$d_index)]
  cv <- cv[, c("lambda", "d", "rmse", "mse", "mae",
               "convergence_rate")]
  best_index <- which.min(cv[[metric]])
  best <- cv[best_index, , drop = FALSE]

  selected_fit <- NULL
  if (isTRUE(refit)) {
    selected_fit <- gwr_sl_fit(
      x = x,
      y = y,
      neighbors = full_neighbors,
      kernel = kernel$name,
      bandwidth = if (bandwidth$adaptive) NULL else bandwidth$fixed,
      lambda = best$lambda,
      alpha = alpha,
      d = best$d,
      standardize = standardize,
      control = control
    )
  }

  object <- list(
    call = call,
    metric = metric,
    best = best,
    cv = cv,
    fold_results = do.call(rbind, fold_results),
    fold = fold,
    fold_labels = fold_labels,
    lambda = lambda,
    lambda_max = lambda_max,
    d = d,
    alpha = alpha,
    fit = selected_fit,
    full_neighbors = full_neighbors
  )
  class(object) <- "gwrs_cv"
  object
}

#' @export
print.gwrs_cv <- function(x, ...) {
  cat(x$method %||% "Spatial cross-validation for GWR-SL", "\n")
  cat("  folds:", length(x$fold_labels), "\n")
  cat("  candidates:", nrow(x$cv), "\n")
  cat("  selected lambda:", format(x$best$lambda), "\n")
  if (inherits(x, "gwrs_nonconvex_cv")) {
    cat("  penalty:", toupper(x$penalty),
        "  gamma:", format(x$gamma), "\n")
  } else {
    cat("  selected d:", format(x$best$d), "\n")
  }
  cat("  ", x$metric, ": ", format(x$best[[x$metric]]), "\n", sep = "")
  invisible(x)
}
