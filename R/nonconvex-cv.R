run_cross_gwr_nonconvex_path <- function(x_train,
                                         y_train,
                                         x_target,
                                         neighbors,
                                         kernel,
                                         bandwidth,
                                         lambda,
                                         penalty,
                                         gamma,
                                         screening,
                                         keep_coefficients,
                                         control) {
  validate_cross_neighbors(neighbors, nrow(x_train), nrow(x_target))
  specification <- resolve_nonconvex_penalty(penalty, gamma)
  cpp_gwr_nonconvex_path_predict(
    x_train,
    y_train,
    x_target,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    lambda,
    specification$code,
    specification$gamma,
    control$tolerance,
    control$max_iterations,
    isTRUE(screening),
    isTRUE(keep_coefficients),
    FALSE,
    control$n_threads,
    control$grain_size,
    nonconvex_guarded(control)
  )
}

#' Spatial cross-validation for GWR-SCAD and GWR-MCP
#'
#' Evaluates a common decreasing lambda path using held-out spatial blocks.
#' Standardization and nearest-neighbor construction are repeated inside each
#' training fold, preventing transformation and neighborhood leakage.
#'
#' With the guarded-block solver, a lambda is eligible only if all held-out
#' observations have finite, converged predictions. Otherwise all its error
#' scores are `NA`; no observation is removed. `cv$observation_count` retains
#' the planned full denominator and `cv$valid` reports eligibility. If none
#' are valid, the object has zero rows in `best`, a `NULL` fit and a warning.
#' Refit uses the existing single-lambda initialization, not the CV path's
#' warm start, and is independently checked for convergence.
#'
#' @inheritParams cv_gwr_sl
#' @inheritParams gwr_nonconvex_fit
#' @param cv_cache Optional read-only preparation from [gwr_nonconvex_cv_cache()].
#'   Inputs must match exactly. Penalty paths and responses are not shared.
#' @param ... Additional arguments passed by [cv_gwr_scad()] or
#'   [cv_gwr_mcp()] to `cv_gwr_nonconvex()`.
#'
#' @return An object inheriting from `gwrs_nonconvex_cv` and `gwrs_cv`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fold <- spatial_folds(coords, n_folds = 3, seed = 42)
#' cv <- cv_gwr_mcp(
#'   x, quake$mag, coords, k = 8, fold = fold,
#'   lambda = c(0.12, 0.04), refit = FALSE,
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#' cv$best
#'
#' @export
cv_gwr_nonconvex <- function(x,
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
                             penalty = c("scad", "mcp"),
                             gamma = NULL,
                             standardize = TRUE,
                             screening = TRUE,
                             metric = c("rmse", "mse", "mae"),
                             refit = TRUE,
                             control = gwrs_control(),
                             cv_cache = NULL) {
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
  specification <- resolve_nonconvex_penalty(penalty, gamma)
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

  if (is.null(cv_cache)) {
    full_neighbors <- gwr_neighbors(coords, k = k)
    full_prepared <- standardize_design(x, standardize)
  } else {
    validate_nonconvex_cv_cache(cv_cache, x, coords, k, fold, standardize)
    full_neighbors <- cv_cache$full_neighbors
    full_prepared <- cv_cache$full_prepared
  }
  lambda_max <- NULL
  if (is.null(lambda)) {
    lambda_max <- cpp_gwr_nonconvex_lambda_max(
      full_prepared$x,
      y,
      full_neighbors$index,
      full_neighbors$distance,
      kernel$code,
      bandwidth$adaptive,
      bandwidth$fixed,
      specification$code,
      specification$gamma,
      control$n_threads,
      control$grain_size
    )
    lambda <- make_lambda_sequence(lambda_max, n_lambda, lambda_min_ratio)
  } else {
    lambda <- validate_lambda_path(lambda)
  }

  n_lambda <- length(lambda)
  squared_error <- numeric(n_lambda)
  absolute_error <- numeric(n_lambda)
  observation_count <- integer(n_lambda)
  convergence_count <- numeric(n_lambda)
  fold_results <- vector("list", nlevels(fold_factor))
  guarded <- nonconvex_guarded(control)
  eligible <- rep(TRUE, n_lambda)

  for (fold_index in seq_len(nlevels(fold_factor))) {
    test <- fold_id == fold_index
    train <- !test
    if (is.null(cv_cache)) {
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
    } else {
      prepared <- cv_cache$folds[[fold_index]]$prepared
      x_test <- cv_cache$folds[[fold_index]]$x_test
      cross_neighbors <- cv_cache$folds[[fold_index]]$neighbors
    }
    raw <- run_cross_gwr_nonconvex_path(
      prepared$x,
      y[train],
      x_test,
      cross_neighbors,
      kernel,
      bandwidth,
      lambda,
      specification$name,
      specification$gamma,
      screening,
      FALSE,
      control
    )
    if (any(!is.finite(raw$predictions)) && !guarded) {
      stop(
        sprintf(
          "Fold %s contains a target with zero or invalid local weights; increase `k` or `bandwidth`.",
          fold_labels[fold_index]
        ),
        call. = FALSE
      )
    }
    errors <- sweep(raw$predictions, 1L, y[test], FUN = "-")
    fold_valid <- colSums(nonconvex_valid(raw)) == sum(test)
    if (guarded) {
      eligible <- eligible & fold_valid
      # Invalidate the complete lambda score, never a subset of observations.
      errors[, !fold_valid] <- NA_real_
    }
    squared_error <- squared_error + colSums(errors^2)
    absolute_error <- absolute_error + colSums(abs(errors))
    observation_count <- observation_count + nrow(errors)
    convergence_count <- convergence_count + colSums(raw$converged != 0L)
    fold_results[[fold_index]] <- data.frame(
      fold = fold_labels[fold_index],
      lambda = lambda,
      rmse = sqrt(colMeans(errors^2)),
      mae = colMeans(abs(errors)),
      convergence_rate = colMeans(raw$converged != 0L),
      max_stationarity_violation = if (guarded) nonconvex_column_max(raw$stationarity) else apply(
        raw$stationarity, 2L, max, na.rm = TRUE
      ),
      max_coordinate_gap = if (guarded) nonconvex_column_max(raw$coordinate_gap) else apply(
        raw$coordinate_gap, 2L, max, na.rm = TRUE
      )
    )
    if (guarded) {
      fold_results[[fold_index]]$max_stationarity_violation <-
        nonconvex_column_max(raw$stationarity)
      fold_results[[fold_index]]$max_coordinate_gap <- nonconvex_column_max(raw$coordinate_gap)
      fold_results[[fold_index]] <- cbind(fold_results[[fold_index]], nonconvex_solver_summary(raw))
      fold_results[[fold_index]]$observation_count <- sum(test)
    }
  }

  mse <- squared_error / observation_count
  cv <- data.frame(
    lambda = lambda,
    rmse = sqrt(mse),
    mse = mse,
    mae = absolute_error / observation_count,
    convergence_rate = convergence_count / observation_count
  )
  if (guarded) {
    cv$valid <- eligible & observation_count == nrow(x)
    cv$observation_count <- observation_count
    cv[!cv$valid, c("rmse", "mse", "mae")] <- NA_real_
  }
  best_index <- which.min(cv[[metric]])
  best <- cv[best_index, , drop = FALSE]

  selected_fit <- NULL
  if (guarded && length(best_index) == 0L) {
    warning("No lambda converged for every held-out observation; best and fit are unavailable.", call. = FALSE)
  }
  if (isTRUE(refit) && length(best_index) > 0L) {
    selected_fit <- gwr_nonconvex_fit(
      x = x,
      y = y,
      neighbors = full_neighbors,
      kernel = kernel$name,
      bandwidth = if (bandwidth$adaptive) NULL else bandwidth$fixed,
      lambda = best$lambda,
      penalty = specification$name,
      gamma = specification$gamma,
      standardize = standardize,
      control = control
    )
  }

  object <- list(
    call = call,
    method = paste("Spatial cross-validation for", specification$method),
    metric = metric,
    penalty = specification$name,
    gamma = specification$gamma,
    best = best,
    cv = cv,
    fold_results = do.call(rbind, fold_results),
    fold = fold,
    fold_labels = fold_labels,
    lambda = lambda,
    lambda_max = lambda_max,
    fit = selected_fit,
    full_neighbors = full_neighbors
  )
  if (guarded) object$control <- control
  class(object) <- c("gwrs_nonconvex_cv", "gwrs_cv")
  object
}

#' @rdname cv_gwr_nonconvex
#' @export
cv_gwr_scad <- function(x, y, coords, k, fold, ..., gamma = 3.7) {
  cv_gwr_nonconvex(
    x, y, coords, k, fold, ..., penalty = "scad", gamma = gamma
  )
}

#' @rdname cv_gwr_nonconvex
#' @export
cv_gwr_mcp <- function(x, y, coords, k, fold, ..., gamma = 3) {
  cv_gwr_nonconvex(
    x, y, coords, k, fold, ..., penalty = "mcp", gamma = gamma
  )
}
