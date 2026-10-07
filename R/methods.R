#' @export
print.gwrs_fit <- function(x, ...) {
  cat(x$method, "fit\n")
  cat("  locations:", nrow(x$coefficients), "\n")
  cat("  predictors:", ncol(x$coefficients) - 1L, "\n")
  if (!is.null(x$lambda)) {
    if (inherits(x, "gwr_nonconvex_fit")) {
      cat("  lambda:", format(x$lambda),
          "  penalty:", toupper(x$penalty),
          "  gamma:", format(x$gamma), "\n")
    } else {
      cat("  lambda:", format(x$lambda),
          "  alpha:", format(x$alpha), "  d:", format(x$d), "\n")
    }
  }
  if (!is.null(x$lambda_diffusion)) {
    cat("  lambda_diffusion:", format(x$lambda_diffusion),
        "  gamma:", format(x$gamma), "\n")
    if (!is.null(x$solver)) {
      cat("  coefficient solver:", x$solver, "\n")
    }
  }
  convergence <- if (is.data.frame(x$diagnostics)) {
    mean(x$diagnostics$converged)
  } else if (is.list(x$diagnostics) &&
             !is.null(x$diagnostics$converged)) {
    as.numeric(x$diagnostics$converged)
  } else {
    NA_real_
  }
  if (is.finite(convergence)) {
    cat("  convergence:", sprintf("%.1f%%", 100 * convergence), "\n")
  }
  invisible(x)
}

#' @export
coef.gwrs_fit <- function(object, ...) {
  object$coefficients
}

#' @export
fitted.gwrs_fit <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.gwrs_fit <- function(object, ...) {
  object$residuals
}

#' Summarize a fitted gwrs model
#'
#' @param object A fitted `gwrs_fit` object.
#' @param ... Unused.
#'
#' @return An object of class `summary.gwrs_fit`.
#'
#' @examples
#' quake <- datasets::quakes[1:32, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' fit <- gwr_lasso_fit(
#'   x, quake$mag, coords = coords, k = 20, lambda = 0.05,
#'   control = gwrs_control(n_threads = 1)
#' )
#' summary(fit)
#'
#' @export
summary.gwrs_fit <- function(object, ...) {
  coefficient_summary <- t(vapply(
    seq_len(ncol(object$coefficients)),
    function(column) {
      values <- object$coefficients[, column]
      c(
        mean = mean(values),
        sd = stats::sd(values),
        min = min(values),
        q25 = unname(stats::quantile(values, 0.25)),
        median = stats::median(values),
        q75 = unname(stats::quantile(values, 0.75)),
        max = max(values),
        zero_fraction = mean(abs(values) <= 1e-10)
      )
    },
    numeric(8L)
  ))
  rownames(coefficient_summary) <- colnames(object$coefficients)
  errors <- object$residuals
  result <- list(
    call = object$call,
    method = object$method,
    coefficient_summary = coefficient_summary,
    performance = c(
      rmse = sqrt(mean(errors^2)),
      mae = mean(abs(errors)),
      bias = mean(errors)
    ),
    diagnostics = object$diagnostics,
    model_diagnostics = object$model_diagnostics
  )
  class(result) <- "summary.gwrs_fit"
  result
}

#' @export
print.summary.gwrs_fit <- function(x, ...) {
  cat(x$method, "summary\n\n")
  print(round(x$coefficient_summary, 5L))
  cat("\nIn-sample residual measures\n")
  print(round(x$performance, 6L))
  if (!is.null(x$model_diagnostics)) {
    global <- x$model_diagnostics$global
    cat("\nStatistical diagnostics\n")
    print(round(unlist(global[c(
      "r_squared", "adjusted_r_squared", "trS", "trStS", "enp", "edf",
      "aic", "aicc", "bic", "gcv"
    )]), 6L))
    cat("Basis:", x$model_diagnostics$validity$basis, "\n")
  }
  invisible(x)
}

new_design_from_fit <- function(object, newdata) {
  if (!is.null(object$terms)) {
    design <- stats::model.matrix(
      stats::delete.response(object$terms),
      data = newdata,
      contrasts.arg = object$contrasts,
      xlev = object$xlevels
    )
    intercept <- match("(Intercept)", colnames(design), nomatch = 0L)
    if (intercept > 0L) design <- design[, -intercept, drop = FALSE]
  } else {
    design <- as_design_matrix(newdata)
  }
  if (ncol(design) != length(object$predictor_names)) {
    stop("`newdata` does not reproduce the training design columns.",
         call. = FALSE)
  }
  storage.mode(design) <- "double"
  design
}

#' Predict from a fitted geographically weighted model
#'
#' For Ridge, Lasso, Elastic Net, SCAD and MCP, new target locations are fitted
#' from retained training data using only their training neighbors. `newdata`
#' supplies predictor values at targets, while `newcoords` controls the local
#' training samples.
#'
#' @param object A fitted `gwrs_fit`.
#' @param newdata Optional target design matrix or formula data frame.
#' @param newcoords Coordinates corresponding to `newdata`.
#' @param neighbors Optional cross-neighbor object constructed with
#'   `gwr_neighbors(training_coords, query_coords = newcoords, ...)`.
#' @param k Number of training neighbors when `neighbors` is not supplied.
#' @param type Return predictions, local coefficients, or both.
#' @param ... Unused.
#'
#' @return A numeric prediction vector, coefficient matrix, or a list of both.
#'
#' @examples
#' train <- datasets::quakes[1:30, ]
#' target <- datasets::quakes[31:33, ]
#' train_coords <- as.matrix(train[, c("long", "lat")])
#' train_x <- as.matrix(train[, c("depth", "stations")])
#' fit <- gwr_en_fit(
#'   train_x, train$mag, coords = train_coords, k = 18,
#'   lambda = 0.05, alpha = 0.7,
#'   control = gwrs_control(n_threads = 1)
#' )
#' prediction <- predict(
#'   fit,
#'   newdata = as.matrix(target[, c("depth", "stations")]),
#'   newcoords = as.matrix(target[, c("long", "lat")]),
#'   k = 18,
#'   type = "both"
#' )
#' prediction$fit
#' prediction$coefficients
#'
#' @export
predict.gwrs_fit <- function(object,
                             newdata = NULL,
                             newcoords = NULL,
                             neighbors = NULL,
                             k = NULL,
                             type = c("response", "coefficients", "both"),
                             ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    if (type == "response") return(object$fitted.values)
    if (type == "coefficients") return(object$coefficients)
    return(list(fit = object$fitted.values,
                coefficients = object$coefficients))
  }
  if (!inherits(object, "gwr_sl_fit") &&
      !inherits(object, "gwr_nonconvex_fit")) {
    stop("New-location prediction is currently defined for Gaussian penalized GWR fits.",
         call. = FALSE)
  }
  if (is.null(object$x) || is.null(object$y)) {
    stop("Refit with `gwrs_control(keep_data = TRUE)` for new predictions.",
         call. = FALSE)
  }
  target_x <- new_design_from_fit(object, newdata)
  if (is.null(neighbors)) {
    if (is.null(newcoords)) {
      stop("Supply `newcoords` or a precomputed cross-neighbor object.",
           call. = FALSE)
    }
    if (is.null(k)) k <- object$neighbors$k
    neighbors <- gwr_neighbors(
      object$neighbors$coords,
      query_coords = newcoords,
      k = k,
      include_self = FALSE
    )
  }
  validate_cross_neighbors(neighbors, nrow(object$x), nrow(target_x))
  x_train <- apply_standardization(
    object$x, object$x_center, object$x_scale
  )
  x_target <- apply_standardization(
    target_x, object$x_center, object$x_scale
  )
  kernel <- kernel_code(object$kernel)
  bandwidth <- resolve_bandwidth(object$bandwidth)
  raw <- if (inherits(object, "gwr_nonconvex_fit")) {
    run_cross_gwr_nonconvex_path(
      x_train,
      object$y,
      x_target,
      neighbors,
      kernel,
      bandwidth,
      object$lambda,
      object$penalty,
      object$gamma,
      FALSE,
      TRUE,
      object$control
    )
  } else {
    run_cross_gwr_sl_path(
      x_train,
      object$y,
      x_target,
      neighbors,
      kernel,
      bandwidth,
      object$lambda,
      object$alpha,
      object$d,
      FALSE,
      TRUE,
      object$control
    )
  }
  if (inherits(object, "gwr_nonconvex_fit")) {
    nonconvex_require_fit(raw, validate_control(object$control))
  }
  if (any(!is.finite(raw$predictions))) {
    stop("At least one target has no valid local weighted fit; increase `k` or `bandwidth`.",
         call. = FALSE)
  }
  coefficient_array <- backtransform_coefficient_path(
    raw$coefficients, object$x_center, object$x_scale
  )
  coefficients <- coefficient_array[, , 1L, drop = FALSE][, , 1L]
  colnames(coefficients) <- c("(Intercept)", object$predictor_names)
  predictions <- as.double(raw$predictions[, 1L])
  if (type == "response") return(predictions)
  if (type == "coefficients") return(coefficients)
  list(fit = predictions, coefficients = coefficients)
}

#' @export
coef.gwrs_cv <- function(object, ...) {
  if (is.null(object$fit)) {
    stop("This cross-validation object was created with `refit = FALSE`.",
         call. = FALSE)
  }
  stats::coef(object$fit, ...)
}

#' @export
predict.gwrs_cv <- function(object, ...) {
  if (is.null(object$fit)) {
    stop("This cross-validation object was created with `refit = FALSE`.",
         call. = FALSE)
  }
  stats::predict(object$fit, ...)
}
