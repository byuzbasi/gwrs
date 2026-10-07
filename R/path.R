#' Compute the loss-based GWR-SL lambda anchor
#'
#' Computes
#' \deqn{\lambda_{\max}=\max_{i,m}|x_{im}^T\bar W_i
#' (y-\bar y_i)|/\alpha,}
#' after the requested global predictor standardization. This is the exact
#' all-zero slope anchor when `d = 0`. For a nonzero GWR center (`d > 0`) it is
#' used only as a reproducible scale for constructing a path; the largest-path
#' solution need not be identically zero.
#'
#' @inheritParams gwr_sl_fit
#' @param alpha Mixing parameter in `(0, 1]`.
#'
#' @return One nonnegative number.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' neighbors <- gwr_neighbors(coords, k = 18)
#' gwr_sl_lambda_max(
#'   x, quake$mag, neighbors = neighbors, alpha = 0.7,
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#'
#' @export
gwr_sl_lambda_max <- function(x,
                              y,
                              coords = NULL,
                              neighbors = NULL,
                              k = NULL,
                              kernel = c("bisquare", "gaussian",
                                         "exponential", "tricube", "boxcar"),
                              bandwidth = NULL,
                              alpha = 0.5,
                              standardize = TRUE,
                              control = gwrs_control()) {
  x <- as_design_matrix(x)
  y <- as_response(y, nrow(x))
  control <- validate_control(control)
  neighbors <- resolve_neighbors(x, coords, neighbors, k)
  kernel <- kernel_code(kernel)
  bandwidth <- resolve_bandwidth(bandwidth)
  validate_penalty_parameters(alpha)
  if (alpha <= 0) {
    stop("`alpha` must be positive when computing `lambda_max`.",
         call. = FALSE)
  }
  prepared <- standardize_design(x, standardize)
  cpp_gwr_sl_lambda_max(
    prepared$x,
    y,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    as.double(alpha),
    control$n_threads,
    control$grain_size
  )
}

#' Fit a warm-started GWR-SL regularization path
#'
#' Local sufficient statistics and nearest neighbors are computed once per
#' target. Solutions then proceed from the largest to the smallest `lambda`
#' with warm starts, sequential screening, and a full KKT correction scan.
#' This avoids materializing a dense spatial weight matrix.
#'
#' @inheritParams gwr_sl_fit
#' @param lambda Decreasing vector of nonnegative regularization strengths.
#'   When `NULL`, a logarithmic sequence is generated from
#'   [gwr_sl_lambda_max()].
#' @param n_lambda Number of automatically generated values.
#' @param lambda_min_ratio Ratio of the smallest to largest generated value.
#' @param screening Whether to use sequential screening. A full KKT scan is
#'   always performed, so screened variables that violate optimality are
#'   restored before a solution is accepted.
#' @param keep_coefficients Whether to retain the full
#'   locations-by-coefficients-by-lambda array. The default avoids a potentially
#'   large allocation.
#' @param return_diagnostics Whether to compute conditional smoother traces,
#'   information criteria, and post-selection refit diagnostics along the path.
#'   These use `O(n * length(lambda))` memory and do not require retaining the
#'   coefficient cube.
#' @param ebic_gamma EBIC model-space multiplier in `[0, 1]`.
#'
#' @return An object of class `gwrs_path`.
#'
#' @examples
#' quake <- datasets::quakes[1:30, ]
#' coords <- as.matrix(quake[, c("long", "lat")])
#' x <- as.matrix(quake[, c("depth", "stations")])
#' path <- gwr_sl_path(
#'   x, quake$mag, coords = coords, k = 18,
#'   n_lambda = 5, lambda_min_ratio = 0.05,
#'   alpha = 0.7, d = 0.5, return_diagnostics = TRUE,
#'   control = gwrs_control(n_threads = 1, diagnostics = "none")
#' )
#' path
#' path$path_summary[, c("lambda", "mean_nonzero", "aicc")]
#'
#' @export
gwr_sl_path <- function(x,
                        y,
                        coords = NULL,
                        neighbors = NULL,
                        k = NULL,
                        kernel = c("bisquare", "gaussian", "exponential",
                                   "tricube", "boxcar"),
                        bandwidth = NULL,
                        lambda = NULL,
                        n_lambda = 50L,
                        lambda_min_ratio = 1e-3,
                        alpha = 0.5,
                        d = 0.5,
                        standardize = TRUE,
                        screening = TRUE,
                        keep_coefficients = FALSE,
                        return_diagnostics = FALSE,
                        ebic_gamma = 0.5,
                        control = gwrs_control()) {
  call <- match.call()
  x <- as_design_matrix(x)
  y <- as_response(y, nrow(x))
  control <- validate_control(control)
  neighbors <- resolve_neighbors(x, coords, neighbors, k)
  kernel <- kernel_code(kernel)
  bandwidth <- resolve_bandwidth(bandwidth)
  validate_penalty_parameters(alpha, d)
  if (length(ebic_gamma) != 1L || !is.finite(ebic_gamma) ||
      ebic_gamma < 0 || ebic_gamma > 1) {
    stop("`ebic_gamma` must be one number in [0, 1].", call. = FALSE)
  }
  if (length(d) != 1L) {
    stop("`d` must be one number in [0, 1].", call. = FALSE)
  }
  prepared <- standardize_design(x, standardize)

  lambda_max <- NULL
  if (is.null(lambda)) {
    if (alpha <= 0) {
      stop("Supply `lambda` explicitly when `alpha = 0`.", call. = FALSE)
    }
    lambda_max <- cpp_gwr_sl_lambda_max(
      prepared$x, y, neighbors$index, neighbors$distance, kernel$code,
      bandwidth$adaptive, bandwidth$fixed, as.double(alpha),
      control$n_threads, control$grain_size
    )
    lambda <- make_lambda_sequence(lambda_max, n_lambda, lambda_min_ratio)
  } else {
    lambda <- validate_lambda_path(lambda)
  }

  raw <- cpp_gwr_sl_path_predict(
    prepared$x,
    y,
    prepared$x,
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
    isTRUE(return_diagnostics),
    control$n_threads,
    control$grain_size
  )
  if (any(!is.finite(raw$predictions))) {
    stop("At least one target has no valid local weighted fit.", call. = FALSE)
  }

  coefficient_path <- NULL
  predictor_names <- colnames(x) %||% paste0("x", seq_len(ncol(x)))
  if (isTRUE(keep_coefficients)) {
    coefficient_path <- backtransform_coefficient_path(
      raw$coefficients, prepared$center, prepared$scale
    )
    dimnames(coefficient_path) <- list(
      location = seq_len(nrow(x)),
      coefficient = c("(Intercept)", predictor_names),
      lambda = format(lambda, digits = 7L)
    )
  }

  path_summary <- data.frame(
    lambda = lambda,
    mean_nonzero = colMeans(raw$nonzero),
    mean_iterations = colMeans(raw$iterations),
    convergence_rate = colMeans(raw$converged != 0L),
    max_kkt_violation = apply(raw$kkt, 2L, max, na.rm = TRUE)
  )
  path_diagnostics <- NULL
  if (isTRUE(return_diagnostics)) {
    path_diagnostics <- assemble_gwr_path_diagnostics(
      y = y,
      fitted_path = raw$predictions,
      hat_diag_path = raw$hat_diag_path,
      s2_diag_path = raw$s2_diag_path,
      post_fitted_path = raw$post_fitted_path,
      post_hat_diag_path = raw$post_hat_diag_path,
      post_s2_diag_path = raw$post_s2_diag_path,
      surface_active_count = raw$surface_active_count,
      ebic_gamma = ebic_gamma
    )
    path_summary <- cbind(path_summary, path_diagnostics)
  }
  object <- list(
    call = call,
    method = "GWR-SL path",
    lambda = lambda,
    lambda_max = lambda_max,
    alpha = as.double(alpha),
    d = as.double(d),
    fitted.values = raw$predictions,
    coefficients = coefficient_path,
    path_summary = path_summary,
    path_diagnostics = path_diagnostics,
    hat_diag_path = raw$hat_diag_path,
    s2_diag_path = raw$s2_diag_path,
    post_fitted_path = raw$post_fitted_path,
    post_hat_diag_path = raw$post_hat_diag_path,
    post_s2_diag_path = raw$post_s2_diag_path,
    surface_active_count = raw$surface_active_count,
    ebic_gamma = ebic_gamma,
    status = status_labels(raw$status),
    weight_sum = raw$weight_sum,
    local_bandwidth = raw$bandwidth,
    kernel = kernel$name,
    adaptive = bandwidth$adaptive,
    bandwidth = if (bandwidth$adaptive) NULL else bandwidth$fixed,
    neighbors = neighbors,
    standardize = isTRUE(standardize),
    x_center = prepared$center,
    x_scale = prepared$scale,
    predictor_names = predictor_names,
    x = if (control$keep_data) x else NULL,
    y = if (control$keep_data) y else NULL,
    control = control
  )
  class(object) <- "gwrs_path"
  object
}

assemble_gwr_path_diagnostics <- function(y,
                                          fitted_path,
                                          hat_diag_path,
                                          s2_diag_path,
                                          post_fitted_path,
                                          post_hat_diag_path,
                                          post_s2_diag_path,
                                          surface_active_count,
                                          ebic_gamma) {
  n <- length(y)
  residual_path <- sweep(fitted_path, 1L, y, FUN = "-")
  rss <- colSums(residual_path^2)
  trS <- colSums(hat_diag_path)
  trStS <- colSums(s2_diag_path)
  enp <- 2 * trS - trStS
  residual_df <- n - enp
  model_df <- trS + 1
  deviance <- n * (
    log(pmax(rss / n, .Machine$double.xmin)) + log(2 * pi) + 1
  )
  aic <- deviance + 2 * model_df
  aicc_denominator <- n - trS - 2
  aicc <- ifelse(
    aicc_denominator > 0,
    deviance + 2 * n * model_df / aicc_denominator,
    Inf
  )
  bic <- deviance + log(n) * model_df
  gcv_denominator <- 1 - trS / n
  gcv <- ifelse(
    gcv_denominator > 0,
    (rss / n) / gcv_denominator^2,
    Inf
  )
  post_residual_path <- sweep(post_fitted_path, 1L, y, FUN = "-")
  post_rss <- colSums(post_residual_path^2)
  post_trS <- colSums(post_hat_diag_path)
  post_trStS <- colSums(post_s2_diag_path)
  post_enp <- 2 * post_trS - post_trStS
  post_residual_df <- n - post_enp
  post_model_df <- post_trS + 1
  post_deviance <- n * (
    log(pmax(post_rss / n, .Machine$double.xmin)) + log(2 * pi) + 1
  )
  active_count <- as.matrix(surface_active_count)
  surface_size <- colSums(active_count > 0)
  p_search <- nrow(active_count)
  model_space <- vapply(
    surface_size,
    function(size) {
      if (size < 0 || size > p_search) return(Inf)
      if (p_search == 0L) return(0)
      lchoose(p_search, size)
    },
    numeric(1L)
  )
  ebic <- bic + 2 * ebic_gamma * model_space
  ebic_refit <- post_deviance + log(n) * post_model_df +
    2 * ebic_gamma * model_space
  ebic_refit[!is.finite(post_residual_df) | post_residual_df <= 1] <- Inf
  data.frame(
    rss = rss,
    deviance = deviance,
    log_likelihood = -0.5 * deviance,
    trS = trS,
    trStS = trStS,
    enp = enp,
    residual_df = residual_df,
    aic = aic,
    aicc = aicc,
    bic = bic,
    gcv = gcv,
    ebic = ebic,
    surface_size = surface_size,
    post_rss = post_rss,
    post_trS = post_trS,
    post_trStS = post_trStS,
    post_enp = post_enp,
    post_residual_df = post_residual_df,
    ebic_refit = ebic_refit
  )
}

#' @export
print.gwrs_path <- function(x, ...) {
  cat(if (inherits(x, "gwrs_nonconvex_path")) {
    paste("Fast", x$method)
  } else {
    "Fast GWR-SL regularization path"
  }, "\n")
  cat("  locations:", nrow(x$fitted.values), "\n")
  cat("  lambda values:", length(x$lambda), "\n")
  if (inherits(x, "gwrs_nonconvex_path")) {
    cat("  penalty:", toupper(x$penalty),
        "  gamma:", format(x$gamma), "\n")
  } else {
    cat("  alpha:", format(x$alpha), "  d:", format(x$d), "\n")
  }
  cat("  convergence range:",
      sprintf("%.1f%%--%.1f%%",
              100 * min(x$path_summary$convergence_rate),
              100 * max(x$path_summary$convergence_rate)), "\n")
  invisible(x)
}
