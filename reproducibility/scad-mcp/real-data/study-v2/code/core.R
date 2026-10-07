# Bandwidth-fraction and selected-fit diagnostic adapter for study-v2.

study_v1_app <- Sys.getenv("GWRS_STUDY_V1_APP")
stopifnot(nzchar(study_v1_app), dir.exists(study_v1_app))
sys.source(file.path(study_v1_app, "code/core.R"), envir = environment())

usa_study_v2_q_id <- function(q) {
  paste0("q", sprintf("%010.0f", as.double(q) * 1e9))
}

usa_study_v2_realized_k <- function(q, n_train, p) {
  stopifnot(
    length(q) == 1L, is.finite(q), q > 0, q <= 1,
    length(n_train) == 1L, n_train >= 1L,
    length(p) == 1L, p >= 1L
  )
  as.integer(min(n_train, max(2L * (p + 1L), ceiling(q * n_train))))
}

usa_grid <- function(cfg) {
  stopifnot(
    length(cfg$q) >= 1L, all(is.finite(cfg$q)),
    all(cfg$q > 0 & cfg$q <= 1), !is.unsorted(cfg$q),
    !anyDuplicated(usa_study_v2_q_id(cfg$q))
  )
  rows <- list()
  add <- function(where, method, q, alpha) {
    reference_k <- if (where == "local") {
      usa_study_v2_realized_k(q, cfg$n, cfg$p)
    } else 0L
    q_label <- if (where == "local") usa_study_v2_q_id(q) else "qglobal"
    id <- paste(where, method, q_label, paste0("a", alpha), sep = "-")
    rows[[length(rows) + 1L]] <<- data.frame(
      id = id, where = where, method = method,
      k = as.integer(reference_k),
      q = if (where == "local") as.double(q) else NA_real_,
      alpha = as.double(alpha), stringsAsFactors = FALSE
    )
  }
  for (method in c("gwr", "ridge", "lasso", "scad", "mcp", "en")) {
    alpha <- if (method == "en") cfg$en_alpha else
      if (method %in% c("ridge", "gwr")) 0 else 1
    for (q in cfg$q) for (value in alpha) add("local", method, q, value)
  }
  for (method in c("ols", "ridge", "lasso", "scad", "mcp", "en")) {
    alpha <- if (method == "en") cfg$en_alpha else
      if (method %in% c("ridge", "ols")) 0 else 1
    for (value in alpha) add("global", method, NA_real_, value)
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

usa_study_v2_resolve_candidate <- function(candidate, train) {
  value <- candidate
  if (identical(value$where, "local")) {
    stopifnot("q" %in% names(value), is.finite(value$q))
    value$k <- usa_study_v2_realized_k(
      value$q, nrow(train$x), ncol(train$x)
    )
  }
  value
}

usa_study_v2_geometry <- function(ctx, train, candidate) {
  neighbors <- usa_neighbors(ctx, train, candidate$k)
  bandwidth <- apply(neighbors$distance, 2L, max)
  safe_bandwidth <- bandwidth
  safe_bandwidth[safe_bandwidth <= 1e-12] <- 1
  weights <- pmax(
    1 - sweep(
      neighbors$distance, 2L, safe_bandwidth * (1 + 1e-10), "/"
    )^2,
    0
  )^2
  effective <- colSums(weights)^2 / colSums(weights^2)
  list(
    requested_q = as.double(candidate$q),
    training_n = nrow(train$x),
    realized_k = as.integer(candidate$k),
    realized_fraction = as.double(candidate$k / nrow(train$x)),
    bandwidth_m = c(
      minimum = min(bandwidth), median = stats::median(bandwidth),
      maximum = max(bandwidth)
    ),
    effective_weight_count = stats::setNames(
      c(min(effective), stats::median(effective), max(effective)),
      c("minimum", "median", "maximum")
    ),
    weight_sum = c(
      minimum = min(colSums(weights)), median = stats::median(colSums(weights)),
      maximum = max(colSums(weights))
    )
  )
}

usa_study_v2_base_path <- usa_path
usa_path <- function(ctx, train, candidate, until = NULL, keep = FALSE,
                     full_diagnostics = FALSE) {
  ctx$native_capture$selected_model_diagnostics <- NULL
  resolved <- usa_study_v2_resolve_candidate(candidate, train)
  result <- usa_study_v2_base_path(
    ctx, train, resolved, until, keep, full_diagnostics
  )
  if (identical(candidate$where, "local")) {
    result$geometry <- usa_study_v2_geometry(ctx, train, resolved)
  }
  if (isTRUE(keep) && isTRUE(train$same_locations) &&
      identical(candidate$where, "local") &&
      candidate$method != "gwr") {
    stopifnot(
      !is.null(result$coefficients_working),
      nrow(result$coefficients_working) == nrow(train$x)
    )
    raw_coefficients <- t(result$coefficients_working)
    coefficients <- ctx$g$backtransform_coefficients(
      raw_coefficients, train$center, train$scale
    )
    colnames(coefficients) <- c("(Intercept)", colnames(train$raw_x))
    fitted <- as.double(result$predictions[, ncol(result$predictions)])
    residuals <- train$y - fitted
    neighbors <- usa_neighbors(ctx, train, resolved$k)
    control <- ctx$g$gwrs_control(
      tolerance = ctx$cfg$tolerance,
      max_iterations = as.integer(ctx$cfg$retry_iterations),
      n_threads = 1L, grain_size = 16L, keep_data = FALSE,
      diagnostics = "standard",
      nonconvex_solver = "guarded_block",
      convex_path_solver = "guarded_working_set"
    )
    common <- list(
      prepared_x = train$x, y = train$y, fitted = fitted,
      residuals = residuals, raw_coefficients = raw_coefficients,
      coefficients = coefficients, neighbors = neighbors,
      kernel = list(code = 1L, name = "bisquare"),
      bandwidth = list(adaptive = TRUE, fixed = 0),
      lambda = tail(result$lambda, 1L),
      x_center = train$center, x_scale = train$scale,
      control = control
    )
    diagnostic <- if (candidate$method %in% c("scad", "mcp")) {
      do.call(ctx$g$compute_gwr_nonconvex_model_diagnostics, c(
        common,
        list(
          penalty = candidate$method,
          gamma = if (candidate$method == "scad") {
            ctx$cfg$scad_gamma
          } else ctx$cfg$mcp_gamma
        )
      ))
    } else {
      do.call(ctx$g$compute_gwr_model_diagnostics, c(
        common,
        list(alpha = candidate$alpha, d = 0)
      ))
    }
    diagnostic$selected_q <- as.double(candidate$q)
    diagnostic$selected_k <- as.integer(resolved$k)
    ctx$native_capture$selected_model_diagnostics <- diagnostic
  }
  result
}
