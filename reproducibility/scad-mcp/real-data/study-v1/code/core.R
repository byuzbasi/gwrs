# Scientific and checkpoint primitives for the USA county application.
# Estimation is delegated to the byte-frozen gwrs 0.4.0.9006 native cores.

usa_hash <- function(value) {
  digest::digest(value, algo = "sha256", serializeVersion = 3)
}

usa_files <- function(files, base) {
  files <- sort(normalizePath(files, mustWork = TRUE))
  base <- normalizePath(base, mustWork = TRUE)
  stopifnot(all(file.exists(files)), !any(nzchar(Sys.readlink(files))))
  data.frame(
    file = substring(files, nchar(base) + 2L),
    bytes = unname(file.info(files)$size),
    sha256 = vapply(
      files, digest::digest, "", algo = "sha256", file = TRUE
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

usa_put <- function(ctx, object, path) {
  if (file.exists(path)) {
    stopifnot(identical(readRDS(path), object))
  } else {
    ctx$io$gwrs_sim_atomic_rds(object, path)
  }
  invisible(path)
}

usa_receipt <- function(ctx, path) {
  receipt <- usa_files(path, dirname(path))
  usa_put(ctx, receipt, paste0(path, ".sha256.rds"))
}

usa_read_shard <- function(path, id, signature) {
  stopifnot(file.exists(path), file.exists(paste0(path, ".sha256.rds")))
  envelope <- readRDS(path)
  stopifnot(
    identical(envelope$id, id),
    identical(envelope$signature, signature),
    identical(envelope$payload_sha256, usa_hash(envelope$payload)),
    identical(
      usa_files(path, dirname(path)),
      readRDS(paste0(path, ".sha256.rds"))
    )
  )
  envelope$payload
}

usa_task <- function(ctx, out, id, signature, fun) {
  stopifnot(grepl("^[A-Za-z0-9.-]+$", id))
  path <- file.path(out, "tasks", paste0(id, ".rds"))
  if (file.exists(path)) {
    if (!file.exists(paste0(path, ".sha256.rds"))) {
      envelope <- readRDS(path)
      stopifnot(
        identical(envelope$id, id),
        identical(envelope$signature, signature),
        identical(
          envelope$payload_sha256, usa_hash(envelope$payload)
        )
      )
      usa_receipt(ctx, path)
    }
    cat("RESUME", id, "\n")
    return(usa_read_shard(path, id, signature))
  }
  cat("START", id, "\n")
  flush.console()
  value <- fun()
  envelope <- list(
    id = id,
    signature = signature,
    payload_sha256 = usa_hash(value),
    payload = value
  )
  usa_put(ctx, envelope, path)
  usa_receipt(ctx, path)
  cat("DONE", id, "\n")
  flush.console()
  usa_read_shard(path, id, signature)
}

usa_grid <- function(cfg) {
  rows <- list()
  add <- function(where, method, k, alpha) {
    id <- paste(
      where, method, paste0("k", k), paste0("a", alpha), sep = "-"
    )
    rows[[length(rows) + 1L]] <<- data.frame(
      id = id, where = where, method = method,
      k = as.integer(k), alpha = as.double(alpha),
      stringsAsFactors = FALSE
    )
  }
  for (method in c("gwr", "ridge", "lasso", "scad", "mcp", "en")) {
    for (k in cfg$k) {
      alpha <- if (method == "en") cfg$en_alpha else
        if (method %in% c("ridge", "gwr")) 0 else 1
      for (value in alpha) add("local", method, k, value)
    }
  }
  for (method in c("ols", "ridge", "lasso", "scad", "mcp", "en")) {
    alpha <- if (method == "en") cfg$en_alpha else
      if (method %in% c("ridge", "ols")) 0 else 1
    for (value in alpha) add("global", method, 0L, value)
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

usa_methods <- function(grid) {
  unique(paste(grid$where, grid$method, sep = "-"))
}

usa_task_plan <- function(grid, folds, mode, calibration = NULL) {
  stopifnot(mode %in% c("smoke", "calibration", "full"))
  if (mode == "calibration") {
    selected <- grid[
      grid$where == "global" | grid$k == calibration$local_k,
      , drop = FALSE
    ]
    return(data.frame(
      id = paste0("cal-", selected$id),
      stage = "calibration",
      outer = calibration$outer,
      inner = calibration$inner,
      validation_fold = NA_integer_,
      candidate_id = selected$id,
      method_key = paste(selected$where, selected$method, sep = "-"),
      cost_class = paste(selected$where, selected$method, sep = "-"),
      stringsAsFactors = FALSE
    ))
  }
  rows <- list()
  append_rows <- function(value) {
    rows[[length(rows) + 1L]] <<- value
  }
  for (outer in folds) {
    for (inner in setdiff(folds, outer)) {
      append_rows(data.frame(
        id = paste0("o", outer, "-i", inner, "-", grid$id),
        stage = "nested_inner", outer = outer, inner = inner,
        validation_fold = NA_integer_, candidate_id = grid$id,
        method_key = paste(grid$where, grid$method, sep = "-"),
        cost_class = paste(grid$where, grid$method, sep = "-"),
        stringsAsFactors = FALSE
      ))
    }
    methods <- usa_methods(grid)
    append_rows(data.frame(
      id = paste0("o", outer, "-refit-", methods),
      stage = "outer_refit", outer = outer, inner = NA_integer_,
      validation_fold = NA_integer_, candidate_id = NA_character_,
      method_key = methods, cost_class = methods,
      stringsAsFactors = FALSE
    ))
  }
  for (validation_fold in folds) {
    append_rows(data.frame(
      id = paste0("f", validation_fold, "-", grid$id),
      stage = "full_tune", outer = NA_integer_, inner = NA_integer_,
      validation_fold = validation_fold, candidate_id = grid$id,
      method_key = paste(grid$where, grid$method, sep = "-"),
      cost_class = paste(grid$where, grid$method, sep = "-"),
      stringsAsFactors = FALSE
    ))
  }
  methods <- usa_methods(grid)
  append_rows(data.frame(
    id = paste0("full-refit-", methods),
    stage = "full_refit", outer = NA_integer_, inner = NA_integer_,
    validation_fold = NA_integer_, candidate_id = NA_character_,
    method_key = methods, cost_class = methods,
    stringsAsFactors = FALSE
  ))
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

usa_nested_split <- function(data, outer, inner) {
  stopifnot(outer != inner)
  fold <- data$fold
  list(
    train = which(fold != outer & fold != inner),
    test = which(fold == inner),
    kind = "nested_inner"
  )
}

usa_outer_split <- function(data, outer) {
  list(
    train = which(data$fold != outer),
    test = which(data$fold == outer),
    kind = "outer_refit"
  )
}

usa_full_tune_split <- function(data, validation_fold) {
  list(
    train = which(data$fold != validation_fold),
    test = which(data$fold == validation_fold),
    kind = "full_tune"
  )
}

usa_full_refit_split <- function(data) {
  index <- seq_len(nrow(data$x))
  list(train = index, test = index, kind = "full_refit")
}

usa_context <- function(ctx, x, y, coords, target_x, target_coords,
                        same_locations = FALSE) {
  stopifnot(
    nrow(x) == length(y), nrow(coords) == nrow(x),
    nrow(target_x) == nrow(target_coords)
  )
  prepared <- ctx$g$standardize_design(x, TRUE)
  list(
    raw_x = x,
    x = prepared$x,
    y = y,
    coords = coords,
    target_raw_x = target_x,
    target_x = ctx$g$apply_standardization(
      target_x, prepared$center, prepared$scale
    ),
    target_coords = target_coords,
    center = prepared$center,
    scale = prepared$scale,
    same_locations = isTRUE(same_locations),
    neighbors = new.env(parent = emptyenv())
  )
}

usa_training <- function(ctx, data, split) {
  stopifnot(
    length(split$train) > 1L, length(split$test) > 0L,
    split$kind == "full_refit" ||
      !length(intersect(split$train, split$test))
  )
  same <- identical(split$kind, "full_refit")
  usa_context(
    ctx,
    data$x[split$train, , drop = FALSE], data$y[split$train],
    data$coords[split$train, , drop = FALSE],
    data$x[split$test, , drop = FALSE],
    data$coords[split$test, , drop = FALSE],
    same_locations = same
  )
}

usa_neighbors <- function(ctx, train, k) {
  key <- paste0("k", k, "-same", as.integer(train$same_locations))
  if (!exists(key, train$neighbors, inherits = FALSE)) {
    value <- if (train$same_locations) {
      ctx$g$gwr_neighbors(train$coords, k = k)
    } else {
      ctx$g$gwr_neighbors(
        train$coords, query_coords = train$target_coords,
        k = k, include_self = FALSE
      )
    }
    assign(key, value, train$neighbors)
  }
  get(key, train$neighbors, inherits = FALSE)
}

usa_anchor_neighbors <- function(ctx, train, k) {
  key <- paste0("anchor-k", k)
  if (!exists(key, train$neighbors, inherits = FALSE)) {
    value <- ctx$g$gwr_neighbors(train$coords, k = k)
    assign(key, value, train$neighbors)
  }
  get(key, train$neighbors, inherits = FALSE)
}

usa_anchor <- function(ctx, train, candidate) {
  method <- candidate$method
  stopifnot(method %in% c("scad", "mcp", "lasso", "en"))
  if (candidate$where == "global") {
    anchor <- max(abs(drop(
      crossprod(train$x, train$y - mean(train$y)) / nrow(train$x)
    )))
    if (method == "en") anchor <- anchor / candidate$alpha
  } else {
    # The path anchor is defined from the training sample. Prediction
    # neighborhoods instead have one column per held-out target and are used
    # only when evaluating that path.
    neighbors <- usa_anchor_neighbors(ctx, train, candidate$k)
    if (method %in% c("scad", "mcp")) {
      anchor <- ctx$g$cpp_gwr_nonconvex_lambda_max(
        train$x, train$y, neighbors$index, neighbors$distance,
        1L, TRUE, 0,
        ctx$g$resolve_nonconvex_penalty(method)$code,
        if (method == "scad") ctx$cfg$scad_gamma else ctx$cfg$mcp_gamma,
        ctx$threads, 16L
      )
    } else {
      anchor <- ctx$g$cpp_gwr_sl_lambda_max(
        train$x, train$y, neighbors$index, neighbors$distance,
        1L, TRUE, 0, candidate$alpha, ctx$threads, 16L
      )
    }
  }
  if (!is.finite(anchor) || anchor <= 0) {
    stop("Nonpositive training lambda anchor requires scientific review")
  }
  as.double(anchor)
}

usa_full_gwr <- function(ctx, train, candidate) {
  stopifnot(
    candidate$where == "local", candidate$method == "gwr",
    train$same_locations
  )
  neighbors <- usa_neighbors(ctx, train, candidate$k)
  attempts <- list()
  fit <- NULL
  for (cap in unique(c(ctx$cfg$max_iterations, ctx$cfg$retry_iterations))) {
    control <- ctx$g$gwrs_control(
      tolerance = ctx$cfg$tolerance,
      max_iterations = as.integer(cap),
      n_threads = 1L,
      grain_size = 16L,
      keep_data = TRUE,
      diagnostics = "full",
      nonconvex_solver = "guarded_block",
      convex_path_solver = "guarded_working_set"
    )
    fit <- ctx$g$gwr_fit(
      train$raw_x, train$y, neighbors = neighbors,
      kernel = "bisquare", standardize = TRUE, control = control
    )
    valid <- all(fit$diagnostics$converged) &&
      all(is.finite(fit$fitted.values)) &&
      all(is.finite(fit$coefficients))
    attempts[[length(attempts) + 1L]] <- list(
      iterations_cap = as.integer(cap),
      valid = valid,
      convergence_rate = mean(fit$diagnostics$converged),
      max_iterations = max(fit$diagnostics$iterations),
      max_stationarity = max(fit$diagnostics$max_kkt_violation)
    )
    if (valid) break
  }
  inference <- ctx$g$gwr_local_inference(
    fit, adjust = "BH", adjust_scope = "coefficient",
    estimator = "penalized"
  )
  tests <- ctx$g$gwr_diagnostic_tests(
    fit,
    block_size = as.integer(ctx$cfg$full_gwr_diagnostics$f_test_block_size),
    memory_limit_mb = ctx$cfg$full_gwr_diagnostics$f_test_memory_limit_mb
  )
  list(
    predictions = matrix(fit$fitted.values, ncol = 1L),
    valid = valid,
    parameter = 0,
    lambda = 0,
    anchor = NA_real_,
    attempts = attempts,
    nonzero = ncol(train$raw_x),
    stationarity = max(fit$diagnostics$max_kkt_violation),
    gap = NA_real_,
    iterations = mean(fit$diagnostics$iterations),
    coefficients_working = NULL,
    coefficients_original = fit$coefficients,
    geometry = list(
      weight_sum = range(fit$diagnostics$weight_sum),
      effective_weight_count = range(
        fit$model_diagnostics$local$effective_sample_size
      ),
      bandwidth_m = range(fit$diagnostics$bandwidth)
    ),
    scaling = list(center = fit$x_center, scale = fit$x_scale),
    full_gwr = list(
      model_diagnostics = fit$model_diagnostics,
      local_inference = inference,
      diagnostic_tests = tests,
      conditional_on_selected_k = TRUE,
      selected_k = as.integer(candidate$k)
    )
  )
}

usa_path <- function(ctx, train, candidate, until = NULL, keep = FALSE,
                     full_diagnostics = FALSE) {
  if (isTRUE(full_diagnostics)) {
    return(usa_full_gwr(ctx, train, candidate))
  }
  method <- candidate$method
  nonconvex <- method %in% c("scad", "mcp")
  global <- candidate$where == "global"
  unpenalized <- method %in% c("gwr", "ols")
  if (unpenalized) {
    anchor <- NA_real_
    lambda <- 0
    parameter <- 0
  } else if (method == "ridge") {
    anchor <- NA_real_
    lambda <- exp(seq(
      log(ctx$cfg$ridge_lambda_max), log(ctx$cfg$ridge_lambda_min),
      length.out = ctx$cfg$n_lambda
    ))
    parameter <- lambda
  } else {
    anchor <- usa_anchor(ctx, train, candidate)
    parameter <- exp(seq(
      0, log(ctx$cfg$lambda_ratio_min), length.out = ctx$cfg$n_lambda
    ))
    lambda <- anchor * parameter
  }
  if (!is.null(until)) {
    stopifnot(until %in% seq_along(lambda))
    lambda <- lambda[seq_len(until)]
    parameter <- parameter[seq_len(until)]
  }
  if (global) {
    neighbors <- ctx$g$gwr_neighbors(
      train$coords,
      query_coords = train$coords[1L, , drop = FALSE],
      k = nrow(train$x), include_self = FALSE
    )
    target <- train$target_x[1L, , drop = FALSE]
    kernel <- 4L
  } else {
    neighbors <- usa_neighbors(ctx, train, candidate$k)
    target <- train$target_x
    kernel <- 1L
  }
  run <- function(cap) {
    arguments <- list(
      x_train = train$x, y_train = train$y, x_target = target,
      neighbor_index = neighbors$index,
      neighbor_distance = neighbors$distance,
      kernel = kernel, adaptive = TRUE, fixed_bandwidth = 0,
      lambda = lambda, tolerance = ctx$cfg$tolerance,
      max_iterations = as.integer(cap), screening = TRUE,
      keep_coefficients = keep || global,
      return_diagnostics = FALSE,
      n_threads = ctx$threads, grain_size = 16L
    )
    if (nonconvex) {
      arguments$penalty <-
        ctx$g$resolve_nonconvex_penalty(method)$code
      arguments$gamma <- if (method == "scad") {
        ctx$cfg$scad_gamma
      } else {
        ctx$cfg$mcp_gamma
      }
      do.call(ctx$g$cpp_gwr_nonconvex_path_predict, arguments)
    } else {
      arguments$alpha <- candidate$alpha
      arguments$d <- 0
      do.call(ctx$g$cpp_gwr_sl_path_predict, arguments)
    }
  }
  check <- function(raw) {
    apply(
      raw$converged == 1L & is.finite(raw$predictions), 2L, all
    )
  }
  attempts <- list()
  for (cap in unique(c(ctx$cfg$max_iterations, ctx$cfg$retry_iterations))) {
    raw <- run(cap)
    valid <- check(raw)
    attempts[[length(attempts) + 1L]] <- list(
      iterations_cap = as.integer(cap),
      valid = valid,
      convergence_rate = colMeans(raw$converged == 1L),
      max_iterations = apply(raw$iterations, 2L, max),
      max_stationarity = apply(
        if (nonconvex) raw$stationarity else raw$kkt, 2L, max
      ),
      max_coordinate_gap = if (nonconvex) {
        apply(raw$coordinate_gap, 2L, max)
      } else {
        NULL
      }
    )
    if (all(valid)) break
  }
  prediction <- raw$predictions
  if (global) {
    coefficients <- matrix(
      raw$coefficients,
      nrow = ncol(train$x) + 1L, ncol = length(lambda)
    )
    prediction <- sweep(
      train$target_x %*% coefficients[-1L, , drop = FALSE],
      2L, coefficients[1L, ], "+"
    )
  }
  valid <- valid & apply(is.finite(prediction), 2L, all)
  bandwidth <- apply(neighbors$distance, 2L, max)
  bandwidth[bandwidth <= 1e-12] <- 1
  weights <- if (global) {
    matrix(1, nrow(neighbors$distance), ncol(neighbors$distance))
  } else {
    pmax(
      1 - sweep(
        neighbors$distance, 2L, bandwidth * (1 + 1e-10), "/"
      )^2,
      0
    )^2
  }
  effective <- colSums(weights)^2 / colSums(weights^2)
  retained <- if (keep) {
    t(matrix(
      raw$coefficients[, , length(lambda), drop = FALSE],
      nrow = ncol(train$x) + 1L
    ))
  } else {
    NULL
  }
  list(
    predictions = prediction,
    valid = valid,
    parameter = parameter,
    lambda = lambda,
    anchor = anchor,
    attempts = attempts,
    nonzero = colMeans(raw$nonzero),
    stationarity = apply(
      if (nonconvex) raw$stationarity else raw$kkt, 2L, max
    ),
    gap = if (nonconvex) {
      apply(raw$coordinate_gap, 2L, max)
    } else {
      rep(NA_real_, length(lambda))
    },
    iterations = colMeans(raw$iterations),
    coefficients_working = retained,
    coefficients_original = NULL,
    geometry = list(
      weight_sum = range(raw$weight_sum),
      effective_weight_count = range(effective),
      bandwidth_m = range(raw$bandwidth)
    ),
    scaling = list(center = train$center, scale = train$scale),
    full_gwr = NULL
  )
}

usa_score <- function(fit, y) {
  stopifnot(nrow(fit$predictions) == length(y))
  errors <- sweep(fit$predictions, 1L, y, "-")
  data.frame(
    index = seq_len(ncol(errors)),
    parameter = fit$parameter,
    lambda = fit$lambda,
    valid = fit$valid,
    n = length(y),
    sse = colSums(errors^2),
    sae = colSums(abs(errors)),
    error_sum = colSums(errors),
    mean_nonzero = fit$nonzero,
    mean_iterations = fit$iterations,
    max_stationarity = fit$stationarity,
    max_coordinate_gap = fit$gap,
    stringsAsFactors = FALSE
  )
}

usa_choose <- function(scores) {
  stopifnot(
    length(scores) >= 1L,
    all(vapply(scores, nrow, 0L) == nrow(scores[[1L]]))
  )
  valid <- Reduce(
    `&`, lapply(scores, function(value) {
      value$valid & is.finite(value$sse)
    })
  )
  total_n <- Reduce(`+`, lapply(scores, `[[`, "n"))
  mse <- Reduce(`+`, lapply(scores, `[[`, "sse")) / total_n
  data.frame(
    index = seq_along(mse), mse = mse, n = total_n, valid = valid
  )
}

usa_common_graph <- function(ctx, ids) {
  path <- file.path(
    ctx$app, "frozen/input/spatial/tables",
    "common-neighbor-graph-directed-v1.csv"
  )
  graph <- utils::read.csv(
    path, stringsAsFactors = FALSE,
    colClasses = c(from_fips = "character", to_fips = "character")
  )
  n <- length(ids)
  stopifnot(
    n == 3107L, !anyDuplicated(ids), nrow(graph) == 18126L,
    all(graph$from_index %in% seq_len(n)),
    all(graph$to_index %in% seq_len(n)),
    all(graph$from_index != graph$to_index),
    identical(ids[graph$from_index], graph$from_fips),
    identical(ids[graph$to_index], graph$to_fips),
    all(graph$binary_weight == 1),
    !anyDuplicated(paste(graph$from_index, graph$to_index, sep = "-"))
  )
  degree <- tabulate(graph$from_index, nbins = n)
  stopifnot(
    min(degree) == 1L, max(degree) == 14L,
    all(abs(graph$row_standard_weight - 1 / degree[graph$from_index]) <
          1e-12)
  )
  k <- max(degree)
  index <- matrix(rep(seq_len(n), each = k), nrow = k, ncol = n)
  distance <- matrix(1, nrow = k, ncol = n)
  position <- integer(n)
  for (edge in seq_len(nrow(graph))) {
    from <- graph$from_index[edge]
    position[from] <- position[from] + 1L
    index[position[from], from] <- graph$to_index[edge]
    distance[position[from], from] <- 0
  }
  stopifnot(identical(position, degree))
  list(
    index = index, distance = distance, degree = degree,
    source = path, source_sha256 = usa_files(path, dirname(path))$sha256
  )
}

usa_moran_one <- function(ctx, residuals, graph, permutations, seed) {
  stopifnot(
    length(residuals) == ncol(graph$index), all(is.finite(residuals)),
    length(permutations) == 1L, permutations >= 0L,
    length(seed) == 1L, seed >= 0L
  )
  ctx$g$cpp_gwr_moran(
    as.double(residuals), graph$index, graph$distance,
    1L, FALSE, 1, TRUE, as.integer(permutations), as.double(seed),
    as.integer(ctx$post_threads), 16L
  )
}

usa_residual_spatial <- function(ctx, predictions, data) {
  cfg <- ctx$cfg$postanalysis
  methods <- unique(predictions$method)
  graph <- usa_common_graph(ctx, data$ids)
  raw <- vector("list", length(methods))
  names(raw) <- methods
  global <- vector("list", length(methods))
  local <- vector("list", length(methods))
  for (method_index in seq_along(methods)) {
    method <- methods[method_index]
    value <- predictions[predictions$method == method, , drop = FALSE]
    value <- value[match(data$ids, value$sample_id), , drop = FALSE]
    stopifnot(
      nrow(value) == length(data$ids), identical(value$sample_id, data$ids),
      identical(value$observed, data$y)
    )
    residual <- value$observed - value$prediction
    answer <- usa_moran_one(
      ctx, residual, graph,
      cfg$moran_lisa_permutations, cfg$moran_lisa_seed
    )
    stopifnot(
      is.finite(answer$I), is.finite(answer$expected_I),
      is.finite(answer$p_value), length(answer$local_I) == length(data$ids),
      length(answer$local_p) == length(data$ids),
      all(is.finite(answer$local_I)), all(is.finite(answer$local_p))
    )
    adjusted <- stats::p.adjust(answer$local_p, method = "BH")
    significant <- is.finite(adjusted) & adjusted <= 0.05
    cluster <- rep("Not Significant", length(residual))
    cluster[significant & answer$centered_residual > 0 &
              answer$spatial_lag > 0] <- "High-High"
    cluster[significant & answer$centered_residual < 0 &
              answer$spatial_lag < 0] <- "Low-Low"
    cluster[significant & answer$centered_residual > 0 &
              answer$spatial_lag < 0] <- "High-Low"
    cluster[significant & answer$centered_residual < 0 &
              answer$spatial_lag > 0] <- "Low-High"
    global[[method_index]] <- data.frame(
      method = method, n = length(residual), I = answer$I,
      expected_I = answer$expected_I, p_value = answer$p_value,
      permutations = as.integer(cfg$moran_lisa_permutations),
      seed = as.integer(cfg$moran_lisa_seed),
      stringsAsFactors = FALSE
    )
    local[[method_index]] <- data.frame(
      sample_id = data$ids, method = method,
      residual = residual,
      residual_centered = answer$centered_residual,
      spatial_lag = answer$spatial_lag,
      local_I = answer$local_I, z_score = answer$local_z,
      p_value = answer$local_p, p_BH_within_method = adjusted,
      significant_BH_005 = significant, cluster_BH_005 = cluster,
      stringsAsFactors = FALSE
    )
    raw[[method_index]] <- answer
  }
  global <- do.call(rbind, global)
  global$p_BH_across_methods <- stats::p.adjust(global$p_value, method = "BH")
  local <- do.call(rbind, local)
  rownames(global) <- NULL
  rownames(local) <- NULL
  list(
    schema = "gwrs-usa-counties-acs2024-residual-spatial-v1",
    graph = list(
      file = basename(graph$source), sha256 = graph$source_sha256,
      directed_edges = sum(graph$degree), minimum_degree = min(graph$degree),
      maximum_degree = max(graph$degree), row_standardized = TRUE
    ),
    permutations = as.integer(cfg$moran_lisa_permutations),
    seed = as.integer(cfg$moran_lisa_seed),
    global = global, local = local, raw = raw
  )
}

usa_evaluate <- function(ctx, train, data, split, candidate,
                         until = NULL, keep = FALSE,
                         full_diagnostics = FALSE) {
  tryCatch({
    warnings <- character()
    fit <- withCallingHandlers(
      usa_path(
        ctx, train, candidate, until, keep,
        full_diagnostics = full_diagnostics
      ),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    )
    score <- usa_score(fit, data$y[split$test])
    if (length(warnings)) score$valid[] <- FALSE
    payload <- list(
      candidate = candidate,
      score = score,
      attempts = fit$attempts,
      anchor = fit$anchor,
      train_ids = data$ids[split$train],
      test_ids = data$ids[split$test],
      scaling_sha256 = usa_hash(fit$scaling),
      geometry = fit$geometry,
      warnings = warnings
    )
    if (keep) {
      selected_index <- ncol(fit$predictions)
      payload$predictions <- fit$predictions[, selected_index]
      payload$observed <- data$y[split$test]
      payload$training_mean <- mean(data$y[split$train])
      if (!is.null(fit$coefficients_original)) {
        coefficients <- fit$coefficients_original
        working <- NULL
      } else {
        working <- fit$coefficients_working
        stopifnot(!is.null(working), all(is.finite(working)))
        coefficients <- ctx$g$backtransform_coefficients(
          t(working), train$center, train$scale
        )
      }
      colnames(coefficients) <- c("(Intercept)", colnames(data$x))
      payload$coefficients <- coefficients
      payload["selected"] <- list(if (
        candidate$method %in% c("scad", "mcp", "lasso", "en")
      ) {
        abs(working[, -1L, drop = FALSE]) > ctx$cfg$selection_threshold
      } else {
        NULL
      })
      payload$selected_index <- selected_index
      payload$full_gwr <- fit$full_gwr
    }
    payload
  }, error = function(error) {
    list(
      candidate = candidate,
      error = conditionMessage(error),
      train_ids = data$ids[split$train],
      test_ids = data$ids[split$test]
    )
  })
}

usa_check_task <- function(value, data, split, candidate, refit = FALSE,
                           cfg, full_refit = FALSE) {
  stopifnot(
    identical(value$train_ids, data$ids[split$train]),
    identical(value$test_ids, data$ids[split$test]),
    identical(value$candidate, candidate)
  )
  if (!is.null(value$error)) {
    if (refit) stop(value$error)
    return(invisible(FALSE))
  }
  expected <- if (candidate$method %in% c("gwr", "ols")) {
    1L
  } else {
    cfg$n_lambda
  }
  if (!refit) stopifnot(nrow(value$score) == expected)
  stopifnot(
    all(value$score$n == length(split$test)),
    all(value$score$index == seq_len(nrow(value$score)))
  )
  if (refit) {
    stopifnot(
      isTRUE(tail(value$score$valid, 1L)),
      length(value$predictions) == length(split$test),
      identical(value$observed, data$y[split$test]),
      all(is.finite(value$predictions)),
      ncol(value$coefficients) == ncol(data$x) + 1L,
      all(is.finite(value$coefficients)),
      nrow(value$coefficients) == if (candidate$where == "local") {
        length(split$test)
      } else {
        1L
      }
    )
    observed_sse <- sum((value$observed - value$predictions)^2)
    expected_sse <- tail(value$score$sse, 1L)
    stopifnot(
      abs(observed_sse - expected_sse) <
        1e-8 * (1 + expected_sse)
    )
    selector <- candidate$method %in% c("scad", "mcp", "lasso", "en")
    if (selector) {
      stopifnot(
        is.matrix(value$selected), is.logical(value$selected),
        !anyNA(value$selected),
        identical(
          dim(value$selected),
          c(
            if (candidate$where == "local") length(split$test) else 1L,
            ncol(data$x)
          )
        )
      )
    } else {
      stopifnot(is.null(value[["selected", exact = TRUE]]))
    }
  }
  if (full_refit && candidate$where == "local" &&
      candidate$method == "gwr") {
    diagnostics <- value$full_gwr
    stopifnot(
      is.list(diagnostics),
      isTRUE(diagnostics$conditional_on_selected_k),
      diagnostics$selected_k == candidate$k,
      nrow(diagnostics$local_inference$local) ==
        nrow(data$x) * (ncol(data$x) + 1L),
      identical(
        diagnostics$diagnostic_tests$global_tests$test,
        c("F1", "F2")
      ),
      nrow(diagnostics$diagnostic_tests$coefficient_tests) ==
        ncol(data$x) + 1L
    )
  }
  invisible(TRUE)
}
