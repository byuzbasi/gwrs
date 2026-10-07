# Scientific helpers for the paper-level penalized GWR selection study.

gwrs_penalized_study_k_candidates <- function(n, p) {
  n <- as.integer(n)
  p <- as.integer(p)
  if (length(n) != 1L || is.na(n) || n < 2L ||
      length(p) != 1L || is.na(p) || p < 4L) {
    stop("`n` and `p` must be valid scalar dimensions.", call. = FALSE)
  }
  candidates <- c(
    max(2L * (p + 1L), ceiling(0.025 * n)),
    max(4L * (p + 1L), ceiling(0.050 * n)),
    max(8L * (p + 1L), ceiling(0.100 * n))
  )
  candidates <- sort(unique(as.integer(candidates)))
  if (any(candidates >= n)) {
    stop("The approved k rule requires every candidate to be below n.",
         call. = FALSE)
  }
  candidates
}

gwrs_penalized_study_cells <- function(mode) {
  mode <- match.arg(
    mode, c("smoke", "pilot", "max_cell_calibration", "production")
  )
  if (identical(mode, "smoke")) {
    return(data.frame(
      section = rep("smoke", 2L),
      scenario_id = c(
        "smooth_multiscale_local_sparse",
        "regional_sparse_discontinuous"
      ),
      design_variant = "baseline",
      sampling_design = "uniform",
      n = 120L,
      p = 6L,
      predictor_rho = 0.2,
      spatial_rho = 0.3,
      local_correlation = 0,
      snr = 3,
      target_grid_side = 8L,
      replications = 1L,
      stringsAsFactors = FALSE
    ))
  }
  if (identical(mode, "pilot")) {
    return(data.frame(
      section = "pilot",
      scenario_id = "regional_sparse_discontinuous",
      design_variant = "baseline",
      sampling_design = "uniform",
      n = 2000L,
      p = 30L,
      predictor_rho = 0.7,
      spatial_rho = 0.7,
      local_correlation = 0,
      snr = 3,
      target_grid_side = 16L,
      replications = 1L,
      stringsAsFactors = FALSE
    ))
  }
  if (identical(mode, "max_cell_calibration")) {
    return(data.frame(
      section = "max_cell_calibration",
      scenario_id = "regional_sparse_discontinuous",
      design_variant = "p100",
      sampling_design = "uniform",
      n = 10000L,
      p = 100L,
      predictor_rho = 0.7,
      spatial_rho = 0.7,
      local_correlation = 0,
      snr = 3,
      target_grid_side = 32L,
      replications = 1L,
      stringsAsFactors = FALSE
    ))
  }

  core <- expand.grid(
    scenario_id = c(
      "smooth_multiscale_local_sparse",
      "regional_sparse_discontinuous"
    ),
    n = c(2000L, 10000L),
    predictor_rho = c(0, 0.7),
    spatial_rho = c(0, 0.7),
    snr = c(1, 3),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  core$section <- "main"
  core$design_variant <- "baseline"
  core$sampling_design <- "uniform"
  core$p <- 30L
  core$local_correlation <- 0
  core$target_grid_side <- 32L
  core$replications <- 50L
  reference <- data.frame(
    section = "reference",
    scenario_id = c("stationary_sparse", "null_slopes"),
    design_variant = "baseline",
    sampling_design = "uniform",
    n = 5000L,
    p = 30L,
    predictor_rho = 0.7,
    spatial_rho = 0.7,
    local_correlation = 0,
    snr = 3,
    target_grid_side = 32L,
    replications = 100L,
    stringsAsFactors = FALSE
  )
  stress <- data.frame(
    section = "stress",
    scenario_id = c(
      "smooth_multiscale_local_sparse",
      "regional_sparse_discontinuous",
      "local_collinearity_sparse",
      "smooth_multiscale_local_sparse",
      "regional_sparse_discontinuous"
    ),
    design_variant = c(
      "p100", "p100", "local_collinearity", "clustered_sampling",
      "l_domain_boundary"
    ),
    sampling_design = c("uniform", "uniform", "uniform", "clustered", "l_shaped"),
    n = 10000L,
    p = c(100L, 100L, 30L, 30L, 30L),
    predictor_rho = c(0.7, 0.7, 0, 0.7, 0.7),
    spatial_rho = 0.7,
    local_correlation = c(0, 0, 0.9, 0, 0),
    snr = 3,
    target_grid_side = 32L,
    replications = 20L,
    stringsAsFactors = FALSE
  )
  cells <- rbind(
    core[c(
      "section", "scenario_id", "design_variant", "sampling_design",
      "n", "p", "predictor_rho", "spatial_rho", "local_correlation",
      "snr", "target_grid_side", "replications"
    )],
    reference,
    stress
  )
  rownames(cells) <- NULL
  cells
}

gwrs_penalized_study_plan <- function(mode, master_seed) {
  mode <- match.arg(
    mode, c("smoke", "pilot", "max_cell_calibration", "production")
  )
  if (identical(mode, "max_cell_calibration")) {
    production <- gwrs_penalized_study_plan("production", master_seed)
    selected <- production[
      production$section == "stress" &
        production$design_variant == "p100" &
        production$scenario_id == "regional_sparse_discontinuous" &
        production$replication == 1L,
      , drop = FALSE
    ]
    if (nrow(selected) != 1L) {
      stop("Could not resolve the frozen production max-cell task.",
           call. = FALSE)
    }
    selected$task_id <- paste0("max-cell-calibration-", selected$task_id)
    selected$cell_id <- paste0("max-cell-calibration-", selected$cell_id)
    selected$section <- "max_cell_calibration"
    rownames(selected) <- NULL
    return(selected)
  }
  cells <- gwrs_penalized_study_cells(mode)
  pieces <- lapply(seq_len(nrow(cells)), function(index) {
    cell <- cells[index, , drop = FALSE]
    cross_code <- sprintf("%03d", round(100 * cell$predictor_rho))
    spatial_code <- sprintf("%03d", round(100 * cell$spatial_rho))
    snr_code <- sprintf("%02d", round(10 * cell$snr))
    cell_id <- sprintf(
      "%s-%s-%s-n%05d-p%03d-rx%s-rs%s-s%s",
      cell$section, cell$scenario_id, cell$design_variant, cell$n, cell$p,
      cross_code, spatial_code, snr_code
    )
    data.frame(
      cell_id = cell_id,
      section = cell$section,
      scenario_id = cell$scenario_id,
      design_variant = cell$design_variant,
      sampling_design = cell$sampling_design,
      n = as.integer(cell$n),
      p = as.integer(cell$p),
      predictor_rho = as.double(cell$predictor_rho),
      spatial_rho = as.double(cell$spatial_rho),
      local_correlation = as.double(cell$local_correlation),
      snr = as.integer(cell$snr),
      target_grid_side = as.integer(cell$target_grid_side),
      replication = seq_len(cell$replications),
      stringsAsFactors = FALSE
    )
  })
  plan <- do.call(rbind, pieces)
  plan$task_id <- sprintf(
    "%s-rep-%05d", plan$cell_id, plan$replication
  )
  task_index <- seq_len(nrow(plan))
  plan$seed <- as.integer(
    ((master_seed + task_index * 104729) %% 2147483646) + 1
  )
  plan$data_seed <- as.integer(
    ((master_seed + task_index * 130363) %% 2147483646) + 1
  )
  plan$fold_seed <- as.integer(
    ((master_seed + task_index * 169087) %% 2147483646) + 1
  )
  plan$k_candidates <- vapply(
    seq_len(nrow(plan)),
    function(index) paste(
      gwrs_penalized_study_k_candidates(plan$n[index], plan$p[index]),
      collapse = ":"
    ),
    character(1L)
  )
  plan <- plan[c(
    "task_id", "cell_id", "section", "scenario_id", "design_variant",
    "sampling_design", "n", "p", "predictor_rho", "spatial_rho",
    "local_correlation", "snr", "target_grid_side", "replication",
    "k_candidates", "data_seed", "fold_seed", "seed"
  )]
  rownames(plan) <- NULL
  plan
}

gwrs_penalized_study_standardize_dgp <- function(x) {
  centered <- sweep(x, 2L, colMeans(x), FUN = "-")
  scale <- sqrt(colMeans(centered^2))
  if (any(!is.finite(scale)) || any(scale <= 1e-12)) {
    stop("The generated predictor matrix is numerically singular.",
         call. = FALSE)
  }
  sweep(centered, 2L, scale, FUN = "/")
}

gwrs_penalized_study_apply_training_standardization <- function(x, n_train) {
  training <- x[seq_len(n_train), , drop = FALSE]
  center <- colMeans(training)
  centered <- sweep(training, 2L, center, FUN = "-")
  scale <- sqrt(colMeans(centered^2))
  if (any(!is.finite(scale)) || any(scale <= 1e-12)) {
    stop("The generated predictor matrix is numerically singular.",
         call. = FALSE)
  }
  list(
    train = sweep(centered, 2L, scale, FUN = "/"),
    target = sweep(
      sweep(x[-seq_len(n_train), , drop = FALSE], 2L, center, FUN = "-"),
      2L, scale, FUN = "/"
    ),
    center = center,
    scale = scale
  )
}

gwrs_penalized_study_seed <- function(seed, offset) {
  as.integer(((as.double(seed) + as.double(offset) * 104729) %% 2147483646) + 1)
}

gwrs_penalized_study_sample_coords <- function(n, sampling_design) {
  n <- as.integer(n)
  if (identical(sampling_design, "uniform")) {
    return(cbind(east = runif(n), north = runif(n)))
  }
  if (identical(sampling_design, "clustered")) {
    centers <- rbind(c(0.20, 0.22), c(0.30, 0.78), c(0.72, 0.30), c(0.78, 0.76))
    membership <- sample.int(4L, n, replace = TRUE, prob = c(0.30, 0.20, 0.30, 0.20))
    coords <- matrix(NA_real_, nrow = n, ncol = 2L)
    pending <- seq_len(n)
    while (length(pending)) {
      candidate <- centers[membership[pending], , drop = FALSE] +
        matrix(rnorm(2L * length(pending), sd = 0.075), ncol = 2L)
      valid <- rowSums(candidate >= 0 & candidate <= 1) == 2L
      coords[pending[valid], ] <- candidate[valid, , drop = FALSE]
      pending <- pending[!valid]
    }
    colnames(coords) <- c("east", "north")
    return(coords)
  }
  if (identical(sampling_design, "l_shaped")) {
    coords <- matrix(NA_real_, nrow = n, ncol = 2L)
    filled <- 0L
    while (filled < n) {
      needed <- n - filled
      candidate <- cbind(runif(2L * needed), runif(2L * needed))
      valid <- !(candidate[, 1L] > 0.55 & candidate[, 2L] > 0.55)
      accepted <- candidate[valid, , drop = FALSE]
      take <- min(needed, nrow(accepted))
      if (take > 0L) {
        rows <- filled + seq_len(take)
        coords[rows, ] <- accepted[seq_len(take), , drop = FALSE]
        filled <- filled + take
      }
    }
    colnames(coords) <- c("east", "north")
    return(coords)
  }
  stop("Unknown sampling design.", call. = FALSE)
}

gwrs_penalized_study_target_grid <- function(side, sampling_design) {
  midpoint <- (seq_len(side) - 0.5) / side
  coords <- as.matrix(expand.grid(east = midpoint, north = midpoint))
  if (identical(sampling_design, "l_shaped")) {
    coords <- coords[!(coords[, 1L] > 0.55 & coords[, 2L] > 0.55), , drop = FALSE]
  }
  colnames(coords) <- c("east", "north")
  coords
}

gwrs_penalized_study_domain_distance <- function(coords, sampling_design) {
  east <- coords[, 1L]
  north <- coords[, 2L]
  outer <- pmin(east, 1 - east, north, 1 - north)
  if (!identical(sampling_design, "l_shaped")) return(outer)
  internal <- ifelse(
    east <= 0.55 & north >= 0.55,
    0.55 - east,
    ifelse(
      north <= 0.55 & east >= 0.55,
      0.55 - north,
      sqrt(pmax(0, 0.55 - east)^2 + pmax(0, 0.55 - north)^2)
    )
  )
  pmin(outer, internal)
}

gwrs_penalized_study_plateau_wendland <- function(coords, center, radius) {
  distance <- sqrt(rowSums(sweep(coords, 2L, center, FUN = "-")^2))
  value <- numeric(length(distance))
  value[distance <= radius / 2] <- 1
  taper <- distance > radius / 2 & distance < radius
  scaled <- (distance[taper] / radius - 0.5) / 0.5
  value[taper] <- (1 - scaled)^4 * (1 + 4 * scaled)
  list(
    value = value,
    support = distance < radius,
    boundary = distance >= 0.75 * radius & distance <= 1.25 * radius
  )
}

gwrs_penalized_study_surfaces <- function(coords, scenario) {
  n <- nrow(coords)
  east <- coords[, 1L]
  north <- coords[, 2L]
  beta <- matrix(0, nrow = n, ncol = 4L)
  support <- matrix(FALSE, nrow = n, ncol = 4L)
  support_boundary <- matrix(FALSE, nrow = n, ncol = 4L)
  intercept <- 0.7 + 0.25 * cos(2 * pi * east) -
    0.20 * sin(2 * pi * north)

  if (identical(scenario, "stationary_sparse")) {
    beta[] <- rep(c(1.20, -1.00, 0.65, 0.45), each = n)
    support[] <- TRUE
  } else if (identical(scenario, "smooth_multiscale_local_sparse")) {
    broad <- 1.00 + 0.35 * sin(2 * pi * east)
    hotspot_2 <- gwrs_penalized_study_plateau_wendland(
      coords, c(0.30, 0.70), 0.40
    )
    hotspot_31 <- gwrs_penalized_study_plateau_wendland(
      coords, c(0.72, 0.25), 0.25
    )
    hotspot_32 <- gwrs_penalized_study_plateau_wendland(
      coords, c(0.74, 0.76), 0.20
    )
    hotspot_4 <- gwrs_penalized_study_plateau_wendland(
      coords, c(0.22, 0.25), 0.16
    )
    beta[, 1L] <- broad
    beta[, 2L] <- -1.00 * hotspot_2$value
    beta[, 3L] <- 0.80 * hotspot_31$value + 0.65 * hotspot_32$value
    beta[, 4L] <- 0.50 * hotspot_4$value
    support[, 1L] <- TRUE
    support[, 2L] <- hotspot_2$support
    support[, 3L] <- hotspot_31$support | hotspot_32$support
    support[, 4L] <- hotspot_4$support
    support_boundary[, 2L] <- hotspot_2$boundary
    support_boundary[, 3L] <- hotspot_31$boundary | hotspot_32$boundary
    support_boundary[, 4L] <- hotspot_4$boundary
  } else if (identical(scenario, "regional_sparse_discontinuous")) {
    support[, 1L] <- east <= 0.55
    support[, 2L] <- east > 0.35 & north > 0.45
    support[, 3L] <- north <= 0.50
    support[, 4L] <- east + north > 1.20
    beta[, 1L] <- 1.20 * support[, 1L]
    beta[, 2L] <- -1.00 * support[, 2L]
    beta[, 3L] <- 0.80 * support[, 3L]
    beta[, 4L] <- 0.60 * support[, 4L]
    support_boundary[, 1L] <- abs(east - 0.55) <= 0.05
    support_boundary[, 2L] <-
      (abs(east - 0.35) <= 0.05 & north >= 0.40) |
      (abs(north - 0.45) <= 0.05 & east >= 0.30)
    support_boundary[, 3L] <- abs(north - 0.50) <= 0.05
    support_boundary[, 4L] <- abs(east + north - 1.20) <= 0.07
  } else if (identical(scenario, "local_collinearity_sparse")) {
    hotspot <- gwrs_penalized_study_plateau_wendland(
      coords, c(0.25, 0.75), 0.22
    )
    beta[, 1L] <- 1.20
    beta[, 3L] <- 0.70 + 0.20 * sin(2 * pi * north)
    beta[, 4L] <- 0.60 * hotspot$value
    support[, 1L] <- TRUE
    support[, 3L] <- TRUE
    support[, 4L] <- hotspot$support
    support_boundary[, 4L] <- hotspot$boundary
  } else if (!identical(scenario, "null_slopes")) {
    stop("Unknown paper-level simulation scenario.", call. = FALSE)
  }
  list(
    intercept = intercept,
    beta = beta,
    support = support,
    support_boundary = support_boundary
  )
}

gwrs_penalized_study_predictors <- function(coords,
                                             n_train,
                                             p,
                                             predictor_rho,
                                             spatial_rho,
                                             local_correlation) {
  n_total <- nrow(coords)
  innovation <- matrix(rnorm(n_total * p), nrow = n_total, ncol = p)
  correlated <- innovation
  if (p > 1L && predictor_rho != 0) {
    innovation_scale <- sqrt(1 - predictor_rho^2)
    for (column in 2:p) {
      correlated[, column] <- predictor_rho * correlated[, column - 1L] +
        innovation_scale * innovation[, column]
    }
  }

  graph <- gwr_neighbors(coords, k = min(9L, n_total))
  external_index <- graph$index[-1L, , drop = FALSE]
  graph_k <- nrow(external_index)
  if (spatial_rho != 0) {
    adjacency <- Matrix::sparseMatrix(
      i = rep(seq_len(n_total), each = graph_k),
      j = as.integer(external_index),
      x = rep.int(1 / graph_k, n_total * graph_k),
      dims = c(n_total, n_total)
    )
    operator <- Matrix::Diagonal(n_total) - spatial_rho * adjacency
    correlated <- as.matrix(Matrix::solve(operator, correlated))
  }

  local_region <- sqrt((coords[, 1L] - 0.70)^2 +
                         (coords[, 2L] - 0.30)^2) < 0.25
  if (p >= 2L && local_correlation > 0) {
    residual <- correlated[, 2L]
    residual <- residual -
      sum(residual[local_region] * correlated[local_region, 1L]) /
      sum(correlated[local_region, 1L]^2) * correlated[, 1L]
    residual_scale <- stats::sd(residual[local_region])
    if (!is.finite(residual_scale) || residual_scale <= 1e-12) {
      stop("Local-collinearity residual is numerically singular.",
           call. = FALSE)
    }
    residual <- residual / residual_scale
    correlated[local_region, 2L] <-
      local_correlation * correlated[local_region, 1L] +
      sqrt(1 - local_correlation^2) * residual[local_region]
  }

  standardized <- gwrs_penalized_study_apply_training_standardization(
    correlated, n_train
  )
  x_all <- rbind(standardized$train, standardized$target)
  train_region <- local_region[seq_len(n_train)]
  external_train <- external_index[, seq_len(n_train), drop = FALSE]
  spatial_product <- mean(
    standardized$train[, 1L][rep(seq_len(n_train), each = graph_k)] *
      x_all[, 1L][as.integer(external_train)]
  )
  safe_cor <- function(a, b) {
    if (length(a) < 3L || stats::sd(a) <= 1e-12 || stats::sd(b) <= 1e-12) {
      return(NA_real_)
    }
    stats::cor(a, b)
  }
  list(
    x = standardized$train,
    x_target = standardized$target,
    local_region = train_region,
    target_local_region = local_region[-seq_len(n_train)],
    realized_predictor_correlation = if (p >= 2L) {
      safe_cor(standardized$train[, 1L], standardized$train[, 2L])
    } else NA_real_,
    realized_local_correlation = if (p >= 2L) {
      safe_cor(
        standardized$train[train_region, 1L],
        standardized$train[train_region, 2L]
      )
    } else NA_real_,
    realized_outside_correlation = if (p >= 2L) {
      safe_cor(
        standardized$train[!train_region, 1L],
        standardized$train[!train_region, 2L]
      )
    } else NA_real_,
    realized_spatial_neighbor_product = spatial_product
  )
}

gwrs_penalized_study_simulate <- function(task, config) {
  n <- as.integer(task$n[[1L]])
  p <- as.integer(task$p[[1L]])
  rho <- as.double(task$predictor_rho[[1L]])
  spatial_rho <- as.double(task$spatial_rho[[1L]])
  local_correlation <- as.double(task$local_correlation[[1L]])
  snr <- as.double(task$snr[[1L]])
  scenario <- task$scenario_id[[1L]]
  sampling_design <- task$sampling_design[[1L]]

  set.seed(gwrs_penalized_study_seed(task$data_seed[[1L]], 1L))
  coords <- gwrs_penalized_study_sample_coords(n, sampling_design)
  target_coords <- gwrs_penalized_study_target_grid(
    task$target_grid_side[[1L]], sampling_design
  )
  combined_coords <- rbind(coords, target_coords)

  set.seed(gwrs_penalized_study_seed(task$data_seed[[1L]], 2L))
  predictor <- gwrs_penalized_study_predictors(
    combined_coords, n, p, rho, spatial_rho, local_correlation
  )
  x <- predictor$x
  x_target <- predictor$x_target
  colnames(x) <- paste0("x", seq_len(p))
  colnames(x_target) <- colnames(x)

  surfaces <- gwrs_penalized_study_surfaces(combined_coords, scenario)
  beta <- matrix(0, nrow = n, ncol = p, dimnames = list(NULL, colnames(x)))
  beta_target <- matrix(
    0, nrow = nrow(target_coords), ncol = p,
    dimnames = list(NULL, colnames(x))
  )
  support <- matrix(FALSE, nrow = n, ncol = p)
  target_support <- matrix(FALSE, nrow = nrow(target_coords), ncol = p)
  support_boundary <- matrix(FALSE, nrow = n, ncol = p)
  beta[, seq_len(4L)] <- surfaces$beta[seq_len(n), , drop = FALSE]
  beta_target[, seq_len(4L)] <- surfaces$beta[-seq_len(n), , drop = FALSE]
  support[, seq_len(4L)] <- surfaces$support[seq_len(n), , drop = FALSE]
  target_support[, seq_len(4L)] <- surfaces$support[-seq_len(n), , drop = FALSE]
  support_boundary[, seq_len(4L)] <-
    surfaces$support_boundary[seq_len(n), , drop = FALSE]
  intercept <- surfaces$intercept[seq_len(n)]
  target_intercept <- surfaces$intercept[-seq_len(n)]

  signal <- intercept + rowSums(x * beta)
  target_signal <- target_intercept + rowSums(x_target * beta_target)
  signal_sd <- stats::sd(signal)
  if (!is.finite(signal_sd) || signal_sd <= .Machine$double.eps) {
    stop("The generated signal has no finite variation.", call. = FALSE)
  }
  noise_sd <- signal_sd / sqrt(snr)
  set.seed(gwrs_penalized_study_seed(task$data_seed[[1L]], 3L))
  y_train <- signal + rnorm(n, sd = noise_sd)
  set.seed(gwrs_penalized_study_seed(task$data_seed[[1L]], 4L))
  y_test <- signal + rnorm(n, sd = noise_sd)
  set.seed(gwrs_penalized_study_seed(task$data_seed[[1L]], 5L))
  y_target <- target_signal + rnorm(nrow(target_coords), sd = noise_sd)
  list(
    coords = coords,
    target_coords = target_coords,
    x = x,
    x_target = x_target,
    y_train = y_train,
    y_test = y_test,
    y_target = y_target,
    signal = signal,
    target_signal = target_signal,
    signal_sd = signal_sd,
    noise_sd = noise_sd,
    intercept = intercept,
    target_intercept = target_intercept,
    beta = beta,
    beta_target = beta_target,
    support = support,
    target_support = target_support,
    support_boundary = support_boundary,
    domain_edge = gwrs_penalized_study_domain_distance(
      coords, sampling_design
    ) <= 0.08,
    target_domain_edge = gwrs_penalized_study_domain_distance(
      target_coords, sampling_design
    ) <= 0.08,
    local_region = predictor$local_region,
    target_local_region = predictor$target_local_region,
    realized_predictor_correlation = predictor$realized_predictor_correlation,
    realized_local_correlation = predictor$realized_local_correlation,
    realized_outside_correlation = predictor$realized_outside_correlation,
    realized_spatial_neighbor_product =
      predictor$realized_spatial_neighbor_product
  )
}

gwrs_penalized_study_select_k <- function(data,
                                           k_candidates,
                                           fold,
                                           config,
                                           control,
                                           heartbeat = function(...) NULL) {
  results <- lapply(k_candidates, function(k) {
    heartbeat(paste0("bandwidth-cv-k", k))
    timing <- system.time({
      cv <- cv_gwr_sl(
        data$x,
        data$y_train,
        data$coords,
        k = k,
        fold = fold,
        kernel = config$kernel,
        lambda = 0,
        alpha = 1,
        d = 0,
        standardize = TRUE,
        metric = "rmse",
        refit = FALSE,
        control = control
      )
    })
    data.frame(
      k = as.integer(k),
      rmse = as.double(cv$best$rmse),
      convergence_rate = as.double(cv$best$convergence_rate),
      elapsed_seconds = unname(timing[["elapsed"]]),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, results)
  best <- order(table$rmse, table$k)[[1L]]
  list(k = table$k[[best]], best = table[best, , drop = FALSE], table = table)
}

gwrs_penalized_study_oracle_union <- function(data,
                                               neighbors,
                                               config,
                                               control) {
  active <- which(colSums(data$support) > 0L)
  n <- nrow(data$x)
  p <- ncol(data$x)
  if (length(active)) {
    raw <- gwr_fit(
      data$x[, active, drop = FALSE],
      data$y_train,
      neighbors = neighbors,
      kernel = config$kernel,
      standardize = TRUE,
      control = control
    )
    coefficients <- matrix(0, nrow = n, ncol = p + 1L)
    coefficients[, 1L] <- raw$coefficients[, 1L]
    coefficients[, active + 1L] <- raw$coefficients[, -1L, drop = FALSE]
  } else {
    raw <- gwr_fit(
      matrix(0, nrow = n, ncol = 1L),
      data$y_train,
      neighbors = neighbors,
      kernel = config$kernel,
      standardize = FALSE,
      control = control
    )
    coefficients <- matrix(0, nrow = n, ncol = p + 1L)
    coefficients[, 1L] <- raw$coefficients[, 1L]
  }
  colnames(coefficients) <- c("(Intercept)", colnames(data$x))
  list(
    method = "Oracle-union GWR",
    coefficients = coefficients,
    fitted.values = raw$fitted.values,
    residuals = raw$residuals,
    diagnostics = raw$diagnostics,
    lambda = 0,
    alpha = NA_real_,
    gamma = NA_real_,
    oracle_active_union = active,
    prediction_fit = raw
  )
}

gwrs_penalized_study_subset_neighbors <- function(neighbors, columns) {
  columns <- as.integer(columns)
  output <- neighbors
  output$index <- neighbors$index[, columns, drop = FALSE]
  output$distance <- neighbors$distance[, columns, drop = FALSE]
  output$query_coords <- neighbors$query_coords[columns, , drop = FALSE]
  output$n_target <- length(columns)
  output$include_self <- FALSE
  class(output) <- "gwrs_neighbors"
  output
}

gwrs_penalized_study_cross_oracle <- function(x_train,
                                               y_train,
                                               x_target,
                                               neighbors,
                                               active,
                                               kernel,
                                               control) {
  active <- as.integer(active)
  if (length(active)) {
    train_design <- x_train[, active, drop = FALSE]
    target_design <- x_target[, active, drop = FALSE]
  } else {
    train_design <- matrix(0, nrow = nrow(x_train), ncol = 1L)
    target_design <- matrix(0, nrow = nrow(x_target), ncol = 1L)
  }
  run_cross <- getFromNamespace("run_cross_gwr_sl_path", "gwrs")
  kernel_code <- getFromNamespace("kernel_code", "gwrs")
  backtransform <- getFromNamespace("backtransform_coefficient_path", "gwrs")
  raw <- run_cross(
    train_design,
    y_train,
    target_design,
    neighbors,
    kernel_code(kernel),
    list(adaptive = TRUE, fixed = 1),
    lambda = 0,
    alpha = 1,
    d = 0,
    screening = FALSE,
    keep_coefficients = TRUE,
    control = control
  )
  coefficient_array <- backtransform(
    raw$coefficients,
    rep.int(0, ncol(train_design)),
    rep.int(1, ncol(train_design))
  )
  coefficients <- matrix(
    coefficient_array[, , 1L],
    nrow = nrow(x_target),
    ncol = ncol(train_design) + 1L
  )
  list(
    predictions = as.double(matrix(
      raw$predictions, nrow = nrow(x_target)
    )[, 1L]),
    coefficients = coefficients,
    converged = as.logical(matrix(
      raw$converged, nrow = nrow(x_target)
    )[, 1L]),
    iterations = as.integer(matrix(
      raw$iterations, nrow = nrow(x_target)
    )[, 1L]),
    kkt = as.double(matrix(raw$kkt, nrow = nrow(x_target))[, 1L])
  )
}

gwrs_penalized_study_oracle_local_set <- function(data,
                                                   x_target,
                                                   support,
                                                   neighbors,
                                                   config,
                                                   control) {
  n_target <- nrow(x_target)
  p <- ncol(data$x)
  coefficients <- matrix(0, nrow = n_target, ncol = p + 1L)
  predictions <- numeric(n_target)
  converged <- logical(n_target)
  iterations <- integer(n_target)
  kkt <- numeric(n_target)
  keys <- apply(support, 1L, function(value) {
    paste0(as.integer(value), collapse = "")
  })
  groups <- split(seq_len(n_target), factor(keys, levels = unique(keys)))
  for (indices in groups) {
    active <- which(support[indices[[1L]], ])
    result <- gwrs_penalized_study_cross_oracle(
      data$x,
      data$y_train,
      x_target[indices, , drop = FALSE],
      gwrs_penalized_study_subset_neighbors(neighbors, indices),
      active,
      config$kernel,
      control
    )
    predictions[indices] <- result$predictions
    coefficients[indices, 1L] <- result$coefficients[, 1L]
    if (length(active)) {
      coefficients[indices, active + 1L] <-
        result$coefficients[, -1L, drop = FALSE]
    }
    converged[indices] <- result$converged
    iterations[indices] <- result$iterations
    kkt[indices] <- result$kkt
  }
  colnames(coefficients) <- c("(Intercept)", colnames(data$x))
  list(
    fit = predictions,
    coefficients = coefficients,
    converged = converged,
    iterations = iterations,
    kkt = kkt,
    support_patterns = length(groups)
  )
}

gwrs_penalized_study_oracle_local <- function(data,
                                               neighbors,
                                               target_neighbors,
                                               config,
                                               control) {
  training <- gwrs_penalized_study_oracle_local_set(
    data, data$x, data$support, neighbors, config, control
  )
  target <- gwrs_penalized_study_oracle_local_set(
    data, data$x_target, data$target_support, target_neighbors, config, control
  )
  diagnostics <- data.frame(
    location = seq_len(nrow(data$x)),
    converged = training$converged,
    iterations = training$iterations,
    max_kkt_violation = training$kkt,
    stringsAsFactors = FALSE
  )
  list(
    method = "Oracle-local GWR",
    coefficients = training$coefficients,
    fitted.values = training$fit,
    residuals = data$y_train - training$fit,
    diagnostics = diagnostics,
    lambda = 0,
    alpha = NA_real_,
    gamma = NA_real_,
    target_prediction = target,
    oracle_support_patterns = training$support_patterns,
    oracle_target_support_patterns = target$support_patterns
  )
}

gwrs_penalized_study_target_prediction <- function(fit,
                                                    method,
                                                    data,
                                                    target_neighbors) {
  if (identical(method, "oracle_local_gwr")) {
    return(list(
      fit = fit$target_prediction$fit,
      coefficients = fit$target_prediction$coefficients,
      convergence_rate = mean(fit$target_prediction$converged)
    ))
  }
  if (identical(method, "oracle_union_gwr")) {
    active <- fit$oracle_active_union
    newdata <- if (length(active)) {
      data$x_target[, active, drop = FALSE]
    } else {
      matrix(0, nrow = nrow(data$x_target), ncol = 1L)
    }
    raw <- predict(
      fit$prediction_fit,
      newdata = newdata,
      neighbors = target_neighbors,
      type = "both"
    )
    coefficients <- matrix(
      0, nrow = nrow(data$x_target), ncol = ncol(data$x) + 1L
    )
    coefficients[, 1L] <- raw$coefficients[, 1L]
    if (length(active)) {
      coefficients[, active + 1L] <- raw$coefficients[, -1L, drop = FALSE]
    }
    colnames(coefficients) <- c("(Intercept)", colnames(data$x))
    return(list(fit = raw$fit, coefficients = coefficients,
                convergence_rate = NA_real_))
  }
  raw <- predict(
    fit,
    newdata = data$x_target,
    neighbors = target_neighbors,
    type = "both"
  )
  list(fit = raw$fit, coefficients = raw$coefficients,
       convergence_rate = NA_real_)
}

gwrs_penalized_study_rmse <- function(error) {
  if (!length(error)) return(NA_real_)
  sqrt(mean(error^2))
}

gwrs_penalized_study_scalar <- function(value, default = NA_real_) {
  if (is.null(value) || !length(value)) return(default)
  as.double(value[[1L]])
}

gwrs_penalized_study_fit_method <- function(method,
                                             data,
                                             neighbors,
                                             target_neighbors,
                                             chosen_k,
                                             fold,
                                             config,
                                             control,
                                             heartbeat = function(...) NULL) {
  heartbeat(paste0("tune-fit-", method))
  tuning_seconds <- 0
  fit_seconds <- 0
  tuning <- NULL
  cv_rmse <- NA_real_
  cv_convergence_rate <- NA_real_
  lambda_max <- NA_real_

  if (identical(method, "gwr")) {
    timing <- system.time({
      fit <- gwr_fit(
        data$x, data$y_train, neighbors = neighbors,
        kernel = config$kernel, standardize = TRUE, control = control
      )
    })
    fit_seconds <- unname(timing[["elapsed"]])
  } else if (identical(method, "oracle_union_gwr")) {
    timing <- system.time({
      fit <- gwrs_penalized_study_oracle_union(
        data, neighbors, config, control
      )
    })
    fit_seconds <- unname(timing[["elapsed"]])
  } else if (identical(method, "oracle_local_gwr")) {
    timing <- system.time({
      fit <- gwrs_penalized_study_oracle_local(
        data, neighbors, target_neighbors, config, control
      )
    })
    fit_seconds <- unname(timing[["elapsed"]])
  } else {
    cv_arguments <- list(
      x = data$x,
      y = data$y_train,
      coords = data$coords,
      k = chosen_k,
      fold = fold,
      kernel = config$kernel,
      standardize = TRUE,
      metric = "rmse",
      refit = FALSE,
      control = control,
      penalty = method
    )
    if (identical(method, "ridge")) {
      cv_arguments$lambda <- config$ridge_lambda
    } else {
      cv_arguments$n_lambda <- config$n_lambda
      cv_arguments$lambda_min_ratio <- config$lambda_min_ratio
    }
    if (identical(method, "elastic_net")) {
      cv_arguments$alpha <- config$en_alpha
    }
    if (identical(method, "scad")) {
      cv_arguments$gamma <- config$scad_gamma
    }
    if (identical(method, "mcp")) {
      cv_arguments$gamma <- config$mcp_gamma
    }
    timing <- system.time({
      cv <- do.call(cv_gwr_penalized, cv_arguments)
    })
    tuning_seconds <- unname(timing[["elapsed"]])
    cv_rmse <- as.double(cv$best$rmse)
    cv_convergence_rate <- as.double(cv$best$convergence_rate)
    lambda_max <- gwrs_penalized_study_scalar(cv$lambda_max)
    fit_arguments <- list(
      x = data$x,
      y = data$y_train,
      neighbors = neighbors,
      kernel = config$kernel,
      lambda = as.double(cv$best$lambda),
      penalty = method,
      standardize = TRUE,
      control = control
    )
    if (identical(method, "elastic_net")) {
      fit_arguments$alpha <- config$en_alpha
    }
    if (identical(method, "scad")) {
      fit_arguments$gamma <- config$scad_gamma
    }
    if (identical(method, "mcp")) {
      fit_arguments$gamma <- config$mcp_gamma
    }
    fit_timing <- system.time({
      fit <- do.call(gwr_penalized_fit, fit_arguments)
    })
    fit_seconds <- unname(fit_timing[["elapsed"]])
    tuning <- list(
      best = cv$best,
      cv = cv$cv,
      fold_results = cv$fold_results
    )
  }
  list(
    fit = fit,
    tuning = tuning,
    tuning_seconds = tuning_seconds,
    fit_seconds = fit_seconds,
    cv_rmse = cv_rmse,
    cv_convergence_rate = cv_convergence_rate,
    lambda_max = lambda_max
  )
}

gwrs_penalized_study_summary <- function(method,
                                          fitted,
                                          data,
                                          task,
                                          config,
                                          chosen_k,
                                          bandwidth_cv,
                                          neighbor_seconds,
                                          target_neighbors,
                                          target_neighbor_seconds) {
  fit <- fitted$fit
  coefficients <- fit$coefficients
  estimated_intercept <- coefficients[, 1L]
  estimated_slopes <- coefficients[, -1L, drop = FALSE]
  active_truth <- data$support
  active_error <- estimated_slopes[active_truth] - data$beta[active_truth]
  inactive_error <- estimated_slopes[!active_truth] - data$beta[!active_truth]
  core_active <- active_truth & !data$support_boundary
  domain_edge <- matrix(data$domain_edge, nrow = nrow(data$x),
                        ncol = ncol(data$x))
  domain_edge_active <- active_truth & domain_edge
  domain_interior_active <- active_truth & !domain_edge
  selection <- gwrs_selection_metrics(
    estimated_slopes,
    data$beta,
    selection_tolerance = config$selection_tolerance,
    truth_tolerance = config$truth_tolerance,
    truth_support = data$support,
    evaluate = method %in% config$selection_methods
  )
  prediction_timing <- system.time({
    target_prediction <- gwrs_penalized_study_target_prediction(
      fit, method, data, target_neighbors
    )
  })
  target_coefficients <- target_prediction$coefficients
  target_slopes <- target_coefficients[, -1L, drop = FALSE]
  target_active_error <- target_slopes[data$target_support] -
    data$beta_target[data$target_support]
  target_inactive_error <- target_slopes[!data$target_support] -
    data$beta_target[!data$target_support]
  converged <- as.logical(fit$diagnostics$converged)
  stationarity <- if (
    "max_stationarity_violation" %in% names(fit$diagnostics)
  ) max(fit$diagnostics$max_stationarity_violation) else NA_real_
  coordinate_gap <- if (
    "max_coordinate_gap" %in% names(fit$diagnostics)
  ) max(fit$diagnostics$max_coordinate_gap) else NA_real_
  bandwidth_selected_k <- if (!is.null(bandwidth_cv$k)) {
    bandwidth_cv$k[[1L]]
  } else {
    bandwidth_cv$best$k[[1L]]
  }
  train_noise_variance <- stats::var(data$y_train - data$signal)
  test_noise_variance <- stats::var(data$y_test - data$signal)
  target_noise_variance <- stats::var(data$y_target - data$target_signal)
  fitted_signal_error <- fit$fitted.values - data$signal

  data.frame(
    cell_id = task$cell_id[[1L]],
    section = task$section[[1L]],
    scenario_id = task$scenario_id[[1L]],
    design_variant = task$design_variant[[1L]],
    sampling_design = task$sampling_design[[1L]],
    replication = task$replication[[1L]],
    method = method,
    data_seed = task$data_seed[[1L]],
    fold_seed = task$fold_seed[[1L]],
    n = task$n[[1L]],
    p = task$p[[1L]],
    predictor_rho = task$predictor_rho[[1L]],
    spatial_rho = task$spatial_rho[[1L]],
    local_correlation = task$local_correlation[[1L]],
    target_grid_side = task$target_grid_side[[1L]],
    n_target = nrow(data$x_target),
    target_snr = task$snr[[1L]],
    realized_signal_variance = data$signal_sd^2,
    noise_sd = data$noise_sd,
    design_snr = data$signal_sd^2 / data$noise_sd^2,
    realized_train_noise_variance = train_noise_variance,
    realized_test_noise_variance = test_noise_variance,
    realized_target_noise_variance = target_noise_variance,
    realized_train_snr = data$signal_sd^2 / train_noise_variance,
    realized_test_snr = data$signal_sd^2 / test_noise_variance,
    realized_target_snr = stats::var(data$target_signal) /
      target_noise_variance,
    realized_predictor_correlation = data$realized_predictor_correlation,
    realized_local_correlation = data$realized_local_correlation,
    realized_outside_correlation = data$realized_outside_correlation,
    realized_spatial_neighbor_product =
      data$realized_spatial_neighbor_product,
    bandwidth_selected_k = bandwidth_selected_k,
    chosen_k = chosen_k,
    bandwidth_cv_rmse = bandwidth_cv$best$rmse,
    lambda = gwrs_penalized_study_scalar(fit$lambda, default = 0),
    lambda_max = fitted$lambda_max,
    alpha = if (identical(method, "elastic_net")) config$en_alpha else
      if (identical(method, "lasso")) 1 else
        if (identical(method, "ridge")) 0 else NA_real_,
    gamma = if (identical(method, "scad")) config$scad_gamma else
      if (identical(method, "mcp")) config$mcp_gamma else NA_real_,
    cv_rmse = fitted$cv_rmse,
    cv_convergence_rate = fitted$cv_convergence_rate,
    all_converged = all(converged),
    convergence_rate = mean(converged),
    train_response_rmse = gwrs_penalized_study_rmse(fit$residuals),
    test_response_rmse = gwrs_penalized_study_rmse(
      fit$fitted.values - data$y_test
    ),
    signal_rmse = gwrs_penalized_study_rmse(
      fitted_signal_error
    ),
    new_location_response_rmse = gwrs_penalized_study_rmse(
      target_prediction$fit - data$y_target
    ),
    new_location_signal_rmse = gwrs_penalized_study_rmse(
      target_prediction$fit - data$target_signal
    ),
    intercept_rmse = gwrs_penalized_study_rmse(
      estimated_intercept - data$intercept
    ),
    coefficient_rmse = gwrs_penalized_study_rmse(
      estimated_slopes - data$beta
    ),
    active_coefficient_rmse = gwrs_penalized_study_rmse(active_error),
    inactive_coefficient_rmse = gwrs_penalized_study_rmse(inactive_error),
    core_active_coefficient_rmse = gwrs_penalized_study_rmse(
      estimated_slopes[core_active] - data$beta[core_active]
    ),
    support_boundary_coefficient_rmse = gwrs_penalized_study_rmse(
      estimated_slopes[data$support_boundary] -
        data$beta[data$support_boundary]
    ),
    active_attenuation_bias = if (length(active_error)) {
      mean(active_error * sign(data$beta[active_truth]))
    } else NA_real_,
    domain_edge_active_coefficient_rmse = gwrs_penalized_study_rmse(
      estimated_slopes[domain_edge_active] - data$beta[domain_edge_active]
    ),
    domain_interior_active_coefficient_rmse = gwrs_penalized_study_rmse(
      estimated_slopes[domain_interior_active] -
        data$beta[domain_interior_active]
    ),
    domain_edge_signal_rmse = gwrs_penalized_study_rmse(
      fitted_signal_error[data$domain_edge]
    ),
    domain_interior_signal_rmse = gwrs_penalized_study_rmse(
      fitted_signal_error[!data$domain_edge]
    ),
    new_location_coefficient_rmse = gwrs_penalized_study_rmse(
      target_slopes - data$beta_target
    ),
    new_location_active_coefficient_rmse = gwrs_penalized_study_rmse(
      target_active_error
    ),
    new_location_inactive_coefficient_rmse = gwrs_penalized_study_rmse(
      target_inactive_error
    ),
    selection_evaluable = selection$selection_evaluable,
    true_positive = selection$true_positive,
    false_positive = selection$false_positive,
    false_negative = selection$false_negative,
    true_negative = selection$true_negative,
    true_positive_rate = selection$true_positive_rate,
    false_positive_rate = selection$false_positive_rate,
    specificity = selection$specificity,
    precision = selection$precision,
    false_discovery_rate = selection$false_discovery_rate,
    f1 = selection$f1,
    mcc = selection$mcc,
    support_iou = selection$support_iou,
    selection_accuracy = selection$selection_accuracy,
    mean_selected_predictors = selection$mean_selected_predictors,
    maximum_stationarity_violation = stationarity,
    maximum_coordinate_gap = coordinate_gap,
    target_convergence_rate = target_prediction$convergence_rate,
    bandwidth_elapsed_seconds = sum(bandwidth_cv$table$elapsed_seconds),
    neighbor_elapsed_seconds = neighbor_seconds,
    target_neighbor_elapsed_seconds = target_neighbor_seconds,
    tuning_elapsed_seconds = fitted$tuning_seconds,
    fit_elapsed_seconds = fitted$fit_seconds,
    prediction_elapsed_seconds = unname(prediction_timing[["elapsed"]]),
    retained_fit_mib = as.numeric(object.size(fit)) / 1024^2,
    signal_sum = sum(data$signal),
    signal_sum_squares = sum(data$signal^2),
    x_sum = sum(data$x),
    y_train_sum = sum(data$y_train),
    y_test_sum = sum(data$y_test),
    target_signal_sum = sum(data$target_signal),
    target_signal_sum_squares = sum(data$target_signal^2),
    x_target_sum = sum(data$x_target),
    y_target_sum = sum(data$y_target),
    support_count = sum(data$support),
    target_support_count = sum(data$target_support),
    oracle_scope = if (identical(method, "oracle_union_gwr")) {
      "global_union_truth"
    } else if (identical(method, "oracle_local_gwr")) {
      "location_specific_truth"
    } else NA_character_,
    stringsAsFactors = FALSE
  )
}
