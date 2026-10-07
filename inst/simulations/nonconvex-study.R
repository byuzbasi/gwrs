# Resumable recovery study for GWR, GWR-Lasso, GWR-SCAD, and GWR-MCP.
# This validates software and estimator behavior at an explicit lambda-anchor
# fraction. It is not a substitute for an approved paper-level tuning design.

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
if (length(file_argument) != 1L) {
  stop("Run this study with Rscript.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", file_argument))
source(file.path(dirname(script_path), "resume-utils.R"))
source(file.path(dirname(script_path), "selection-metrics.R"))

read_integer <- function(name, default, minimum = 1L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(raw))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop("`", name, "` must be an integer >= ", minimum, ".",
         call. = FALSE)
  }
  value
}

read_number <- function(name, default, lower = -Inf, upper = Inf) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.double(raw))
  if (length(value) != 1L || !is.finite(value) ||
      value < lower || value > upper) {
    stop(
      "`", name, "` must be finite and in [", lower, ", ", upper, "].",
      call. = FALSE
    )
  }
  value
}

read_flag <- function(name, default = FALSE) {
  raw <- Sys.getenv(name, unset = if (default) "YES" else "NO")
  if (!raw %in% c("YES", "NO")) {
    stop("`", name, "` must be YES or NO.", call. = FALSE)
  }
  identical(raw, "YES")
}

read_choices <- function(name, default, choices) {
  raw <- Sys.getenv(name, unset = paste(default, collapse = ","))
  values <- unique(trimws(strsplit(raw, ",", fixed = TRUE)[[1L]]))
  if (!length(values) || any(!nzchar(values)) ||
      any(!values %in% choices)) {
    stop(
      "`", name, "` must be a comma-separated subset of: ",
      paste(choices, collapse = ", "), ".", call. = FALSE
    )
  }
  values
}

output_directory <- Sys.getenv("GWRS_NONCONVEX_SIM_OUTPUT")
if (!nzchar(output_directory)) {
  stop(
    "Set `GWRS_NONCONVEX_SIM_OUTPUT` to a new, versioned directory.",
    call. = FALSE
  )
}
action <- match.arg(
  Sys.getenv("GWRS_NONCONVEX_SIM_ACTION", unset = "run"),
  c("run", "status")
)
if (identical(action, "status")) {
  gwrs_sim_status(output_directory)
  quit(save = "no", status = 0L)
}

library(gwrs)

mode <- match.arg(
  Sys.getenv("GWRS_NONCONVEX_SIM_MODE", unset = "smoke"),
  c("smoke", "standard", "production")
)
defaults <- switch(
  mode,
  smoke = list(n = 240L, p = 6L, k = 40L, replications = 2L),
  standard = list(n = 1500L, p = 12L, k = 60L, replications = 10L),
  production = list(n = 10000L, p = 30L, k = 100L, replications = 50L)
)
n <- read_integer("GWRS_NONCONVEX_SIM_N", defaults$n)
p <- read_integer("GWRS_NONCONVEX_SIM_P", defaults$p, minimum = 5L)
k <- read_integer("GWRS_NONCONVEX_SIM_K", defaults$k)
replications <- read_integer(
  "GWRS_NONCONVEX_SIM_REPLICATIONS", defaults$replications
)
if (p >= k) {
  stop("Use `k > p` so local rank does not dominate method recovery.",
       call. = FALSE)
}
if (k > n) {
  stop("`GWRS_NONCONVEX_SIM_K` cannot exceed N.", call. = FALSE)
}

scenarios <- read_choices(
  "GWRS_NONCONVEX_SIM_SCENARIOS",
  c("smooth_sparse", "regional_sparse"),
  c("smooth_sparse", "regional_sparse")
)
methods <- read_choices(
  "GWRS_NONCONVEX_SIM_METHODS",
  c("gwr", "lasso", "scad", "mcp"),
  c("gwr", "lasso", "scad", "mcp")
)
kernel <- match.arg(
  Sys.getenv("GWRS_NONCONVEX_SIM_KERNEL", unset = "gaussian"),
  c("gaussian", "bisquare", "exponential", "tricube", "boxcar")
)
master_seed <- read_integer(
  "GWRS_NONCONVEX_SIM_SEED", 20260826L, minimum = 0L
)
n_threads <- read_integer("GWRS_NONCONVEX_SIM_THREADS", 1L)
lambda_fraction <- read_number(
  "GWRS_NONCONVEX_SIM_LAMBDA_FRACTION", 0.20,
  lower = .Machine$double.eps, upper = 1
)
scad_gamma <- read_number(
  "GWRS_NONCONVEX_SIM_SCAD_GAMMA", 3.7,
  lower = 2 + .Machine$double.eps
)
mcp_gamma <- read_number(
  "GWRS_NONCONVEX_SIM_MCP_GAMMA", 3,
  lower = 1 + .Machine$double.eps
)
noise_sd <- read_number(
  "GWRS_NONCONVEX_SIM_NOISE_SD", 0.35,
  lower = .Machine$double.eps
)
predictor_rho <- read_number(
  "GWRS_NONCONVEX_SIM_PREDICTOR_RHO", 0.35,
  lower = -0.95, upper = 0.95
)
selection_tolerance <- read_number(
  "GWRS_NONCONVEX_SIM_SELECTION_TOLERANCE", 1e-8,
  lower = 0
)
truth_tolerance <- read_number(
  "GWRS_NONCONVEX_SIM_TRUTH_TOLERANCE", 1e-12,
  lower = 0
)
solver_tolerance <- read_number(
  "GWRS_NONCONVEX_SIM_SOLVER_TOLERANCE", 1e-7,
  lower = .Machine$double.eps
)
max_iterations <- read_integer(
  "GWRS_NONCONVEX_SIM_MAX_ITERATIONS", 2000L
)
grain_size <- read_integer("GWRS_NONCONVEX_SIM_GRAIN_SIZE", 16L)
keep_task_data <- read_flag("GWRS_NONCONVEX_SIM_KEEP_TASK_DATA")
confirm_large <- read_flag("GWRS_NONCONVEX_SIM_CONFIRM_LARGE")
retry_failed <- read_flag(
  "GWRS_NONCONVEX_SIM_RETRY_FAILED", default = TRUE
)
recover_lock <- read_flag("GWRS_NONCONVEX_SIM_RECOVER_LOCK")
max_tasks_raw <- Sys.getenv("GWRS_NONCONVEX_SIM_MAX_TASKS")
max_tasks <- if (nzchar(max_tasks_raw)) {
  read_integer("GWRS_NONCONVEX_SIM_MAX_TASKS", 1L)
} else {
  Inf
}
if (identical(mode, "production") && !confirm_large) {
  stop(
    "Production mode requires `GWRS_NONCONVEX_SIM_CONFIRM_LARGE=YES`.",
    call. = FALSE
  )
}

config <- list(
  schema = gwrs_sim_schema,
  study = "gwr-nonconvex-recovery-validation-v1",
  scope = "fixed-anchor-fraction-software-validation",
  mode = mode,
  package_version = as.character(utils::packageVersion("gwrs")),
  n = n,
  p = p,
  k = k,
  replications = replications,
  scenarios = scenarios,
  methods = methods,
  master_seed = master_seed,
  kernel = kernel,
  lambda_fraction = lambda_fraction,
  scad_gamma = scad_gamma,
  mcp_gamma = mcp_gamma,
  noise_sd = noise_sd,
  predictor_rho = predictor_rho,
  selection_tolerance = selection_tolerance,
  truth_tolerance = truth_tolerance,
  solver_tolerance = solver_tolerance,
  max_iterations = max_iterations,
  n_threads = n_threads,
  grain_size = grain_size,
  keep_task_data = keep_task_data
)

plan <- expand.grid(
  scenario_id = scenarios,
  replication = seq_len(replications),
  method = methods,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
plan$task_id <- sprintf(
  "%s-rep-%05d-%s", plan$scenario_id, plan$replication, plan$method
)
task_index <- seq_len(nrow(plan))
plan$seed <- as.integer(
  ((master_seed + task_index * 104729) %% 2147483646) + 1
)
data_key <- match(
  paste(plan$scenario_id, plan$replication, sep = ":"),
  unique(paste(plan$scenario_id, plan$replication, sep = ":"))
)
plan$data_seed <- as.integer(
  ((master_seed + data_key * 130363) %% 2147483646) + 1
)
plan <- plan[c(
  "task_id", "scenario_id", "replication", "method", "data_seed", "seed"
)]

runtime_signature <- gwrs_sim_runtime_signature("gwrs")
workflow_files <- c(
  resume_utils = file.path(dirname(script_path), "resume-utils.R"),
  selection_metrics = file.path(dirname(script_path), "selection-metrics.R"),
  study_script = script_path
)
runtime_signature <- rbind(
  runtime_signature,
  data.frame(
    field = paste0(names(workflow_files), "_sha256"),
    value = gwrs_sim_sha256(unname(workflow_files)),
    stringsAsFactors = FALSE
  )
)

simulate_spatial_data <- function(config, scenario, data_seed) {
  set.seed(data_seed)
  coords <- cbind(east = runif(config$n), north = runif(config$n))
  x <- matrix(0, nrow = config$n, ncol = config$p)
  x[, 1L] <- rnorm(config$n)
  innovation_scale <- sqrt(1 - config$predictor_rho^2)
  if (config$p > 1L) {
    for (column in 2:config$p) {
      x[, column] <- config$predictor_rho * x[, column - 1L] +
        innovation_scale * rnorm(config$n)
    }
  }
  colnames(x) <- paste0("x", seq_len(config$p))
  beta <- matrix(
    0, nrow = config$n, ncol = config$p,
    dimnames = list(NULL, colnames(x))
  )
  east <- coords[, 1L]
  north <- coords[, 2L]
  intercept <- 0.7 + 0.25 * cos(2 * pi * east) -
    0.20 * sin(2 * pi * north)

  if (identical(scenario, "smooth_sparse")) {
    beta[, 1L] <- 1.20 + 0.50 * sin(2 * pi * east)
    beta[, 2L] <- -1.00 + 0.40 * cos(2 * pi * north)
    beta[, 3L] <- 0.65 + 0.25 * sin(2 * pi * (east + north))
    beta[, 4L] <- 0.45 + 0.15 * (east - north)
  } else if (identical(scenario, "regional_sparse")) {
    beta[, 1L] <- 1.20 * (east <= 0.55)
    beta[, 2L] <- -1.00 * (east > 0.35 & north > 0.45)
    beta[, 3L] <- 0.80 * (north <= 0.50)
    beta[, 4L] <- 0.60 * ((east + north) > 1.20)
  } else {
    stop("Unknown simulation scenario.", call. = FALSE)
  }

  signal <- intercept + rowSums(x * beta)
  y <- signal + rnorm(config$n, sd = config$noise_sd)
  list(
    coords = coords,
    x = x,
    y = y,
    signal = signal,
    intercept = intercept,
    beta = beta
  )
}

task_function <- function(task, config, heartbeat) {
  scenario <- task$scenario_id[[1L]]
  method <- task$method[[1L]]
  data_seed <- task$data_seed[[1L]]
  heartbeat("simulate-data")
  data <- simulate_spatial_data(config, scenario, data_seed)

  heartbeat("build-neighbours")
  neighbor_time <- system.time({
    neighbors <- gwr_neighbors(data$coords, k = config$k)
  })[["elapsed"]]
  control <- gwrs_control(
    tolerance = config$solver_tolerance,
    max_iterations = config$max_iterations,
    n_threads = config$n_threads,
    grain_size = config$grain_size,
    keep_data = FALSE,
    diagnostics = "none"
  )

  lambda_max <- 0
  lambda <- 0
  anchor_time <- 0
  if (!identical(method, "gwr")) {
    heartbeat("compute-lambda-anchor")
    anchor_timing <- system.time({
      lambda_max <- switch(
        method,
        lasso = gwr_sl_lambda_max(
          data$x, data$y, neighbors = neighbors,
          kernel = config$kernel, alpha = 1, control = control
        ),
        scad = gwr_nonconvex_lambda_max(
          data$x, data$y, neighbors = neighbors,
          kernel = config$kernel, penalty = "scad",
          gamma = config$scad_gamma, control = control
        ),
        mcp = gwr_nonconvex_lambda_max(
          data$x, data$y, neighbors = neighbors,
          kernel = config$kernel, penalty = "mcp",
          gamma = config$mcp_gamma, control = control
        )
      )
    })
    anchor_time <- unname(anchor_timing[["elapsed"]])
    lambda <- config$lambda_fraction * lambda_max
  }

  heartbeat(paste0("fit-", method))
  fit_timing <- system.time({
    fit <- switch(
      method,
      gwr = gwr_fit(
        data$x, data$y, neighbors = neighbors,
        kernel = config$kernel, control = control
      ),
      lasso = gwr_lasso_fit(
        data$x, data$y, neighbors = neighbors,
        kernel = config$kernel, lambda = lambda, control = control
      ),
      scad = gwr_scad_fit(
        data$x, data$y, neighbors = neighbors,
        kernel = config$kernel, lambda = lambda,
        gamma = config$scad_gamma, control = control
      ),
      mcp = gwr_mcp_fit(
        data$x, data$y, neighbors = neighbors,
        kernel = config$kernel, lambda = lambda,
        gamma = config$mcp_gamma, control = control
      )
    )
  })
  fit_time <- unname(fit_timing[["elapsed"]])

  heartbeat("summarize-recovery")
  estimated_intercept <- fit$coefficients[, 1L]
  estimated_slopes <- fit$coefficients[, -1L, drop = FALSE]
  active_truth <- abs(data$beta) > config$truth_tolerance
  selection <- gwrs_selection_metrics(
    estimated_slopes,
    data$beta,
    selection_tolerance = config$selection_tolerance,
    truth_tolerance = config$truth_tolerance,
    evaluate = method %in% c("lasso", "scad", "mcp")
  )
  active_error <- estimated_slopes[active_truth] - data$beta[active_truth]
  inactive_error <- estimated_slopes[!active_truth] - data$beta[!active_truth]
  stationarity <- if (
    "max_stationarity_violation" %in% names(fit$diagnostics)
  ) {
    max(fit$diagnostics$max_stationarity_violation)
  } else {
    NA_real_
  }
  coordinate_gap <- if ("max_coordinate_gap" %in% names(fit$diagnostics)) {
    max(fit$diagnostics$max_coordinate_gap)
  } else {
    NA_real_
  }

  summary <- data.frame(
    scenario_id = scenario,
    replication = task$replication[[1L]],
    method = method,
    data_seed = data_seed,
    n = config$n,
    p = config$p,
    k = config$k,
    lambda_max = lambda_max,
    lambda = lambda,
    lambda_fraction = if (identical(method, "gwr")) NA_real_ else
      config$lambda_fraction,
    gamma = switch(
      method, scad = config$scad_gamma, mcp = config$mcp_gamma, NA_real_
    ),
    all_converged = all(fit$diagnostics$converged),
    response_rmse = sqrt(mean(fit$residuals^2)),
    signal_rmse = sqrt(mean((fit$fitted.values - data$signal)^2)),
    intercept_rmse = sqrt(mean((estimated_intercept - data$intercept)^2)),
    coefficient_rmse = sqrt(mean((estimated_slopes - data$beta)^2)),
    active_coefficient_rmse = sqrt(mean(active_error^2)),
    inactive_coefficient_rmse = sqrt(mean(inactive_error^2)),
    selection_evaluable = selection$selection_evaluable,
    true_positive = selection$true_positive,
    false_positive = selection$false_positive,
    false_negative = selection$false_negative,
    true_negative = selection$true_negative,
    true_positive_rate = selection$true_positive_rate,
    false_positive_rate = selection$false_positive_rate,
    precision = selection$precision,
    f1 = selection$f1,
    mcc = selection$mcc,
    support_iou = selection$support_iou,
    selection_accuracy = selection$selection_accuracy,
    mean_selected_predictors = selection$mean_selected_predictors,
    maximum_stationarity_violation = stationarity,
    maximum_coordinate_gap = coordinate_gap,
    neighbor_elapsed_seconds = unname(neighbor_time),
    anchor_elapsed_seconds = anchor_time,
    fit_elapsed_seconds = fit_time,
    retained_fit_mib = as.numeric(object.size(fit)) / 1024^2,
    signal_sum = sum(data$signal),
    signal_sum_squares = sum(data$signal^2),
    stringsAsFactors = FALSE
  )
  result <- list(summary = summary)
  if (isTRUE(config$keep_task_data)) {
    result$task_data <- list(
      coords = data$coords,
      x = data$x,
      y = data$y,
      signal = data$signal,
      intercept = data$intercept,
      beta = data$beta,
      coefficients = fit$coefficients,
      fitted.values = fit$fitted.values
    )
  }
  result
}

validate_result <- function(result, task) {
  summary <- result$summary
  core_metrics <- c(
    "response_rmse", "signal_rmse", "intercept_rmse",
    "coefficient_rmse", "active_coefficient_rmse",
    "inactive_coefficient_rmse", "fit_elapsed_seconds",
    "signal_sum", "signal_sum_squares"
  )
  valid <- is.list(result) && is.data.frame(summary) && nrow(summary) == 1L &&
    identical(summary$scenario_id, task$scenario_id) &&
    identical(summary$replication, task$replication) &&
    identical(summary$method, task$method) &&
    identical(summary$data_seed, task$data_seed) &&
    isTRUE(summary$all_converged) &&
    all(is.finite(unlist(summary[core_metrics])))
  if (isTRUE(valid)) {
    selection_method <- task$method[[1L]] %in% c("lasso", "scad", "mcp")
    selection_counts <- c(
      "true_positive", "false_positive", "false_negative", "true_negative",
      "selection_accuracy", "mean_selected_predictors"
    )
    selection_ratios <- c(
      "true_positive_rate", "false_positive_rate", "precision", "f1",
      "mcc", "support_iou"
    )
    valid <- identical(summary$selection_evaluable, selection_method) &&
      if (selection_method) {
        all(is.finite(unlist(summary[selection_counts]))) &&
          all(vapply(
            summary[selection_ratios],
            function(value) is.finite(value) || is.na(value),
            logical(1L)
          ))
      } else {
        all(is.na(unlist(summary[c(selection_counts, selection_ratios)])))
      }
  }
  if (isTRUE(valid) && task$method[[1L]] %in% c("scad", "mcp")) {
    valid <- is.finite(summary$maximum_stationarity_violation) &&
      is.finite(summary$maximum_coordinate_gap)
  }
  isTRUE(valid)
}

status <- gwrs_sim_run(
  output_directory = output_directory,
  config = config,
  plan = plan,
  task_function = task_function,
  validate_result = validate_result,
  runtime_signature = runtime_signature,
  max_tasks = max_tasks,
  retry_failed = retry_failed,
  recover_lock = recover_lock,
  confirm_large = confirm_large,
  limits = list(n = 20000L, p = 100L, k = 200L, tasks = 500L)
)
print(status$progress, row.names = FALSE)
cat("Output:", status$output_directory, "\n")
