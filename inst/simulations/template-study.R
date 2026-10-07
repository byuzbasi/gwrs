# A small, resumable gwrs simulation workflow template.
# This is a software example, not a scientific Monte Carlo design.

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
if (length(file_argument) != 1L) {
  stop("Run this template with Rscript.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", file_argument))
source(file.path(dirname(script_path), "resume-utils.R"))

read_integer <- function(name, default, minimum = 1L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(raw))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop("`", name, "` must be an integer >= ", minimum, ".",
         call. = FALSE)
  }
  value
}

read_number <- function(name, default) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.double(raw))
  if (length(value) != 1L || !is.finite(value)) {
    stop("`", name, "` must be one finite number.", call. = FALSE)
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

output_directory <- Sys.getenv("GWRS_SIM_OUTPUT")
if (!nzchar(output_directory)) {
  stop("Set `GWRS_SIM_OUTPUT` to a new, versioned directory.",
       call. = FALSE)
}
action <- match.arg(
  Sys.getenv("GWRS_SIM_ACTION", unset = "run"),
  c("run", "status")
)
if (identical(action, "status")) {
  gwrs_sim_status(output_directory)
  quit(save = "no", status = 0L)
}

library(gwrs)

mode <- match.arg(
  Sys.getenv("GWRS_SIM_MODE", unset = "smoke"),
  c("smoke", "standard")
)
defaults <- if (identical(mode, "smoke")) {
  list(n = 300L, p = 6L, k = 40L, replications = 3L)
} else {
  list(n = 2500L, p = 12L, k = 60L, replications = 20L)
}
n <- read_integer("GWRS_SIM_N", defaults$n)
p <- read_integer("GWRS_SIM_P", defaults$p)
k <- read_integer("GWRS_SIM_K", defaults$k)
replications <- read_integer(
  "GWRS_SIM_REPLICATIONS", defaults$replications
)
if (k > n) stop("`GWRS_SIM_K` cannot exceed `GWRS_SIM_N`.", call. = FALSE)

master_seed <- read_integer("GWRS_SIM_SEED", 20260826L, minimum = 0L)
n_threads <- read_integer("GWRS_SIM_THREADS", 1L)
lambda <- read_number("GWRS_SIM_LAMBDA", 0.06)
alpha <- read_number("GWRS_SIM_ALPHA", 0.7)
if (lambda < 0 || alpha < 0 || alpha > 1) {
  stop("Lambda must be nonnegative and alpha must be in [0, 1].",
       call. = FALSE)
}
keep_full_fit <- read_flag("GWRS_SIM_KEEP_FULL_FIT")
confirm_large <- read_flag("GWRS_SIM_CONFIRM_LARGE")
retry_failed <- read_flag("GWRS_SIM_RETRY_FAILED", default = TRUE)
recover_lock <- read_flag("GWRS_SIM_RECOVER_LOCK")
max_tasks_raw <- Sys.getenv("GWRS_SIM_MAX_TASKS")
max_tasks <- if (nzchar(max_tasks_raw)) {
  read_integer("GWRS_SIM_MAX_TASKS", 1L)
} else {
  Inf
}

config <- list(
  schema = gwrs_sim_schema,
  study = "elastic-net-workflow-template",
  mode = mode,
  package_version = as.character(utils::packageVersion("gwrs")),
  n = n,
  p = p,
  k = k,
  replications = replications,
  master_seed = master_seed,
  penalty = "elastic_net",
  lambda = lambda,
  alpha = alpha,
  d = 0.5,
  noise_sd = 0.25,
  n_threads = n_threads,
  keep_full_fit = keep_full_fit
)
plan <- gwrs_sim_make_plan(
  replications = replications,
  master_seed = master_seed,
  scenario_id = "elastic-net"
)
runtime_signature <- gwrs_sim_runtime_signature("gwrs")
workflow_files <- c(
  resume_utils = file.path(dirname(script_path), "resume-utils.R"),
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

task_function <- function(task, config, heartbeat) {
  heartbeat("simulate-data")
  coords <- cbind(east = runif(config$n), north = runif(config$n))
  x <- matrix(
    rnorm(config$n * config$p), nrow = config$n, ncol = config$p,
    dimnames = list(NULL, paste0("x", seq_len(config$p)))
  )
  beta <- numeric(config$p)
  beta[seq_len(min(3L, config$p))] <- c(1.2, -0.8, 0.45)[
    seq_len(min(3L, config$p))
  ]
  signal <- 0.7 + as.double(x %*% beta)
  y <- signal + rnorm(config$n, sd = config$noise_sd)

  heartbeat("build-neighbours")
  neighbors <- gwr_neighbors(coords, k = config$k)
  heartbeat("fit-model")
  fit <- gwr_penalized_fit(
    x, y,
    neighbors = neighbors,
    kernel = "gaussian",
    lambda = config$lambda,
    penalty = config$penalty,
    alpha = config$alpha,
    d = config$d,
    control = gwrs_control(
      n_threads = config$n_threads,
      keep_data = FALSE,
      diagnostics = "none"
    )
  )

  heartbeat("summarize")
  estimated_slopes <- fit$coefficients[, -1L, drop = FALSE]
  truth <- matrix(
    beta, nrow = config$n, ncol = config$p, byrow = TRUE
  )
  selected <- abs(estimated_slopes) > 1e-8
  active_truth <- matrix(
    abs(beta) > 0, nrow = config$n, ncol = config$p, byrow = TRUE
  )
  summary <- data.frame(
    scenario_id = task$scenario_id[[1L]],
    replication = task$replication[[1L]],
    n = config$n,
    p = config$p,
    k = config$k,
    all_converged = all(fit$diagnostics$converged),
    response_rmse = sqrt(mean(fit$residuals^2)),
    signal_rmse = sqrt(mean((fit$fitted.values - signal)^2)),
    coefficient_rmse = sqrt(mean((estimated_slopes - truth)^2)),
    selection_accuracy = mean(selected == active_truth),
    mean_selected_predictors = mean(rowSums(selected)),
    retained_fit_mib = as.numeric(object.size(fit)) / 1024^2,
    stringsAsFactors = FALSE
  )
  result <- list(summary = summary)
  if (isTRUE(config$keep_full_fit)) result$fit <- fit
  result
}

validate_result <- function(result, task) {
  is.list(result) && is.data.frame(result$summary) &&
    nrow(result$summary) == 1L &&
    identical(result$summary$replication, task$replication) &&
    isTRUE(result$summary$all_converged) &&
    all(is.finite(unlist(result$summary[c(
      "response_rmse", "signal_rmse", "coefficient_rmse",
      "selection_accuracy", "mean_selected_predictors"
    )])))
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
  confirm_large = confirm_large
)
print(status$progress, row.names = FALSE)
cat("Output:", status$output_directory, "\n")
