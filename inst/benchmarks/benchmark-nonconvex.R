# Reproducible and resumable benchmark for the GWR-SCAD/MCP path core.
# The smoke default validates the protocol; production execution is external.

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
if (length(file_argument) != 1L) {
  stop("Run this benchmark with Rscript.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", file_argument))
source(file.path(dirname(script_path), "..", "simulations", "resume-utils.R"))

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

output_directory <- Sys.getenv("GWRS_NONCONVEX_BENCH_OUTPUT")
if (!nzchar(output_directory)) {
  stop(
    "Set `GWRS_NONCONVEX_BENCH_OUTPUT` to a new, versioned directory.",
    call. = FALSE
  )
}
action <- match.arg(
  Sys.getenv("GWRS_NONCONVEX_BENCH_ACTION", unset = "run"),
  c("run", "status")
)
if (identical(action, "status")) {
  gwrs_sim_status(output_directory)
  quit(save = "no", status = 0L)
}

library(gwrs)

mode <- match.arg(
  Sys.getenv("GWRS_NONCONVEX_BENCH_MODE", unset = "smoke"),
  c("smoke", "standard", "production")
)
defaults <- switch(
  mode,
  smoke = list(n = 240L, p = 6L, k = 40L, n_lambda = 4L),
  standard = list(n = 5000L, p = 20L, k = 80L, n_lambda = 15L),
  production = list(n = 25000L, p = 30L, k = 120L, n_lambda = 25L)
)
n <- read_integer("GWRS_NONCONVEX_BENCH_N", defaults$n)
p <- read_integer("GWRS_NONCONVEX_BENCH_P", defaults$p)
k <- read_integer("GWRS_NONCONVEX_BENCH_K", defaults$k)
n_lambda <- read_integer(
  "GWRS_NONCONVEX_BENCH_N_LAMBDA", defaults$n_lambda, minimum = 2L
)
if (p >= k) {
  stop("Use `k > p` so local rank does not dominate the benchmark.",
       call. = FALSE)
}
if (k > n) {
  stop("`GWRS_NONCONVEX_BENCH_K` cannot exceed N.", call. = FALSE)
}

master_seed <- read_integer(
  "GWRS_NONCONVEX_BENCH_SEED", 20260826L, minimum = 0L
)
lambda_min_ratio <- read_number(
  "GWRS_NONCONVEX_BENCH_LAMBDA_MIN_RATIO", 0.05,
  lower = .Machine$double.eps, upper = 1
)
tolerance <- read_number(
  "GWRS_NONCONVEX_BENCH_TOLERANCE", 1e-7,
  lower = .Machine$double.eps
)
max_iterations <- read_integer(
  "GWRS_NONCONVEX_BENCH_MAX_ITERATIONS", 2000L
)
grain_size <- read_integer("GWRS_NONCONVEX_BENCH_GRAIN_SIZE", 16L)
confirm_large <- read_flag("GWRS_NONCONVEX_BENCH_CONFIRM_LARGE")
recover_lock <- read_flag("GWRS_NONCONVEX_BENCH_RECOVER_LOCK")
retry_failed <- read_flag(
  "GWRS_NONCONVEX_BENCH_RETRY_FAILED", default = TRUE
)
max_tasks_raw <- Sys.getenv("GWRS_NONCONVEX_BENCH_MAX_TASKS")
max_tasks <- if (nzchar(max_tasks_raw)) {
  read_integer("GWRS_NONCONVEX_BENCH_MAX_TASKS", 1L)
} else {
  Inf
}
if (identical(mode, "production") && !confirm_large) {
  stop(
    "Production mode requires `GWRS_NONCONVEX_BENCH_CONFIRM_LARGE=YES`.",
    call. = FALSE
  )
}

config <- list(
  schema = gwrs_sim_schema,
  study = "gwr-nonconvex-core-benchmark-v1",
  mode = mode,
  package_version = as.character(utils::packageVersion("gwrs")),
  n = n,
  p = p,
  k = k,
  n_lambda = n_lambda,
  lambda_min_ratio = lambda_min_ratio,
  master_seed = master_seed,
  data_seed = master_seed,
  penalties = c("scad", "mcp"),
  scad_gamma = 3.7,
  mcp_gamma = 3,
  tolerance = tolerance,
  max_iterations = max_iterations,
  grain_size = grain_size
)
plan <- gwrs_sim_make_plan(
  replications = 2L,
  master_seed = master_seed,
  scenario_id = "nonconvex-benchmark"
)
plan$penalty <- config$penalties
plan$task_id <- paste0("nonconvex-benchmark-", plan$penalty)

runtime_signature <- gwrs_sim_runtime_signature("gwrs")
workflow_files <- c(
  resume_utils = file.path(
    dirname(script_path), "..", "simulations", "resume-utils.R"
  ),
  benchmark_script = script_path
)
runtime_signature <- rbind(
  runtime_signature,
  data.frame(
    field = paste0(names(workflow_files), "_sha256"),
    value = gwrs_sim_sha256(unname(workflow_files)),
    stringsAsFactors = FALSE
  )
)

make_benchmark_data <- function(config) {
  set.seed(config$data_seed)
  coords <- cbind(east = runif(config$n), north = runif(config$n))
  x <- matrix(rnorm(config$n * config$p), nrow = config$n)
  colnames(x) <- paste0("x", seq_len(config$p))
  active <- seq_len(min(5L, config$p))
  beta <- numeric(config$p)
  beta[active] <- c(1.4, -1.1, 0.8, -0.55, 0.35)[active]
  spatial_signal <-
    0.5 * sin(2 * pi * coords[, 1L]) +
    0.35 * cos(2 * pi * coords[, 2L])
  y <- 0.7 + as.double(x %*% beta) + spatial_signal +
    rnorm(config$n, sd = 0.45)
  list(coords = coords, x = x, y = y)
}

task_function <- function(task, config, heartbeat) {
  penalty <- task$penalty[[1L]]
  gamma <- if (identical(penalty, "scad")) {
    config$scad_gamma
  } else {
    config$mcp_gamma
  }
  heartbeat("simulate-data")
  data <- make_benchmark_data(config)

  heartbeat("build-neighbours")
  neighbor_time <- system.time({
    neighbors <- gwr_neighbors(data$coords, k = config$k)
  })[["elapsed"]]

  base_control <- function(n_threads) {
    gwrs_control(
      tolerance = config$tolerance,
      max_iterations = config$max_iterations,
      n_threads = n_threads,
      grain_size = config$grain_size,
      keep_data = FALSE,
      diagnostics = "none"
    )
  }
  cases <- data.frame(
    case = c(
      "serial_unscreened", "serial_screened",
      "parallel_unscreened", "parallel_screened"
    ),
    n_threads = c(1L, 1L, -1L, -1L),
    screening = c(FALSE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  paths <- vector("list", nrow(cases))
  elapsed <- numeric(nrow(cases))

  for (index in seq_len(nrow(cases))) {
    heartbeat(paste0("fit-", cases$case[[index]]))
    timing <- system.time({
      paths[[index]] <- gwr_nonconvex_path(
        data$x,
        data$y,
        neighbors = neighbors,
        penalty = penalty,
        gamma = gamma,
        n_lambda = config$n_lambda,
        lambda_min_ratio = config$lambda_min_ratio,
        screening = cases$screening[[index]],
        keep_coefficients = FALSE,
        return_diagnostics = FALSE,
        control = base_control(cases$n_threads[[index]])
      )
    })
    elapsed[[index]] <- unname(timing[["elapsed"]])
  }

  heartbeat("compare-results")
  baseline <- paths[[1L]]$fitted.values
  prediction_differences <- vapply(
    paths,
    function(path) max(abs(path$fitted.values - baseline)),
    numeric(1L)
  )
  convergence <- vapply(
    paths,
    function(path) min(path$path_summary$convergence_rate),
    numeric(1L)
  )
  stationarity <- vapply(
    paths,
    function(path) max(path$path_summary$max_stationarity_violation),
    numeric(1L)
  )
  coordinate_gap <- vapply(
    paths,
    function(path) max(path$path_summary$max_coordinate_gap),
    numeric(1L)
  )
  names(elapsed) <- cases$case

  summary <- data.frame(
    penalty = penalty,
    gamma = gamma,
    n = config$n,
    p = config$p,
    k = config$k,
    n_lambda = config$n_lambda,
    available_threads = RcppParallel::defaultNumThreads(),
    neighbor_elapsed_seconds = unname(neighbor_time),
    serial_unscreened_seconds = elapsed[["serial_unscreened"]],
    serial_screened_seconds = elapsed[["serial_screened"]],
    parallel_unscreened_seconds = elapsed[["parallel_unscreened"]],
    parallel_screened_seconds = elapsed[["parallel_screened"]],
    serial_screening_speedup =
      elapsed[["serial_unscreened"]] / elapsed[["serial_screened"]],
    screened_parallel_speedup =
      elapsed[["serial_screened"]] / elapsed[["parallel_screened"]],
    maximum_prediction_difference = max(prediction_differences),
    minimum_convergence_rate = min(convergence),
    maximum_stationarity_violation = max(stationarity),
    maximum_coordinate_gap = max(coordinate_gap),
    lambda_max = paths[[1L]]$lambda[[1L]],
    lambda_min = paths[[1L]]$lambda[[config$n_lambda]],
    retained_paths_mib = sum(vapply(
      paths, function(path) as.numeric(object.size(path)) / 1024^2,
      numeric(1L)
    )),
    dense_weight_matrix_gib_avoided = 8 * config$n^2 / 1024^3,
    stringsAsFactors = FALSE
  )
  list(summary = summary)
}

validate_result <- function(result, task) {
  summary <- result$summary
  is.list(result) && is.data.frame(summary) && nrow(summary) == 1L &&
    identical(summary$penalty, task$penalty) &&
    isTRUE(summary$minimum_convergence_rate == 1) &&
    is.finite(summary$maximum_prediction_difference) &&
    summary$maximum_prediction_difference <= 1e-6 &&
    is.finite(summary$maximum_stationarity_violation) &&
    is.finite(summary$maximum_coordinate_gap) &&
    all(is.finite(unlist(summary[c(
      "serial_unscreened_seconds", "serial_screened_seconds",
      "parallel_unscreened_seconds", "parallel_screened_seconds"
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
  confirm_large = confirm_large,
  limits = list(n = 20000L, p = 100L, k = 200L, tasks = 500L)
)
print(status$progress, row.names = FALSE)
cat("Output:", status$output_directory, "\n")
