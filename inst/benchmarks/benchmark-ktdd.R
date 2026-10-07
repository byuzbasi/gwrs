# Reproducible, resumable GWR-KTDD benchmark; run after installing gwrs.
#
# Required:
#   GWRS_KTDD_OUTPUT=/new/versioned/output/directory
#
# Production runs additionally require:
#   GWRS_KTDD_MODE=production
#   GWRS_KTDD_CONFIRM_PRODUCTION=YES

library(gwrs)

benchmark_schema <- "gwrs-ktdd-benchmark-v0.3.2"

read_integer <- function(name, default, minimum = 1L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(raw))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(sprintf("`%s` must be an integer >= %d.", name, minimum),
         call. = FALSE)
  }
  value
}

read_double <- function(name, default, allow_infinite = FALSE) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.double(raw))
  valid <- length(value) == 1L && !is.na(value) && value > 0 &&
    (is.finite(value) || (allow_infinite && is.infinite(value)))
  if (!valid) {
    stop(sprintf("`%s` must be a positive number.", name), call. = FALSE)
  }
  value
}

read_csv_values <- function(name, default) {
  raw <- Sys.getenv(name, unset = default)
  values <- trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])
  values[nzchar(values)]
}

read_threads <- function() {
  raw <- read_csv_values("GWRS_KTDD_THREADS", "1,-1")
  value <- suppressWarnings(as.integer(raw))
  if (anyNA(value) || any(value == 0L) || any(value < -1L)) {
    stop("`GWRS_KTDD_THREADS` must contain -1 or positive integers.",
         call. = FALSE)
  }
  unique(value)
}

atomic_save_rds <- function(object, path) {
  if (file.exists(path)) {
    stop(sprintf("Refusing to overwrite `%s`.", path), call. = FALSE)
  }
  temporary <- tempfile(pattern = ".partial-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3, compress = FALSE)
  if (!file.rename(temporary, path)) {
    stop(sprintf("Could not atomically create `%s`.", path), call. = FALSE)
  }
  invisible(path)
}

atomic_write_lines <- function(lines, path) {
  if (file.exists(path)) {
    stop(sprintf("Refusing to overwrite `%s`.", path), call. = FALSE)
  }
  temporary <- tempfile(pattern = ".partial-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) {
    stop(sprintf("Could not atomically create `%s`.", path), call. = FALSE)
  }
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  if (file.exists(path)) {
    stop(sprintf("Refusing to overwrite `%s`.", path), call. = FALSE)
  }
  temporary <- tempfile(pattern = ".partial-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    stop(sprintf("Could not atomically create `%s`.", path), call. = FALSE)
  }
  invisible(path)
}

sha256sum_portable <- function(files) {
  if ("sha256sum" %in% getNamespaceExports("tools")) {
    return(unname(as.character(tools::sha256sum(files))))
  }
  executable <- Sys.which("shasum")
  if (nzchar(executable)) {
    return(vapply(files, function(path) {
      output <- system2(executable, c("-a", "256", shQuote(path)),
                        stdout = TRUE, stderr = TRUE)
      if (!length(output)) stop("`shasum` returned no output.", call. = FALSE)
      strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
    }, character(1)))
  }
  executable <- Sys.which("sha256sum")
  if (nzchar(executable)) {
    return(vapply(files, function(path) {
      output <- system2(executable, shQuote(path), stdout = TRUE, stderr = TRUE)
      if (!length(output)) {
        stop("`sha256sum` returned no output.", call. = FALSE)
      }
      strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
    }, character(1)))
  }
  stop("No SHA-256 implementation is available.", call. = FALSE)
}

mode <- match.arg(
  Sys.getenv("GWRS_KTDD_MODE", unset = "smoke"),
  c("smoke", "production")
)
if (identical(mode, "production") &&
    !identical(Sys.getenv("GWRS_KTDD_CONFIRM_PRODUCTION"), "YES")) {
  stop(
    paste0(
      "Production mode requires `GWRS_KTDD_CONFIRM_PRODUCTION=YES`. ",
      "This guard prevents an accidental long run."
    ),
    call. = FALSE
  )
}

output_directory <- Sys.getenv("GWRS_KTDD_OUTPUT")
if (!nzchar(output_directory)) {
  stop("Set `GWRS_KTDD_OUTPUT` to a new, versioned output directory.",
       call. = FALSE)
}
output_directory <- normalizePath(output_directory, mustWork = FALSE)

defaults <- if (identical(mode, "production")) {
  list(n = 100000L, p = 30L, k = 100L)
} else {
  list(n = 300L, p = 6L, k = 40L)
}
n <- read_integer("GWRS_KTDD_N", defaults$n)
p <- read_integer("GWRS_KTDD_P", defaults$p)
k <- min(read_integer("GWRS_KTDD_K", defaults$k), n)
seed <- read_integer("GWRS_KTDD_SEED", 20260825L, minimum = 0L)
solvers <- read_csv_values("GWRS_KTDD_SOLVERS", "auto")
if (!length(solvers) || any(!solvers %in% c(
  "auto", "factorized", "matrix_free"
))) {
  stop(
    "`GWRS_KTDD_SOLVERS` must contain auto, factorized, or matrix_free.",
    call. = FALSE
  )
}
solvers <- unique(solvers)
threads <- read_threads()
factor_cache_limit_mib <- read_double(
  "GWRS_KTDD_FACTOR_CACHE_LIMIT_MIB", 512, allow_infinite = TRUE
)
outer_tolerance <- read_double("GWRS_KTDD_OUTER_TOLERANCE", 1e-6)
outer_max_iterations <- read_integer(
  "GWRS_KTDD_OUTER_MAX_ITERATIONS", 500L
)
coefficient_tolerance <- read_double(
  "GWRS_KTDD_COEFFICIENT_TOLERANCE", 1e-8
)
coefficient_max_iterations <- read_integer(
  "GWRS_KTDD_COEFFICIENT_MAX_ITERATIONS", 2000L
)
process_tolerance <- read_double("GWRS_KTDD_PROCESS_TOLERANCE", 1e-8)
process_max_iterations <- read_integer(
  "GWRS_KTDD_PROCESS_MAX_ITERATIONS", 2000L
)

configuration <- list(
  schema = benchmark_schema,
  mode = mode,
  package = "gwrs",
  package_version = as.character(utils::packageVersion("gwrs")),
  seed = seed,
  n = n,
  p = p,
  k = k,
  solvers = solvers,
  threads = threads,
  factor_cache_limit_mib = factor_cache_limit_mib,
  outer_tolerance = outer_tolerance,
  outer_max_iterations = outer_max_iterations,
  coefficient_tolerance = coefficient_tolerance,
  coefficient_max_iterations = coefficient_max_iterations,
  process_tolerance = process_tolerance,
  process_max_iterations = process_max_iterations
)

configuration_file <- file.path(output_directory, "run-config.rds")
if (dir.exists(output_directory)) {
  if (!file.exists(configuration_file)) {
    stop(
      paste0(
        "The output directory exists but has no run-config.rds; refusing ",
        "to write into a possibly unrelated directory."
      ),
      call. = FALSE
    )
  }
  recorded_configuration <- readRDS(configuration_file)
  if (!identical(recorded_configuration, configuration)) {
    stop(
      "The requested configuration differs from the recorded run; use a new output directory.",
      call. = FALSE
    )
  }
  message("Resuming the matching benchmark in: ", output_directory)
} else {
  if (!dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create the output directory.", call. = FALSE)
  }
  atomic_save_rds(configuration, configuration_file)
  flattened_configuration <- data.frame(
    field = names(configuration),
    value = vapply(configuration, function(value) {
      paste(value, collapse = ",")
    }, character(1)),
    row.names = NULL
  )
  atomic_write_csv(
    flattened_configuration,
    file.path(output_directory, "run-config.csv")
  )
  atomic_write_lines(
    capture.output(utils::sessionInfo()),
    file.path(output_directory, "session-info.txt")
  )
}

arma_uword_bits <- gwrs:::cpp_arma_uword_bits()
if (n >= 65536L && arma_uword_bits < 64L) {
  stop(
    paste0(
      "This benchmark requires a gwrs build with 64-bit Armadillo sparse ",
      "indices. Reinstall the current source package before resuming."
    ),
    call. = FALSE
  )
}
dll_directory <- system.file("libs", package = "gwrs")
dll_candidates <- list.files(
  dll_directory,
  pattern = paste0("[.]", sub("^[.]", "", .Platform$dynlib.ext), "$"),
  full.names = TRUE,
  recursive = TRUE
)
if (length(dll_candidates) != 1L) {
  stop("Could not identify the installed gwrs shared library.", call. = FALSE)
}
runtime_signature <- data.frame(
  field = c(
    "package_version", "arma_uword_bits", "shared_library_sha256",
    "R_version", "platform", "RcppArmadillo_version",
    "RcppParallel_version", "available_threads"
  ),
  value = c(
    as.character(utils::packageVersion("gwrs")),
    as.character(arma_uword_bits),
    sha256sum_portable(dll_candidates),
    R.version.string,
    R.version$platform,
    as.character(utils::packageVersion("RcppArmadillo")),
    as.character(utils::packageVersion("RcppParallel")),
    as.character(RcppParallel::defaultNumThreads())
  ),
  stringsAsFactors = FALSE
)
runtime_signature_file <- file.path(
  output_directory, "runtime-signature.csv"
)
if (file.exists(runtime_signature_file)) {
  recorded_runtime_signature <- utils::read.csv(
    runtime_signature_file, stringsAsFactors = FALSE
  )
  if (!identical(recorded_runtime_signature, runtime_signature)) {
    stop(
      "The installed numerical runtime differs from the recorded run; use a new output directory.",
      call. = FALSE
    )
  }
} else {
  atomic_write_csv(runtime_signature, runtime_signature_file)
}

inputs_file <- file.path(output_directory, "inputs.rds")
if (file.exists(inputs_file)) {
  inputs <- readRDS(inputs_file)
  message("Loaded input checkpoint.")
} else {
  set.seed(seed)
  coords <- cbind(east = runif(n), north = runif(n))
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  process_truth <-
    0.7 * sin(2 * pi * coords[, 1L]) * cos(pi * coords[, 2L]) +
    0.25 * cos(2 * pi * coords[, 2L])
  process_truth <- process_truth - mean(process_truth)
  coefficient <- numeric(p)
  coefficient[seq_len(min(3L, p))] <- c(1.4, -0.8, 0.45)[
    seq_len(min(3L, p))
  ]
  y <- 0.7 + as.double(x %*% coefficient) + process_truth +
    rnorm(n, sd = 0.2)

  neighbor_time <- system.time({
    neighbors <- gwr_neighbors(coords, k = k)
  })
  laplacian_time <- system.time({
    laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
  })
  source <- as.double(laplacian %*% process_truth)
  inputs <- list(
    schema = benchmark_schema,
    configuration = configuration,
    coords = coords,
    x = x,
    y = y,
    process_truth = process_truth,
    neighbors = neighbors,
    laplacian = laplacian,
    source = source,
    neighbor_time = neighbor_time,
    laplacian_time = laplacian_time
  )
  atomic_save_rds(inputs, inputs_file)
  message("Created input checkpoint.")
}

if (!identical(inputs$configuration, configuration)) {
  stop("Input checkpoint configuration mismatch.", call. = FALSE)
}
if (!identical(dim(inputs$x), c(n, p)) || length(inputs$y) != n ||
    !identical(dim(inputs$neighbors$index), c(k, n)) ||
    !identical(dim(inputs$laplacian), c(n, n)) ||
    length(inputs$source) != n) {
  stop("Input checkpoint dimensions are incomplete or inconsistent.",
       call. = FALSE)
}

case_grid <- expand.grid(
  solver = solvers,
  n_threads = threads,
  stringsAsFactors = FALSE
)
case_grid$case_id <- sprintf(
  "%s-threads-%s",
  case_grid$solver,
  ifelse(case_grid$n_threads == -1L, "auto", case_grid$n_threads)
)

run_case <- function(solver, n_threads, case_id) {
  path <- file.path(output_directory, paste0("case-", case_id, ".rds"))
  if (file.exists(path)) {
    result <- readRDS(path)
    if (!identical(result$configuration, configuration) ||
        !identical(result$case_id, case_id)) {
      stop(sprintf("Checkpoint mismatch for `%s`.", case_id), call. = FALSE)
    }
    message("Loaded completed case: ", case_id)
    return(result)
  }

  message("Running case: ", case_id)
  elapsed <- system.time({
    fit <- gwr_ktdd_fit(
      inputs$x, inputs$y,
      neighbors = inputs$neighbors,
      kernel = "gaussian",
      diffusion_operator = inputs$laplacian,
      source = inputs$source,
      coefficient_laplacian = inputs$laplacian,
      lambda_diffusion = 1,
      gamma = 0.03,
      solver = solver,
      factor_cache_limit_mib = factor_cache_limit_mib,
      coefficient_tolerance = coefficient_tolerance,
      coefficient_max_iterations = coefficient_max_iterations,
      linear_tolerance = process_tolerance,
      linear_max_iterations = process_max_iterations,
      control = gwrs_control(
        tolerance = outer_tolerance,
        max_iterations = outer_max_iterations,
        n_threads = n_threads,
        keep_data = FALSE,
        diagnostics = "none"
      )
    )
  })

  result <- list(
    schema = benchmark_schema,
    configuration = configuration,
    case_id = case_id,
    requested_solver = solver,
    resolved_solver = fit$solver,
    arma_uword_bits = arma_uword_bits,
    n_threads = n_threads,
    elapsed = elapsed,
    coefficients = fit$coefficients,
    process = fit$process,
    objective = fit$objective,
    loss_components = fit$loss_components,
    diagnostics = fit$diagnostics,
    workspace = fit$workspace,
    retained_fit_mib = as.numeric(object.size(fit)) / 1024^2
  )
  atomic_save_rds(result, path)
  result
}

results <- lapply(seq_len(nrow(case_grid)), function(index) {
  run_case(
    solver = case_grid$solver[[index]],
    n_threads = case_grid$n_threads[[index]],
    case_id = case_grid$case_id[[index]]
  )
})

summary_table <- do.call(rbind, lapply(results, function(result) {
  data.frame(
    case_id = result$case_id,
    requested_solver = result$requested_solver,
    resolved_solver = result$resolved_solver,
    arma_uword_bits = result$arma_uword_bits,
    n_threads = result$n_threads,
    elapsed_seconds = unname(result$elapsed[["elapsed"]]),
    outer_converged = result$diagnostics$converged,
    outer_iterations = result$diagnostics$iterations,
    coefficient_pcg_converged =
      result$diagnostics$coefficient_linear_converged,
    coefficient_pcg_iterations =
      result$diagnostics$coefficient_linear_iterations,
    process_pcg_converged = result$diagnostics$linear_converged,
    process_pcg_iterations = result$diagnostics$linear_iterations,
    objective = result$objective,
    objective_finite = is.finite(result$objective),
    coefficient_rows = nrow(result$coefficients),
    coefficient_columns = ncol(result$coefficients),
    process_length = length(result$process),
    factor_cache_mib = result$workspace$factor_cache_mib,
    approximate_working_mib = result$workspace$approximate_working_mib,
    retained_fit_mib = result$retained_fit_mib,
    stringsAsFactors = FALSE
  )
}))

summary_file <- file.path(output_directory, "case-summary.csv")
if (!file.exists(summary_file)) atomic_write_csv(summary_table, summary_file)

agreement_rows <- list()
for (solver in unique(case_grid$solver)) {
  selected <- results[vapply(results, function(result) {
    identical(result$requested_solver, solver)
  }, logical(1))]
  if (length(selected) > 1L) {
    reference <- selected[[1L]]
    for (candidate in selected[-1L]) {
      agreement_rows[[length(agreement_rows) + 1L]] <- data.frame(
        requested_solver = solver,
        reference_case = reference$case_id,
        candidate_case = candidate$case_id,
        coefficient_max_abs = max(abs(
          reference$coefficients - candidate$coefficients
        )),
        process_max_abs = max(abs(reference$process - candidate$process)),
        objective_abs = abs(reference$objective - candidate$objective),
        stringsAsFactors = FALSE
      )
    }
  }
}
if (length(agreement_rows)) {
  agreement_file <- file.path(output_directory, "numerical-agreement.csv")
  if (!file.exists(agreement_file)) {
    atomic_write_csv(do.call(rbind, agreement_rows), agreement_file)
  }
}

complete <- with(
  summary_table,
  outer_converged & coefficient_pcg_converged & process_pcg_converged &
    objective_finite & coefficient_rows == n &
    coefficient_columns == p + 1L & process_length == n
)
if (!all(complete)) {
  stop(
    paste0(
      "At least one case failed a convergence requirement. Check ",
      "case-summary.csv; the manifest was not finalized."
    ),
    call. = FALSE
  )
}

manifest_file <- file.path(output_directory, "manifest-sha256.csv")
if (!file.exists(manifest_file)) {
  files <- list.files(output_directory, full.names = TRUE, recursive = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[basename(files) != basename(manifest_file)]
  manifest <- data.frame(
    file = basename(files),
    bytes = unname(file.info(files)$size),
    sha256 = sha256sum_portable(files),
    stringsAsFactors = FALSE
  )
  manifest <- manifest[order(manifest$file), , drop = FALSE]
  atomic_write_csv(manifest, manifest_file)
}

print(summary_table, row.names = FALSE)
message("Benchmark complete. SHA-256 manifest: ", manifest_file)
