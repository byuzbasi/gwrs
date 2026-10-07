# Resumable local simulation utilities for gwrs research workflows.
# These functions orchestrate files and tasks; numerical fitting remains in
# the package's compiled cores.

gwrs_sim_schema <- "gwrs-resumable-simulation-v1"

gwrs_sim_task_rng_kind <- function() {
  c(
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
}

gwrs_sim_utc <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

gwrs_sim_sha256 <- function(files) {
  if (!length(files)) return(character())
  if ("sha256sum" %in% getNamespaceExports("tools")) {
    return(unname(as.character(tools::sha256sum(files))))
  }
  executable <- Sys.which("shasum")
  arguments <- "-a"
  if (!nzchar(executable)) {
    executable <- Sys.which("sha256sum")
    arguments <- character()
  }
  if (!nzchar(executable)) {
    stop("No SHA-256 implementation is available.", call. = FALSE)
  }
  vapply(files, function(path) {
    command_arguments <- c(
      arguments, if (length(arguments)) "256", shQuote(path)
    )
    output <- system2(
      executable, command_arguments, stdout = TRUE, stderr = TRUE
    )
    if (!length(output)) {
      stop("The SHA-256 command returned no output.", call. = FALSE)
    }
    strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  }, character(1L), USE.NAMES = FALSE)
}

gwrs_sim_atomic_rds <- function(object, path, replace = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !replace) {
    stop("Refusing to overwrite `", path, "`.", call. = FALSE)
  }
  temporary <- tempfile(".partial-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3, compress = FALSE)
  if (!file.rename(temporary, path)) {
    if (!replace || !file.exists(path)) {
      stop("Could not atomically create `", path, "`.", call. = FALSE)
    }
    unlink(path)
    if (!file.rename(temporary, path)) {
      stop("Could not replace derived state `", path, "`.", call. = FALSE)
    }
  }
  invisible(path)
}

gwrs_sim_atomic_csv <- function(object, path, replace = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !replace) {
    stop("Refusing to overwrite `", path, "`.", call. = FALSE)
  }
  temporary <- tempfile(".partial-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    if (!replace || !file.exists(path)) {
      stop("Could not atomically create `", path, "`.", call. = FALSE)
    }
    unlink(path)
    if (!file.rename(temporary, path)) {
      stop("Could not replace derived state `", path, "`.", call. = FALSE)
    }
  }
  invisible(path)
}

gwrs_sim_atomic_text <- function(text, path, replace = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !replace) {
    stop("Refusing to overwrite `", path, "`.", call. = FALSE)
  }
  temporary <- tempfile(".partial-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(enc2utf8(as.character(text)), temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) {
    if (!replace || !file.exists(path)) {
      stop("Could not atomically create `", path, "`.", call. = FALSE)
    }
    unlink(path)
    if (!file.rename(temporary, path)) {
      stop("Could not replace derived state `", path, "`.", call. = FALSE)
    }
  }
  invisible(path)
}

gwrs_sim_flatten_config <- function(config) {
  data.frame(
    field = names(config),
    value = vapply(config, function(value) {
      paste(value, collapse = ",")
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}

gwrs_sim_make_plan <- function(replications,
                               master_seed,
                               scenario_id = "base") {
  replications <- as.integer(replications)
  master_seed <- as.double(master_seed)
  if (length(replications) != 1L || is.na(replications) ||
      replications < 1L) {
    stop("`replications` must be a positive integer.", call. = FALSE)
  }
  if (length(master_seed) != 1L || !is.finite(master_seed) ||
      master_seed < 0 || master_seed > 2147483646) {
    stop("`master_seed` must be between 0 and 2147483646.", call. = FALSE)
  }
  if (!is.character(scenario_id) || length(scenario_id) != 1L ||
      !grepl("^[A-Za-z0-9_-]+$", scenario_id)) {
    stop("`scenario_id` must contain only letters, digits, _ or -.",
         call. = FALSE)
  }
  replication <- seq_len(replications)
  seed <- ((master_seed + replication * 104729) %% 2147483646) + 1
  data.frame(
    task_id = sprintf("%s-rep-%05d", scenario_id, replication),
    scenario_id = scenario_id,
    replication = replication,
    seed = as.integer(seed),
    stringsAsFactors = FALSE
  )
}

gwrs_sim_validate_plan <- function(plan) {
  required <- c("task_id", "seed")
  if (!is.data.frame(plan) || !nrow(plan) ||
      !all(required %in% names(plan))) {
    stop("The task plan must contain `task_id` and `seed` columns.",
         call. = FALSE)
  }
  if (anyNA(plan$task_id) || any(!grepl(
    "^[A-Za-z0-9_-]+$", plan$task_id
  )) || anyDuplicated(plan$task_id)) {
    stop("Task identifiers must be unique, path-safe strings.",
         call. = FALSE)
  }
  seeds <- suppressWarnings(as.integer(plan$seed))
  if (anyNA(seeds) || any(seeds < 1L) || anyDuplicated(seeds)) {
    stop("Task seeds must be unique positive integers.", call. = FALSE)
  }
  plan$task_id <- as.character(plan$task_id)
  plan$seed <- seeds
  rownames(plan) <- NULL
  plan
}

gwrs_sim_preflight <- function(config,
                               plan,
                               confirm_large = FALSE,
                               limits = list(
                                 n = 20000L,
                                 p = 100L,
                                 k = 200L,
                                 tasks = 500L
                               )) {
  values <- c(
    n = as.double(if (is.null(config$n)) NA_real_ else config$n),
    p = as.double(if (is.null(config$p)) NA_real_ else config$p),
    k = as.double(if (is.null(config$k)) NA_real_ else config$k),
    tasks = nrow(plan)
  )
  limits_vector <- as.double(unlist(limits[names(values)]))
  exceeded <- is.finite(values) & values > limits_vector
  if (any(exceeded) && !isTRUE(confirm_large)) {
    detail <- paste(
      paste0(names(values)[exceeded], "=", values[exceeded],
             " (limit ", limits_vector[exceeded], ")"),
      collapse = ", "
    )
    stop(
      "Large-run guard triggered: ", detail,
      ". Reduce the design or set explicit large-run confirmation.",
      call. = FALSE
    )
  }
  n <- values[["n"]]
  p <- values[["p"]]
  k <- values[["k"]]
  data.frame(
    field = c(
      "n", "p", "k", "tasks", "estimated_design_mib",
      "estimated_neighbor_mib", "large_run_confirmed"
    ),
    value = c(
      n, p, k, nrow(plan),
      if (is.finite(n) && is.finite(p)) 8 * n * p / 1024^2 else NA,
      if (is.finite(n) && is.finite(k)) 16 * n * k / 1024^2 else NA,
      isTRUE(confirm_large)
    ),
    stringsAsFactors = FALSE
  )
}

gwrs_sim_runtime_signature <- function(package = "gwrs") {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package `", package, "` is not installed.", call. = FALSE)
  }
  dll_directory <- system.file("libs", package = package)
  dll_pattern <- paste0(
    "[.]", sub("^[.]", "", .Platform$dynlib.ext), "$"
  )
  dll <- list.files(
    dll_directory, pattern = dll_pattern, full.names = TRUE, recursive = TRUE
  )
  dll_hash <- if (length(dll) == 1L) gwrs_sim_sha256(dll) else NA_character_
  package_directory <- system.file(package = package)
  runtime_files <- c(
    file.path(package_directory, c("DESCRIPTION", "NAMESPACE")),
    list.files(
      file.path(package_directory, "R"), full.names = TRUE,
      recursive = TRUE, include.dirs = FALSE
    ),
    list.files(
      file.path(package_directory, "simulations"), full.names = TRUE,
      recursive = TRUE, include.dirs = FALSE
    )
  )
  runtime_files <- sort(unique(runtime_files[file.exists(runtime_files)]))
  if (!length(runtime_files)) {
    stop("The installed package runtime files are unavailable.", call. = FALSE)
  }
  runtime_relative <- substring(
    runtime_files, nchar(package_directory) + 2L
  )
  runtime_payload <- paste(
    runtime_relative,
    as.double(unname(file.info(runtime_files)$size)),
    unname(gwrs_sim_sha256(runtime_files)),
    sep = "\t"
  )
  runtime_payload_path <- tempfile("gwrs-installed-runtime-")
  on.exit(unlink(runtime_payload_path), add = TRUE)
  writeLines(runtime_payload, runtime_payload_path, useBytes = TRUE)
  installed_runtime_hash <- gwrs_sim_sha256(runtime_payload_path)
  offline_root <- Sys.getenv("GWRS_OFFLINE_REPO", unset = "")
  offline_manifest <- if (nzchar(offline_root)) {
    file.path(offline_root, "offline-manifest-sha256.csv")
  } else {
    ""
  }
  offline_manifest_hash <- if (nzchar(offline_manifest) &&
                               file.exists(offline_manifest)) {
    gwrs_sim_sha256(offline_manifest)
  } else {
    ""
  }
  index_bits <- tryCatch(
    get("cpp_arma_uword_bits", envir = asNamespace(package))(),
    error = function(error) NA_integer_
  )
  dependency_packages <- c(
    "FNN", "Matrix", "Rcpp", "RcppArmadillo", "RcppParallel"
  )
  dependency_versions <- vapply(dependency_packages, function(dependency) {
    if (requireNamespace(dependency, quietly = TRUE)) {
      as.character(utils::packageVersion(dependency))
    } else {
      NA_character_
    }
  }, character(1L))
  workflow_raw <- Sys.getenv("GWRS_WORKFLOW_FILES", unset = "")
  workflow_files <- if (nzchar(workflow_raw)) {
    unique(strsplit(workflow_raw, ":", fixed = TRUE)[[1L]])
  } else {
    character()
  }
  if (length(workflow_files) && any(!file.exists(workflow_files))) {
    stop("A file named in GWRS_WORKFLOW_FILES does not exist.",
         call. = FALSE)
  }
  workflow_fields <- if (length(workflow_files)) {
    paste0(
      "workflow_file_",
      gsub("[^A-Za-z0-9]+", "_", basename(workflow_files)),
      "_sha256"
    )
  } else {
    character()
  }
  fields <- c(
    "package", "package_version", "shared_library_sha256",
    "installed_package_runtime_sha256", "offline_manifest_sha256",
    "armadillo_index_bits", "R_version", "platform", "RNG_kind",
    "library_paths", "RCPP_PARALLEL_NUM_THREADS", "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
    "workflow_version", "workflow_stage",
    paste0("dependency_", dependency_packages, "_version"),
    workflow_fields
  )
  values <- c(
    package, as.character(utils::packageVersion(package)), dll_hash,
    installed_runtime_hash, offline_manifest_hash,
    index_bits, R.version.string, R.version$platform,
    paste(gwrs_sim_task_rng_kind(), collapse = "/"),
    paste(.libPaths(), collapse = "|"),
    Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = ""),
    Sys.getenv("OMP_NUM_THREADS", unset = ""),
    Sys.getenv("OPENBLAS_NUM_THREADS", unset = ""),
    Sys.getenv("MKL_NUM_THREADS", unset = ""),
    Sys.getenv("BLIS_NUM_THREADS", unset = ""),
    Sys.getenv("GWRS_WORKFLOW_VERSION", unset = ""),
    Sys.getenv("GWRS_SLURM_STAGE", unset = ""),
    unname(dependency_versions),
    if (length(workflow_files)) gwrs_sim_sha256(workflow_files) else character()
  )
  data.frame(
    field = fields,
    value = as.character(values),
    stringsAsFactors = FALSE
  )
}

gwrs_sim_session_info <- function() {
  environment_names <- c(
    "SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_JOB_PARTITION",
    "R_LIBS", "R_LIBS_USER", "RCPP_PARALLEL_NUM_THREADS",
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "TMPDIR"
  )
  environment <- Sys.getenv(environment_names, unset = "")
  c(
    paste0("recorded_utc: ", gwrs_sim_utc()),
    paste0("working_directory: ", normalizePath(getwd())),
    "environment:",
    paste0("  ", environment_names, "=", environment),
    "session_info:",
    paste0("  ", capture.output(utils::sessionInfo()))
  )
}

gwrs_sim_scientific_signature <- function(output_directory) {
  relative <- c(
    "run-config.rds", "task-plan.csv", "runtime-signature.csv"
  )
  paths <- file.path(output_directory, relative)
  if (!all(file.exists(paths))) {
    stop("Cannot create the scientific signature before immutable inputs.",
         call. = FALSE)
  }
  data.frame(
    artifact = relative,
    bytes = unname(file.info(paths)$size),
    sha256 = gwrs_sim_sha256(paths),
    stringsAsFactors = FALSE
  )
}

gwrs_sim_verify_scientific_signature <- function(output_directory) {
  signature_path <- file.path(output_directory, "scientific-signature.csv")
  if (!file.exists(signature_path)) return(FALSE)
  signature <- tryCatch(
    utils::read.csv(signature_path, stringsAsFactors = FALSE),
    error = function(error) NULL
  )
  if (is.null(signature) || !identical(
    names(signature), c("artifact", "bytes", "sha256")
  )) return(FALSE)
  files <- file.path(output_directory, signature$artifact)
  if (!all(file.exists(files))) return(FALSE)
  identical(
    as.double(unname(file.info(files)$size)), as.double(signature$bytes)
  ) && identical(
    unname(gwrs_sim_sha256(files)), as.character(signature$sha256)
  )
}

gwrs_sim_checkpoint_path <- function(output_directory, task_id) {
  file.path(output_directory, "tasks", paste0(task_id, ".rds"))
}

gwrs_sim_error_files <- function(output_directory, task_id) {
  list.files(
    file.path(output_directory, "errors"),
    pattern = paste0("^", task_id, "-attempt-[0-9]+[.]rds$"),
    full.names = TRUE
  )
}

gwrs_sim_read_checkpoint <- function(path) {
  tryCatch(readRDS(path), error = function(error) NULL)
}

gwrs_sim_checkpoint_valid <- function(checkpoint,
                                      task,
                                      validate_result = NULL) {
  if (!is.list(checkpoint) ||
      !identical(checkpoint$schema, gwrs_sim_schema) ||
      !identical(checkpoint$task_id, as.character(task$task_id[[1L]])) ||
      !identical(checkpoint$seed, as.integer(task$seed[[1L]])) ||
      !identical(checkpoint$status, "completed") ||
      is.null(checkpoint$result)) {
    return(FALSE)
  }
  if (is.null(validate_result)) return(TRUE)
  isTRUE(tryCatch(
    validate_result(checkpoint$result, task),
    error = function(error) FALSE
  ))
}

gwrs_sim_quarantine <- function(path, output_directory, reason) {
  if (!file.exists(path)) return(invisible(NULL))
  directory <- file.path(output_directory, "quarantine")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
  destination <- file.path(
    directory,
    paste0(basename(path), ".", reason, ".", stamp, ".", Sys.getpid())
  )
  if (!file.rename(path, destination)) {
    stop("Could not quarantine invalid checkpoint `", path, "`.",
         call. = FALSE)
  }
  invisible(destination)
}

gwrs_sim_scan <- function(output_directory,
                          plan,
                          validate_result = NULL) {
  rows <- lapply(seq_len(nrow(plan)), function(index) {
    task <- plan[index, , drop = FALSE]
    path <- gwrs_sim_checkpoint_path(output_directory, task$task_id[[1L]])
    checkpoint <- if (file.exists(path)) gwrs_sim_read_checkpoint(path) else NULL
    valid <- gwrs_sim_checkpoint_valid(
      checkpoint, task, validate_result = validate_result
    )
    error_files <- gwrs_sim_error_files(
      output_directory, task$task_id[[1L]]
    )
    status <- if (valid) {
      "completed"
    } else if (file.exists(path)) {
      "invalid"
    } else if (length(error_files)) {
      "failed"
    } else {
      "pending"
    }
    data.frame(
      task_id = task$task_id[[1L]],
      seed = task$seed[[1L]],
      status = status,
      attempts = length(error_files) + as.integer(valid),
      elapsed_seconds = if (valid) checkpoint$elapsed_seconds else NA_real_,
      checkpoint = if (valid) basename(path) else "",
      bytes = if (valid) unname(file.info(path)$size) else NA_real_,
      sha256 = if (valid) gwrs_sim_sha256(path) else "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gwrs_sim_progress <- function(task_status,
                              state = "idle",
                              current_task = "",
                              phase = "") {
  if (!is.data.frame(task_status) || !nrow(task_status) ||
      !all(c("status", "elapsed_seconds") %in% names(task_status))) {
    stop("`task_status` must be a nonempty simulation status table.",
         call. = FALSE)
  }
  completed <- sum(task_status$status == "completed")
  failed <- sum(task_status$status %in% c("failed", "invalid"))
  elapsed <- sum(task_status$elapsed_seconds, na.rm = TRUE)
  mean_elapsed <- if (completed) elapsed / completed else NA_real_
  pending <- nrow(task_status) - completed
  progress <- data.frame(
    schema = gwrs_sim_schema,
    updated_utc = gwrs_sim_utc(),
    state = state,
    phase = phase,
    current_task_id = current_task,
    total_tasks = nrow(task_status),
    completed_tasks = completed,
    failed_or_invalid_tasks = failed,
    pending_tasks = pending,
    percent_complete = 100 * completed / nrow(task_status),
    elapsed_completed_seconds = elapsed,
    estimated_remaining_seconds = if (is.finite(mean_elapsed)) {
      mean_elapsed * pending
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

gwrs_sim_write_progress <- function(output_directory,
                                    task_status,
                                    state = "idle",
                                    current_task = "",
                                    phase = "") {
  progress <- gwrs_sim_progress(
    task_status, state = state, current_task = current_task, phase = phase
  )
  gwrs_sim_atomic_csv(
    progress, file.path(output_directory, "progress.csv"), replace = TRUE
  )
  invisible(progress)
}

gwrs_sim_write_cached_status <- function(output_directory,
                                         task_status,
                                         state = "idle",
                                         current_task = "",
                                         phase = "") {
  gwrs_sim_atomic_csv(
    task_status, file.path(output_directory, "task-status.csv"),
    replace = TRUE
  )
  progress <- gwrs_sim_write_progress(
    output_directory, task_status, state = state,
    current_task = current_task, phase = phase
  )
  invisible(list(progress = progress, tasks = task_status))
}

gwrs_sim_write_status <- function(output_directory,
                                  plan,
                                  validate_result = NULL,
                                  state = "idle",
                                  current_task = "",
                                  phase = "") {
  task_status <- gwrs_sim_scan(
    output_directory, plan, validate_result = validate_result
  )
  gwrs_sim_write_cached_status(
    output_directory, task_status, state = state,
    current_task = current_task, phase = phase
  )
}

gwrs_sim_make_heartbeat <- function(output_directory,
                                    task_status,
                                    current_task,
                                    interval_seconds = getOption(
                                      "gwrs.sim.heartbeat_interval_seconds",
                                      5
                                    )) {
  interval_seconds <- suppressWarnings(as.double(interval_seconds))
  if (length(interval_seconds) != 1L || !is.finite(interval_seconds) ||
      interval_seconds < 0) {
    stop("The heartbeat interval must be one nonnegative number.",
         call. = FALSE)
  }
  last_write <- -Inf
  function(value = "running") {
    phase <- as.character(value)
    if (length(phase) != 1L || is.na(phase)) phase <- "running"
    now <- unname(proc.time()[["elapsed"]])
    if (now - last_write >= interval_seconds) {
      gwrs_sim_write_progress(
        output_directory, task_status, state = "running",
        current_task = current_task, phase = phase
      )
      last_write <<- now
    }
    invisible(NULL)
  }
}

gwrs_sim_lock_owner_valid <- function(owner) {
  scalar_text <- function(value, allow_empty = FALSE) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      (allow_empty || nzchar(value))
  }
  if (!is.list(owner) ||
      !identical(owner$schema, gwrs_sim_schema) ||
      !scalar_text(owner$token) ||
      !scalar_text(owner$host) ||
      !scalar_text(owner$slurm_job_id, allow_empty = TRUE) ||
      !scalar_text(owner$started_utc)) {
    return(FALSE)
  }
  pid <- suppressWarnings(as.double(owner$pid))
  if (length(pid) != 1L || !is.finite(pid) || pid < 1 ||
      pid != floor(pid) || pid > .Machine$integer.max) {
    return(FALSE)
  }
  !nzchar(owner$slurm_job_id) || grepl("^[0-9]+$", owner$slurm_job_id)
}

gwrs_sim_acquire_lock <- function(output_directory, recover_lock = FALSE) {
  lock_path <- file.path(output_directory, "run.lock")
  if (dir.exists(lock_path)) {
    if (!isTRUE(recover_lock)) {
      stop(
        "The run directory is locked. Use status mode, wait for the active ",
        "process, or explicitly recover a stale lock.", call. = FALSE
      )
    }
    owner_path <- file.path(lock_path, "owner.rds")
    previous_owner <- if (file.exists(owner_path)) {
      gwrs_sim_read_checkpoint(owner_path)
    } else {
      NULL
    }
    if (!gwrs_sim_lock_owner_valid(previous_owner)) {
      stop(
        "The existing run lock has missing or invalid owner metadata; ",
        "refusing stale-lock recovery.", call. = FALSE
      )
    }
    previous_job <- previous_owner$slurm_job_id
    if (nzchar(previous_job)) {
      squeue <- Sys.which("squeue")
      if (!nzchar(squeue)) {
        stop(
          "Cannot verify SLURM job ", previous_job,
          " because `squeue` is unavailable; refusing stale-lock recovery.",
          call. = FALSE
        )
      }
      queue_state <- tryCatch(
        suppressWarnings(system2(
          squeue,
          c("--noheader", "--jobs", previous_job, "--format", "%T"),
          stdout = TRUE, stderr = TRUE
        )),
        error = function(error) error
      )
      if (inherits(queue_state, "error")) {
        stop(
          "The `squeue` query for SLURM job ", previous_job,
          " failed; refusing stale-lock recovery.", call. = FALSE
        )
      }
      queue_status <- attr(queue_state, "status")
      queue_succeeded <- is.null(queue_status) || (
        length(queue_status) == 1L && !is.na(queue_status) &&
          as.integer(queue_status) == 0L
      )
      queue_state <- trimws(as.character(queue_state))
      queue_state <- queue_state[nzchar(queue_state)]
      normalized_state <- tolower(sub("[.]$", "", queue_state))
      invalid_job_id <- !is.null(queue_status) &&
        length(queue_status) == 1L && !is.na(queue_status) &&
        as.integer(queue_status) == 1L && length(normalized_state) == 1L &&
        normalized_state %in% c(
          "invalid job id specified",
          "squeue: error: invalid job id specified",
          "slurm_load_jobs error: invalid job id specified"
        )
      # TRUBA debug jobs separately audit the cluster's exact diagnostic.
      if (!queue_succeeded && !invalid_job_id) {
        stop(
          "The `squeue` query for SLURM job ", previous_job,
          " failed with an unverified diagnostic; refusing stale-lock ",
          "recovery. Audit the exact TRUBA diagnostic in the debug stage.",
          call. = FALSE
        )
      }
      if (queue_succeeded && length(queue_state)) {
        stop(
          "SLURM job ", previous_job,
          " is still present in `squeue`; refusing stale-lock recovery.",
          call. = FALSE
        )
      }
    } else {
      current_host <- as.character(Sys.info()[["nodename"]])
      if (!identical(previous_owner$host, current_host)) {
        stop(
          "The existing local lock belongs to another or unknown host; ",
          "manual intervention is required.", call. = FALSE
        )
      }
      previous_pid <- as.integer(previous_owner$pid)
      process_alive <- tryCatch(
        suppressWarnings(tools::pskill(previous_pid, signal = 0L)),
        error = function(error) NA
      )
      if (length(process_alive) != 1L || is.na(process_alive)) {
        stop(
          "Could not verify whether the existing local lock owner is ",
          "inactive; refusing stale-lock recovery.", call. = FALSE
        )
      }
      if (isTRUE(process_alive)) {
        stop("The existing local run lock belongs to an active process.",
             call. = FALSE)
      }
    }
    quarantine <- file.path(output_directory, "quarantine")
    dir.create(quarantine, recursive = TRUE, showWarnings = FALSE)
    destination <- file.path(
      quarantine,
      paste0("run.lock.", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"),
             ".", Sys.getpid(), ".", basename(tempfile("archive-")))
    )
    if (!file.rename(lock_path, destination)) {
      stop("Could not archive the existing run lock.", call. = FALSE)
    }
  }
  if (!dir.create(lock_path, showWarnings = FALSE)) {
    stop("Could not acquire the run lock.", call. = FALSE)
  }
  token <- basename(tempfile("gwrs-sim-lock-"))
  owner <- list(
    schema = gwrs_sim_schema,
    token = token,
    pid = Sys.getpid(),
    host = Sys.info()[["nodename"]],
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = ""),
    slurm_array_job_id = Sys.getenv("SLURM_ARRAY_JOB_ID", unset = ""),
    slurm_array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""),
    slurm_submit_dir = Sys.getenv("SLURM_SUBMIT_DIR", unset = ""),
    started_utc = gwrs_sim_utc()
  )
  gwrs_sim_atomic_rds(owner, file.path(lock_path, "owner.rds"))
  list(path = lock_path, token = token)
}

gwrs_sim_release_lock <- function(lock) {
  if (is.null(lock) || !dir.exists(lock$path)) return(invisible(FALSE))
  owner_path <- file.path(lock$path, "owner.rds")
  owner <- if (file.exists(owner_path)) gwrs_sim_read_checkpoint(owner_path) else NULL
  if (!is.list(owner) || !identical(owner$token, lock$token)) {
    stop("Refusing to release a lock owned by another process.",
         call. = FALSE)
  }
  unlink(lock$path, recursive = TRUE)
  invisible(!dir.exists(lock$path))
}

gwrs_sim_initialize <- function(output_directory,
                                config,
                                plan,
                                runtime_signature,
                                confirm_large = FALSE,
                                limits = list(
                                  n = 20000L,
                                  p = 100L,
                                  k = 200L,
                                  tasks = 500L
                                )) {
  if (!is.character(output_directory) || length(output_directory) != 1L ||
      !nzchar(output_directory)) {
    stop("`output_directory` must be one nonempty path.", call. = FALSE)
  }
  if (!is.list(config) || !length(config) || is.null(names(config)) ||
      any(!nzchar(names(config)))) {
    stop("`config` must be a named nonempty list.", call. = FALSE)
  }
  plan <- gwrs_sim_validate_plan(plan)
  runtime_signature <- data.frame(
    field = as.character(runtime_signature$field),
    value = as.character(runtime_signature$value),
    stringsAsFactors = FALSE
  )
  preflight <- gwrs_sim_preflight(
    config, plan, confirm_large = confirm_large, limits = limits
  )
  output_directory <- normalizePath(output_directory, mustWork = FALSE)
  config_path <- file.path(output_directory, "run-config.rds")
  session_path <- file.path(output_directory, "session-info.txt")
  signature_path <- file.path(
    output_directory, "scientific-signature.csv"
  )
  if (dir.exists(output_directory)) {
    if (!file.exists(config_path)) {
      stop("Existing output directory has no run-config.rds.",
           call. = FALSE)
    }
    if (!identical(readRDS(config_path), config)) {
      stop("The requested configuration differs from the recorded run.",
           call. = FALSE)
    }
    recorded_plan <- utils::read.csv(
      file.path(output_directory, "task-plan.csv"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    for (field in intersect(names(plan), names(recorded_plan))) {
      if (is.integer(plan[[field]])) {
        recorded_plan[[field]] <- as.integer(recorded_plan[[field]])
      } else if (is.double(plan[[field]])) {
        recorded_plan[[field]] <- as.double(recorded_plan[[field]])
      } else if (is.logical(plan[[field]])) {
        recorded_plan[[field]] <- as.logical(recorded_plan[[field]])
      } else if (is.character(plan[[field]])) {
        recorded_plan[[field]] <- as.character(recorded_plan[[field]])
      }
    }
    if (!identical(recorded_plan, plan)) {
      stop("The requested task plan differs from the recorded run.",
           call. = FALSE)
    }
    recorded_runtime <- utils::read.csv(
      file.path(output_directory, "runtime-signature.csv"),
      stringsAsFactors = FALSE
    )
    if (!identical(recorded_runtime, runtime_signature)) {
      stop("The numerical runtime differs from the recorded run.",
           call. = FALSE)
    }
    if (file.exists(signature_path) &&
        !gwrs_sim_verify_scientific_signature(output_directory)) {
      stop("The recorded scientific signature is invalid.",
           call. = FALSE)
    }
    if (!file.exists(file.path(output_directory, "manifest-sha256.csv"))) {
      if (!file.exists(session_path)) {
        gwrs_sim_atomic_text(gwrs_sim_session_info(), session_path)
      }
      if (!file.exists(signature_path)) {
        gwrs_sim_atomic_csv(
          gwrs_sim_scientific_signature(output_directory), signature_path
        )
      }
    }
  } else {
    if (!dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create the output directory.", call. = FALSE)
    }
    dir.create(file.path(output_directory, "tasks"))
    dir.create(file.path(output_directory, "errors"))
    dir.create(file.path(output_directory, "quarantine"))
    gwrs_sim_atomic_rds(config, config_path)
    gwrs_sim_atomic_csv(
      gwrs_sim_flatten_config(config),
      file.path(output_directory, "run-config.csv")
    )
    gwrs_sim_atomic_csv(plan, file.path(output_directory, "task-plan.csv"))
    gwrs_sim_atomic_csv(
      runtime_signature,
      file.path(output_directory, "runtime-signature.csv")
    )
    gwrs_sim_atomic_csv(
      preflight, file.path(output_directory, "preflight.csv")
    )
    gwrs_sim_atomic_text(gwrs_sim_session_info(), session_path)
    gwrs_sim_atomic_csv(
      gwrs_sim_scientific_signature(output_directory), signature_path
    )
  }
  list(
    output_directory = output_directory,
    config = config,
    plan = plan,
    runtime_signature = runtime_signature,
    preflight = preflight
  )
}

gwrs_sim_save_error <- function(output_directory, task, error, elapsed) {
  existing <- gwrs_sim_error_files(output_directory, task$task_id[[1L]])
  attempt <- length(existing) + 1L
  path <- file.path(
    output_directory, "errors",
    sprintf("%s-attempt-%03d.rds", task$task_id[[1L]], attempt)
  )
  record <- list(
    schema = gwrs_sim_schema,
    task_id = task$task_id[[1L]],
    seed = task$seed[[1L]],
    attempt = attempt,
    failed_utc = gwrs_sim_utc(),
    elapsed_seconds = elapsed,
    message = conditionMessage(error),
    call = paste(deparse(conditionCall(error)), collapse = " ")
  )
  gwrs_sim_atomic_rds(record, path)
  invisible(path)
}

gwrs_sim_collect_summaries <- function(output_directory, plan) {
  summaries <- lapply(seq_len(nrow(plan)), function(index) {
    task <- plan[index, , drop = FALSE]
    checkpoint <- readRDS(gwrs_sim_checkpoint_path(
      output_directory, task$task_id[[1L]]
    ))
    summary <- checkpoint$result$summary
    if (!is.data.frame(summary) || nrow(summary) < 1L) {
      stop("Each completed result must contain a nonempty `summary`.",
           call. = FALSE)
    }
    summary_rows <- nrow(summary)
    cbind(
      data.frame(
        task_id = rep(task$task_id[[1L]], summary_rows),
        seed = rep(task$seed[[1L]], summary_rows),
        elapsed_seconds = rep(checkpoint$elapsed_seconds, summary_rows),
        stringsAsFactors = FALSE
      ),
      summary
    )
  })
  do.call(rbind, summaries)
}

gwrs_sim_task_manifest <- function(output_directory,
                                   plan,
                                   validate_result = NULL) {
  status <- gwrs_sim_scan(
    output_directory, plan, validate_result = validate_result
  )
  if (!all(status$status == "completed")) {
    stop("Cannot create a task manifest for an incomplete run.",
         call. = FALSE)
  }
  completed_utc <- vapply(seq_len(nrow(plan)), function(index) {
    checkpoint <- readRDS(gwrs_sim_checkpoint_path(
      output_directory, plan$task_id[[index]]
    ))
    as.character(checkpoint$completed_utc)
  }, character(1L))
  completed_slurm_job_id <- vapply(seq_len(nrow(plan)), function(index) {
    checkpoint <- readRDS(gwrs_sim_checkpoint_path(
      output_directory, plan$task_id[[index]]
    ))
    if (is.null(checkpoint$slurm_job_id)) "" else
      as.character(checkpoint$slurm_job_id)
  }, character(1L))
  completed_slurm_step_id <- vapply(seq_len(nrow(plan)), function(index) {
    checkpoint <- readRDS(gwrs_sim_checkpoint_path(
      output_directory, plan$task_id[[index]]
    ))
    if (is.null(checkpoint$slurm_step_id)) "" else
      as.character(checkpoint$slurm_step_id)
  }, character(1L))
  cbind(
    plan,
    data.frame(
      status = status$status,
      attempts = status$attempts,
      elapsed_seconds = status$elapsed_seconds,
      completed_utc = completed_utc,
      slurm_job_id = completed_slurm_job_id,
      slurm_step_id = completed_slurm_step_id,
      checkpoint = file.path("tasks", status$checkpoint),
      bytes = as.double(status$bytes),
      sha256 = status$sha256,
      stringsAsFactors = FALSE
    )
  )
}

gwrs_sim_completion_path <- function(output_directory) {
  file.path(output_directory, "COMPLETED")
}

gwrs_sim_read_completion <- function(output_directory) {
  path <- gwrs_sim_completion_path(output_directory)
  if (!file.exists(path)) return(NULL)
  record <- tryCatch(
    as.list(as.data.frame(read.dcf(path), stringsAsFactors = FALSE)[1L, ]),
    error = function(error) NULL
  )
  if (is.null(record)) return(NULL)
  lapply(record, as.character)
}

gwrs_sim_verify_completion <- function(output_directory) {
  record <- gwrs_sim_read_completion(output_directory)
  required <- c(
    "Schema", "Completed-UTC", "Total-Tasks", "Manifest-SHA256",
    "Scientific-Signature-SHA256"
  )
  if (is.null(record) || !all(required %in% names(record)) ||
      !identical(record$Schema, gwrs_sim_schema)) return(FALSE)
  manifest_path <- file.path(output_directory, "manifest-sha256.csv")
  signature_path <- file.path(output_directory, "scientific-signature.csv")
  task_manifest_path <- file.path(output_directory, "task-manifest.csv")
  if (!all(file.exists(c(
    manifest_path, signature_path, task_manifest_path
  )))) return(FALSE)
  tasks <- tryCatch(
    utils::read.csv(task_manifest_path, stringsAsFactors = FALSE),
    error = function(error) NULL
  )
  if (is.null(tasks) || !nrow(tasks) ||
      !all(tasks$status == "completed")) return(FALSE)
  identical(record$`Manifest-SHA256`, gwrs_sim_sha256(manifest_path)) &&
    identical(
      record$`Scientific-Signature-SHA256`,
      gwrs_sim_sha256(signature_path)
    ) && identical(
      suppressWarnings(as.integer(record$`Total-Tasks`)),
      as.integer(nrow(tasks))
    ) && gwrs_sim_verify_scientific_signature(output_directory) &&
    gwrs_sim_verify_manifest(output_directory)
}

gwrs_sim_finalize <- function(output_directory, plan, validate_result = NULL) {
  status <- gwrs_sim_write_status(
    output_directory, plan, validate_result = validate_result,
    state = "finalizing", phase = "validate-all-tasks"
  )
  if (!all(status$tasks$status == "completed")) {
    stop("Cannot finalize an incomplete simulation run.", call. = FALSE)
  }
  summaries <- gwrs_sim_collect_summaries(output_directory, plan)
  gwrs_sim_atomic_csv(
    summaries, file.path(output_directory, "task-summary.csv"),
    replace = TRUE
  )
  task_manifest <- gwrs_sim_task_manifest(
    output_directory, plan, validate_result = validate_result
  )
  gwrs_sim_atomic_csv(
    task_manifest, file.path(output_directory, "task-manifest.csv"),
    replace = TRUE
  )
  if (!gwrs_sim_verify_scientific_signature(output_directory)) {
    stop("Cannot finalize an invalid scientific signature.", call. = FALSE)
  }
  manifest_path <- file.path(output_directory, "manifest-sha256.csv")
  relative <- list.files(
    output_directory, recursive = TRUE, all.files = TRUE,
    no.. = TRUE, include.dirs = FALSE
  )
  relative <- relative[
    !relative %in% c(
      "manifest-sha256.csv", "COMPLETED", "progress.csv",
      "progress.json", "progress.tsv", "task-status.csv"
    ) &
      !grepl("(^|/)run[.]lock/", relative) &
      !grepl("(^|/)quarantine/", relative) &
      !grepl("(^|/)[.]partial-", relative)
  ]
  paths <- file.path(output_directory, relative)
  manifest <- data.frame(
    file = relative,
    bytes = as.double(unname(file.info(paths)$size)),
    sha256 = gwrs_sim_sha256(paths),
    stringsAsFactors = FALSE
  )
  manifest <- manifest[order(manifest$file), , drop = FALSE]
  if (file.exists(manifest_path)) {
    recorded <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    recorded$bytes <- as.double(recorded$bytes)
    if (!identical(recorded, manifest)) {
      stop("The existing final manifest does not match the run.",
           call. = FALSE)
    }
  } else {
    gwrs_sim_atomic_csv(manifest, manifest_path)
  }
  if (!gwrs_sim_verify_manifest(output_directory)) {
    stop("The final size/checksum manifest failed verification.",
         call. = FALSE)
  }
  completion_path <- gwrs_sim_completion_path(output_directory)
  if (file.exists(completion_path)) {
    if (!gwrs_sim_verify_completion(output_directory)) {
      stop("The existing completion marker is invalid.", call. = FALSE)
    }
  } else {
    completion <- c(
      paste0("Schema: ", gwrs_sim_schema),
      paste0("Completed-UTC: ", gwrs_sim_utc()),
      paste0("Total-Tasks: ", nrow(task_manifest)),
      paste0("Manifest-SHA256: ", gwrs_sim_sha256(manifest_path)),
      paste0(
        "Scientific-Signature-SHA256: ",
        gwrs_sim_sha256(file.path(
          output_directory, "scientific-signature.csv"
        ))
      )
    )
    gwrs_sim_atomic_text(completion, completion_path)
    if (!gwrs_sim_verify_completion(output_directory)) {
      stop("The completion marker failed verification.", call. = FALSE)
    }
  }
  gwrs_sim_write_status(
    output_directory, plan, validate_result = validate_result,
    state = "complete", phase = "validated-completion"
  )
  invisible(manifest)
}

gwrs_sim_verify_manifest <- function(output_directory) {
  path <- file.path(output_directory, "manifest-sha256.csv")
  if (!file.exists(path)) return(FALSE)
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!identical(names(manifest), c("file", "bytes", "sha256")) ||
      anyDuplicated(manifest$file)) return(FALSE)
  files <- file.path(output_directory, manifest$file)
  all(file.exists(files)) && identical(
    as.double(unname(file.info(files)$size)), as.double(manifest$bytes)
  ) && identical(
    unname(gwrs_sim_sha256(files)), as.character(manifest$sha256)
  )
}

gwrs_sim_evaluate_task <- function(task,
                                   config,
                                   task_function,
                                   validate_result,
                                   heartbeat) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(
      ".Random.seed", envir = .GlobalEnv, inherits = FALSE
    )) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  do.call(RNGkind, as.list(gwrs_sim_task_rng_kind()))
  set.seed(task$seed[[1L]])
  tryCatch(
    {
      result <- task_function(task, config, heartbeat)
      if (!isTRUE(validate_result(result, task))) {
        stop("Task result failed scientific validation.", call. = FALSE)
      }
      list(ok = TRUE, result = result)
    },
    error = function(error) list(ok = FALSE, error = error)
  )
}

gwrs_sim_run <- function(output_directory,
                         config,
                         plan,
                         task_function,
                         validate_result,
                         runtime_signature,
                         max_tasks = Inf,
                         retry_failed = TRUE,
                         recover_lock = FALSE,
                         confirm_large = FALSE,
                         limits = list(
                           n = 20000L,
                           p = 100L,
                           k = 200L,
                           tasks = 500L
                         )) {
  if (!is.function(task_function) || !is.function(validate_result)) {
    stop("Task and validation callbacks must be functions.", call. = FALSE)
  }
  if (length(max_tasks) != 1L || is.na(max_tasks) || max_tasks <= 0) {
    stop("`max_tasks` must be positive or Inf.", call. = FALSE)
  }
  run <- gwrs_sim_initialize(
    output_directory, config, plan, runtime_signature,
    confirm_large = confirm_large, limits = limits
  )
  manifest_path <- file.path(run$output_directory, "manifest-sha256.csv")
  if (file.exists(manifest_path)) {
    if (!gwrs_sim_verify_manifest(run$output_directory)) {
      stop("The finalized run manifest is invalid.", call. = FALSE)
    }
    completion_path <- gwrs_sim_completion_path(run$output_directory)
    if (file.exists(completion_path)) {
      if (!gwrs_sim_verify_completion(run$output_directory)) {
        stop("The completion marker is invalid.", call. = FALSE)
      }
      if (dir.exists(file.path(run$output_directory, "run.lock"))) {
        completed_lock <- gwrs_sim_acquire_lock(
          run$output_directory, recover_lock = recover_lock
        )
        gwrs_sim_release_lock(completed_lock)
      }
      recorded_status <- gwrs_sim_status(
        run$output_directory, quiet = TRUE
      )
      progress_complete <- !is.null(recorded_status$progress) &&
        nrow(recorded_status$progress) == 1L &&
        identical(recorded_status$progress$state[[1L]], "complete")
      if (!progress_complete) {
        gwrs_sim_write_status(
          run$output_directory, run$plan,
          validate_result = validate_result,
          state = "complete", phase = "validated-completion"
        )
      }
      return(invisible(gwrs_sim_status(run$output_directory, quiet = TRUE)))
    }
  }

  lock <- gwrs_sim_acquire_lock(
    run$output_directory, recover_lock = recover_lock
  )
  on.exit(gwrs_sim_release_lock(lock), add = TRUE)
  attempted <- 0L
  status_cache <- gwrs_sim_write_status(
    run$output_directory, run$plan, validate_result = validate_result,
    state = "running", phase = "resume-scan"
  )

  for (index in seq_len(nrow(run$plan))) {
    task <- run$plan[index, , drop = FALSE]
    task_index <- match(task$task_id[[1L]], status_cache$tasks$task_id)
    if (is.na(task_index)) {
      stop("A planned task is absent from the cached task status.",
           call. = FALSE)
    }
    if (identical(
      status_cache$tasks$status[[task_index]], "completed"
    )) next
    path <- gwrs_sim_checkpoint_path(
      run$output_directory, task$task_id[[1L]]
    )
    if (file.exists(path)) {
      gwrs_sim_quarantine(path, run$output_directory, "invalid")
    }
    failed_before <- length(gwrs_sim_error_files(
      run$output_directory, task$task_id[[1L]]
    )) > 0L
    if (failed_before && !isTRUE(retry_failed)) next
    if (attempted >= max_tasks) break

    heartbeat <- gwrs_sim_make_heartbeat(
      run$output_directory, status_cache$tasks,
      current_task = task$task_id[[1L]]
    )
    heartbeat("starting-task")
    started <- proc.time()[["elapsed"]]
    outcome <- gwrs_sim_evaluate_task(
      task = task,
      config = run$config,
      task_function = task_function,
      validate_result = validate_result,
      heartbeat = heartbeat
    )
    elapsed <- proc.time()[["elapsed"]] - started

    if (isTRUE(outcome$ok)) {
      checkpoint <- list(
        schema = gwrs_sim_schema,
        task_id = task$task_id[[1L]],
        seed = task$seed[[1L]],
        status = "completed",
        completed_utc = gwrs_sim_utc(),
        slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = ""),
        slurm_step_id = Sys.getenv("SLURM_STEP_ID", unset = ""),
        elapsed_seconds = unname(elapsed),
        result = outcome$result
      )
      gwrs_sim_atomic_rds(checkpoint, path)
    } else {
      gwrs_sim_save_error(
        run$output_directory, task, outcome$error, unname(elapsed)
      )
    }
    attempted <- attempted + 1L
    updated_task <- gwrs_sim_scan(
      run$output_directory, task, validate_result = validate_result
    )
    if (!identical(updated_task$task_id, task$task_id[[1L]])) {
      stop("The task-level status update returned an unexpected task.",
           call. = FALSE)
    }
    status_cache$tasks[task_index, ] <- updated_task
    status_cache <- gwrs_sim_write_cached_status(
      run$output_directory, status_cache$tasks, state = "running",
      phase = if (isTRUE(outcome$ok)) "checkpoint-saved" else "task-failed"
    )
  }

  status <- gwrs_sim_write_status(
    run$output_directory, run$plan, validate_result = validate_result,
    state = "incomplete", phase = "waiting-for-resume"
  )
  if (all(status$tasks$status == "completed")) {
    gwrs_sim_finalize(
      run$output_directory, run$plan, validate_result = validate_result
    )
  }
  gwrs_sim_release_lock(lock)
  lock <- NULL
  invisible(gwrs_sim_status(run$output_directory, quiet = TRUE))
}

gwrs_sim_status <- function(output_directory, quiet = FALSE) {
  output_directory <- normalizePath(output_directory, mustWork = TRUE)
  plan_path <- file.path(output_directory, "task-plan.csv")
  if (!file.exists(plan_path)) {
    stop("The directory has no task plan.", call. = FALSE)
  }
  plan <- utils::read.csv(
    plan_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  plan$seed <- as.integer(plan$seed)
  tasks <- gwrs_sim_scan(output_directory, plan)
  progress_path <- file.path(output_directory, "progress.csv")
  progress <- if (file.exists(progress_path)) {
    utils::read.csv(progress_path, stringsAsFactors = FALSE)
  } else {
    NULL
  }
  result <- list(
    output_directory = output_directory,
    progress = progress,
    tasks = tasks,
    locked = dir.exists(file.path(output_directory, "run.lock")),
    finalized = file.exists(file.path(
      output_directory, "manifest-sha256.csv"
    )),
    manifest_valid = gwrs_sim_verify_manifest(output_directory),
    completed = file.exists(gwrs_sim_completion_path(output_directory)),
    completion_valid = gwrs_sim_verify_completion(output_directory)
  )
  if (!quiet) {
    if (!is.null(progress)) print(progress, row.names = FALSE)
    print(tasks, row.names = FALSE)
    cat("locked:", result$locked, " finalized:", result$finalized,
        " manifest_valid:", result$manifest_valid,
        " completed:", result$completed,
        " completion_valid:", result$completion_valid, "\n")
  }
  result
}
