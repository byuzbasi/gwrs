# Reuse the installed workflow's atomic writers, SHA-256 and owned locks.
selection_io_utils <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      path <- system.file("simulations", "resume-utils.R", package = "gwrs")
      if (!nzchar(path)) stop("Installed checkpoint utilities are missing.", call. = FALSE)
      cached <<- new.env(parent = baseenv())
      sys.source(path, envir = cached)
    }
    cached
  }
})

selection_hash <- function(object) {
  path <- tempfile("gwrs-signature-")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 3, compress = FALSE)
  if (requireNamespace("digest", quietly = TRUE))
    return(digest::digest(file = path, algo = "sha256"))
  selection_io_utils()$gwrs_sim_sha256(path)
}

selection_provenance <- function(prep, settings, weights) {
  functions <- c("gwr_select", "selection_prepare", "selection_grid",
    "selection_evaluate", "selection_fit", "selection_order", "selection_result",
    "gwr_fit", "gwr_sl_fit", "gwr_local_inference",
    "gaussian_smoother_diagnostics", "compute_gwrs_weighted_moran")
  code <- lapply(functions, function(name) {
    f <- get(name, envir = asNamespace("gwrs"))
    list(formals(f), body(f))
  })
  dll <- getLoadedDLLs()[["gwrs"]][["path"]]
  package_files <- system.file(c("R/gwrs.rdb", "DESCRIPTION"), package = "gwrs")
  runtime <- list(R = R.version.string, platform = R.version$platform,
    versions = vapply(c("gwrs", "Rcpp", "RcppArmadillo", "RcppParallel", "FNN", "Matrix"),
      function(p) as.character(utils::packageVersion(p)), ""),
    code_sha256 = selection_hash(code),
    dll_sha256 = selection_io_utils()$gwrs_sim_sha256(dll),
    installed_files_sha256 = selection_io_utils()$gwrs_sim_sha256(package_files[nzchar(package_files)]),
    math_libraries = extSoftVersion()[intersect(c("BLAS", "LAPACK"), names(extSoftVersion()))])
  inputs <- list(x = prep$full$x, y = prep$full$y, coords = prep$coords,
    row_names = rownames(prep$full$model), labels = prep$labels,
    assignment = prep$assignment, parents = prep$parents)
  signature <- list(schema = "gwrs-joint-selection-v1", data_sha256 = selection_hash(inputs),
    settings_sha256 = selection_hash(settings), weights_sha256 = selection_hash(weights),
    runtime = runtime)
  list(signature = signature, sha256 = selection_hash(signature),
    sample_n = length(prep$full$y), sample_rows = rownames(prep$full$model),
    formula = paste(deparse(stats::formula(prep$full$terms)), collapse = " "),
    moran = if (is.null(weights)) list(status = "not_calculated: weights_not_supplied") else
      list(status = "fixed_user_weights", sha256 = signature$weights_sha256,
        row_standardize = settings$moran_row_standardize,
        zero_policy = settings$moran_zero_policy,
        islands = which(Matrix::rowSums(weights) == 0)),
    preprocessing = list(standardize = settings$standardize, missing = "fail",
      sample = "same full-scope observations for every candidate"))
}

selection_manifest_write <- function(io) {
  payload <- io$manifest
  payload$sha256 <- NULL
  io$manifest$sha256 <- selection_hash(payload)
  io$utils$gwrs_sim_atomic_rds(io$manifest, file.path(io$path, "manifest.rds"), replace = TRUE)
}

selection_checkpoint_read <- function(path) {
  utils <- selection_io_utils()
  manifest <- readRDS(file.path(path, "manifest.rds"))
  payload <- manifest
  payload$sha256 <- NULL
  if (!identical(manifest$schema, "gwrs-joint-selection-v1") ||
      !identical(manifest$sha256, selection_hash(payload)))
    stop("Checkpoint manifest signature is invalid.", call. = FALSE)
  records <- list()
  for (entry in manifest$files) {
    if (!is.character(entry$file) || length(entry$file) != 1L ||
        basename(entry$file) != entry$file || !grepl("^(pack|result|ranking)-[0-9]+[.](rds|csv)$", entry$file))
      stop("Invalid checkpoint file entry.", call. = FALSE)
    file <- file.path(path, entry$file)
    if (!file.exists(file) || !identical(as.double(file.info(file)$size), entry$bytes) ||
        !identical(utils$gwrs_sim_sha256(file), entry$sha256))
      stop("Checkpoint size/checksum mismatch: ", entry$file, call. = FALSE)
    if (entry$kind == "pack") {
      pack <- readRDS(file)
      if (!identical(pack$signature, manifest$provenance$sha256) ||
          !is.list(pack$records) || !length(pack$records) ||
          anyDuplicated(names(pack$records)) || any(names(pack$records) %in% names(records)))
        stop("Invalid or duplicate checkpoint records.", call. = FALSE)
      for (id in names(pack$records)) {
        row <- pack$records[[id]]$row
        if (!is.data.frame(row) || nrow(row) != 1L ||
            !identical(row$candidate_id, id) || !is.logical(row$eligible) ||
            length(row$eligible) != 1L || is.na(row$eligible) ||
            !all(selection_criteria() %in% names(row)))
          stop("Invalid candidate checkpoint schema.", call. = FALSE)
      }
      records <- c(records, pack$records)
    }
  }
  if (length(records) != manifest$durable_count)
    stop("Checkpoint candidate count is inconsistent.", call. = FALSE)
  list(manifest = manifest, records = records)
}

selection_checkpoint_open <- function(path, resume, provenance) {
  if (is.null(path)) {
    if (resume) stop("resume requires a checkpoint directory.", call. = FALSE)
    return(NULL)
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("checkpoint must be a directory path.", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Checkpoint progress requires the suggested package jsonlite.", call. = FALSE)
  if (file.exists(path) && !resume)
    stop("Checkpoint path already exists; use a new path or resume = TRUE.", call. = FALSE)
  if (resume && !file.exists(file.path(path, "manifest.rds")))
    stop("No checkpoint manifest to resume.", call. = FALSE)
  if (!resume && !dir.create(path, recursive = TRUE, showWarnings = FALSE))
    stop("Could not create checkpoint directory.", call. = FALSE)
  io <- new.env(parent = emptyenv())
  io$path <- normalizePath(path, mustWork = TRUE)
  io$utils <- selection_io_utils()
  io$lock <- io$utils$gwrs_sim_acquire_lock(io$path)
  success <- FALSE
  on.exit(if (!success) selection_checkpoint_unlock(io), add = TRUE)
  if (resume) {
    saved <- selection_checkpoint_read(io$path)
    io$manifest <- saved$manifest
    io$records <- saved$records
    if (!identical(io$manifest$provenance$sha256, provenance$sha256))
      stop("Checkpoint data, scientific settings or software signature differs.", call. = FALSE)
  } else {
    io$records <- list()
    io$manifest <- list(schema = "gwrs-joint-selection-v1",
      run_id = basename(io$path), created_utc = io$utils$gwrs_sim_utc(),
      provenance = provenance, durable_count = 0L, files = list(), search_status = NULL)
    selection_manifest_write(io)
  }
  success <- TRUE
  io
}

selection_checkpoint_unlock <- function(io) io$utils$gwrs_sim_release_lock(io$lock)

selection_checkpoint_next <- function(io, prefix, extension) {
  # Uncommitted orphan packs are retained, ignored and never overwritten.
  files <- list.files(io$path, pattern = paste0("^", prefix, "-[0-9]+[.]", extension, "$"))
  ids <- as.integer(sub(paste0("^", prefix, "-([0-9]+)[.].*$"), "\\1", files))
  paste0(prefix, "-", sprintf("%06d", if (length(ids)) max(ids) + 1L else 1L), ".", extension)
}

selection_checkpoint_entry <- function(io, file, kind) {
  path <- file.path(io$path, file)
  io$manifest$files[[length(io$manifest$files) + 1L]] <- list(file = file,
    kind = kind, bytes = as.double(file.info(path)$size), sha256 = io$utils$gwrs_sim_sha256(path))
}

selection_checkpoint_commit <- function(io, records) {
  file <- selection_checkpoint_next(io, "pack", "rds")
  io$utils$gwrs_sim_atomic_rds(list(signature = io$manifest$provenance$sha256,
    records = records), file.path(io$path, file))
  selection_checkpoint_entry(io, file, "pack")
  io$manifest$durable_count <- io$manifest$durable_count + length(records)
  io$manifest$search_status <- NULL
  selection_manifest_write(io)
}

selection_checkpoint_finish <- function(io, result) {
  file <- selection_checkpoint_next(io, "result", "rds")
  io$utils$gwrs_sim_atomic_rds(result, file.path(io$path, file))
  selection_checkpoint_entry(io, file, "result")
  file <- selection_checkpoint_next(io, "ranking", "csv")
  io$utils$gwrs_sim_atomic_csv(result$ranking, file.path(io$path, file))
  selection_checkpoint_entry(io, file, "ranking")
  io$manifest$search_status <- result$search_status
  selection_manifest_write(io)
}

selection_progress <- function(io, records, expected, evaluated, uncommitted,
    elapsed, phase, verbose, max_models) {
  tab <- selection_record_table(records)
  count <- length(records)
  throughput <- if (evaluated > 0 && elapsed > 0) evaluated / elapsed else NA_real_
  pending <- if (is.na(expected)) NA_real_ else max(0, expected - count)
  eta <- if (is.finite(throughput)) pending / throughput else NA_real_
  progress <- list(run_id = if (is.null(io)) "in_memory" else io$manifest$run_id,
    phase = phase, work_unit = "one subset/kernel/bandwidth candidate including diagnostics",
    total = expected, calculated = count,
    eligible = sum(tab$eligible), failed = sum(tab$status == "failed"),
    invalid = sum(tab$status == "invalid"), criterion_unavailable = sum(tab$status == "criterion_unavailable"),
    running = 0L, pending = pending,
    durable = if (is.null(io)) 0L else io$manifest$durable_count,
    uncommitted = if (is.null(io)) count else uncommitted,
    percent = if (is.finite(expected) && expected > 0) 100 * count / expected else NA_real_,
    elapsed_seconds = elapsed, new_evaluations = evaluated, candidates_per_second = throughput,
    search_eta_seconds = eta, eta_scope = "remaining candidate search; excludes final refits/output validation",
    estimated_search_finish_utc = if (is.finite(eta)) format(Sys.time() + eta, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC") else NA_character_,
    heartbeat_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    last_candidate = if (count) names(records)[count] else NA_character_, maximum_candidates = max_models)
  if (!is.null(io)) {
    io$utils$gwrs_sim_atomic_text(jsonlite::toJSON(progress, auto_unbox = TRUE, na = "null", pretty = TRUE),
      file.path(io$path, "progress.json"), replace = TRUE)
    path <- file.path(io$path, "progress.tsv")
    utils::write.table(as.data.frame(progress), path, sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = !file.exists(path), append = file.exists(path))
  }
  if (verbose) {
    cat(sprintf("[%s] %s: %d candidates, %d eligible, %d failed; durable %d; search ETA %s\n",
      progress$heartbeat_utc, phase, count, progress$eligible, progress$failed, progress$durable,
      if (is.finite(eta)) paste(round(eta, 1), "s") else "unknown"))
    utils::flush.console()
  }
  invisible(progress)
}

#' Verify a joint GWR selection checkpoint
#' @param checkpoint Directory created by [gwr_select()].
#' @return Verification summary with candidate counts, search status and
#'   provenance. All recorded packs/results are checked for size and SHA-256;
#'   errors stop verification. Uncommitted orphan files are not counted.
#' @details This is read-only and does not fit models. Completion and numerical
#'   eligibility are separate: a complete search can contain failed candidates.
#'   A stale run.lock is never removed automatically. After confirming its owner
#'   is no longer running, archive that lock directory before resuming.
#' @export
gwr_select_verify <- function(checkpoint) {
  saved <- selection_checkpoint_read(normalizePath(checkpoint, mustWork = TRUE))
  tab <- selection_record_table(saved$records)
  status <- saved$manifest$search_status
  if (!is.null(status) && (status$evaluated != nrow(tab) ||
      (isTRUE(status$complete) && is.finite(status$expected) && status$expected != nrow(tab))))
    stop("Final search status does not match verified candidate counts.", call. = FALSE)
  list(valid = TRUE, durable_candidates = nrow(tab), eligible = sum(tab$eligible),
    failed = sum(tab$status == "failed"), search_status = status,
    locked = dir.exists(file.path(checkpoint, "run.lock")), provenance = saved$manifest$provenance)
}
