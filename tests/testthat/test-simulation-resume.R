simulation_candidates <- c(
  file.path("inst", "simulations", "resume-utils.R"),
  file.path("..", "..", "inst", "simulations", "resume-utils.R")
)
# Direct source-tree tests must exercise the checkout, not an older installed
# copy. Package checks fall back to the installed inst/ payload.
simulation_utils <- if (any(file.exists(simulation_candidates))) {
  simulation_candidates[file.exists(simulation_candidates)][[1L]]
} else {
  system.file(
    "simulations", "resume-utils.R", package = "gwrs", mustWork = TRUE
  )
}
simulation_utils <- normalizePath(simulation_utils, mustWork = TRUE)
sys.source(simulation_utils, envir = environment())

simulation_hash_environment <- function(backend) {
  isolated <- new.env(parent = environment(gwrs_sim_run))
  sys.source(simulation_utils, envir = isolated)
  if (identical(backend, "tools")) {
    skip_if_not("sha256sum" %in% getNamespaceExports("tools"))
  } else {
    executable <- Sys.which(backend)
    skip_if_not(nzchar(executable))
    isolated$getNamespaceExports <- function(namespace) {
      if (identical(namespace, "tools")) character() else
        base::getNamespaceExports(namespace)
    }
    isolated$Sys.which <- function(command) {
      if (identical(command, backend)) executable else ""
    }
  }
  isolated
}

for (hash_backend in c("tools", "shasum", "sha256sum")) {
  local({
    backend <- hash_backend
    test_that(paste("SHA-256 returns unnamed values using", backend), {
      isolated <- simulation_hash_environment(backend)
      directory <- tempfile("gwrs-hash-names-")
      dir.create(directory)
      on.exit(unlink(directory, recursive = TRUE), add = TRUE)
      paths <- file.path(directory, c("abc with spaces.txt", "empty.txt"))
      writeBin(charToRaw("abc"), paths[[1L]])
      writeBin(raw(), paths[[2L]])
      expected <- c(
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
      )
      expect_identical(isolated$gwrs_sim_sha256(character()), character())
      expect_identical(isolated$gwrs_sim_sha256(paths[[1L]]), expected[[1L]])
      expect_identical(isolated$gwrs_sim_sha256(paths), expected)
      expect_identical(
        isolated$gwrs_sim_sha256(setNames(paths, c("first", "second"))),
        expected
      )
      manifest <- data.frame(file = paths, sha256 = isolated$gwrs_sim_sha256(paths))
      expect_identical(rownames(manifest), c("1", "2"))
    })
  })
}

simulation_fixture <- function(replications = 4L) {
  config <- list(
    schema = gwrs_sim_schema,
    study = "unit-test",
    n = 40L,
    p = 3L,
    k = 10L,
    replications = replications
  )
  plan <- gwrs_sim_make_plan(replications, master_seed = 73L)
  runtime <- data.frame(
    field = "engine", value = "unit-test", stringsAsFactors = FALSE
  )
  task_function <- function(task, config, heartbeat) {
    heartbeat("unit-computation")
    values <- rnorm(config$p)
    list(summary = data.frame(
      replication = task$replication[[1L]],
      draw_sum = sum(values),
      stringsAsFactors = FALSE
    ))
  }
  validate_result <- function(result, task) {
    is.list(result) && is.data.frame(result$summary) &&
      nrow(result$summary) == 1L &&
      identical(result$summary$replication, task$replication) &&
      is.finite(result$summary$draw_sum)
  }
  list(
    config = config,
    plan = plan,
    runtime = runtime,
    task_function = task_function,
    validate_result = validate_result
  )
}

for (hash_backend in c("shasum", "sha256sum")) {
  local({
    backend <- hash_backend
    test_that(paste("external hashing supports finalization and resume using", backend), {
      isolated <- simulation_hash_environment(backend)
      fixture <- simulation_fixture(replications = 2L)
      directory <- tempfile("gwrs-external-resume-")
      on.exit(unlink(directory, recursive = TRUE), add = TRUE)
      invoke <- function(task_function = fixture$task_function, max_tasks = Inf) {
        isolated$gwrs_sim_run(
          directory, fixture$config, fixture$plan, task_function,
          fixture$validate_result, fixture$runtime, max_tasks = max_tasks
        )
      }
      first <- invoke(max_tasks = 1L)
      expect_identical(first$tasks$status, c("completed", "pending"))
      first_path <- isolated$gwrs_sim_checkpoint_path(directory, fixture$plan$task_id[[1L]])
      first_hash <- isolated$gwrs_sim_sha256(first_path)
      first_bytes <- file.info(first_path)$size
      resumed <- invoke()
      expect_true(resumed$completion_valid)
      expect_true(resumed$manifest_valid)
      expect_identical(resumed$tasks$attempts, c(1L, 1L))
      expect_identical(isolated$gwrs_sim_sha256(first_path), first_hash)
      expect_identical(file.info(first_path)$size, first_bytes)

      manifest <- utils::read.csv(file.path(directory, "manifest-sha256.csv"))
      protected <- file.path(directory, manifest$file)
      protected_hashes <- isolated$gwrs_sim_sha256(protected)
      forbidden <- function(...) stop("A valid shard must not be recomputed")
      repeated <- invoke(task_function = forbidden)
      expect_true(repeated$completion_valid)
      expect_identical(isolated$gwrs_sim_sha256(protected), protected_hashes)

      # Test-fixture-only interruption: exercise final manifest row-name equality
      # when the marker must be recreated. Never used on an actual study.
      unlink(isolated$gwrs_sim_completion_path(directory))
      finalized <- invoke(task_function = forbidden)
      expect_true(finalized$completion_valid)
      expect_identical(isolated$gwrs_sim_sha256(protected), protected_hashes)

      writeBin(charToRaw("tampered fixture"), first_path)
      expect_false(isolated$gwrs_sim_verify_manifest(directory))
      expect_false(isolated$gwrs_sim_verify_completion(directory))
      expect_error(invoke(task_function = forbidden), "manifest")
    })
  })
}

test_that("simulation plans are deterministic and guarded", {
  first <- gwrs_sim_make_plan(5L, master_seed = 19L, scenario_id = "small")
  second <- gwrs_sim_make_plan(5L, master_seed = 19L,
                               scenario_id = "small")
  expect_identical(first, second)
  expect_identical(anyDuplicated(first$task_id), 0L)
  expect_identical(anyDuplicated(first$seed), 0L)

  fixture <- simulation_fixture()
  large <- fixture$config
  large$n <- 20001L
  expect_error(
    gwrs_sim_preflight(large, fixture$plan),
    "Large-run guard triggered"
  )
  expect_s3_class(
    gwrs_sim_preflight(large, fixture$plan, confirm_large = TRUE),
    "data.frame"
  )
})

test_that("runtime audit reports the RNG kind used inside each task", {
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
  RNGkind("Mersenne-Twister", "Box-Muller", "Rejection")

  expected <- paste(gwrs_sim_task_rng_kind(), collapse = "/")
  runtime <- gwrs_sim_runtime_signature("gwrs")
  task <- data.frame(seed = 1729L)
  outcome <- gwrs_sim_evaluate_task(
    task = task,
    config = list(),
    task_function = function(task, config, heartbeat) {
      list(rng_kind = paste(RNGkind(), collapse = "/"))
    },
    validate_result = function(result, task) TRUE,
    heartbeat = function(...) invisible(NULL)
  )

  expect_identical(
    runtime$value[runtime$field == "RNG_kind"], expected
  )
  expect_true(outcome$ok)
  expect_identical(outcome$result$rng_kind, expected)
  expect_identical(RNGkind(), c(
    "Mersenne-Twister", "Box-Muller", "Rejection"
  ))
})

test_that("heartbeat writes cached progress atomically and is throttled", {
  output_directory <- tempfile("gwrs-sim-heartbeat-")
  dir.create(output_directory)
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)
  task_status <- data.frame(
    task_id = c("task-1", "task-2"),
    seed = c(1L, 2L),
    status = c("completed", "pending"),
    attempts = c(1L, 0L),
    elapsed_seconds = c(2, NA_real_),
    checkpoint = c("task-1.rds", ""),
    bytes = c(100, NA_real_),
    sha256 = c("cached-hash", ""),
    stringsAsFactors = FALSE
  )

  heartbeat <- gwrs_sim_make_heartbeat(
    output_directory, task_status, current_task = "task-2",
    interval_seconds = 3600
  )
  heartbeat("first-phase")
  first <- utils::read.csv(
    file.path(output_directory, "progress.csv"), stringsAsFactors = FALSE
  )
  heartbeat("throttled-phase")
  throttled <- utils::read.csv(
    file.path(output_directory, "progress.csv"), stringsAsFactors = FALSE
  )

  expect_identical(first$phase, "first-phase")
  expect_identical(throttled$phase, "first-phase")
  expect_identical(first$completed_tasks, 1L)
  expect_identical(first$current_task_id, "task-2")
  expect_false(file.exists(file.path(output_directory, "task-status.csv")))
  expect_false(any(grepl("^[.]partial-", list.files(
    output_directory, all.files = TRUE
  ))))

  unthrottled <- gwrs_sim_make_heartbeat(
    output_directory, task_status, current_task = "task-2",
    interval_seconds = 0
  )
  unthrottled("fresh-phase")
  fresh <- utils::read.csv(
    file.path(output_directory, "progress.csv"), stringsAsFactors = FALSE
  )
  expect_identical(fresh$phase, "fresh-phase")
  expect_error(
    gwrs_sim_make_heartbeat(
      output_directory, task_status, "task-2", interval_seconds = -1
    ),
    "nonnegative"
  )
})

test_that("frequent heartbeats do not trigger full-plan scans", {
  fixture <- simulation_fixture(replications = 3L)
  output_directory <- tempfile("gwrs-sim-heartbeat-scans-")
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)
  original_scan <- gwrs_sim_scan
  scan_sizes <- integer()
  simulation_environment <- environment(gwrs_sim_run)
  assign(
    "gwrs_sim_scan",
    function(output_directory, plan, validate_result = NULL) {
      scan_sizes <<- c(scan_sizes, nrow(plan))
      original_scan(
        output_directory, plan, validate_result = validate_result
      )
    },
    envir = simulation_environment
  )
  on.exit(assign(
    "gwrs_sim_scan", original_scan, envir = simulation_environment
  ), add = TRUE)
  fixture$task_function <- function(task, config, heartbeat) {
    for (index in seq_len(100L)) {
      heartbeat(paste0("inner-phase-", index))
    }
    values <- rnorm(config$p)
    list(summary = data.frame(
      replication = task$replication[[1L]],
      draw_sum = sum(values),
      stringsAsFactors = FALSE
    ))
  }

  result <- gwrs_sim_run(
    output_directory,
    fixture$config,
    fixture$plan,
    fixture$task_function,
    fixture$validate_result,
    fixture$runtime,
    max_tasks = 1L
  )

  expect_false(result$finalized)
  expect_identical(scan_sizes, c(3L, 1L, 3L, 3L))
  expect_equal(sum(result$tasks$status == "completed"), 1L)
})

test_that("interrupted simulation resumes to deterministic results", {
  fixture <- simulation_fixture()
  # CSV readers infer an all-zero numeric column as integer. The resume
  # contract must restore the type declared by the current frozen plan.
  fixture$plan$zero_double <- rep(0, nrow(fixture$plan))
  expect_type(fixture$plan$zero_double, "double")
  resumed_directory <- tempfile("gwrs-sim-resumed-")
  fresh_directory <- tempfile("gwrs-sim-fresh-")
  on.exit(unlink(resumed_directory, recursive = TRUE), add = TRUE)
  on.exit(unlink(fresh_directory, recursive = TRUE), add = TRUE)

  set.seed(991L)
  rng_kind <- RNGkind()
  rng_state <- .Random.seed
  interrupted <- gwrs_sim_run(
    resumed_directory,
    fixture$config,
    fixture$plan,
    fixture$task_function,
    fixture$validate_result,
    fixture$runtime,
    max_tasks = 2L
  )
  expect_identical(RNGkind(), rng_kind)
  expect_identical(.Random.seed, rng_state)
  expect_equal(sum(interrupted$tasks$status == "completed"), 2L)
  expect_false(interrupted$finalized)
  expect_false(interrupted$completed)
  expect_false(interrupted$completion_valid)
  expect_false(interrupted$locked)

  resumed <- gwrs_sim_run(
    resumed_directory,
    fixture$config,
    fixture$plan,
    fixture$task_function,
    fixture$validate_result,
    fixture$runtime
  )
  fresh <- gwrs_sim_run(
    fresh_directory,
    fixture$config,
    fixture$plan,
    fixture$task_function,
    fixture$validate_result,
    fixture$runtime
  )
  expect_true(resumed$finalized)
  expect_true(resumed$manifest_valid)
  expect_true(resumed$completed)
  expect_true(resumed$completion_valid)
  expect_true(fresh$manifest_valid)
  expect_true(fresh$completion_valid)
  expected_metadata <- c(
    "session-info.txt", "scientific-signature.csv", "task-manifest.csv",
    "manifest-sha256.csv", "COMPLETED"
  )
  expect_true(all(file.exists(file.path(
    resumed_directory, expected_metadata
  ))))
  task_manifest <- utils::read.csv(
    file.path(resumed_directory, "task-manifest.csv"),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(task_manifest), nrow(fixture$plan))
  expect_true(all(task_manifest$status == "completed"))
  expect_true(all(task_manifest$bytes > 0))
  expect_true(all(nzchar(task_manifest$sha256)))
  final_manifest <- utils::read.csv(
    file.path(resumed_directory, "manifest-sha256.csv"),
    stringsAsFactors = FALSE
  )
  expect_false(any(final_manifest$file %in% c(
    "COMPLETED", "progress.csv", "task-status.csv"
  )))
  resumed_summary <- utils::read.csv(
    file.path(resumed_directory, "task-summary.csv")
  )
  fresh_summary <- utils::read.csv(
    file.path(fresh_directory, "task-summary.csv")
  )
  comparison_columns <- c("task_id", "seed", "replication", "draw_sum")
  expect_equal(
    resumed_summary[comparison_columns],
    fresh_summary[comparison_columns],
    tolerance = 0
  )
  expect_true(gwrs_sim_verify_manifest(resumed_directory))

  changed <- fixture$config
  changed$n <- changed$n + 1L
  expect_error(
    gwrs_sim_run(
      resumed_directory, changed, fixture$plan,
      fixture$task_function, fixture$validate_result, fixture$runtime
    ),
    "configuration differs"
  )
  changed_runtime <- fixture$runtime
  changed_runtime$value <- "changed-engine"
  expect_error(
    gwrs_sim_run(
      resumed_directory, fixture$config, fixture$plan,
      fixture$task_function, fixture$validate_result, changed_runtime
    ),
    "runtime differs"
  )
})

test_that("summary collection supports one or more rows per task", {
  one_row <- simulation_fixture(replications = 1L)
  one_row_directory <- tempfile("gwrs-sim-one-row-")
  multi_row_directory <- tempfile("gwrs-sim-multi-row-")
  on.exit(unlink(one_row_directory, recursive = TRUE), add = TRUE)
  on.exit(unlink(multi_row_directory, recursive = TRUE), add = TRUE)

  one_row_result <- gwrs_sim_run(
    one_row_directory,
    one_row$config,
    one_row$plan,
    one_row$task_function,
    one_row$validate_result,
    one_row$runtime
  )
  one_row_summary <- utils::read.csv(
    file.path(one_row_directory, "task-summary.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(one_row_result$completion_valid)
  expect_equal(nrow(one_row_summary), 1L)
  expect_identical(one_row_summary$task_id, one_row$plan$task_id)
  expect_identical(as.integer(one_row_summary$seed), one_row$plan$seed)

  multi_row <- simulation_fixture(replications = 2L)
  multi_row$config$study <- "multi-row-unit-test"
  multi_row$task_function <- function(task, config, heartbeat) {
    heartbeat("multi-row-computation")
    list(summary = data.frame(
      replication = rep(task$replication[[1L]], 2L),
      method = c("first", "second"),
      estimate = rnorm(2L),
      stringsAsFactors = FALSE
    ))
  }
  multi_row$validate_result <- function(result, task) {
    is.list(result) && is.data.frame(result$summary) &&
      nrow(result$summary) == 2L &&
      identical(result$summary$method, c("first", "second")) &&
      all(result$summary$replication == task$replication[[1L]]) &&
      all(is.finite(result$summary$estimate))
  }
  multi_row_result <- gwrs_sim_run(
    multi_row_directory,
    multi_row$config,
    multi_row$plan,
    multi_row$task_function,
    multi_row$validate_result,
    multi_row$runtime
  )
  multi_row_summary <- utils::read.csv(
    file.path(multi_row_directory, "task-summary.csv"),
    stringsAsFactors = FALSE
  )
  task_manifest <- utils::read.csv(
    file.path(multi_row_directory, "task-manifest.csv"),
    stringsAsFactors = FALSE
  )

  expect_true(multi_row_result$manifest_valid)
  expect_true(multi_row_result$completion_valid)
  expect_equal(nrow(multi_row_summary), 2L * nrow(multi_row$plan))
  expect_identical(
    multi_row_summary$task_id,
    rep(multi_row$plan$task_id, each = 2L)
  )
  expect_identical(
    as.integer(multi_row_summary$seed),
    rep(multi_row$plan$seed, each = 2L)
  )
  expect_identical(
    multi_row_summary$method,
    rep(c("first", "second"), nrow(multi_row$plan))
  )
  expect_equal(
    multi_row_summary$elapsed_seconds,
    rep(task_manifest$elapsed_seconds, each = 2L),
    tolerance = 0
  )
})

test_that("invalid checkpoints are quarantined and recomputed", {
  fixture <- simulation_fixture(replications = 2L)
  output_directory <- tempfile("gwrs-sim-corrupt-")
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)

  gwrs_sim_run(
    output_directory,
    fixture$config,
    fixture$plan,
    fixture$task_function,
    fixture$validate_result,
    fixture$runtime,
    max_tasks = 1L
  )
  checkpoint <- gwrs_sim_checkpoint_path(
    output_directory, fixture$plan$task_id[[1L]]
  )
  writeLines("not an RDS checkpoint", checkpoint)

  result <- gwrs_sim_run(
    output_directory,
    fixture$config,
    fixture$plan,
    fixture$task_function,
    fixture$validate_result,
    fixture$runtime
  )
  expect_true(result$manifest_valid)
  quarantined <- list.files(
    file.path(output_directory, "quarantine"),
    pattern = "[.]invalid[.]", full.names = TRUE
  )
  expect_length(quarantined, 1L)
  expect_true(file.exists(checkpoint))
})

test_that("failed tasks retain attempts and retry on resume", {
  fixture <- simulation_fixture(replications = 2L)
  output_directory <- tempfile("gwrs-sim-retry-")
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)
  attempts <- new.env(parent = emptyenv())
  transient_task <- function(task, config, heartbeat) {
    heartbeat("transient-test")
    task_id <- task$task_id[[1L]]
    if (!exists(task_id, envir = attempts, inherits = FALSE)) {
      assign(task_id, TRUE, envir = attempts)
      stop("intentional transient failure", call. = FALSE)
    }
    fixture$task_function(task, config, heartbeat)
  }

  first <- gwrs_sim_run(
    output_directory,
    fixture$config,
    fixture$plan,
    transient_task,
    fixture$validate_result,
    fixture$runtime
  )
  expect_false(first$finalized)
  expect_true(all(first$tasks$status == "failed"))
  expect_length(list.files(file.path(output_directory, "errors")), 2L)

  resumed <- gwrs_sim_run(
    output_directory,
    fixture$config,
    fixture$plan,
    transient_task,
    fixture$validate_result,
    fixture$runtime
  )
  expect_true(resumed$manifest_valid)
  expect_true(all(resumed$tasks$attempts == 2L))
})

test_that("run lock prevents concurrent writers", {
  fixture <- simulation_fixture(replications = 1L)
  output_directory <- tempfile("gwrs-sim-lock-")
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)
  gwrs_sim_initialize(
    output_directory, fixture$config, fixture$plan, fixture$runtime
  )
  lock <- gwrs_sim_acquire_lock(output_directory)
  on.exit(gwrs_sim_release_lock(lock), add = TRUE)
  owner <- readRDS(file.path(output_directory, "run.lock", "owner.rds"))
  expect_true(all(c(
    "slurm_job_id", "slurm_array_job_id", "slurm_array_task_id",
    "slurm_submit_dir"
  ) %in% names(owner)))

  status <- gwrs_sim_status(output_directory, quiet = TRUE)
  expect_true(status$locked)
  expect_error(
    gwrs_sim_run(
      output_directory, fixture$config, fixture$plan,
      fixture$task_function, fixture$validate_result, fixture$runtime
    ),
    "directory is locked"
  )
  expect_error(
    gwrs_sim_run(
      output_directory, fixture$config, fixture$plan,
      fixture$task_function, fixture$validate_result, fixture$runtime,
      recover_lock = TRUE
    ),
    "active process"
  )
})

test_that("lock recovery fails closed during owner publication races", {
  output_directory <- tempfile("gwrs-sim-lock-race-")
  dir.create(output_directory)
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)
  lock_path <- file.path(output_directory, "run.lock")
  dir.create(lock_path)

  expect_error(
    gwrs_sim_acquire_lock(output_directory, recover_lock = TRUE),
    "missing or invalid owner metadata"
  )
  expect_true(dir.exists(lock_path))
  expect_false(dir.exists(file.path(output_directory, "quarantine")))

  writeLines("incomplete owner publication", file.path(lock_path, "owner.rds"))
  expect_error(
    gwrs_sim_acquire_lock(output_directory, recover_lock = TRUE),
    "missing or invalid owner metadata"
  )
  expect_true(dir.exists(lock_path))
})

test_that("SLURM lock recovery requires a conclusive squeue result", {
  testthat::skip_on_os("windows")
  output_directory <- tempfile("gwrs-sim-slurm-lock-")
  bin_directory <- tempfile("gwrs-sim-empty-path-")
  dir.create(output_directory)
  dir.create(bin_directory)
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)
  on.exit(unlink(bin_directory, recursive = TRUE), add = TRUE)
  previous_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = previous_path), add = TRUE)
  lock_path <- file.path(output_directory, "run.lock")
  dir.create(lock_path)
  stale_owner <- list(
    schema = gwrs_sim_schema,
    token = "slurm-lock-test",
    pid = Sys.getpid(),
    host = Sys.info()[["nodename"]],
    slurm_job_id = "12345",
    started_utc = gwrs_sim_utc()
  )
  gwrs_sim_atomic_rds(
    stale_owner,
    file.path(lock_path, "owner.rds")
  )

  Sys.setenv(PATH = bin_directory)
  expect_error(
    gwrs_sim_acquire_lock(output_directory, recover_lock = TRUE),
    "squeue.*unavailable"
  )
  expect_true(dir.exists(lock_path))

  squeue <- file.path(bin_directory, "squeue")
  writeLines(c("#!/bin/sh", "exit 2"), squeue)
  Sys.chmod(squeue, mode = "0755")
  expect_error(
    gwrs_sim_acquire_lock(output_directory, recover_lock = TRUE),
    "squeue.*failed"
  )
  expect_true(dir.exists(lock_path))

  writeLines(c("#!/bin/sh", "echo RUNNING", "exit 0"), squeue)
  Sys.chmod(squeue, mode = "0755")
  expect_error(
    gwrs_sim_acquire_lock(output_directory, recover_lock = TRUE),
    "still present"
  )
  expect_true(dir.exists(lock_path))

  writeLines(c(
    "#!/bin/sh",
    "echo 'slurm_load_jobs error: Invalid job id specified' >&2",
    "exit 1"
  ), squeue)
  Sys.chmod(squeue, mode = "0755")
  invalid_id_recovered <- gwrs_sim_acquire_lock(
    output_directory, recover_lock = TRUE
  )
  expect_true(dir.exists(invalid_id_recovered$path))
  expect_true(gwrs_sim_release_lock(invalid_id_recovered))

  dir.create(lock_path)
  gwrs_sim_atomic_rds(stale_owner, file.path(lock_path, "owner.rds"))
  writeLines(c("#!/bin/sh", "exit 0"), squeue)
  Sys.chmod(squeue, mode = "0755")
  recovered <- gwrs_sim_acquire_lock(
    output_directory, recover_lock = TRUE
  )
  expect_true(dir.exists(recovered$path))
  expect_true(gwrs_sim_release_lock(recovered))
  expect_length(list.dirs(
    file.path(output_directory, "quarantine"),
    recursive = FALSE, full.names = TRUE
  ), 2L)
})

test_that("completion is bound to the final checksum manifest", {
  fixture <- simulation_fixture(replications = 1L)
  output_directory <- tempfile("gwrs-sim-completion-")
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)
  result <- gwrs_sim_run(
    output_directory,
    fixture$config,
    fixture$plan,
    fixture$task_function,
    fixture$validate_result,
    fixture$runtime
  )
  expect_true(result$completion_valid)
  summary_path <- file.path(output_directory, "task-summary.csv")
  write("tampered", file = summary_path, append = TRUE)
  expect_false(gwrs_sim_verify_manifest(output_directory))
  expect_false(gwrs_sim_verify_completion(output_directory))
})

test_that("a completed run safely recovers finalization tracking state", {
  fixture <- simulation_fixture(replications = 1L)
  output_directory <- tempfile("gwrs-sim-finalize-recovery-")
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)
  gwrs_sim_run(
    output_directory,
    fixture$config,
    fixture$plan,
    fixture$task_function,
    fixture$validate_result,
    fixture$runtime
  )
  unlink(gwrs_sim_completion_path(output_directory))

  lock_path <- file.path(output_directory, "run.lock")
  dir.create(lock_path)
  gwrs_sim_atomic_rds(
    list(
      schema = gwrs_sim_schema,
      token = "stale-finalize-lock",
      pid = 2147483647L,
      host = Sys.info()[["nodename"]],
      slurm_job_id = "",
      started_utc = gwrs_sim_utc()
    ),
    file.path(lock_path, "owner.rds")
  )
  progress_path <- file.path(output_directory, "progress.csv")
  progress <- utils::read.csv(progress_path, stringsAsFactors = FALSE)
  progress$state <- "finalizing"
  progress$phase <- "interrupted-before-completion-marker"
  gwrs_sim_atomic_csv(progress, progress_path, replace = TRUE)

  recovered <- gwrs_sim_run(
    output_directory,
    fixture$config,
    fixture$plan,
    fixture$task_function,
    fixture$validate_result,
    fixture$runtime,
    recover_lock = TRUE
  )
  expect_false(recovered$locked)
  expect_true(recovered$manifest_valid)
  expect_true(recovered$completion_valid)
  expect_identical(recovered$progress$state[[1L]], "complete")
  expect_length(list.dirs(
    file.path(output_directory, "quarantine"),
    recursive = FALSE, full.names = TRUE
  ), 1L)
})
