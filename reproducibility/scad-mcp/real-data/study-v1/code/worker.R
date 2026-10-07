arguments <- commandArgs(TRUE)
stopifnot(length(arguments) == 7L)
action <- arguments[[1L]]
app <- normalizePath(arguments[[2L]], mustWork = TRUE)
run <- normalizePath(arguments[[3L]], mustWork = FALSE)
mode <- arguments[[4L]]
subject <- arguments[[5L]]
implementation <- arguments[[6L]]
event_path <- arguments[[7L]]
stopifnot(
  action %in% c(
    "preflight", "initialize", "task", "scan", "select",
    "finalize", "verify"
  ),
  mode %in% c("smoke", "calibration", "full"),
  identical(
    readLines(
      file.path(app, "config/implementation-signature.txt"),
      warn = FALSE
    ),
    implementation
  )
)

source(file.path(app, "code/bootstrap.R"))
study <- usa_study_v1_load(app, mode)
e <- study$e
ctx <- study$ctx
data <- study$data
policy <- study$policy
grid <- e$usa_grid(ctx$cfg)
folds <- if (mode == "smoke") {
  as.integer(policy$smoke$outer_folds)
} else {
  as.integer(policy$full$outer_folds)
}
plan <- e$usa_task_plan(grid, folds, mode, policy$calibration)

if (mode == "full") {
  stopifnot(
    nrow(grid) == policy$full$grid_paths,
    length(e$usa_methods(grid)) == policy$full$method_families,
    sum(plan$stage == "nested_inner") == policy$full$nested_inner_tasks,
    sum(plan$stage == "outer_refit") == policy$full$outer_refit_tasks,
    sum(plan$stage == "full_tune") == policy$full$full_tuning_tasks,
    sum(plan$stage == "full_refit") == policy$full$full_refit_tasks,
    nrow(plan) == policy$full$total_tasks
  )
} else {
  stopifnot(
    nrow(plan) == if (mode == "smoke") {
      policy$smoke$tasks
    } else {
      policy$calibration$tasks
    }
  )
}
stopifnot(!anyDuplicated(plan$id))

binding <- jsonlite::fromJSON(file.path(app, "config/binding.json"))
data_hash <- e$usa_hash(data[c("x", "y", "coords", "ids")])
fold_hash <- e$usa_hash(data$fold)
scientific_identity <- list(
  schema = "gwrs-usa-counties-acs2024-study-v1-scientific-v1",
  implementation = implementation,
  mode = mode,
  dataset = if (mode == "smoke") "synthetic" else "ACS2024-counties",
  smoke_seed = if (mode == "smoke") policy$smoke$seed else NA_integer_,
  nonconvex_solver = policy$nonconvex_solver,
  convex_path_solver = policy$convex_path_solver,
  config = ctx$cfg,
  grid = grid,
  plan = plan,
  ids = data$ids,
  predictor_names = colnames(data$x),
  data_sha256 = data_hash,
  folds_sha256 = fold_hash,
  frozen_files = binding$files
)
scientific_signature <- e$usa_hash(scientific_identity)

usa_study_v1_hash_file <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

usa_study_v1_runtime_identity <- function() {
  packages <- c(
    "gwrs", "Rcpp", "RcppArmadillo", "RcppParallel", "FNN",
    "Matrix", "digest", "jsonlite"
  )
  dll_directory <- system.file("libs", package = "gwrs")
  dll <- list.files(
    dll_directory,
    pattern = paste0("[.]", sub("^[.]", "", .Platform$dynlib.ext), "$"),
    full.names = TRUE, recursive = TRUE
  )
  stopifnot(length(dll) == 1L)
  thread_names <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    "RCPP_PARALLEL_NUM_THREADS"
  )
  list(
    schema = "gwrs-usa-counties-acs2024-study-v1-runtime-v1",
    R = R.version.string,
    platform = R.version$platform,
    libraries = .libPaths(),
    packages = vapply(
      packages,
      function(name) as.character(utils::packageVersion(name)),
      ""
    ),
    gwrs_shared_library_sha256 = usa_study_v1_hash_file(dll),
    armadillo_index_bits = gwrs:::cpp_arma_uword_bits(),
    BLAS = unname(extSoftVersion()[["BLAS"]]),
    threads = as.list(Sys.getenv(thread_names))
  )
}

runtime_identity <- usa_study_v1_runtime_identity()
stopifnot(runtime_identity$armadillo_index_bits == 64L)
runtime_signature <- e$usa_hash(runtime_identity)
numeric <- file.path(run, "numeric")
task_path <- function(id) {
  file.path(numeric, "tasks", paste0(id, ".rds"))
}

usa_study_v1_publish_text <- function(text, path) {
  text <- enc2utf8(as.character(text))
  if (file.exists(path)) {
    stopifnot(identical(readLines(path, warn = FALSE), text))
  } else {
    ctx$io$gwrs_sim_atomic_text(text, path)
  }
  invisible(path)
}

usa_study_v1_publish_json <- function(value, path) {
  text <- unlist(strsplit(
    jsonlite::toJSON(
      value, auto_unbox = TRUE, pretty = TRUE, null = "null",
      na = "null", digits = NA
    ),
    "\n", fixed = TRUE
  ), use.names = FALSE)
  usa_study_v1_publish_text(text, path)
}

usa_study_v1_publish_csv <- function(value, path) {
  text <- capture.output(utils::write.csv(
    value, row.names = FALSE, na = ""
  ))
  usa_study_v1_publish_text(text, path)
}

usa_study_v1_receipted_object <- function(path) {
  stopifnot(file.exists(path), file.exists(paste0(path, ".sha256.rds")))
  stopifnot(identical(
    e$usa_files(path, dirname(path)),
    readRDS(paste0(path, ".sha256.rds"))
  ))
  readRDS(path)
}

outer_selection_path <- function(outer) {
  file.path(numeric, paste0("selection-o", outer, ".rds"))
}

full_selection_path <- function() {
  file.path(numeric, "selection-full.rds")
}

read_outer_selection <- function(outer) {
  usa_study_v1_receipted_object(outer_selection_path(outer))
}

read_full_selection <- function() {
  usa_study_v1_receipted_object(full_selection_path())
}

usa_study_v1_info <- function(id) {
  row <- plan[match(id, plan$id), , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  if (row$stage == "calibration") {
    split <- e$usa_nested_split(data, row$outer, row$inner)
    split$test <- head(split$test, policy$calibration$first_validation_targets)
    candidate <- grid[match(row$candidate_id, grid$id), , drop = FALSE]
    return(list(
      row = row, split = split, candidate = candidate,
      refit = FALSE, selected_index = NULL, full_refit = FALSE
    ))
  }
  if (row$stage == "nested_inner") {
    split <- e$usa_nested_split(data, row$outer, row$inner)
    candidate <- grid[match(row$candidate_id, grid$id), , drop = FALSE]
    return(list(
      row = row, split = split, candidate = candidate,
      refit = FALSE, selected_index = NULL, full_refit = FALSE
    ))
  }
  if (row$stage == "full_tune") {
    split <- e$usa_full_tune_split(data, row$validation_fold)
    candidate <- grid[match(row$candidate_id, grid$id), , drop = FALSE]
    return(list(
      row = row, split = split, candidate = candidate,
      refit = FALSE, selected_index = NULL, full_refit = FALSE
    ))
  }
  selection <- if (row$stage == "outer_refit") {
    read_outer_selection(row$outer)
  } else {
    stopifnot(row$stage == "full_refit")
    read_full_selection()
  }
  choice <- selection$best[
    paste(selection$best$where, selection$best$method, sep = "-") ==
      row$method_key,
    , drop = FALSE
  ]
  stopifnot(nrow(choice) == 1L)
  candidate <- grid[match(choice$id, grid$id), , drop = FALSE]
  full_refit <- row$stage == "full_refit"
  list(
    row = row,
    split = if (full_refit) {
      e$usa_full_refit_split(data)
    } else {
      e$usa_outer_split(data, row$outer)
    },
    candidate = candidate, refit = TRUE,
    selected_index = choice$index, full_refit = full_refit
  )
}

usa_study_v1_validate <- function(value, id) {
  info <- usa_study_v1_info(id)
  stopifnot(
    is.list(value$timing),
    length(value$timing$cpu_seconds) == 1L,
    length(value$timing$utc_elapsed_seconds) == 1L,
    is.finite(value$timing$cpu_seconds), value$timing$cpu_seconds >= 0,
    is.finite(value$timing$utc_elapsed_seconds),
    value$timing$utc_elapsed_seconds >= 0
  )
  e$usa_check_task(
    value, data, info$split, info$candidate,
    refit = info$refit, cfg = ctx$cfg, full_refit = info$full_refit
  )
  if (!is.null(value$error)) {
    stopifnot(
      !info$refit, is.character(value$error), length(value$error) == 1L
    )
    return(invisible(value))
  }
  score <- value$score
  stopifnot(
    !anyNA(score$valid),
    all(score$n == length(info$split$test)),
    all(is.finite(score$sse[score$valid])),
    all(score$sse[score$valid] >= 0),
    all(is.na(score$sse[!score$valid])),
    all(is.na(score$sae[!score$valid])),
    all(is.na(score$error_sum[!score$valid])),
    all(!score$valid | tail(value$attempts, 1L)[[1L]]$valid),
    !length(value$warnings) || !any(score$valid)
  )
  if (info$candidate$method %in% c("scad", "mcp")) {
    stopifnot(length(value$native_states) == length(value$attempts))
    for (index in seq_along(value$native_states)) {
      state <- value$native_states[[index]]
      stopifnot(
        ncol(state$path_state) == nrow(score),
        nrow(state$path_state) == if (info$candidate$where == "global") {
          1L
        } else {
          length(info$split$test)
        },
        identical(state$path_state == 1L, state$converged == 1L),
        all(state$path_state %in% 0:2),
        all(state$iterations <= value$attempts[[index]]$iterations_cap)
      )
    }
    stopifnot(all(
      !score$valid |
        apply(tail(value$native_states, 1L)[[1L]]$path_state == 1L,
              2L, all)
    ))
  }
  if (info$candidate$method %in% c("lasso", "en", "ridge")) {
    stopifnot(length(value$convex_native_states) == length(value$attempts))
    for (index in seq_along(value$convex_native_states)) {
      state <- value$convex_native_states[[index]]
      stopifnot(
        ncol(state$converged) == nrow(score),
        nrow(state$converged) == if (info$candidate$where == "global") {
          1L
        } else {
          length(info$split$test)
        },
        all(state$converged %in% 0:1),
        all(state$iterations <= value$attempts[[index]]$iterations_cap)
      )
    }
    stopifnot(all(
      !score$valid |
        apply(tail(value$convex_native_states, 1L)[[1L]]$converged == 1L,
              2L, all)
    ))
  }
  if (info$candidate$method %in% c("lasso", "en")) {
    diagnostics <- value$working_diagnostics
    stopifnot(is.data.frame(diagnostics), nrow(diagnostics) > 0L)
    final <- diagnostics[
      diagnostics$cap == max(diagnostics$cap), , drop = FALSE
    ]
    expected_targets <- if (info$candidate$where == "global") {
      1L
    } else {
      length(info$split$test)
    }
    stopifnot(
      nrow(final) == expected_targets * nrow(score),
      all(final$native_joint_certificate_pass[final$converged == 1L]),
      all(final$independent_kkt[final$converged == 1L] <=
            policy$independent_kkt_review_threshold),
      all(final$relative_duality_gap[final$converged == 1L] <=
            policy$duality_gap_review_relative_threshold),
      all(final$screen_fit_batches == final$working_attempts),
      all(final$screen_fit_batches <=
            policy$screen_fit_max_batches * final$screen_fit_events),
      all(final$final_polish_attempts == 0),
      all(final$final_polish_accepted == 0)
    )
  } else {
    stopifnot(is.null(value$working_diagnostics))
  }
  if (mode %in% c("calibration", "smoke")) stopifnot(any(score$valid))
  if (info$refit) {
    stopifnot(
      value$selected_index == info$selected_index,
      identical(value$training_mean, mean(data$y[info$split$train]))
    )
    local <- info$candidate$where == "local"
    coefficients <- value$coefficients
    products <- if (local) {
      data$x[info$split$test, , drop = FALSE] *
        coefficients[, -1L, drop = FALSE]
    } else {
      sweep(
        data$x[info$split$test, , drop = FALSE],
        2L, coefficients[1L, -1L], "*"
      )
    }
    reconstructed <- rowSums(products) + if (local) {
      coefficients[, 1L]
    } else {
      coefficients[1L, 1L]
    }
    stopifnot(all(
      abs(reconstructed - value$predictions) <=
        1e-8 * (1 + abs(reconstructed) + rowSums(abs(products)))
    ))
  }
  invisible(value)
}

usa_study_v1_read_task <- function(id, repair_receipt = TRUE) {
  path <- task_path(id)
  if (!file.exists(path)) stop("Missing task shard: ", id)
  if (repair_receipt && !file.exists(paste0(path, ".sha256.rds"))) {
    envelope <- readRDS(path)
    stopifnot(
      identical(envelope$id, id),
      identical(envelope$signature, scientific_signature),
      identical(
        envelope$payload_sha256, e$usa_hash(envelope$payload)
      )
    )
    e$usa_receipt(ctx, path)
  }
  value <- e$usa_read_shard(path, id, scientific_signature)
  usa_study_v1_validate(value, id)
  value
}

usa_study_v1_choose <- function(scope, outer = NULL) {
  stopifnot(scope %in% c("outer", "full"))
  validation_folds <- if (scope == "outer") {
    setdiff(folds, outer)
  } else {
    folds
  }
  values_by_candidate <- vector("list", nrow(grid))
  for (validation_fold in validation_folds) {
    for (candidate_index in seq_len(nrow(grid))) {
      id <- if (scope == "outer") {
        paste0(
          "o", outer, "-i", validation_fold, "-",
          grid$id[candidate_index]
        )
      } else {
        paste0("f", validation_fold, "-", grid$id[candidate_index])
      }
      value <- usa_study_v1_read_task(id)
      values_by_candidate[[candidate_index]][[
        length(values_by_candidate[[candidate_index]]) + 1L
      ]] <- value
    }
  }
  candidates <- lapply(seq_len(nrow(grid)), function(candidate_index) {
    values <- values_by_candidate[[candidate_index]]
    if (any(vapply(
      values, function(value) !is.null(value$error), TRUE
    ))) return(NULL)
    score <- e$usa_choose(lapply(values, `[[`, "score"))
    cbind(
      grid[candidate_index, , drop = FALSE][rep(1L, nrow(score)), ],
      score, row.names = NULL
    )
  })
  available <- Filter(Negate(is.null), candidates)
  if (!length(available)) stop("No complete tuning candidates for ", scope)
  candidates <- do.call(rbind, available)
  best <- lapply(e$usa_methods(grid), function(method_key) {
    candidate <- candidates[
      paste(candidates$where, candidates$method, sep = "-") ==
        method_key & candidates$valid & is.finite(candidates$mse),
      , drop = FALSE
    ]
    if (!nrow(candidate)) {
      stop("No valid tuning candidate for ", method_key, " in ", scope)
    }
    selected <- which(candidate$mse == min(candidate$mse))
    if (length(selected) != 1L) {
      stop("Exact tuning tie requires author review: ", method_key)
    }
    answer <- candidate[selected, , drop = FALSE]
    last <- if (answer$method %in% c("gwr", "ols")) {
      1L
    } else {
      ctx$cfg$n_lambda
    }
    answer$path_boundary <- answer$index %in% c(1L, last)
    answer
  })
  list(
    scope = scope,
    outer = if (scope == "outer") outer else NA_integer_,
    validation_folds = validation_folds,
    best = do.call(rbind, best),
    candidates = candidates
  )
}

usa_study_v1_stable_files <- function() {
  files <- list.files(
    run, recursive = TRUE, full.names = TRUE, all.files = TRUE,
    include.dirs = FALSE, no.. = TRUE
  )
  relative <- substring(files, nchar(run) + 2L)
  excluded <- relative %in% c(
    "controller.lock", "progress.json", "progress.tsv",
    "manifest-sha256.csv", "COMPLETED"
  ) | grepl("^(invocations|controller-logs)/", relative)
  files[!excluded]
}

usa_study_v1_completion <- function() {
  path <- file.path(run, "COMPLETED")
  if (!file.exists(path)) return(NULL)
  tryCatch(
    as.list(as.data.frame(read.dcf(path), stringsAsFactors = FALSE)[1L, ]),
    error = function(error) NULL
  )
}

usa_study_v1_verify_manifest <- function() {
  manifest_path <- file.path(run, "manifest-sha256.csv")
  completion <- usa_study_v1_completion()
  stopifnot(file.exists(manifest_path), !is.null(completion))
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  manifest$bytes <- as.numeric(manifest$bytes)
  current <- e$usa_files(usa_study_v1_stable_files(), run)
  stopifnot(
    identical(manifest, current),
    identical(
      completion$Schema,
      "gwrs-usa-counties-acs2024-study-v1-run-v1"
    ),
    identical(completion$`Scientific-Signature`, scientific_signature),
    identical(completion$`Implementation-Signature`, implementation),
    identical(completion$`Runtime-Signature`, runtime_signature),
    as.integer(completion$`Total-Tasks`) == nrow(plan),
    identical(
      completion$`Manifest-SHA256`,
      usa_study_v1_hash_file(manifest_path)
    )
  )
  invisible(TRUE)
}

usa_study_v1_finalize <- function(write_new = TRUE) {
  required <- c(
    file.path(numeric, "task-summary.rds"),
    file.path(numeric, "task-summary.csv"),
    file.path(numeric, "task-manifest.csv"),
    file.path(numeric, "review-status.json"),
    file.path(numeric, "REPORT.md"),
    file.path(run, "manifest-sha256.csv"),
    file.path(run, "COMPLETED")
  )
  if (mode %in% c("smoke", "full")) {
    outer_paths <- file.path(
      numeric, paste0("selection-o", folds, ".rds")
    )
    required <- c(
      required, outer_paths, paste0(outer_paths, ".sha256.rds"),
      full_selection_path(), paste0(full_selection_path(), ".sha256.rds"),
      file.path(numeric, "outer-selection-summary.rds"),
      file.path(numeric, "outer-selection-summary.csv"),
      file.path(numeric, "full-selection-summary.rds"),
      file.path(numeric, "full-selection-summary.csv"),
      file.path(numeric, "oof-predictions.rds"),
      file.path(numeric, "oof-predictions.csv"),
      file.path(numeric, "oof-metrics.rds"),
      file.path(numeric, "oof-metrics.csv"),
      file.path(numeric, "full-fit-index.rds"),
      file.path(numeric, "full-fit-index.csv"),
      file.path(numeric, "gwr-full-model-diagnostics.rds"),
      file.path(numeric, "gwr-full-local-inference.rds"),
      file.path(numeric, "gwr-full-f-tests.rds"),
      file.path(numeric, "gwr-full-local-inference.csv"),
      file.path(numeric, "gwr-full-f1-f2.csv"),
      file.path(numeric, "gwr-full-f3.csv")
    )
    if (mode == "full") {
      required <- c(
        required,
        file.path(numeric, "oof-residual-spatial.rds"),
        file.path(numeric, "oof-moran-global.csv"),
        file.path(numeric, "oof-lisa-local.csv")
      )
    }
  }
  if (!write_new) {
    missing <- required[!file.exists(required)]
    if (length(missing)) {
      stop(
        "Read-only verification requires existing finalized outputs: ",
        paste(basename(missing), collapse = ", ")
      )
    }
  }

  values <- lapply(plan$id, usa_study_v1_read_task)
  names(values) <- plan$id
  task_summary <- do.call(rbind, lapply(seq_along(values), function(index) {
    value <- values[[index]]
    data.frame(
      task = plan$id[index], stage = plan$stage[index],
      outer = plan$outer[index], inner = plan$inner[index],
      validation_fold = plan$validation_fold[index],
      method_key = plan$method_key[index],
      cpu_seconds = value$timing$cpu_seconds,
      utc_seconds = value$timing$utc_elapsed_seconds,
      task_error = !is.null(value$error),
      scored_parameters = if (is.null(value$score)) 0L else nrow(value$score),
      valid_parameters = if (is.null(value$score)) 0L else
        sum(value$score$valid),
      invalid_parameters = if (is.null(value$score)) 0L else
        sum(!value$score$valid),
      stringsAsFactors = FALSE
    )
  }))
  task_manifest <- do.call(rbind, lapply(seq_len(nrow(plan)), function(index) {
    path <- task_path(plan$id[index])
    envelope <- readRDS(path)
    data.frame(
      plan[index, , drop = FALSE], status = "completed",
      task_error = !is.null(envelope$payload$error),
      completed_utc = envelope$payload$completed_utc,
      slurm_job_id = envelope$payload$slurm_job_id,
      bytes = unname(file.info(path)$size),
      sha256 = usa_study_v1_hash_file(path),
      receipt = basename(paste0(path, ".sha256.rds")),
      stringsAsFactors = FALSE
    )
  }))
  review <- list(
    schema = "gwrs-usa-counties-acs2024-study-v1-review-v1",
    mode = mode,
    tasks = nrow(plan),
    valid_parameters = sum(task_summary$valid_parameters),
    invalid_parameters = sum(task_summary$invalid_parameters),
    scored_parameters = sum(task_summary$scored_parameters),
    task_errors = sum(task_summary$task_error),
    all_parameters_converged =
      sum(task_summary$invalid_parameters) == 0L &&
      sum(task_summary$task_error) == 0L,
    interpretation = if (mode == "calibration") {
      "Thirty-two-target Linux resource calibration only; not a CV result."
    } else if (mode == "smoke") {
      "Synthetic workflow validation only; not an empirical result."
    } else {
      paste(
        "Closed nested spatial OOF, common-graph Moran/LISA, and descriptive",
        "full-refit analysis; publication rendering remains separate."
      )
    }
  )

  if (mode %in% c("smoke", "full")) {
    outer_selection <- do.call(rbind, lapply(folds, function(outer) {
      cbind(outer = outer, read_outer_selection(outer)$best, row.names = NULL)
    }))
    full_selection <- read_full_selection()$best
    methods <- e$usa_methods(grid)
    predictions <- do.call(rbind, lapply(methods, function(method_key) {
      do.call(rbind, lapply(folds, function(outer) {
        value <- values[[paste0("o", outer, "-refit-", method_key)]]
        data.frame(
          sample_id = value$test_ids, outer_fold = outer,
          method = method_key, observed = value$observed,
          prediction = value$predictions,
          baseline = rep(value$training_mean, length(value$test_ids)),
          stringsAsFactors = FALSE
        )
      }))
    }))
    metrics <- do.call(rbind, lapply(methods, function(method_key) {
      value <- predictions[
        predictions$method == method_key, , drop = FALSE
      ]
      stopifnot(
        nrow(value) == length(data$ids), !anyDuplicated(value$sample_id),
        setequal(value$sample_id, data$ids), all(is.finite(value$prediction))
      )
      error <- value$prediction - value$observed
      sse <- sum(error^2)
      data.frame(
        method = method_key, n = nrow(value),
        rmse = sqrt(sse / nrow(value)), mae = mean(abs(error)),
        bias = mean(error),
        predictive_R2 = 1 - sse /
          sum((value$observed - value$baseline)^2),
        stringsAsFactors = FALSE
      )
    }))
    full_fit_index <- do.call(rbind, lapply(methods, function(method_key) {
      id <- paste0("full-refit-", method_key)
      value <- values[[id]]
      choice <- full_selection[
        paste(full_selection$where, full_selection$method, sep = "-") ==
          method_key,
        , drop = FALSE
      ]
      stopifnot(nrow(choice) == 1L, value$selected_index == choice$index)
      data.frame(
        task = id, method = method_key, candidate_id = choice$id,
        k = choice$k, alpha = choice$alpha,
        selected_index = choice$index,
        selected_lambda = tail(value$score$lambda, 1L),
        coefficient_rows = nrow(value$coefficients),
        coefficient_columns = ncol(value$coefficients),
        selected_fraction = if (is.null(value$selected)) NA_real_ else
          mean(value$selected),
        stringsAsFactors = FALSE
      )
    }))
    expected_n <- if (mode == "full") {
      policy$full$expected_oof_predictions_per_method
    } else {
      policy$smoke$n
    }
    stopifnot(all(metrics$n == expected_n))
    review$methods <- nrow(metrics)
    review$oof_predictions_per_method <-
      as.integer(table(predictions$method))
    review$outer_boundary_selected <- sum(outer_selection$path_boundary)
    review$full_boundary_selected <- sum(full_selection$path_boundary)
    e$usa_put(
      ctx, outer_selection,
      file.path(numeric, "outer-selection-summary.rds")
    )
    e$usa_put(
      ctx, full_selection,
      file.path(numeric, "full-selection-summary.rds")
    )
    e$usa_put(ctx, predictions, file.path(numeric, "oof-predictions.rds"))
    e$usa_put(ctx, metrics, file.path(numeric, "oof-metrics.rds"))
    e$usa_put(ctx, full_fit_index, file.path(numeric, "full-fit-index.rds"))
    usa_study_v1_publish_csv(
      outer_selection,
      file.path(numeric, "outer-selection-summary.csv")
    )
    usa_study_v1_publish_csv(
      full_selection,
      file.path(numeric, "full-selection-summary.csv")
    )
    usa_study_v1_publish_csv(
      predictions, file.path(numeric, "oof-predictions.csv")
    )
    usa_study_v1_publish_csv(
      metrics, file.path(numeric, "oof-metrics.csv")
    )
    usa_study_v1_publish_csv(
      full_fit_index, file.path(numeric, "full-fit-index.csv")
    )

    gwr_value <- values[["full-refit-local-gwr"]]
    gwr_diagnostics <- gwr_value$full_gwr
    stopifnot(is.list(gwr_diagnostics))
    e$usa_put(
      ctx, gwr_diagnostics$model_diagnostics,
      file.path(numeric, "gwr-full-model-diagnostics.rds")
    )
    e$usa_put(
      ctx, gwr_diagnostics$local_inference,
      file.path(numeric, "gwr-full-local-inference.rds")
    )
    e$usa_put(
      ctx, gwr_diagnostics$diagnostic_tests,
      file.path(numeric, "gwr-full-f-tests.rds")
    )
    usa_study_v1_publish_csv(
      gwr_diagnostics$local_inference$local,
      file.path(numeric, "gwr-full-local-inference.csv")
    )
    usa_study_v1_publish_csv(
      gwr_diagnostics$diagnostic_tests$global_tests,
      file.path(numeric, "gwr-full-f1-f2.csv")
    )
    usa_study_v1_publish_csv(
      gwr_diagnostics$diagnostic_tests$coefficient_tests,
      file.path(numeric, "gwr-full-f3.csv")
    )
    if (mode == "full") {
      residual_spatial <- e$usa_residual_spatial(ctx, predictions, data)
      e$usa_put(
        ctx, residual_spatial,
        file.path(numeric, "oof-residual-spatial.rds")
      )
      usa_study_v1_publish_csv(
        residual_spatial$global,
        file.path(numeric, "oof-moran-global.csv")
      )
      usa_study_v1_publish_csv(
        residual_spatial$local,
        file.path(numeric, "oof-lisa-local.csv")
      )
      review$residual_spatial_methods <- nrow(residual_spatial$global)
      review$residual_spatial_permutations <- residual_spatial$permutations
      review$residual_spatial_seed <- residual_spatial$seed
    }
  }

  if (mode %in% c("smoke", "calibration")) {
    stopifnot(review$task_errors == 0L, review$invalid_parameters == 0L)
  }
  e$usa_put(ctx, task_summary, file.path(numeric, "task-summary.rds"))
  usa_study_v1_publish_csv(
    task_summary, file.path(numeric, "task-summary.csv")
  )
  usa_study_v1_publish_csv(
    task_manifest, file.path(numeric, "task-manifest.csv")
  )
  usa_study_v1_publish_json(
    review, file.path(numeric, "review-status.json")
  )
  usa_study_v1_publish_text(c(
    "# gwrs USA county study-v1 numerical closure", "",
    paste0("Mode: ", mode),
    paste0("Tasks: ", nrow(plan)),
    paste0(
      "Valid scored parameters: ", review$valid_parameters,
      "/", review$scored_parameters
    ),
    paste0("Retained invalid parameters: ", review$invalid_parameters),
    paste0("Retained task errors: ", review$task_errors), "",
    review$interpretation
  ), file.path(numeric, "REPORT.md"))

  manifest_path <- file.path(run, "manifest-sha256.csv")
  manifest <- e$usa_files(usa_study_v1_stable_files(), run)
  if (file.exists(manifest_path)) {
    recorded <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    recorded$bytes <- as.numeric(recorded$bytes)
    stopifnot(identical(recorded, manifest))
  } else {
    usa_study_v1_publish_csv(manifest, manifest_path)
  }
  completion_path <- file.path(run, "COMPLETED")
  if (!file.exists(completion_path)) {
    usa_study_v1_publish_text(c(
      "Schema: gwrs-usa-counties-acs2024-study-v1-run-v1",
      paste0(
        "Completed-UTC: ",
        format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
      ),
      paste0("Mode: ", mode),
      paste0("Total-Tasks: ", nrow(plan)),
      paste0("Scientific-Signature: ", scientific_signature),
      paste0("Implementation-Signature: ", implementation),
      paste0("Runtime-Signature: ", runtime_signature),
      paste0(
        "Manifest-SHA256: ", usa_study_v1_hash_file(manifest_path)
      ),
      paste0("Invalid-Parameters: ", review$invalid_parameters),
      paste0("Task-Errors: ", review$task_errors)
    ), completion_path)
  }
  usa_study_v1_verify_manifest()
  invisible(review)
}

if (action == "preflight") {
  stopifnot(
    nrow(data$x) == ctx$cfg$n, ncol(data$x) == ctx$cfg$p,
    identical(colnames(data$x), data$predictors)
  )
  split_sizes <- unlist(lapply(folds, function(outer) {
    vapply(setdiff(folds, outer), function(inner) {
      split <- e$usa_nested_split(data, outer, inner)
      stopifnot(
        max(ctx$cfg$k) <= length(split$train),
        qr(cbind(1, data$x[split$train, , drop = FALSE]))$rank ==
          ncol(data$x) + 1L
      )
      length(split$train)
    }, integer(1L))
  }))
  cat(jsonlite::toJSON(list(
    schema = "gwrs-usa-counties-acs2024-study-v1-preflight-v1",
    mode = mode, n = nrow(data$x), p = ncol(data$x),
    tasks = nrow(plan),
    nested_inner_tasks = sum(plan$stage == "nested_inner"),
    outer_refit_tasks = sum(plan$stage == "outer_refit"),
    full_tuning_tasks = sum(plan$stage == "full_tune"),
    full_refit_tasks = sum(plan$stage == "full_refit"),
    minimum_nested_training_n = min(split_sizes),
    methods = e$usa_methods(grid),
    scientific_signature = scientific_signature,
    runtime_signature = runtime_signature,
    external_workers = Sys.getenv("GWRS_STUDY_WORKERS", unset = "1"),
    native_threads_per_worker = 1,
    model_fits = 0
  ), auto_unbox = TRUE, pretty = TRUE), "\n")
  cat("Frozen ACS data, folds, grid, package and runtime preflight: OK. No model fitted.\n")
  quit(save = "no", status = 0L)
}

dir.create(numeric, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(numeric, "tasks"), recursive = TRUE, showWarnings = FALSE)

if (action == "initialize") {
  e$usa_put(ctx, scientific_identity, file.path(numeric, "identity.rds"))
  e$usa_put(ctx, runtime_identity, file.path(numeric, "runtime-identity.rds"))
  e$usa_put(ctx, plan, file.path(numeric, "task-plan.rds"))
  usa_study_v1_publish_csv(plan, file.path(numeric, "task-plan.csv"))
  usa_study_v1_publish_text(
    scientific_signature, file.path(numeric, "scientific-signature.txt")
  )
  usa_study_v1_publish_text(
    runtime_signature, file.path(numeric, "runtime-signature.txt")
  )
  usa_study_v1_publish_csv(data.frame(
    field = c(
      "scientific_signature", "implementation_signature",
      "runtime_signature", "data_sha256", "folds_sha256",
      "grid_sha256", "plan_sha256"
    ),
    value = c(
      scientific_signature, implementation, runtime_signature,
      data_hash, fold_hash, e$usa_hash(grid), e$usa_hash(plan)
    ),
    stringsAsFactors = FALSE
  ), file.path(numeric, "scientific-signature.csv"))
  session <- c(
    paste0(
      "recorded_utc: ",
      format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    ),
    paste0("working_directory: ", normalizePath(getwd())),
    paste0("SLURM_JOB_ID=", Sys.getenv("SLURM_JOB_ID")),
    paste0("SLURM_CPUS_PER_TASK=", Sys.getenv("SLURM_CPUS_PER_TASK")),
    paste0("GWRS_STUDY_WORKERS=", Sys.getenv("GWRS_STUDY_WORKERS")),
    paste0(
      "GWRS_STUDY_POST_THREADS=",
      Sys.getenv("GWRS_STUDY_POST_THREADS", unset = "1")
    ),
    paste0("R_LIBS_USER=", Sys.getenv("R_LIBS_USER")),
    paste0("R_LIBS=", Sys.getenv("R_LIBS")),
    capture.output(sessionInfo()), capture.output(.libPaths())
  )
  session_path <- file.path(numeric, "session-info.txt")
  if (!file.exists(session_path)) {
    usa_study_v1_publish_text(session, session_path)
  }
  cat("Initialized", mode, "run with", nrow(plan), "tasks. No model fitted.\n")
  quit(save = "no", status = 0L)
}

stopifnot(
  identical(readRDS(file.path(numeric, "identity.rds")), scientific_identity),
  identical(
    readRDS(file.path(numeric, "runtime-identity.rds")), runtime_identity
  ),
  identical(readRDS(file.path(numeric, "task-plan.rds")), plan)
)

if (action == "task") {
  info <- usa_study_v1_info(subject)
  started_cpu <- sum(proc.time()[c("user.self", "sys.self")])
  started_utc <- as.numeric(Sys.time())
  train <- e$usa_training(ctx, data, info$split)
  value <- e$usa_task(ctx, numeric, subject, scientific_signature, function() {
    result <- e$usa_evaluate(
      ctx, train, data, info$split, info$candidate,
      until = if (info$refit) info$selected_index else NULL,
      keep = info$refit,
      full_diagnostics = info$full_refit &&
        info$candidate$where == "local" &&
        info$candidate$method == "gwr"
    )
    result$timing <- list(
      cpu_seconds = unname(
        sum(proc.time()[c("user.self", "sys.self")]) - started_cpu
      ),
      utc_elapsed_seconds = as.numeric(Sys.time()) - started_utc,
      scope = paste(
        "one model path including training setup and inline numerical",
        "certificates; publication rendering excluded"
      )
    )
    result$completed_utc <- format(
      Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
    )
    result$slurm_job_id <- Sys.getenv("SLURM_JOB_ID")
    result$slurm_step_id <- Sys.getenv("SLURM_STEP_ID")
    usa_study_v1_validate(result, subject)
    result
  })
  usa_study_v1_validate(value, subject)
  usa_study_v1_publish_json(list(
    schema = "gwrs-usa-counties-acs2024-study-v1-task-event-v1",
    id = subject, stage = info$row$stage,
    cost_class = info$row$cost_class,
    elapsed_seconds = value$timing$utc_elapsed_seconds,
    cpu_seconds = value$timing$cpu_seconds,
    invalid_parameters = if (is.null(value$score)) 0L else
      sum(!value$score$valid),
    task_error = !is.null(value$error),
    scientific_signature = scientific_signature,
    shard_sha256 = usa_study_v1_hash_file(task_path(subject))
  ), event_path)
  quit(save = "no", status = 0L)
}

if (action == "scan") {
  rows <- if (subject == "all") {
    plan
  } else {
    plan[plan$stage == subject, , drop = FALSE]
  }
  stopifnot(nrow(rows) > 0L)
  valid <- list()
  missing <- character()
  invalid <- list()
  for (index in seq_len(nrow(rows))) {
    id <- rows$id[index]
    if (!file.exists(task_path(id))) {
      missing <- c(missing, id)
      next
    }
    result <- tryCatch(
      list(ok = TRUE, value = usa_study_v1_read_task(id)),
      error = function(error) {
        list(ok = FALSE, message = conditionMessage(error))
      }
    )
    if (!result$ok) {
      invalid[[length(invalid) + 1L]] <- list(
        id = id, message = result$message
      )
    } else {
      value <- result$value
      valid[[length(valid) + 1L]] <- list(
        id = id, stage = rows$stage[index],
        cost_class = rows$cost_class[index],
        elapsed_seconds = value$timing$utc_elapsed_seconds,
        invalid_parameters = if (is.null(value$score)) 0L else
          sum(!value$score$valid),
        task_error = !is.null(value$error)
      )
    }
  }
  usa_study_v1_publish_json(list(
    schema = "gwrs-usa-counties-acs2024-study-v1-scan-v1",
    requested_stage = subject, total = nrow(rows), valid = valid,
    missing = unname(missing), invalid = invalid,
    scientific_signature = scientific_signature
  ), event_path)
  quit(save = "no", status = if (length(invalid)) 2L else 0L)
}

if (action == "select") {
  if (identical(subject, "full")) {
    choice <- usa_study_v1_choose("full")
    path <- full_selection_path()
    scope <- "full"
    outer <- NA_integer_
  } else {
    outer <- as.integer(subject)
    stopifnot(!is.na(outer), outer %in% folds)
    choice <- usa_study_v1_choose("outer", outer)
    path <- outer_selection_path(outer)
    scope <- "outer"
  }
  e$usa_put(ctx, choice, path)
  if (!file.exists(paste0(path, ".sha256.rds"))) e$usa_receipt(ctx, path)
  stopifnot(identical(usa_study_v1_receipted_object(path), choice))
  usa_study_v1_publish_json(list(
    schema = "gwrs-usa-counties-acs2024-study-v1-selection-event-v1",
    scope = scope, outer = outer, selected = nrow(choice$best),
    boundary_selected = sum(choice$best$path_boundary),
    scientific_signature = scientific_signature,
    selection_sha256 = usa_study_v1_hash_file(path)
  ), event_path)
  quit(save = "no", status = 0L)
}

if (action == "finalize") {
  review <- usa_study_v1_finalize(TRUE)
  usa_study_v1_publish_json(list(
    schema = "gwrs-usa-counties-acs2024-study-v1-final-event-v1",
    mode = mode, tasks = review$tasks,
    valid_parameters = review$valid_parameters,
    invalid_parameters = review$invalid_parameters,
    task_errors = review$task_errors,
    scientific_signature = scientific_signature,
    manifest_sha256 = usa_study_v1_hash_file(
      file.path(run, "manifest-sha256.csv")
    )
  ), event_path)
  quit(save = "no", status = 0L)
}

if (action == "verify") {
  usa_study_v1_finalize(FALSE)
  usa_study_v1_verify_manifest()
  cat("All tasks, selections, OOF/full-refit coverage, signatures and SHA-256 manifest: OK.\n")
  cat("Read-only verification fitted no model.\n")
  quit(save = "no", status = 0L)
}

stop("Unhandled action")
