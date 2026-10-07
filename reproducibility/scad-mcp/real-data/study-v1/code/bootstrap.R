usa_study_v1_read_data <- function(app) {
  input_root <- Sys.getenv(
    "GWRS_STUDY_INPUT_ROOT", unset = file.path(app, "frozen/input")
  )
  model_path <- file.path(
    input_root, "spatial/data/county-model-data-v1.csv"
  )
  fold_path <- file.path(
    input_root, "spatial/tables/spatial-folds-v1.csv"
  )
  catalog_path <- file.path(
    input_root, "config/variable-catalog-v1.csv"
  )
  model <- utils::read.csv(
    model_path, check.names = FALSE,
    colClasses = c(
      fips = "character", state_fips = "character",
      county_fips = "character"
    )
  )
  folds <- utils::read.csv(
    fold_path, check.names = FALSE,
    colClasses = c(fips = "character", state_fips = "character")
  )
  catalog <- utils::read.csv(catalog_path, check.names = FALSE)
  response <- catalog$alias[catalog$role == "response"]
  predictors <- catalog$alias[catalog$role == "predictor"]
  stopifnot(
    length(response) == 1L,
    length(predictors) == 26L,
    nrow(model) == 3107L,
    identical(model$fips, folds$fips),
    identical(as.integer(model$spatial_fold), as.integer(folds$fold)),
    all(c(response, predictors, "easting_m", "northing_m") %in%
          names(model))
  )
  x <- as.matrix(model[, predictors, drop = FALSE])
  storage.mode(x) <- "double"
  y <- as.double(model[[response]])
  coords <- as.matrix(model[, c("easting_m", "northing_m")])
  storage.mode(coords) <- "double"
  stopifnot(
    all(is.finite(x)), all(is.finite(y)), all(is.finite(coords)),
    ncol(x) == 26L, !anyDuplicated(model$fips),
    setequal(unique(as.integer(folds$fold)), 1:5)
  )
  list(
    x = x, y = y, coords = coords, ids = model$fips,
    fold = as.integer(folds$fold), model = model, catalog = catalog,
    response = response, predictors = predictors
  )
}

usa_study_v1_smoke_data <- function(policy) {
  smoke <- policy$smoke
  set.seed(smoke$seed)
  x <- matrix(rnorm(smoke$n * smoke$p), smoke$n, smoke$p)
  colnames(x) <- paste0("x", seq_len(smoke$p))
  coords <- cbind(runif(smoke$n), runif(smoke$n))
  y <- 1.5 + 2.5 * x[, 1L] - 1.5 * x[, 2L] +
    0.4 * coords[, 1L] * x[, 3L] + rnorm(smoke$n, sd = 0.15)
  fold <- rep(smoke$outer_folds, length.out = smoke$n)
  ids <- sprintf("smoke-%03d", seq_len(smoke$n))
  list(
    x = x, y = y, coords = coords, ids = ids,
    fold = as.integer(fold), model = NULL, catalog = NULL,
    response = "synthetic_y", predictors = colnames(x)
  )
}

usa_study_v1_load <- function(app, mode) {
  app <- normalizePath(app, mustWork = TRUE)
  stopifnot(mode %in% c("smoke", "calibration", "full"))
  policy <- jsonlite::fromJSON(file.path(app, "config/policy.json"))
  cfg <- jsonlite::fromJSON(file.path(app, "config/numerics.json"))
  thread_names <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    "RCPP_PARALLEL_NUM_THREADS"
  )
  stopifnot(all(Sys.getenv(thread_names) == "1"))

  package_library <- Sys.getenv("GWRS_STUDY_R_LIB")
  stopifnot(nzchar(package_library), dir.exists(package_library))
  .libPaths(unique(c(package_library, .libPaths())))
  library(gwrs, lib.loc = package_library, character.only = FALSE)
  stopifnot(
    as.character(utils::packageVersion("gwrs")) ==
      policy$candidate_package_version,
    normalizePath(find.package("gwrs")) ==
      normalizePath(file.path(package_library, "gwrs"))
  )

  environment <- new.env(parent = globalenv())
  sys.source(file.path(app, "code/core.R"), environment)
  sys.source(file.path(app, "code/diagnostics.R"), environment)
  io <- new.env(parent = globalenv())
  sys.source(file.path(app, "frozen/resume-utils.R"), io)

  package_namespace <- asNamespace("gwrs")
  api <- new.env(parent = package_namespace)
  required <- c(
    "standardize_design", "apply_standardization",
    "backtransform_coefficients", "gwr_neighbors",
    "resolve_nonconvex_penalty", "cpp_gwr_nonconvex_lambda_max",
    "cpp_gwr_nonconvex_path_predict", "cpp_gwr_sl_lambda_max",
    "cpp_gwr_sl_path_predict", "cpp_gwr_moran", "gwrs_control", "gwr_fit",
    "gwr_local_inference", "gwr_diagnostic_tests"
  )
  for (name in required) {
    api[[name]] <- get(name, envir = package_namespace)
  }
  native_nonconvex <- api$cpp_gwr_nonconvex_path_predict
  native_convex <- api$cpp_gwr_sl_path_predict
  capture <- new.env(parent = emptyenv())
  capture$nonconvex_calls <- list()
  capture$convex_calls <- list()
  capture$working_metrics <- list()

  api$cpp_gwr_nonconvex_path_predict <- function(...) {
    result <- native_nonconvex(..., guarded_block = TRUE)
    capture$nonconvex_calls[[length(capture$nonconvex_calls) + 1L]] <-
      result[c(
        "path_state", "converged", "iterations", "stationarity",
        "coordinate_gap", "objective", "pair_accepted", "block_checks",
        "block_accepted"
      )]
    result
  }

  api$cpp_gwr_sl_path_predict <- function(...) {
    original <- list(...)
    arguments <- original
    working <- arguments$alpha > 0 && arguments$d == 0 &&
      all(arguments$lambda > 0)
    arguments$guarded_acceleration <- working
    arguments$guarded_block <- FALSE
    arguments$guarded_working_set <- working
    if (working) arguments$keep_coefficients <- TRUE
    result <- do.call(native_convex, arguments)
    fields <- c(
      "converged", "iterations", "kkt", "status",
      "acceleration_accepted", "acceleration_gain"
    )
    capture$convex_calls[[length(capture$convex_calls) + 1L]] <-
      result[intersect(fields, names(result))]
    if (working) {
      policy_value <- result$working_set_policy
      diagnostics <- result$working_set_diagnostics
      certificate_bound <- matrix(
        rep(
          arguments$tolerance *
            (1 + arguments$lambda * arguments$alpha),
          each = nrow(result$converged)
        ),
        nrow(result$converged), length(arguments$lambda)
      )
      stopifnot(
        identical(policy_value$schema, "gwrs-convex-working-set-v3"),
        policy_value$max_size == 128,
        policy_value$max_inner_updates == 8,
        policy_value$max_screen_fit_batches == 16,
        identical(
          dim(result$coefficients),
          c(
            ncol(arguments$x_train) + 1L,
            nrow(arguments$x_target), length(arguments$lambda)
          )
        ),
        all(is.finite(result$coefficients)),
        all(result$kkt[result$converged == 1L] <=
              certificate_bound[result$converged == 1L] * (1 + 1e-10)),
        all(diagnostics$max_selected_size <= policy_value$max_size),
        all(diagnostics$certificate_failures <=
              diagnostics$certificate_checks),
        all(diagnostics$duality_certificate_failures <=
              diagnostics$duality_certificate_checks)
      )
      capture$working_metrics[[length(capture$working_metrics) + 1L]] <-
        environment$usa_study_v1_working_metrics(list(
          args = arguments, raw = result
        ))
      if (!isTRUE(original$keep_coefficients)) {
        result$coefficients <- array(numeric(), c(0L, 0L, 0L))
      }
    }
    result
  }

  context <- list(
    app = app, cfg = cfg, g = api, io = io, threads = 1L,
    post_threads = as.integer(Sys.getenv(
      "GWRS_STUDY_POST_THREADS", unset = "1"
    )),
    native_capture = capture
  )
  stopifnot(
    length(context$post_threads) == 1L,
    !is.na(context$post_threads), context$post_threads >= 1L
  )
  if (mode == "smoke") {
    context$cfg$n <- policy$smoke$n
    context$cfg$p <- policy$smoke$p
    context$cfg$k <- policy$smoke$k
    context$cfg$en_alpha <- policy$smoke$en_alpha
    context$cfg$n_lambda <- policy$smoke$n_lambda
    context$cfg$lambda_ratio_min <- policy$smoke$lambda_ratio_min
    data <- usa_study_v1_smoke_data(policy)
  } else {
    data <- usa_study_v1_read_data(app)
  }
  sys.source(file.path(app, "code/models.R"), environment)
  environment$usa_study_v1_install_models(environment)
  list(e = environment, ctx = context, data = data, policy = policy)
}
