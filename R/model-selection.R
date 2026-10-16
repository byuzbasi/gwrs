#' Joint subset, kernel and bandwidth selection for unpenalized GWR
#'
#' Compares ordinary Gaussian GWR models on one frozen sample. All five
#' criteria are reported; the selected criterion ranks eligible models.
#' No coefficient penalty is fitted.
#' @param formula Formula defining the response and full candidate scope.
#'   An intercept is required. Factors/multi-column terms are selected as groups.
#' @param data Data frame. Missing/nonfinite values in the full scope fail.
#' @param coords Coordinate matrix or two column names. Distances are Euclidean;
#'   use appropriate projected coordinates. No projection is performed.
#' @param kernels Candidate kernels supported by [gwr_fit()].
#' @param bandwidth "adaptive", "fixed", or both.
#' @param k Adaptive neighbor grid including the target. NULL searches all
#'   integers from p + 2 to n per subset (p excludes the intercept).
#' @param bw Fixed-distance grid. NULL uses 30 log-spaced values from the
#'   smallest positive interpoint distance to twice the largest distance.
#' @param criterion "aicc" (default), "aic", "bic", "gcv", or "loo_cv".
#'   "loo" aliases conditional LOO SSE.
#' @param search "exhaustive" or "stepwise". Stepwise is deterministic
#'   best-improvement addition/deletion after full-model and hierarchy-closed
#'   required-plus-single-term seeds. Every subset searches all requested
#'   kernels/bandwidths. It has no global optimum guarantee.
#' @param required Character term labels or a one-sided formula of required
#'   terms. Their interaction parents are also retained.
#' @param include_ols Compare corresponding OLS fits and, when no terms are
#'   required, intercept-only OLS.
#'   best_model is always GWR; best_overall can be OLS.
#' @param standardize Passed to [gwr_fit()]. FALSE preserves the input scale.
#' @param moran_weights Fixed numeric/sparse weight matrix in data row order,
#'   shared by all candidates. NULL reports Moran as not calculated.
#' @param moran_permutations,moran_seed Exploratory Moran permutation count/seed.
#' @param moran_row_standardize Row-standardize the supplied weights.
#' @param moran_zero_policy "error" rejects islands; "keep" retains all rows.
#' @param f_tests Compute conditional F1/F2/F3 diagnostics for valid GWR fits.
#' @param max_models Whole-search candidate cap. Oversized exhaustive designs
#'   fail before fitting. Stepwise stops explicitly incomplete at the cap.
#' @param max_evaluations New candidate cap for this invocation; resume can
#'   continue a provisional result. Inf means no invocation cap.
#' @param memory_limit_mb Workspace guard, not a hard process memory limit.
#' @param checkpoint New output directory for immutable checkpoint packs and a
#'   SHA-256 manifest. Requires the suggested package jsonlite.
#' @param resume Resume verified, data/settings/software-matching packs.
#' @param checkpoint_every Maximum number of new candidate records per pack.
#' @param dry_run Return design/resource bounds without fitting or testing.
#' @param control [gwrs_control()] settings. Selection always retains training
#'   data and standard diagnostics.
#' @param verbose Print progress every 30 seconds and at phase changes.
#' @param local_inference Compute exploratory local coefficient p-values and
#'   BH-adjusted p-values for the selected GWR only. TRUE by default. Candidate
#'   fits retain standard diagnostics; only the final selected GWR uses full
#'   diagnostics. FALSE skips this additional calculation.
#' @details Adaptive Gaussian/exponential kernels retain only k neighbors,
#'   matching existing gwrs semantics; fixed-distance fits retain all n.
#'   Finite fixed grids do not guarantee a continuous-bandwidth optimum.
#'   Exact criterion ties prefer fewer design columns, OLS for overall ranking,
#'   larger bandwidths, then canonical kernel/subset order.
#'
#'   Eligibility requires full local predictor rank, direct converged solves,
#'   finite coefficients/predictions and positive residual effective df.
#'   Each criterion separately requires a finite value. CN > 30, VIF > 10,
#'   effective sample size at most the local parameter count, leverage > 0.99,
#'   or residual effective df < 5 are descriptive warnings, not exclusion gates
#'   or universal inferential cutoffs. All failures remain in the result.
#'
#'   F3 tests spatial variability, not zero coefficients. F references are
#'   moment-matched approximations; Moran uses total residual randomization.
#'   P-values are raw, exploratory and conditional, never ranking objectives
#'   or post-selection/multiplicity corrections. Conditional LOO fixes fitted
#'   weights/neighborhoods. Predictive claims require external/nested spatial
#'   validation repeating selection inside training folds.
#'
#'   Local coefficient p-values concern zero coefficients, not spatial
#'   constancy (F3). The default BH adjustment is separate for each coefficient
#'   across locations, including the intercept. Neither BH nor the local
#'   reference distribution accounts for response-based subset, kernel or
#'   bandwidth selection. These results are exploratory, not calibrated
#'   post-selection inference, and never affect ranking.
#' @return A gwrs_selection with best_model, best_overall, ranking,
#'   best_by_criterion, global_tests, coefficient_tests, candidate_diagnostics,
#'   local_inference, local_inference_status, settings, provenance and
#'   search_status. local_inference is a [gwr_local_inference()] result, or NULL
#'   if disabled/unavailable. No valid GWR gives best_model=NULL
#'   and status "no_valid_gwr". Incomplete searches are provisional.
#'   dry_run returns a gwrs_selection_plan.
#' @references Gollini et al. (2015), \doi{10.18637/jss.v063.i17}.
#'   Hurvich et al. (1998), \doi{10.1111/1467-9868.00125}.
#'   Leung et al. (2000), \doi{10.1068/a3162}.
#' @examples
#' set.seed(1401)
#' d <- data.frame(east = runif(24), north = runif(24),
#'   x1 = rnorm(24), x2 = rnorm(24))
#' d$y <- 1 + (1 + d$east) * d$x1 + rnorm(24, sd = 0.5)
#' s <- gwr_select(y ~ x1 + x2, d, c("east", "north"),
#'   kernels = c("bisquare", "gaussian"), k = c(16, 24),
#'   f_tests = FALSE, control = gwrs_control(n_threads = 1))
#' s$best_by_criterion
#' summary(s)
#' @export
gwr_select <- function(formula, data, coords,
    kernels = c("bisquare", "gaussian", "exponential", "tricube", "boxcar"),
    bandwidth = "adaptive", k = NULL, bw = NULL, criterion = "aicc",
    search = c("exhaustive", "stepwise"), required = character(),
    include_ols = TRUE, standardize = FALSE,
    moran_weights = NULL, moran_permutations = 999L, moran_seed = 1L,
    moran_row_standardize = TRUE, moran_zero_policy = c("error", "keep"),
    f_tests = TRUE, max_models = 50000L, max_evaluations = Inf,
    memory_limit_mb = 512, checkpoint = NULL, resume = FALSE,
    checkpoint_every = 25L, dry_run = FALSE, control = gwrs_control(),
    verbose = interactive(), local_inference = TRUE) {
  call <- match.call()
  search <- match.arg(search)
  moran_zero_policy <- match.arg(moran_zero_policy)
  criterion <- tolower(criterion)
  if (identical(criterion, "loo")) criterion <- "loo_cv"
  if (length(criterion) != 1L || !criterion %in% selection_criteria())
    stop("Unknown selection criterion.", call. = FALSE)
  for (arg in c("include_ols", "standardize", "moran_row_standardize",
                "f_tests", "resume", "dry_run", "verbose", "local_inference")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop(arg, " must be TRUE or FALSE.", call. = FALSE)
  }
  selection_integer(max_models, "max_models", 1)
  selection_integer(checkpoint_every, "checkpoint_every", 1)
  if (!identical(max_evaluations, Inf)) selection_integer(max_evaluations, "max_evaluations", 0)
  selection_integer(moran_permutations, "moran_permutations", 0)
  selection_integer(moran_seed, "moran_seed", 0)
  if (!is.numeric(memory_limit_mb) || length(memory_limit_mb) != 1L ||
      !is.finite(memory_limit_mb) || memory_limit_mb <= 0)
    stop("memory_limit_mb must be positive and finite.", call. = FALSE)
  if (!is.character(kernels) || !length(kernels) || anyNA(kernels) ||
      any(!kernels %in% c("bisquare", "gaussian", "exponential", "tricube", "boxcar")))
    stop("Invalid candidate kernels.", call. = FALSE)
  kernels <- sort(unique(kernels))
  if (!is.character(bandwidth) || !length(bandwidth) || anyNA(bandwidth) ||
      any(!bandwidth %in% c("adaptive", "fixed")))
    stop("bandwidth must contain adaptive and/or fixed.", call. = FALSE)
  bandwidth <- sort(unique(bandwidth))
  if (!is.null(k) && !"adaptive" %in% bandwidth)
    stop("k is only an adaptive search grid.", call. = FALSE)
  if (!is.null(bw) && !"fixed" %in% bandwidth)
    stop("bw is only a fixed search grid.", call. = FALSE)
  control <- validate_control(control)
  control$keep_data <- TRUE
  control$diagnostics <- TRUE
  control$diagnostic_level <- "standard"
  prep <- selection_prepare(formula, data, coords, required)
  n <- length(prep$full$y)
  if (!is.null(k)) {
    if (!is.numeric(k) || !length(k) || any(!is.finite(k)) ||
        any(k != floor(k)) || any(k < 2 | k > n))
      stop("k must contain integers between 2 and n.", call. = FALSE)
    k <- sort(unique(as.integer(k)))
  }
  if (!is.null(bw) && (!is.numeric(bw) || !length(bw) ||
                       any(!is.finite(bw) | bw <= 0)))
    stop("bw must contain positive finite distances.", call. = FALSE)
  if (!is.null(moran_weights)) {
    for (labels in dimnames(moran_weights)) {
      if (!is.null(labels) && !identical(labels, rownames(prep$full$model)))
        stop("Named Moran weights must match data row order exactly.", call. = FALSE)
    }
    moran_weights <- validate_gwrs_moran_weights(moran_weights, n, moran_zero_policy)
  }
  settings <- list(kernels = kernels, bandwidth = bandwidth, k = k, bw = bw,
    criterion = criterion, search = search, required = prep$labels[prep$required],
    include_ols = include_ols, standardize = standardize,
    moran_permutations = moran_permutations, moran_seed = moran_seed,
    moran_row_standardize = moran_row_standardize,
    moran_zero_policy = moran_zero_policy, f_tests = f_tests,
    local_inference = local_inference,
    memory_limit_mb = memory_limit_mb, control = control)
  subsets <- selection_initial_subsets(prep, search, max_models)
  # Candidate count precedes the potentially quadratic neighbor allocation.
  if ("fixed" %in% bandwidth && is.null(bw)) settings$bw <- rep(1, 30L)
  null_ols <- include_ols && !length(prep$required)
  initial_count <- sum(vapply(subsets, selection_grid_count, numeric(1),
    prep = prep, settings = settings)) + as.integer(null_ols)
  if (search == "exhaustive" && initial_count > max_models)
    stop("Exhaustive design exceeds max_models; restrict grids/scope or explicitly choose stepwise.", call. = FALSE)
  max_k <- if ("fixed" %in% bandwidth || is.null(k)) n else max(k)
  neighbor_mb <- 16 * as.double(n) * max_k / 1024^2
  if (neighbor_mb > memory_limit_mb)
    stop("Cached-neighbor workspace exceeds memory_limit_mb; restrict adaptive k or increase the explicit budget.", call. = FALSE)
  neighbors <- gwr_neighbors(prep$coords, k = max_k)
  if ("fixed" %in% bandwidth) {
    if (is.null(bw)) {
      positive <- neighbors$distance[neighbors$distance > 0]
      if (!length(positive)) stop("Automatic fixed bandwidths require positive distances.", call. = FALSE)
      bw <- exp(seq(log(min(positive)), log(max(positive)) + log(2), length.out = 30L))
    }
    settings$bw <- sort(unique(bw))
  }
  grids <- if (search == "exhaustive")
    lapply(subsets, selection_grid, prep = prep, settings = settings) else NULL
  if (!is.null(grids)) names(grids) <- vapply(subsets, selection_subset_key, "")
  null_grid <- if (null_ols) selection_grid(integer(), prep, settings) else NULL
  expected <- if (search == "exhaustive")
    sum(vapply(grids, nrow, integer(1))) + as.integer(null_ols) else NA_integer_
  if (dry_run) {
    design <- if (search == "exhaustive") do.call(rbind, unname(grids)) else NULL
    if (null_ols && !is.null(design)) design <- rbind(null_grid, design)
    plan <- list(scope = prep$labels, required = settings$required,
      n = n, design_columns = ncol(prep$full$x), search = search,
      expected_candidates = expected, initial_candidates = initial_count,
      maximum_candidates = max_models, neighbor_workspace_mb = neighbor_mb,
      settings = settings, grids = design,
      moran = if (is.null(moran_weights)) "not_calculated: weights_not_supplied" else "fixed_user_weights")
    class(plan) <- "gwrs_selection_plan"
    return(plan)
  }
  provenance <- selection_provenance(prep, settings, moran_weights)
  io <- selection_checkpoint_open(checkpoint, resume, provenance)
  if (!is.null(io)) on.exit(selection_checkpoint_unlock(io), add = TRUE)
  records <- if (is.null(io)) list() else io$records
  fresh <- list()
  evaluated <- 0L
  halted <- FALSE
  reason <- "completed"
  started <- proc.time()[["elapsed"]]
  last_progress <- -Inf
  emit <- function(phase, force = FALSE) {
    elapsed <- proc.time()[["elapsed"]] - started
    if (!force && elapsed - last_progress < 30) return(invisible(NULL))
    last_progress <<- elapsed
    selection_progress(io, records, expected, evaluated, length(fresh), elapsed,
      phase, verbose, max_models)
  }
  commit <- function() {
    if (!is.null(io) && length(fresh)) selection_checkpoint_commit(io, fresh)
    fresh <<- list()
  }
  evaluate <- function(grid) {
    for (i in seq_len(nrow(grid))) {
      candidate <- grid[i, , drop = FALSE]
      id <- candidate$candidate_id
      if (!is.null(records[[id]])) next
      if (evaluated >= max_evaluations || length(records) >= max_models) {
        halted <<- TRUE
        reason <<- if (evaluated >= max_evaluations) "invocation_budget" else "search_budget"
        return(invisible(NULL))
      }
      rec <- selection_evaluate(candidate, prep, settings, neighbors, moran_weights)
      records[[id]] <<- rec
      fresh[[id]] <<- rec
      evaluated <<- evaluated + 1L
      if (length(fresh) >= checkpoint_every) commit()
      emit("search")
    }
    invisible(NULL)
  }
  emit("search", TRUE)
  if (null_ols) evaluate(null_grid)
  visited <- character()
  frontier <- subsets
  repeat {
    for (s in frontier) {
      if (halted) break
      grid <- if (search == "exhaustive") grids[[selection_subset_key(s)]] else
        selection_grid(s, prep, settings)
      evaluate(grid)
      if (!halted) visited <- union(visited, selection_subset_key(s))
    }
    if (halted || search == "exhaustive") break
    tab <- selection_record_table(records)
    ranked <- selection_order(tab, criterion, "GWR")
    centers <- if (length(ranked)) list(selection_decode(tab$subset_id[ranked[1L]])) else frontier
    next_sets <- unlist(lapply(centers, selection_adjacent, prep = prep), recursive = FALSE)
    next_sets <- selection_unique_subsets(next_sets)
    frontier <- Filter(function(s) !selection_subset_key(s) %in% visited, next_sets)
    if (!length(frontier)) break
  }
  commit()
  emit("finalizing", TRUE)
  result <- selection_result(records, prep, settings, neighbors, provenance, call,
    complete = !halted, reason = reason, evaluated = evaluated, expected = expected,
    max_models = max_models)
  if (!is.null(io)) selection_checkpoint_finish(io, result)
  emit(if (halted) "incomplete" else "complete", TRUE)
  result
}

selection_criteria <- function() c("aic", "aicc", "bic", "gcv", "loo_cv")

selection_integer <- function(value, name, minimum) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < minimum || value >= .Machine$integer.max)
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
}

selection_prepare <- function(formula, data, coords, required) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("Supply a two-sided model formula.", call. = FALSE)
  full <- formula_model_components(formula, data, coords)
  if (length(attr(full$terms, "offset")))
    stop("Offsets are unsupported in GWR subset selection.", call. = FALSE)
  if (!is.numeric(full$y) || !is.null(dim(full$y)) || any(!is.finite(full$y)) ||
      length(full$y) < 4L || stats::var(full$y) == 0)
    stop("The response must be a finite, nonconstant numeric vector with at least four observations.", call. = FALSE)
  full$x <- as_design_matrix(full$x)
  coords <- as_coordinate_matrix(full$coords)
  if (nrow(coords) != length(full$y))
    stop("Coordinates and data must have identical row counts/order.", call. = FALSE)
  labels <- attr(full$terms, "term.labels")
  if (!length(labels)) stop("Supply at least one candidate predictor term.", call. = FALSE)
  if (inherits(required, "formula")) {
    if (length(required) != 2L) stop("required must be a one-sided formula.", call. = FALSE)
    required <- attr(stats::terms(required), "term.labels")
  }
  if (!is.character(required) || anyNA(required) || any(!required %in% labels))
    stop("Every required term must occur in the full formula.", call. = FALSE)
  design <- stats::model.matrix(full$terms, full$model)
  assignment <- attr(design, "assign")
  assignment <- assignment[assignment != 0L]
  factors <- attr(full$terms, "factors") > 0
  parents <- lapply(seq_along(labels), function(j) {
    vars <- which(factors[, j])
    if (length(vars) <= 1L) return(integer())
    main <- vapply(vars, function(v) any(colSums(factors) == 1L & factors[v, ]), logical(1))
    if (!all(main)) stop("Interaction scope must explicitly include all main effects.", call. = FALSE)
    which(vapply(seq_along(labels), function(i) {
      other <- which(factors[, i])
      length(other) < length(vars) && all(other %in% vars)
    }, logical(1)))
  })
  prep <- list(full = full, coords = coords, labels = labels, assignment = assignment,
    parents = parents, required = match(unique(required), labels))
  prep$required <- selection_close(prep$required, prep)
  prep
}

selection_close <- function(s, prep) {
  repeat {
    next_s <- sort(unique(c(s, unlist(prep$parents[s], use.names = FALSE))))
    if (identical(as.integer(s), as.integer(next_s))) return(as.integer(next_s))
    s <- next_s
  }
}

selection_subset_key <- function(s) if (length(s)) paste(sort(s), collapse = ".") else "0"
selection_decode <- function(key) if (identical(key, "0")) integer() else as.integer(strsplit(key, ".", fixed = TRUE)[[1L]])
selection_unique_subsets <- function(x) x[!duplicated(vapply(x, selection_subset_key, ""))]

selection_initial_subsets <- function(prep, search, max_models) {
  all <- seq_along(prep$labels)
  req <- prep$required
  if (search == "stepwise") {
    x <- c(list(all), if (length(req)) list(req), lapply(all, function(j) selection_close(union(req, j), prep)))
    return(selection_unique_subsets(x))
  }
  optional <- setdiff(all, req)
  if (length(optional) >= 31L || 2^length(optional) - as.integer(!length(req)) > max_models)
    stop("Subset enumeration exceeds max_models; restrict scope or explicitly choose stepwise.", call. = FALSE)
  out <- if (length(req)) list(req) else list()
  for (m in seq_along(optional))
    out <- c(out, utils::combn(optional, m, function(s) sort(union(s, req)), simplify = FALSE))
  Filter(function(s) identical(s, selection_close(s, prep)), out)
}

selection_adjacent <- function(s, prep) {
  add <- lapply(setdiff(seq_along(prep$labels), s), function(j) sort(c(s, j)))
  drop <- lapply(setdiff(s, prep$required), function(j) setdiff(s, j))
  Filter(function(x) length(x) && identical(x, selection_close(x, prep)), c(add, drop))
}

selection_components <- function(s, prep) {
  out <- prep$full
  out$x <- out$x[, prep$assignment %in% s, drop = FALSE]
  remove <- setdiff(seq_along(prep$labels), s)
  if (length(remove)) out$terms <- stats::drop.terms(out$terms, remove, keep.response = TRUE)
  used <- all.vars(stats::delete.response(out$terms))
  out$contrasts <- out$contrasts[intersect(names(out$contrasts), used)]
  out$xlevels <- out$xlevels[intersect(names(out$xlevels), used)]
  out
}

selection_grid_count <- function(s, prep, settings) {
  p <- sum(prep$assignment %in% s)
  n <- length(prep$full$y)
  adaptive <- if ("adaptive" %in% settings$bandwidth)
    if (is.null(settings$k)) max(0, n - p - 1L) else length(settings$k) else 0
  fixed <- if ("fixed" %in% settings$bandwidth) length(settings$bw) else 0
  as.double(settings$include_ols) + length(settings$kernels) * (adaptive + fixed)
}

selection_grid <- function(s, prep, settings) {
  key <- selection_subset_key(s)
  p <- sum(prep$assignment %in% s)
  n <- length(prep$full$y)
  label <- if (length(s)) paste(prep$labels[s], collapse = " + ") else "1"
  rows <- list()
  push <- function(family, kernel, type, k, bw, boundary, support) {
    id <- paste(key, family, kernel, type, format(k, scientific = FALSE, trim = TRUE),
      sprintf("%.17g", bw), sep = "|")
    rows[[length(rows) + 1L]] <<- data.frame(candidate_id = id, subset_id = key,
      terms = label, n_terms = length(s), n_columns = p, family = family,
      kernel = kernel, bandwidth_type = type, k = k, bw = bw,
      boundary = boundary, support = support, row.names = NULL, stringsAsFactors = FALSE)
  }
  if (settings$include_ols) push("OLS", NA_character_, "global", NA_integer_,
    NA_real_, FALSE, "all_observations")
  if (length(s)) {
    for (kernel in settings$kernels) for (type in settings$bandwidth) {
      if (type == "adaptive") {
        values <- settings$k %||% if (p + 2L <= n) seq.int(p + 2L, n) else integer()
        if (length(values)) push("GWR", kernel, type, values, NA_real_,
          values %in% range(values), if (kernel %in% c("gaussian", "exponential"))
            "knn_truncated" else "compact")
      } else {
        if (length(settings$bw)) push("GWR", kernel, type, n, settings$bw,
          settings$bw %in% range(settings$bw), if (kernel %in% c("gaussian", "exponential"))
            "all_observations" else "compact")
      }
    }
  }
  if (!length(rows)) return(data.frame(candidate_id = character(), subset_id = character(),
    terms = character(), n_terms = integer(), n_columns = integer(),
    family = character(), kernel = character(), bandwidth_type = character(),
    k = integer(), bw = double(), boundary = logical(), support = character()))
  do.call(rbind, rows)
}
