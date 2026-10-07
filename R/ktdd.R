as_sparse_numeric_matrix <- function(x, name) {
  if (!inherits(x, "sparseMatrix")) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }
  x <- methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix") * 1
  if (any(!is.finite(x@x))) {
    stop(sprintf("`%s` must contain only finite values.", name),
         call. = FALSE)
  }
  x
}

#' Build a symmetric sparse graph Laplacian
#'
#' Directed kernel weights from a reusable neighbor structure are averaged
#' with their transpose. If `A` is the resulting adjacency matrix, the return
#' value is `diag(rowSums(A)) - A`. Self edges are removed. The adjacency uses
#' raw kernel weights rather than the target-wise normalized weights in the
#' local regression loss, so the numerical meaning of a downstream smoothness
#' parameter depends on the kernel, bandwidth, and graph density. This helper
#' creates a coefficient-smoothing graph; it does not infer a physical PDE,
#' mesh, source, or boundary condition.
#'
#' @param neighbors A same-location object from [gwr_neighbors()].
#' @param kernel Spatial kernel.
#' @param bandwidth Fixed bandwidth, or `NULL` for the target-specific adaptive
#'   bandwidth stored through neighbor distances.
#'
#' @return A symmetric sparse `dgCMatrix` graph Laplacian.
#'
#' @examples
#' coords <- cbind(east = seq_len(20), north = rep(0, 20))
#' neighbors <- gwr_neighbors(coords, k = 10)
#' laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
#' laplacian
#' max(abs(Matrix::rowSums(laplacian)))
#'
#' @export
gwr_graph_laplacian <- function(neighbors,
                                kernel = c("bisquare", "gaussian",
                                           "exponential", "tricube",
                                           "boxcar"),
                                bandwidth = NULL) {
  if (!inherits(neighbors, "gwrs_neighbors")) {
    stop("`neighbors` must be created with `gwr_neighbors()`.", call. = FALSE)
  }
  n_train <- neighbors$n_train %||% neighbors$n
  n_target <- neighbors$n_target %||% neighbors$n
  if (!identical(as.integer(n_train), as.integer(n_target))) {
    stop("A graph Laplacian requires same-location neighbors.", call. = FALSE)
  }
  kernel <- kernel_code(kernel)
  bandwidth <- resolve_bandwidth(bandwidth)
  cpp_graph_laplacian(
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed
  )
}

validate_coefficient_laplacian <- function(laplacian, n) {
  laplacian <- as_sparse_numeric_matrix(laplacian,
                                        "coefficient_laplacian")
  if (!identical(dim(laplacian), c(n, n))) {
    stop("`coefficient_laplacian` must be n by n.", call. = FALSE)
  }
  if (!isTRUE(Matrix::isSymmetric(laplacian, tol = 1e-9))) {
    stop("`coefficient_laplacian` must be symmetric.", call. = FALSE)
  }
  row_error <- max(abs(Matrix::rowSums(laplacian)))
  scale <- max(1, if (length(laplacian@x)) max(abs(laplacian@x)) else 0)
  if (!is.finite(row_error) || row_error > 1e-8 * scale) {
    stop("`coefficient_laplacian` must have zero row sums.", call. = FALSE)
  }
  if (any(Matrix::diag(laplacian) < -1e-10 * scale)) {
    stop("`coefficient_laplacian` has a negative diagonal.", call. = FALSE)
  }
  positive <- which(laplacian@x > 1e-10 * scale)
  if (length(positive)) {
    column <- findInterval(positive - 1L, laplacian@p) - 1L
    if (any(laplacian@i[positive] != column)) {
      stop(
        "`coefficient_laplacian` must have nonpositive off-diagonal entries.",
        call. = FALSE
      )
    }
  }
  laplacian
}

estimate_ktdd_workspace <- function(n, p, k, solver) {
  p1 <- p + 1
  mib <- 1024^2
  factor_cache <- 8 * n * p1^2 / mib
  common <- 8 * n * (2 * k + 4 * p1 + 16) / mib
  solver_specific <- if (identical(solver, "factorized")) {
    factor_cache
  } else {
    8 * n * 7 * p1 / mib
  }
  list(
    factor_cache_mib = factor_cache,
    approximate_working_mib = common + solver_specific,
    armadillo_index_bits = cpp_arma_uword_bits()
  )
}

resolve_ktdd_solver <- function(solver,
                                factor_cache_limit_mib,
                                n,
                                p,
                                k) {
  solver <- match.arg(solver, c("auto", "factorized", "matrix_free"))
  if (!is.numeric(factor_cache_limit_mib) ||
      length(factor_cache_limit_mib) != 1L ||
      is.na(factor_cache_limit_mib) || factor_cache_limit_mib <= 0) {
    stop("`factor_cache_limit_mib` must be positive or Inf.", call. = FALSE)
  }
  factor_cache_mib <- 8 * n * (p + 1)^2 / 1024^2
  resolved <- if (identical(solver, "auto")) {
    if (factor_cache_mib <= factor_cache_limit_mib) {
      "factorized"
    } else {
      "matrix_free"
    }
  } else {
    solver
  }
  if (identical(solver, "factorized") &&
      factor_cache_mib > factor_cache_limit_mib) {
    stop(
      sprintf(
        paste0(
          "The factorized KTDD cache is approximately %.1f MiB, above ",
          "`factor_cache_limit_mib = %.1f`. Use `solver = \"matrix_free\"`, ",
          "`solver = \"auto\"`, or explicitly raise the limit."
        ),
        factor_cache_mib, factor_cache_limit_mib
      ),
      call. = FALSE
    )
  }
  estimates <- estimate_ktdd_workspace(n, p, k, resolved)
  c(
    list(requested = solver, resolved = resolved,
         factor_cache_limit_mib = factor_cache_limit_mib),
    estimates
  )
}

#' Fit PDE-constrained geographically weighted regression (GWR-KTDD)
#'
#' GWR-KTDD jointly estimates geographically varying coefficients and a latent
#' spatial process informed by a user-supplied partial differential equation
#' (PDE). KTDD is the Turkish abbreviation for a partial differential equation;
#' the method can therefore also be described as PDE-GWR. Its dimensionally
#' consistent objective is
#' \deqn{\frac12\sum_i\sum_{j\in N_i}\bar w_{ij}
#' (y_j-z_j^T\beta_i-U_j)^2
#' +\frac{\lambda_D}{2}\|LU-S\|_2^2
#' +\frac{\gamma}{2}\operatorname{tr}(B^T L_W B),}
#' where \eqn{\bar w_{ij}=w_{ij}/\sum_{j\in N_i}w_{ij}}, `z_j` includes the
#' intercept, the `n` rows of `B` are the location-specific coefficient
#' vectors, and `U` has one process value per observation. The supplied
#' operator `L` may be rectangular but must have `n` columns, its source `S`
#' must have one value per operator row, and `L_W` is an `n`-by-`n` graph
#' Laplacian applied separately to every coefficient surface, including the
#' intercept. This avoids applying a location graph to a single coefficient
#' vector of incompatible dimension.
#'
#' The first term measures normalized local data fit, the second penalizes
#' disagreement with the discretized equation `L %*% U = S`, and the third
#' smooths coefficient surfaces over `L_W`. Thus `lambda_diffusion` controls
#' physical consistency and `gamma` controls coefficient smoothness; their
#' numerical scales depend on the units and scaling of the supplied operators.
#' When `standardize = TRUE`, estimation and the smoothness penalty operate on
#' the globally standardized predictor scale, although returned coefficients
#' are transformed back to the original scale.
#'
#' The convex quadratic objective is solved by alternating coefficient and
#' process updates. `solver = "factorized"` uses cached local factorizations
#' inside parallel block-Jacobi coefficient updates. It is usually fastest when
#' its `8 * n * (p + 1)^2`-byte cache is moderate. `solver = "matrix_free"`
#' instead solves the joint coefficient normal equations with diagonally
#' preconditioned conjugate gradients (PCG), applying local Gram blocks and the
#' coefficient graph directly. It does not allocate the factor cube and solves
#' the same objective. `solver = "auto"` uses the projected factor-cache size
#' to select between these numerical routes.
#'
#' The process right-hand side is accumulated in parallel by observation, and
#' its diagonally preconditioned PCG updates use parallel row- and
#' column-oriented sparse matrix-vector products when multiple threads are
#' requested, without forming `crossprod(L)`. No solver forms a dense
#' `n`-by-`n` weight, smoother, or diffusion cross-product matrix. The returned
#' `workspace` values are approximate allocation guides, not measurements of
#' peak resident memory. Matrix-free working storage scales primarily with
#' neighbor edges, coefficient vectors, and sparse operators rather than
#' `n * (p + 1)^2`. On 64-bit platforms the package compiles Armadillo with
#' 64-bit indices, allowing sparse `n`-by-`n` operators whose nominal dimension
#' product exceeds the 32-bit element-count limit. This widens C++ sparse and
#' neighbor-index storage; the R `dgCMatrix` nonzero count must still fit its
#' integer slot representation.
#' If `L` annihilates the constant vector, the decomposition into local
#' intercepts and `U` has an arbitrary constant shift. In that case the solver
#' selects the equivalent mean-zero process representation automatically;
#' fitted values and every objective component are unchanged.
#'
#' A no-flux operator commonly also requires a compatible source (for example,
#' a zero-sum source for a connected Laplacian). An incompatible source is
#' still handled as a penalized least-squares target, but it leaves irreducible
#' diffusion loss. The package does not infer a PDE, boundary condition, source,
#' or universal tuning parameters. Users should inspect both outer and PCG
#' convergence diagnostics. KTDD currently supports fitted locations only and
#' exposes descriptive, loss-compatible diagnostics rather than Gaussian AICc,
#' local coefficient tests, or GWR F-tests.
#'
#' @inheritParams gwr_sl_fit
#' @param diffusion_operator Sparse or dense discretized PDE operator `L` with
#'   `n` columns. Its row count may differ from `n`.
#' @param source Numeric source vector `S`, with one value per row of
#'   `diffusion_operator`. Its scientific definition and boundary compatibility
#'   are the caller's responsibility.
#' @param coefficient_laplacian Optional `n`-by-`n` symmetric graph Laplacian
#'   `L_W`, with zero row sums, nonnegative diagonal, and nonpositive
#'   off-diagonal entries. The default is constructed from the raw kernel
#'   adjacency on `neighbors` using [gwr_graph_laplacian()].
#' @param lambda_diffusion Positive weight on PDE disagreement. Its scale is
#'   specific to the units and mesh scaling of `diffusion_operator`.
#' @param gamma Nonnegative weight smoothing all coefficient surfaces over
#'   `coefficient_laplacian`.
#' @param solver Coefficient-system solver. `"factorized"` caches one local
#'   factorization per location and is fastest when that cache is modest.
#'   `"matrix_free"` uses a diagonally preconditioned joint PCG update without
#'   retaining the `n * (p + 1)^2` factor cube. `"auto"` selects between them
#'   from the projected factor-cache size.
#' @param factor_cache_limit_mib Maximum projected factor cache, in MiB, used
#'   by `solver = "auto"`. An explicitly requested factorized solver also
#'   refuses a larger cache unless this limit is deliberately raised or set to
#'   `Inf`.
#' @param coefficient_tolerance Relative PCG tolerance for a matrix-free
#'   coefficient update.
#' @param coefficient_max_iterations Maximum PCG iterations for each
#'   matrix-free coefficient update.
#' @param linear_tolerance Relative PCG tolerance for the process update.
#' @param linear_max_iterations Maximum PCG iterations per outer update.
#'
#' @return An object of class `gwr_ktdd_fit` and `gwrs_fit`. The `solver`,
#'   `workspace`, and `diagnostics` components record the resolved numerical
#'   route, projected memory, coefficient PCG, process PCG, and outer
#'   convergence information.
#'
#' @examples
#' n <- 20
#' coords <- cbind(east = seq_len(n), north = rep(0, n))
#' x <- cbind(trend = seq(-1, 1, length.out = n), wave = sin(seq_len(n)))
#' neighbors <- gwr_neighbors(coords, k = 12)
#' laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
#' process <- sin(seq_len(n) / 3)
#' process <- process - mean(process)
#' source <- as.double(laplacian %*% process)
#' y <- 0.5 + x[, "trend"] - 0.25 * x[, "wave"] + process
#' fit <- gwr_ktdd_fit(
#'   x, y, neighbors = neighbors, kernel = "gaussian",
#'   diffusion_operator = laplacian, source = source,
#'   coefficient_laplacian = laplacian,
#'   lambda_diffusion = 1, gamma = 0.05,
#'   linear_tolerance = 1e-7, linear_max_iterations = 1000,
#'   control = gwrs_control(
#'     tolerance = 1e-6, max_iterations = 500, n_threads = 1
#'   )
#' )
#' fit
#' fit$loss_components
#' fit$workspace
#'
#' @export
gwr_ktdd_fit <- function(x,
                         y,
                         coords = NULL,
                         neighbors = NULL,
                         k = NULL,
                         kernel = c("bisquare", "gaussian", "exponential",
                                    "tricube", "boxcar"),
                         bandwidth = NULL,
                         diffusion_operator,
                         source,
                         coefficient_laplacian = NULL,
                         lambda_diffusion = 1,
                         gamma = 0.1,
                         solver = c("auto", "factorized", "matrix_free"),
                         factor_cache_limit_mib = 512,
                         coefficient_tolerance = 1e-8,
                         coefficient_max_iterations = 2000L,
                         standardize = TRUE,
                         linear_tolerance = 1e-8,
                         linear_max_iterations = 2000L,
                         control = gwrs_control(max_iterations = 500L)) {
  call <- match.call()
  x <- as_design_matrix(x)
  y <- as_response(y, nrow(x))
  control <- validate_control(control)
  neighbors <- resolve_neighbors(x, coords, neighbors, k)
  kernel <- kernel_code(kernel)
  bandwidth <- resolve_bandwidth(bandwidth)
  n <- nrow(x)

  if (missing(diffusion_operator)) {
    stop("Supply an explicit discretized `diffusion_operator`.",
         call. = FALSE)
  }
  diffusion_operator <- as_sparse_numeric_matrix(diffusion_operator,
                                                  "diffusion_operator")
  if (ncol(diffusion_operator) != n) {
    stop("`diffusion_operator` must have one column per observation.",
         call. = FALSE)
  }
  if (missing(source)) {
    stop("Supply the physical source vector `source` explicitly.",
         call. = FALSE)
  }
  source <- as.double(source)
  if (length(source) != nrow(diffusion_operator) || any(!is.finite(source))) {
    stop("`source` must be finite and match the operator row count.",
         call. = FALSE)
  }
  if (is.null(coefficient_laplacian)) {
    coefficient_laplacian <- gwr_graph_laplacian(
      neighbors, kernel = kernel$name,
      bandwidth = if (bandwidth$adaptive) NULL else bandwidth$fixed
    )
  }
  coefficient_laplacian <- validate_coefficient_laplacian(
    coefficient_laplacian, n
  )
  if (length(lambda_diffusion) != 1L || !is.finite(lambda_diffusion) ||
      lambda_diffusion <= 0) {
    stop("`lambda_diffusion` must be one positive finite number.",
         call. = FALSE)
  }
  if (length(gamma) != 1L || !is.finite(gamma) || gamma < 0) {
    stop("`gamma` must be one nonnegative finite number.", call. = FALSE)
  }
  coefficient_max_iterations <- as.integer(coefficient_max_iterations)
  if (length(coefficient_tolerance) != 1L ||
      !is.finite(coefficient_tolerance) || coefficient_tolerance <= 0 ||
      length(coefficient_max_iterations) != 1L ||
      is.na(coefficient_max_iterations) || coefficient_max_iterations < 1L) {
    stop("Coefficient-solver controls must be positive scalars.",
         call. = FALSE)
  }
  linear_max_iterations <- as.integer(linear_max_iterations)
  if (length(linear_tolerance) != 1L || !is.finite(linear_tolerance) ||
      linear_tolerance <= 0 || length(linear_max_iterations) != 1L ||
      is.na(linear_max_iterations) || linear_max_iterations < 1L) {
    stop("Linear-solver controls must be positive scalars.", call. = FALSE)
  }

  prepared <- standardize_design(x, standardize)
  solver_info <- resolve_ktdd_solver(
    solver = solver,
    factor_cache_limit_mib = factor_cache_limit_mib,
    n = n,
    p = ncol(x),
    k = nrow(neighbors$index)
  )
  constant_image <- as.double(diffusion_operator %*% rep.int(1, n))
  operator_scale <- max(1, sqrt(sum(diffusion_operator@x^2)))
  center_process <- sqrt(sum(constant_image^2)) <= 1e-10 * operator_scale
  raw <- cpp_gwr_ktdd_fit(
    prepared$x,
    y,
    neighbors$index,
    neighbors$distance,
    kernel$code,
    bandwidth$adaptive,
    bandwidth$fixed,
    diffusion_operator,
    source,
    coefficient_laplacian,
    as.double(lambda_diffusion),
    as.double(gamma),
    match(solver_info$resolved, c("factorized", "matrix_free")) - 1L,
    as.double(coefficient_tolerance),
    coefficient_max_iterations,
    control$tolerance,
    control$max_iterations,
    as.double(linear_tolerance),
    linear_max_iterations,
    center_process,
    control$n_threads,
    control$grain_size
  )

  coefficients <- backtransform_coefficients(
    raw$coefficients, prepared$center, prepared$scale
  )
  predictor_names <- colnames(x) %||% paste0("x", seq_len(ncol(x)))
  colnames(coefficients) <- c("(Intercept)", predictor_names)
  diagnostics <- list(
    converged = isTRUE(raw$converged),
    iterations = as.integer(raw$iterations),
    coefficient_change = as.double(raw$coefficient_change),
    process_change = as.double(raw$process_change),
    linear_converged = isTRUE(raw$linear_converged),
    linear_iterations = as.integer(raw$linear_iterations),
    linear_residual = as.double(raw$linear_residual),
    coefficient_linear_converged =
      isTRUE(raw$coefficient_linear_converged),
    coefficient_linear_iterations =
      as.integer(raw$coefficient_linear_iterations),
    coefficient_linear_residual =
      as.double(raw$coefficient_linear_residual),
    coefficient_solver = solver_info$resolved,
    coefficient_solver_requested = solver_info$requested,
    factor_cache_mib = solver_info$factor_cache_mib,
    approximate_working_mib = solver_info$approximate_working_mib,
    armadillo_index_bits = solver_info$armadillo_index_bits,
    process_centered = isTRUE(raw$process_centered),
    local_solver = status_labels(raw$status),
    cached_factor = if (identical(solver_info$resolved, "factorized")) {
      ifelse(raw$factor_type == 0L, "cholesky", "pseudoinverse")
    } else {
      rep.int("not_cached", n)
    },
    weight_sum = as.double(raw$weight_sum),
    bandwidth = as.double(raw$bandwidth),
    observation_weight = as.double(raw$observation_weight)
  )
  model_diagnostics <- compute_descriptive_model_diagnostics(
    y = y,
    fitted = as.double(raw$fitted),
    residuals = as.double(raw$residuals),
    neighbors = neighbors,
    kernel = kernel,
    bandwidth = bandwidth,
    control = control,
    method = "GWR-KTDD",
    extra = list(
      joint_objective = as.double(raw$objective),
      local_loss = as.double(raw$loss_components[[1L]]),
      diffusion_loss = as.double(raw$loss_components[[2L]]),
      coefficient_smoothness_loss = as.double(raw$loss_components[[3L]])
    )
  )

  object <- list(
    call = call,
    method = "GWR-KTDD",
    coefficients = coefficients,
    process = as.double(raw$process),
    fitted.values = as.double(raw$fitted),
    residuals = as.double(raw$residuals),
    objective = as.double(raw$objective),
    loss_components = raw$loss_components,
    lambda_diffusion = as.double(lambda_diffusion),
    gamma = as.double(gamma),
    solver = solver_info$resolved,
    solver_requested = solver_info$requested,
    factor_cache_limit_mib = solver_info$factor_cache_limit_mib,
    workspace = solver_info,
    kernel = kernel$name,
    adaptive = bandwidth$adaptive,
    bandwidth = if (bandwidth$adaptive) NULL else bandwidth$fixed,
    neighbors = neighbors,
    diffusion_operator = if (control$keep_data) diffusion_operator else NULL,
    source = if (control$keep_data) source else NULL,
    coefficient_laplacian = if (control$keep_data) {
      coefficient_laplacian
    } else NULL,
    standardize = isTRUE(standardize),
    x_center = prepared$center,
    x_scale = prepared$scale,
    diagnostics = diagnostics,
    solver_diagnostics = diagnostics,
    model_diagnostics = model_diagnostics,
    local_r2 = model_diagnostics$local$local_r_squared %||% NULL,
    terms = NULL,
    predictor_names = predictor_names,
    x = if (control$keep_data) x else NULL,
    y = if (control$keep_data) y else NULL,
    control = control
  )
  class(object) <- c("gwr_ktdd_fit", "gwrs_fit")
  object
}

#' Formula interface for PDE-constrained GWR (GWR-KTDD)
#'
#' This is the formula interface to [gwr_ktdd_fit()]. The discretized PDE
#' operator and source remain explicit required inputs; no physical model or
#' boundary condition is inferred from the coordinates.
#'
#' @param formula Model formula.
#' @param data Data frame containing model variables.
#' @param coords Coordinate matrix or two column names in `data`.
#' @param ... Additional arguments passed to [gwr_ktdd_fit()].
#'
#' @return An object of class `gwr_ktdd_fit` and `gwrs_fit`.
#'
#' @examples
#' n <- 20
#' dat <- data.frame(
#'   east = seq_len(n),
#'   north = rep(0, n),
#'   trend = seq(-1, 1, length.out = n),
#'   wave = sin(seq_len(n))
#' )
#' neighbors <- gwr_neighbors(as.matrix(dat[, c("east", "north")]), k = 12)
#' laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
#' process <- sin(seq_len(n) / 3)
#' process <- process - mean(process)
#' source <- as.double(laplacian %*% process)
#' dat$y <- 0.5 + dat$trend - 0.25 * dat$wave + process
#' fit <- gwr_ktdd(
#'   y ~ trend + wave, dat, coords = c("east", "north"),
#'   k = 12, kernel = "gaussian",
#'   diffusion_operator = laplacian, source = source,
#'   coefficient_laplacian = laplacian,
#'   lambda_diffusion = 1, gamma = 0.05,
#'   control = gwrs_control(
#'     tolerance = 1e-6, max_iterations = 500, n_threads = 1
#'   )
#' )
#' fit
#'
#' @export
gwr_ktdd <- function(formula, data, coords, ...) {
  call <- match.call()
  components <- formula_model_components(formula, data, coords)
  fit <- gwr_ktdd_fit(
    components$x, components$y, coords = components$coords, ...
  )
  decorate_formula_fit(fit, components, call)
}
