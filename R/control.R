#' Control numerical fitting in `gwrs`
#'
#' @param tolerance Relative coefficient convergence tolerance.
#' @param max_iterations Maximum number of solver iterations per location.
#' @param n_threads Number of parallel worker threads. Use `1` for serial
#'   execution or `-1` to let `RcppParallel` choose.
#' @param grain_size Minimum number of target locations assigned per parallel
#'   task.
#' @param keep_data Whether fitted objects retain the training matrices.
#' @param diagnostics Statistical diagnostic level. `"standard"` computes
#'   compact fit and effective-degrees-of-freedom diagnostics, `"full"` also
#'   retains quantities needed for local coefficient inference, and `"none"`
#'   skips statistical diagnostics. Logical values remain supported for
#'   backward compatibility (`TRUE` maps to `"standard"`). Numerical solver
#'   diagnostics are always retained because they are small and are required
#'   to audit convergence.
#' @param nonconvex_solver Gaussian SCAD/MCP solver: `"coordinate"` retains
#'   the legacy updates and failure behavior (default); `"guarded_block"`
#'   adds objective-decreasing pair and flat-tail QR proposals. In guarded
#'   mode a failed target stops its remaining lambda path. CV scores require
#'   convergence for every held-out observation. Other estimators ignore this
#'   option. No environment variable is needed.
#'
#' @return A list of class `gwrs_control`.
#'
#' @examples
#' control <- gwrs_control(
#'   tolerance = 1e-6,
#'   max_iterations = 500,
#'   n_threads = 1,
#'   diagnostics = "full"
#' )
#' control
#'
#' @export
gwrs_control <- function(tolerance = 1e-7,
                         max_iterations = 1000L,
                         n_threads = -1L,
                         grain_size = 16L,
                         keep_data = TRUE,
                         diagnostics = TRUE,
                         nonconvex_solver = c("coordinate", "guarded_block")) {
  nonconvex_solver <- match.arg(nonconvex_solver)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("`tolerance` must be one positive finite number.", call. = FALSE)
  }
  max_iterations <- as.integer(max_iterations)
  n_threads <- as.integer(n_threads)
  grain_size <- as.integer(grain_size)
  if (is.na(max_iterations) || max_iterations < 1L) {
    stop("`max_iterations` must be a positive integer.", call. = FALSE)
  }
  if (is.na(n_threads) || n_threads == 0L || n_threads < -1L) {
    stop("`n_threads` must be -1 or a positive integer.", call. = FALSE)
  }
  if (is.na(grain_size) || grain_size < 1L) {
    stop("`grain_size` must be a positive integer.", call. = FALSE)
  }
  diagnostic_level <- normalize_diagnostic_level(diagnostics)
  structure(
    list(
      tolerance = tolerance,
      max_iterations = max_iterations,
      n_threads = n_threads,
      grain_size = grain_size,
      keep_data = isTRUE(keep_data),
      diagnostics = diagnostic_level != "none",
      diagnostic_level = diagnostic_level,
      nonconvex_solver = nonconvex_solver
    ),
    class = "gwrs_control"
  )
}

validate_control <- function(control) {
  if (is.null(control)) return(gwrs_control())
  if (!inherits(control, "gwrs_control")) {
    stop("`control` must be created with `gwrs_control()`.", call. = FALSE)
  }
  if (is.null(control$diagnostic_level)) {
    control$diagnostic_level <- normalize_diagnostic_level(
      control$diagnostics %||% TRUE
    )
  }
  control$nonconvex_solver <- match.arg(
    control$nonconvex_solver %||% "coordinate", c("coordinate", "guarded_block")
  )
  control
}

normalize_diagnostic_level <- function(diagnostics) {
  if (is.logical(diagnostics)) {
    if (length(diagnostics) != 1L || is.na(diagnostics)) {
      stop("`diagnostics` must be one logical value or a diagnostic level.",
           call. = FALSE)
    }
    return(if (diagnostics) "standard" else "none")
  }
  if (!is.character(diagnostics) || length(diagnostics) != 1L ||
      is.na(diagnostics)) {
    stop("`diagnostics` must be one of 'none', 'standard', or 'full'.",
         call. = FALSE)
  }
  match.arg(diagnostics, c("standard", "full", "none"))
}

kernel_code <- function(kernel) {
  choices <- c("gaussian", "bisquare", "exponential", "tricube", "boxcar")
  if (length(kernel) > 1L) kernel <- kernel[[1L]]
  kernel <- match.arg(kernel, choices)
  list(name = kernel, code = match(kernel, choices) - 1L)
}
