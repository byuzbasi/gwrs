# R only transports configuration and audits native results; no optimization.
nonconvex_guarded <- function(control) {
  identical(control$nonconvex_solver, "guarded_block")
}

nonconvex_valid <- function(raw) {
  raw$path_state == 1L & is.finite(raw$predictions)
}

nonconvex_state_labels <- function(state) {
  labels <- c("unattempted", "converged", "failed")
  matrix(labels[state + 1L], nrow(state), ncol(state))
}

nonconvex_require_fit <- function(raw, control) {
  if (nonconvex_guarded(control) && !all(nonconvex_valid(raw))) {
    condition <- structure(list(
      message = paste("Guarded SCAD/MCP fit did not converge at every target;",
        "no fitted model or prediction was accepted. Inspect condition$solver_result."),
      call = NULL, solver_result = raw
    ), class = c("gwrs_nonconvex_convergence_error", "error", "condition"))
    stop(condition)
  }
  invisible(raw)
}

nonconvex_column_max <- function(x) {
  apply(x, 2L, function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE))
}

nonconvex_solver_summary <- function(raw) {
  data.frame(
    converged_targets = colSums(raw$path_state == 1L),
    failed_targets = colSums(raw$path_state == 2L),
    unattempted_targets = colSums(raw$path_state == 0L),
    valid = colSums(nonconvex_valid(raw)) == nrow(raw$predictions),
    pair_attempts = colSums(raw$pair_attempts),
    pair_accepted = colSums(raw$pair_accepted),
    block_checks = colSums(raw$block_checks),
    block_accepted = colSums(raw$block_accepted)
  )
}
