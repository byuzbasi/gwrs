selection_example <- function(n = 28L, p = 2L) {
  e <- make_spatial_example(n = n, p = p, seed = 1401)
  d <- data.frame(y = e$y, e$x, east = e$coords[, 1], north = e$coords[, 2])
  names(d)[seq_len(p) + 1L] <- paste0("x", seq_len(p))
  d
}

small_selection <- function(d = selection_example(), ...) {
  gwr_select(y ~ x1 + x2, d, c("east", "north"),
    kernels = "gaussian", k = c(18, nrow(d)), f_tests = FALSE,
    control = gwrs_control(n_threads = 1), verbose = FALSE, ...)
}
