make_spatial_example <- function(n = 60L, p = 3L, seed = 123L) {
  set.seed(seed)
  coords <- cbind(runif(n), runif(n))
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("x", seq_len(p))
  y <- 0.7 + 1.5 * x[, 1L] - 0.8 * x[, 2L] +
    0.4 * coords[, 1L] + rnorm(n, sd = 0.15)
  list(x = x, y = y, coords = coords)
}

adaptive_weights <- function(neighbors, target, kernel = "bisquare") {
  distance <- neighbors$distance[, target]
  bandwidth <- max(distance) * (1 + 1e-10)
  ratio <- distance / bandwidth
  weight <- switch(
    kernel,
    gaussian = exp(-0.5 * ratio^2),
    bisquare = ifelse(ratio < 1, (1 - ratio^2)^2, 0),
    exponential = exp(-ratio),
    tricube = ifelse(ratio < 1, (1 - ratio^3)^3, 0),
    boxcar = as.numeric(ratio <= 1)
  )
  weight / sum(weight)
}

serial_control <- function(tolerance = 1e-8, max_iterations = 3000L) {
  gwrs_control(
    tolerance = tolerance,
    max_iterations = max_iterations,
    n_threads = 1L,
    grain_size = 2L
  )
}

dense_gwr_operators <- function(x, neighbors, kernel = "bisquare") {
  n <- nrow(x)
  design <- cbind(`(Intercept)` = 1, x)
  smoother <- matrix(0, n, n)
  coefficient_sum <- matrix(0, ncol(design), n)
  coefficient_sumsq <- matrix(0, ncol(design), n)
  coefficient_map_ss <- matrix(0, n, ncol(design))
  coefficient_operators <- lapply(seq_len(ncol(design)), function(j) matrix(0, n, n))
  for (target in seq_len(n)) {
    index <- neighbors$index[, target]
    weight <- adaptive_weights(neighbors, target, kernel)
    local_design <- design[index, , drop = FALSE]
    coefficient_map <- solve(
      crossprod(local_design, local_design * weight),
      t(local_design * weight)
    )
    smoother[target, index] <- design[target, ] %*% coefficient_map
    for (j in seq_len(ncol(design))) {
      coefficient_operators[[j]][target, index] <- coefficient_map[j, ]
    }
    coefficient_sum[, index] <- coefficient_sum[, index] + coefficient_map
    coefficient_sumsq[, index] <- coefficient_sumsq[, index] +
      coefficient_map^2
    coefficient_map_ss[target, ] <- rowSums(coefficient_map^2)
  }
  list(
    smoother = smoother,
    coefficient_sum = coefficient_sum,
    coefficient_sumsq = coefficient_sumsq,
    coefficient_map_ss = coefficient_map_ss,
    coefficient_operators = coefficient_operators
  )
}
