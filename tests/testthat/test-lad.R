test_that("LAD ADMM reaches the independently evaluated convex objective", {
  n <- 24L
  x <- matrix(seq(-1, 1, length.out = n), ncol = 1L)
  coords <- cbind(seq_len(n), rep(0, n))
  y <- 0.4 + 1.8 * x[, 1L] + sin(seq_len(n)) * 0.04
  y[5L] <- y[5L] + 3
  neighbors <- gwr_neighbors(coords, k = n)
  lambda <- 0.12
  alpha <- 0.65
  d <- 0.35
  fit <- lad_gwr_sl_fit(
    x, y, neighbors = neighbors, kernel = "gaussian", bandwidth = 1e6,
    lambda = lambda, alpha = alpha, d = d, standardize = FALSE,
    control = serial_control(tolerance = 1e-7, max_iterations = 5000)
  )
  target <- 1L
  index <- neighbors$index[, target]
  distance <- neighbors$distance[, target]
  weight <- exp(-0.5 * (distance / 1e6)^2)
  weight <- weight / sum(weight)
  center <- fit$gwr_center[target, 1L]
  objective <- function(theta) {
    residual <- y[index] - theta[1L] - x[index, 1L] * theta[2L]
    sum(weight * abs(residual)) + lambda * alpha * abs(theta[2L]) +
      0.5 * lambda * (1 - alpha) * (theta[2L] - d * center)^2
  }
  reference <- stats::optim(
    c(0, 0), objective, method = "Nelder-Mead",
    control = list(maxit = 50000, reltol = 1e-12)
  )
  expect_true(fit$diagnostics$converged[target])
  expect_equal(fit$diagnostics$objective[target], reference$value,
               tolerance = 2e-5)
  expect_equal(unname(objective(fit$coefficients[target, ])), reference$value,
               tolerance = 2e-5)
})

test_that("LAD loss is less distorted by a gross response outlier", {
  n <- 35L
  x <- matrix(seq(-2, 2, length.out = n), ncol = 1L)
  coords <- cbind(seq_len(n), rep(0, n))
  y <- 1 + 2 * x[, 1L]
  y[18L] <- y[18L] + 100
  neighbors <- gwr_neighbors(coords, k = n)
  control <- serial_control(tolerance = 1e-6, max_iterations = 5000)
  least_squares <- gwr_sl_fit(
    x, y, neighbors = neighbors, kernel = "gaussian", bandwidth = 1e6,
    lambda = 1e-4, alpha = 1, d = 0, standardize = FALSE,
    control = control
  )
  lad <- lad_gwr_sl_fit(
    x, y, neighbors = neighbors, kernel = "gaussian", bandwidth = 1e6,
    lambda = 1e-4, alpha = 1, d = 0, standardize = FALSE,
    control = control
  )
  expect_lt(abs(lad$coefficients[1L, 1L] - 1),
            abs(least_squares$coefficients[1L, 1L] - 1))
  expect_lt(abs(lad$coefficients[1L, 2L] - 2), 1e-3)
})
