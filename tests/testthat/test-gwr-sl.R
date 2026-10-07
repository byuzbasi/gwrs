test_that("unpenalized core agrees with direct weighted least squares", {
  example <- make_spatial_example(n = 60, p = 3, seed = 10)
  neighbors <- gwr_neighbors(example$coords, k = 35)
  fit <- gwr_fit(
    example$x, example$y, neighbors = neighbors, standardize = FALSE,
    control = serial_control(tolerance = 1e-10)
  )

  for (target in c(1L, 23L, 60L)) {
    index <- neighbors$index[, target]
    weight <- adaptive_weights(neighbors, target, "bisquare")
    direct <- stats::lm.wfit(
      cbind(1, example$x[index, , drop = FALSE]),
      example$y[index],
      w = weight
    )$coefficients
    expect_equal(unname(fit$coefficients[target, ]), unname(direct),
                 tolerance = 1e-8)
  }
})

test_that("GWR-SL coordinate update matches an independent R implementation", {
  example <- make_spatial_example(n = 55, p = 3, seed = 19)
  neighbors <- gwr_neighbors(example$coords, k = 32)
  lambda <- 0.14
  alpha <- 0.65
  d <- 0.4
  fit <- gwr_sl_fit(
    example$x, example$y, neighbors = neighbors, lambda = lambda,
    alpha = alpha, d = d, standardize = FALSE,
    control = serial_control(tolerance = 1e-11)
  )

  target <- 17L
  index <- neighbors$index[, target]
  weight <- adaptive_weights(neighbors, target)
  local_x <- example$x[index, , drop = FALSE]
  local_y <- example$y[index]
  x_mean <- colSums(local_x * weight)
  y_mean <- sum(local_y * weight)
  xc <- sweep(local_x, 2L, x_mean)
  yc <- local_y - y_mean
  gram <- crossprod(xc, xc * weight)
  score <- crossprod(xc, weight * yc)[, 1L]
  center <- solve(gram, score)
  beta <- center
  residual <- as.double(yc - xc %*% beta)
  lambda1 <- lambda * alpha
  lambda2 <- lambda * (1 - alpha)
  for (iteration in seq_len(10000L)) {
    old <- beta
    for (column in seq_along(beta)) {
      partial <- sum(xc[, column] * weight * residual) +
        gram[column, column] * beta[column]
      shifted <- partial + lambda2 * d * center[column]
      updated <- sign(shifted) * max(abs(shifted) - lambda1, 0) /
        (gram[column, column] + lambda2)
      residual <- residual - xc[, column] * (updated - beta[column])
      beta[column] <- updated
    }
    if (max(abs(beta - old)) <= 1e-12 * (1 + max(abs(beta)))) break
  }
  reference <- c(y_mean - sum(x_mean * beta), beta)
  expect_equal(unname(fit$coefficients[target, ]), unname(reference),
               tolerance = 2e-8)
  expect_lt(max(fit$diagnostics$max_kkt_violation), 1e-7)
})

test_that("serial and parallel target evaluation agree", {
  example <- make_spatial_example(n = 75, p = 4, seed = 27)
  neighbors <- gwr_neighbors(example$coords, k = 30)
  serial <- gwr_sl_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.09,
    alpha = 0.7, d = 0.5, control = serial_control()
  )
  parallel <- gwr_sl_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.09,
    alpha = 0.7, d = 0.5,
    control = gwrs_control(tolerance = 1e-8, max_iterations = 3000,
                           n_threads = 2, grain_size = 3)
  )
  expect_equal(parallel$coefficients, serial$coefficients,
               tolerance = 1e-12)
  expect_equal(parallel$fitted.values, serial$fitted.values,
               tolerance = 1e-12)
})

test_that("comparison wrappers are exact GWR-SL special cases", {
  example <- make_spatial_example(n = 45, seed = 31)
  neighbors <- gwr_neighbors(example$coords, k = 25)
  control <- serial_control()
  lasso <- gwr_lasso_fit(example$x, example$y, neighbors = neighbors,
                         lambda = 0.1, control = control)
  direct <- gwr_sl_fit(example$x, example$y, neighbors = neighbors,
                       lambda = 0.1, alpha = 1, d = 0, control = control)
  expect_equal(lasso$coefficients, direct$coefficients)

  ridge <- gwr_ridge_fit(example$x, example$y, neighbors = neighbors,
                         lambda = 0.1, control = control)
  direct_ridge <- gwr_sl_fit(example$x, example$y, neighbors = neighbors,
                             lambda = 0.1, alpha = 0, d = 0,
                             control = control)
  expect_equal(ridge$coefficients, direct_ridge$coefficients)
})
