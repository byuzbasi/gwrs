test_that("local collinearity equals a direct weighted calculation", {
  example <- make_spatial_example(n = 36, p = 2, seed = 121)
  example$x[, 2] <- 0.7 * example$x[, 1] + 0.4 * example$x[, 2]
  neighbors <- gwr_neighbors(example$coords, k = 22)
  fit <- gwr_fit(
    example$x, example$y, neighbors = neighbors,
    control = gwrs_control(n_threads = 1)
  )
  result <- gwr_local_collinearity(fit, return_pairwise = TRUE)
  target <- 11L
  index <- neighbors$index[, target]
  weight <- adaptive_weights(neighbors, target)
  local_x <- example$x[index, , drop = FALSE]
  means <- colSums(local_x * weight)
  centered <- sweep(local_x, 2, means)
  covariance <- crossprod(centered, centered * weight)
  correlation <- covariance / sqrt(outer(diag(covariance), diag(covariance)))
  r <- correlation[1, 2]

  expect_equal(result$vif[target, ], rep(1 / (1 - r^2), 2),
               tolerance = 2e-9, ignore_attr = TRUE)
  expect_equal(
    result$condition_number$condition_number[target],
    sqrt((1 + abs(r)) / (1 - abs(r))), tolerance = 2e-9
  )
  pair <- subset(result$pairwise_correlations, location == target)
  expect_equal(pair$correlation, r, tolerance = 2e-9)

  one <- gwr_fit(
    example$x[, 1, drop = FALSE], example$y, neighbors = neighbors,
    control = gwrs_control(n_threads = 1)
  )
  one_result <- gwr_local_collinearity(one, return_pairwise = TRUE)
  expect_null(one_result$pairwise_correlations)
  expect_equal(one_result$vif[, 1], rep(1, nrow(example$x)))
})

test_that("sparse Moran statistics equal an explicit weight matrix", {
  example <- make_spatial_example(n = 38, p = 2, seed = 127)
  neighbors <- gwr_neighbors(example$coords, k = 16)
  fit <- gwr_fit(
    example$x, example$y, neighbors = neighbors,
    control = gwrs_control(n_threads = 1)
  )
  result <- gwr_moran(fit, permutations = 0)
  n <- length(fit$residuals)
  weights <- matrix(0, n, n)
  for (target in seq_len(n)) {
    index <- neighbors$index[, target]
    weight <- adaptive_weights(neighbors, target)
    weight[index == target] <- 0
    weight <- weight / sum(weight)
    weights[target, index] <- weight
  }
  z <- fit$residuals - mean(fit$residuals)
  expected <- n / sum(weights) * drop(crossprod(z, weights %*% z)) /
    sum(z^2)
  expect_equal(result$I, expected, tolerance = 2e-12)
  expect_true(is.na(result$p_value))

  serial <- gwr_moran(fit, permutations = 129, seed = 44)
  parallel_fit <- fit
  parallel_fit$control$n_threads <- 2L
  parallel <- gwr_moran(parallel_fit, permutations = 129, seed = 44)
  expect_identical(parallel, serial)
  expect_identical(
    gwr_moran(parallel_fit, permutations = 129, seed = 44),
    parallel
  )
  serial_local <- gwr_local_moran(fit, permutations = 129, seed = 45)
  parallel_local <- gwr_local_moran(
    parallel_fit, permutations = 129, seed = 45
  )
  expect_identical(parallel_local, serial_local)
  expect_identical(
    gwr_local_moran(parallel_fit, permutations = 129, seed = 45),
    parallel_local
  )
  lisa <- gwr_local_moran(fit, permutations = 19, seed = 45)
  expect_s3_class(lisa, "gwrs_local_moran")
  expect_equal(nrow(lisa), n)
})

test_that("GWR diagnostic F components equal dense reference calculations", {
  example <- make_spatial_example(n = 34, p = 2, seed = 131)
  neighbors <- gwr_neighbors(example$coords, k = 21)
  fit <- gwr_fit(
    example$x, example$y, neighbors = neighbors, standardize = FALSE,
    control = gwrs_control(diagnostics = "standard", n_threads = 1)
  )
  result <- gwr_diagnostic_tests(fit)
  reference <- dense_gwr_operators(example$x, neighbors)
  smoother <- reference$smoother
  q_diag <- 1 - 2 * diag(smoother) + colSums(smoother^2)
  l_delta1 <- sum(q_diag)
  q <- crossprod(diag(nrow(smoother)) - smoother)
  l_delta2 <- sum(q^2)
  global <- stats::lm.fit(cbind(1, example$x), example$y)
  rss_ols <- sum(global$residuals^2)
  df_ols <- length(example$y) - global$rank
  rss_local <- sum(fit$residuals^2)
  expected_f1 <- (rss_local / l_delta1) / (rss_ols / df_ols)
  expected_f2 <- ((rss_ols - rss_local) / (df_ols - l_delta1)) /
    (rss_ols / df_ols)

  expect_equal(result$diagnostics$L_delta1, l_delta1, tolerance = 2e-8)
  expect_equal(result$diagnostics$L_delta2, l_delta2, tolerance = 2e-8)
  expect_equal(result$global_tests$statistic[1], expected_f1,
               tolerance = 2e-8)
  expect_equal(result$global_tests$statistic[2], expected_f2,
               tolerance = 2e-8)
  expect_equal(result$global_tests$df1[1], l_delta1^2 / l_delta2, tolerance = 2e-8)
  expect_equal(result$global_tests$p_value[1],
               pf(expected_f1, l_delta1^2 / l_delta2, df_ols), tolerance = 2e-8)
  f2_df <- (df_ols - l_delta1)^2 / (df_ols - 2 * l_delta1 + l_delta2)
  expect_equal(result$global_tests$df1[2], f2_df, tolerance = 2e-8)
  expect_equal(result$global_tests$p_value[2],
               pf(expected_f2, f2_df, df_ols, lower.tail = FALSE), tolerance = 2e-8)

  delta1 <- nrow(example$x) - 2 * sum(diag(smoother)) + sum(smoother^2)
  sigma2 <- rss_local / delta1
  for (column in seq_len(ncol(fit$coefficients))) {
    beta <- fit$coefficients[, column]
    variability <- (sum(beta^2) - sum(beta)^2 / length(beta)) /
      length(beta)
    diagonal_bj <- (
      reference$coefficient_sumsq[column, ] -
        reference$coefficient_sum[column, ]^2 / length(beta)
    ) / length(beta)
    gamma1 <- sum(diagonal_bj)
    operator <- reference$coefficient_operators[[column]]
    centered <- sweep(operator, 2, colMeans(operator))
    bj <- crossprod(centered) / nrow(operator)
    gamma2 <- sum(bj^2)
    expected <- (variability / gamma1) / sigma2
    expect_equal(result$coefficient_tests$statistic[column], expected,
                 tolerance = 2e-8)
    expect_equal(result$coefficient_tests$df1[column], gamma1^2 / gamma2,
                 tolerance = 2e-8)
    expect_equal(result$coefficient_tests$p_value[column],
                 pf(expected, gamma1^2 / gamma2, l_delta1^2 / l_delta2,
                    lower.tail = FALSE), tolerance = 2e-8)
  }
})

test_that("penalized F-tests require an explicit approximation opt-in", {
  example <- make_spatial_example(n = 32, p = 2, seed = 137)
  neighbors <- gwr_neighbors(example$coords, k = 20)
  fit <- gwr_lasso_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.06,
    control = gwrs_control(n_threads = 1)
  )
  expect_error(gwr_diagnostic_tests(fit), "allow_approximate")
  approximate <- gwr_diagnostic_tests(fit, allow_approximate = TRUE)
  expect_false(approximate$validity$classical)
  expect_true(approximate$validity$exploratory)
})
