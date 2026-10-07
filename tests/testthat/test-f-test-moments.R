# Independent, deliberately dense operators are restricted to tiny test fixtures.
dense_conditional_f_operators <- function(fit) {
  n <- nrow(fit$x)
  p <- ncol(fit$x)
  x <- sweep(sweep(fit$x, 2, fit$x_center), 2, fit$x_scale, "/")
  beta <- sweep(fit$coefficients[, -1, drop = FALSE], 2, fit$x_scale, "*")
  coefficient <- lapply(seq_len(p + 1), function(j) matrix(0, n, n))
  smoother <- matrix(0, n, n)
  for (i in seq_len(n)) {
    idx <- fit$neighbors$index[, i]
    w <- adaptive_weights(fit$neighbors, i, fit$kernel)
    local <- x[idx, , drop = FALSE]
    mu <- colSums(local * w)
    xc <- sweep(local, 2, mu)
    gram <- crossprod(xc, xc * w)
    rhs <- t(xc * w)
    a <- if (fit$lambda <= 0 || (is.null(fit$penalty) && fit$alpha == 0)) {
      seq_len(p)
    } else which(abs(beta[i, ]) > 1e-10)
    derivative <- matrix(0, p, length(idx))
    if (length(a)) {
      system <- gram[a, a, drop = FALSE]
      right <- rhs[a, , drop = FALSE]
      if (inherits(fit, "gwr_nonconvex_fit")) {
        t <- abs(beta[i, a])
        curvature <- if (fit$penalty == "mcp") {
          ifelse(t < fit$gamma * fit$lambda, -1 / fit$gamma, 0)
        } else {
          ifelse(t > fit$lambda & t < fit$gamma * fit$lambda,
                 -1 / (fit$gamma - 1), 0)
        }
        diag(system) <- diag(system) + curvature
      } else {
        l2 <- fit$lambda * (1 - fit$alpha)
        diag(system) <- diag(system) + l2
        if (l2 * fit$d != 0) right <- right + l2 * fit$d * solve(gram, rhs)[a, , drop = FALSE]
      }
      derivative[a, ] <- solve(system, right)
    }
    map <- rbind(w - drop(mu %*% derivative),
                 sweep(derivative, 1, fit$x_scale, "/"))
    map[1, ] <- map[1, ] - drop(fit$x_center %*% map[-1, , drop = FALSE])
    for (j in seq_len(p + 1)) coefficient[[j]][i, idx] <- map[j, ]
    smoother[i, idx] <- w + drop((x[i, ] - mu) %*% derivative)
  }
  q <- crossprod(diag(n) - smoother)
  gamma <- vapply(coefficient, function(op) {
    b <- crossprod(sweep(op, 2, colMeans(op))) / n
    c(gamma1 = sum(diag(b)), gamma2 = sum(b^2))
  }, numeric(2))
  list(smoother = smoother, coefficient = coefficient, q = q, gamma = gamma)
}

test_that("exact F moments are invariant to blocks and thread count", {
  e <- make_spatial_example(n = 36, p = 2, seed = 701)
  fit <- gwr_fit(e$x, e$y, coords = e$coords, k = 23,
                 control = gwrs_control(n_threads = 1))
  original <- fit
  one <- gwr_diagnostic_tests(fit, block_size = 1)
  fit$control$n_threads <- 2L
  two <- gwr_diagnostic_tests(fit, block_size = 7)
  expect_equal(one$global_tests, two$global_tests, tolerance = 1e-10)
  expect_equal(one$coefficient_tests, two$coefficient_tests, tolerance = 1e-10)
  expect_equal(one$diagnostics$coefficient_moments,
               two$diagnostics$coefficient_moments, tolerance = 1e-12)
  expect_identical(two$diagnostics$moment_method, "exact_neighbor_block")
  expect_true(two$validity$exact_operator_moments)
  expect_match(two$validity$reference_distribution, "approximation")
  expect_identical(fit$coefficients, original$coefficients)
  expect_identical(fit$fitted.values, original$fitted.values)
  large <- gwr_diagnostic_tests(original, block_size = 1000)
  expect_lt(large$diagnostics$block_size, nrow(e$x))
  expect_equal(one$global_tests, large$global_tests, tolerance = 1e-10)
})

test_that("SCAD and MCP derivative moments use working-scale coefficients", {
  e <- make_spatial_example(n = 38, p = 2, seed = 702)
  e$x[, 1] <- 5 * e$x[, 1] + 70
  e$x[, 2] <- 50 * e$x[, 2] - 10
  for (penalty in c("scad", "mcp")) {
    fit <- gwr_nonconvex_fit(e$x, e$y, coords = e$coords, k = 27,
      lambda = 0.1, penalty = penalty, control = gwrs_control(n_threads = 1))
    test <- gwr_diagnostic_tests(fit, allow_approximate = TRUE, block_size = 5)
    ref <- dense_conditional_f_operators(fit)
    expect_equal(test$diagnostics$L_delta1, sum(diag(ref$q)), tolerance = 2e-9)
    expect_equal(test$diagnostics$L_delta2, sum(ref$q^2), tolerance = 2e-9)
    expect_equal(test$diagnostics$coefficient_moments$gamma1,
                 unname(ref$gamma[1, ]), tolerance = 2e-9)
    expect_equal(test$diagnostics$coefficient_moments$gamma2,
                 unname(ref$gamma[2, ]), tolerance = 2e-9)
    h <- tcrossprod(qr.Q(qr(cbind(1, e$x))))
    contrast <- diag(nrow(e$x)) - h - ref$q
    expected_df <- sum(diag(contrast))^2 / sum(contrast^2)
    if (is.finite(test$global_tests$df1[2])) {
      expect_equal(test$global_tests$df1[2], expected_df, tolerance = 2e-8)
    }
    expect_true(test$validity$exploratory)
    expect_false(test$validity$post_selection_validated)
  }
})

test_that("centered convex shrinkage retains its pilot derivative", {
  e <- make_spatial_example(n = 30, p = 2, seed = 703)
  fit <- gwr_sl_fit(e$x, e$y, coords = e$coords, k = 22,
    lambda = 0.12, alpha = 0.4, d = 0.6,
    control = gwrs_control(n_threads = 1))
  test <- gwr_diagnostic_tests(fit, allow_approximate = TRUE, block_size = 4)
  ref <- dense_conditional_f_operators(fit)
  expect_equal(test$diagnostics$L_delta2, sum(ref$q^2), tolerance = 2e-9)
  expect_equal(test$diagnostics$coefficient_moments$gamma2,
               unname(ref$gamma[2, ]), tolerance = 2e-9)
  expect_gt(test$diagnostics$trQH, 0)
})

test_that("HD and zero-variability cases do not manufacture F p-values", {
  e <- make_spatial_example(n = 18, p = 24, seed = 704)
  fit <- gwr_mcp_fit(e$x, e$y, coords = e$coords, k = 15, lambda = 100,
                     control = gwrs_control(n_threads = 1))
  test <- gwr_diagnostic_tests(fit, allow_approximate = TRUE, block_size = 4)
  expect_equal(test$diagnostics$df_ols, 0)
  expect_true(all(is.na(test$global_tests$p_value)))
  expect_true(all(is.na(test$coefficient_tests$p_value[-1])))
  expect_equal(test$diagnostics$coefficient_moments$gamma2[-1], rep(0, 24))
  small <- make_spatial_example(n = 20, p = 2, seed = 705)
  stationary <- gwr_fit(small$x, small$y, coords = small$coords, k = 20,
                        kernel = "boxcar", control = gwrs_control(n_threads = 1))
  result <- gwr_diagnostic_tests(stationary, block_size = 3)
  expect_true(all(is.na(result$coefficient_tests$p_value)))
})

test_that("F workspace and malformed inputs fail before silent fallbacks", {
  e <- make_spatial_example(n = 25, p = 2, seed = 706)
  fit <- gwr_fit(e$x, e$y, coords = e$coords, k = 18,
                 control = gwrs_control(n_threads = 1))
  for (bad in list(0, -1, 1.5, NA_real_, Inf)) {
    expect_error(gwr_diagnostic_tests(fit, block_size = bad), "block_size")
  }
  for (bad in list(0, -1, NA_real_, Inf)) {
    expect_error(gwr_diagnostic_tests(fit, memory_limit_mb = bad), "memory_limit_mb")
  }
  expect_error(gwr_diagnostic_tests(fit, memory_limit_mb = 1e-6), "memory_limit_mb")
  low <- gwr_diagnostic_tests(fit, block_size = 1)
  budget <- low$diagnostics$workspace_estimate_mb * 1.1
  expect_error(gwr_diagnostic_tests(fit, block_size = 24, memory_limit_mb = budget),
               "workspace estimate")
  expect_s3_class(gwr_diagnostic_tests(fit, block_size = 1, memory_limit_mb = budget),
                 "gwrs_diagnostic_tests")
  broken <- fit
  broken$neighbors$index[1, 1] <- broken$neighbors$index[2, 1]
  expect_error(gwr_diagnostic_tests(broken), "duplicate")
  singular <- fit
  singular$x[, 2] <- singular$x[, 1]
  singular$x_center[2] <- singular$x_center[1]
  singular$x_scale[2] <- singular$x_scale[1]
  expect_error(gwr_diagnostic_tests(singular), "local operator unavailable")
})
