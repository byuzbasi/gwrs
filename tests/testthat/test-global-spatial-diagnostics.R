test_that("global CN and VIF agree with direct correlation algebra", {
  e <- make_spatial_example(n = 48, p = 4, seed = 911)
  e$x[, 2] <- .9 * e$x[, 1] + .2 * e$x[, 2]
  result <- gwr_global_collinearity(e$x)
  cr <- cor(e$x)
  ev <- eigen(cr, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(result$correlation, cr, tolerance = 1e-12)
  expect_equal(result$condition_number$condition_number, sqrt(max(ev) / min(ev)),
               tolerance = 1e-10)
  expect_equal(result$vif$vif, unname(diag(solve(cr))), tolerance = 1e-10)
  transformed <- sweep(sweep(e$x, 2, c(-3, 20, .5, 4), "*"),
                        2, c(100, -40, 25, 10), "+")
  other <- gwr_global_collinearity(transformed)
  expect_equal(other$condition_number, result$condition_number, tolerance = 1e-10)
  expect_equal(other$vif, result$vif, tolerance = 1e-10)
  fit <- gwr_fit(e$x, e$y, coords = e$coords, k = nrow(e$x), kernel = "boxcar",
                 control = gwrs_control(n_threads = 1))
  expect_identical(gwr_global_collinearity(fit), result)
  local <- gwr_local_collinearity(fit)
  expect_equal(local$condition_number$condition_number,
               rep(result$condition_number$condition_number, nrow(e$x)),
               tolerance = 1e-10)
  expect_equal(local$vif[1, ], setNames(result$vif$vif, colnames(e$x)),
               tolerance = 1e-10)
})

test_that("singular and constant predictors are not hidden", {
  e <- make_spatial_example(n = 40, p = 2, seed = 912)
  duplicate <- cbind(e$x, copy = e$x[, 1])
  singular <- gwr_global_collinearity(duplicate)
  expect_identical(singular$rank, 2L)
  expect_true(is.infinite(singular$condition_number$condition_number))
  expect_true(all(is.infinite(singular$vif$vif)))
  constant <- gwr_global_collinearity(cbind(e$x, constant = 1))
  expect_identical(constant$constant_columns, "constant")
  expect_true(is.na(constant$rank))
  expect_true(all(is.na(constant$correlation)))
  one <- gwr_global_collinearity(e$x[, 1, drop = FALSE])
  expect_equal(one$condition_number$condition_number, 1)
  expect_equal(one$vif$vif, 1)
  expect_error(gwr_global_collinearity(matrix(1, 1, 3)), "two rows")
  expect_error(gwr_global_collinearity(cbind(e$x, NA)), "finite")
  expect_error(gwr_global_collinearity(data.frame(a = letters[1:3])), "numeric")
  expect_error(gwr_global_collinearity(e$x, threshold_cn = NA), "thresholds")
  fit <- gwr_fit(e$x, e$y, coords = e$coords, k = 15,
                 control = gwrs_control(n_threads = 1))
  idx <- fit$neighbors$index[, 1]
  fit$x[idx, 1] <- 0 # Constant only inside this diagnostic neighborhood.
  for (scale in c(TRUE, FALSE)) {
    local <- gwr_local_collinearity(fit, scale = scale, return_pairwise = TRUE)
    expect_true(is.infinite(local$condition_number$condition_number[1]))
    expect_true(all(is.infinite(local$vif[1, ])))
    expect_true(is.na(local$max_abs_correlation[1]))
  }
})

test_that("combined Moran and LISA reuse exactly the existing outputs", {
  e <- make_spatial_example(n = 40, p = 2, seed = 913)
  fit <- gwr_fit(e$x, e$y, coords = e$coords, k = 22,
                 control = gwrs_control(n_threads = 1))
  for (standardize in c(TRUE, FALSE)) {
    both <- gwr_moran_diagnostics(fit, permutations = 39, seed = 81,
                                  row_standardize = standardize)
    expect_identical(both$global, gwr_moran(fit, permutations = 39, seed = 81,
                                           row_standardize = standardize))
    expect_identical(both$local, gwr_local_moran(fit, permutations = 39, seed = 81,
                                                row_standardize = standardize))
  }
  both <- gwr_moran_diagnostics(fit, permutations = 39, seed = 81)
  fit$control$n_threads <- 2L
  expect_identical(gwr_moran_diagnostics(fit, permutations = 39, seed = 81), both)
  none <- gwr_moran_diagnostics(fit, permutations = 0)
  expect_true(is.na(none$global$p_value))
  expect_true(all(is.na(none$local$p_value)))
  expect_error(gwr_moran_diagnostics(fit, alpha = 2), "alpha")
  # Independent observed-statistic calculation on a small sparse graph.
  w <- matrix(0, nrow(e$x), nrow(e$x))
  for (i in seq_len(nrow(w))) {
    idx <- fit$neighbors$index[, i]
    ww <- adaptive_weights(fit$neighbors, i)
    ww[idx == i] <- 0
    w[i, idx] <- ww / sum(ww)
  }
  z <- residuals(fit) - mean(residuals(fit))
  expect_equal(both$local$local_I, drop(z * (w %*% z)) / mean(z^2),
               tolerance = 1e-12)
  expect_equal(both$global$I, mean(both$local$local_I), tolerance = 1e-12)
})

test_that("reports and maps include global CN, LISA and F1-F3", {
  e <- make_spatial_example(n = 40, p = 2, seed = 914)
  fit <- gwr_fit(e$x, e$y, coords = e$coords, k = 23,
                 control = gwrs_control(n_threads = 1))
  both <- gwr_moran_diagnostics(fit, permutations = 19, seed = 13)
  local <- gwr_local_collinearity(fit)
  global <- gwr_global_collinearity(fit)
  output <- capture.output(report <- gwrs_report(fit, include_inference = FALSE,
    include_f_tests = TRUE, moran = both$global, local_moran = both$local,
    global_collinearity = global, collinearity = local))
  expect_identical(report$local_moran, both$local)
  expect_identical(report$global_collinearity, global)
  expect_equal(report$f_tests$global_tests$test, c("F1", "F2"))
  expect_equal(nrow(report$f_tests$coefficient_tests), 3)
  expect_true(any(grepl("Residual LISA", output)))
  expect_true(any(grepl("global predictor collinearity", output)))
  stats <- summary(local)
  expect_equal(stats$locations, 40)
  expect_equal(stats$p90, unname(quantile(local$condition_number$condition_number,
                                         .9, type = 1)))
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(local, coords = e$coords))
  local$condition_number$condition_number[1] <- Inf
  expect_invisible(plot(local, coords = e$coords))
  expect_invisible(plot(both$local, coords = e$coords))
  expect_error(plot(both$local, coords = e$coords[-1, ]), "one row")
  expect_invisible(plot(gwr_local_moran(fit, permutations = 0), coords = e$coords))
})
