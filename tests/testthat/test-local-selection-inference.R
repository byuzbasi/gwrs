test_that("ordinary GWR inference reproduces independent OLS with global boxcar weights", {
  e <- make_spatial_example(n = 34, p = 2, seed = 1403)
  fit <- gwr_fit(e$x, e$y, coords = e$coords, k = 34, kernel = "boxcar",
    standardize = FALSE,
    control = gwrs_control(diagnostics = "full", n_threads = 1))
  ref <- summary(stats::lm(e$y ~ e$x))$coefficients
  inf <- gwr_local_inference(fit)
  for (j in seq_len(ncol(fit$coefficients))) {
    z <- inf$local[inf$local$coefficient == colnames(fit$coefficients)[j], ]
    expect_equal(z$estimate, rep(ref[j, 1], 34), tolerance = 1e-9)
    expect_equal(z$standard_error, rep(ref[j, 2], 34), tolerance = 1e-9)
    expect_equal(z$statistic, rep(ref[j, 3], 34), tolerance = 1e-9)
    expect_equal(z$p_value, rep(ref[j, 4], 34), tolerance = 1e-9)
    expect_equal(z$p_adjusted, stats::p.adjust(z$p_value, "BH"))
  }
  expect_identical(inf$status, "calculated")
  expect_identical(inf$settings$basis, "unpenalized GWR linear smoother")
  expect_false(inf$settings$reference_calibrated)
  expect_false(inf$settings$selection_conditioned)
  expect_equal(inf$settings$reference_df, 31, tolerance = 1e-9)
})

test_that("selected-model inference agrees with direct fits and does not change ranking", {
  d <- selection_example()
  on <- small_selection(d)
  off <- small_selection(d, local_inference = FALSE)
  clean <- function(z) { z$elapsed_seconds <- NULL; z }
  expect_equal(clean(on$ranking), clean(off$ranking))
  expect_equal(coef(on), coef(off), tolerance = 1e-12)
  expect_identical(off$local_inference_status, "disabled")
  expect_null(off$local_inference)
  expect_identical(gwr_diagnostics(off$best_model)$level, "standard")
  expect_error(gwr_local_inference(off), "local_inference = TRUE")
  expect_identical(on$settings$control$diagnostic_level, "standard")
  expect_identical(gwr_diagnostics(on$best_model)$level, "full")
  expect_identical(on$local_inference_status, "calculated")
  expect_equal(gwr_local_inference(on), on$local_inference)
  expect_equal(summary(on)$local_inference, on$local_inference)
  expect_output(print(summary(on)), "does not establish calibration")
  row <- on$selected
  x <- model.matrix(as.formula(paste("y ~", row$terms)), d)[, -1, drop = FALSE]
  fit <- gwr_fit(x, d$y, coords = as.matrix(d[c("east", "north")]),
    k = row$k, kernel = row$kernel, standardize = FALSE,
    control = gwrs_control(diagnostics = "full", n_threads = 1))
  expect_equal(on$local_inference$local, gwr_local_inference(fit)$local,
    tolerance = 1e-12)
  expect_true(on$local_inference$settings$selection_conditioned)
  expect_equal(on$local_inference$settings$selection$candidate_id, row$candidate_id)
  expect_equal(on$local_inference$settings$selection$sample_rows, rownames(d))
  global <- gwr_local_inference(on, adjust_scope = "global")
  expect_equal(global$local$p_adjusted, stats::p.adjust(global$local$p_value, "BH"))
  for (name in unique(on$local_inference$local$coefficient)) {
    z <- subset(on$local_inference$local, coefficient == name)
    expect_equal(z$p_adjusted, stats::p.adjust(z$p_value, "BH"))
  }
  expect_equal(gwr_local_inference(on, adjust = "none")$local$p_adjusted,
    on$local_inference$local$p_value)
})

test_that("full diagnostics are restricted to the selected GWR final refit", {
  levels <- character()
  actual <- get("selection_fit", envir = asNamespace("gwrs"))
  local_mocked_bindings(selection_fit = function(candidate, prep, settings, neighbors) {
    if (candidate$family == "GWR") levels <<- c(levels, settings$control$diagnostic_level)
    actual(candidate, prep, settings, neighbors)
  }, .package = "gwrs")
  small_selection()
  expect_equal(sum(levels == "full"), 1L)
  expect_gt(sum(levels == "standard"), 1L)
})

test_that("unavailable references and missing tests are explicit", {
  e <- make_spatial_example(n = 24, p = 2, seed = 1404)
  fit <- gwr_fit(e$x, e$y, coords = e$coords, k = 18,
    control = gwrs_control(diagnostics = "full", n_threads = 1))
  bad <- fit
  bad$model_diagnostics$global$edf <- 0
  expect_warning(inf <- gwr_local_inference(bad), "p-values are NA")
  expect_identical(inf$status, "unavailable")
  expect_true(all(is.na(inf$local$p_value)))
  expect_true(all(is.na(inf$local$p_adjusted)))
  expect_true(all(is.na(inf$local$significant)))
  expect_true(all(is.na(inf$summary$min_adjusted_p)))
  expect_true(all(is.na(inf$summary$significant_percent)))
  expect_true(all(inf$summary$available_tests == 0L))
  expect_false(inf$settings$reference_available)
  expect_output(print(inf), "reference: unavailable")
  expect_true(all(is.finite(gwr_local_inference(bad, df = "normal")$local$p_value)))
  bad$model_diagnostics$global$edf <- NA_real_
  expect_warning(gwr_local_inference(bad), "p-values are NA")
  fit$model_diagnostics$inference$statistics[, 2] <- NA_real_
  partial <- gwr_local_inference(fit)
  expect_identical(partial$status, "partially_available")
  expect_true(is.na(partial$summary$min_adjusted_p[2]))
  expect_true(is.na(partial$summary$median_abs_statistic[2]))
  expect_equal(partial$summary$unavailable_tests[2], 24L)
})

test_that("empty and provisional selections retain explicit inference status", {
  empty <- small_selection(max_evaluations = 0)
  expect_null(empty$local_inference)
  expect_identical(empty$local_inference_status, "not_calculated: no_valid_gwr")
  expect_error(gwr_local_inference(empty), "No eligible GWR")
  provisional <- small_selection(max_evaluations = 3)
  expect_false(provisional$search_status$complete)
  expect_false(provisional$local_inference$settings$selection$search_complete)
  expect_true(any(grepl("Provisional", provisional$local_inference$notes)))
  expect_error(small_selection(local_inference = NA), "TRUE or FALSE")
})
