test_that("spatial values resolve coefficient aliases and diagnostics", {
  built <- make_heatmap_fit()
  fit <- built$fit

  named <- gwrs:::extract_gwrs_spatial_values(fit, "income")
  indexed <- gwrs:::extract_gwrs_spatial_values(fit, "beta_1")
  intercept <- gwrs:::extract_gwrs_spatial_values(fit, "beta_0")
  local_r2 <- gwrs:::extract_gwrs_spatial_values(fit, "local_r2")
  residuals <- gwrs:::extract_gwrs_spatial_values(fit, "residuals")
  fitted <- gwrs:::extract_gwrs_spatial_values(fit, "fitted")
  response <- gwrs:::extract_gwrs_spatial_values(fit, "response")

  expect_equal(named$points$value, fit$coefficients[, "income"])
  expect_equal(indexed$points$value, named$points$value)
  expect_equal(intercept$points$value, fit$coefficients[, "(Intercept)"])
  expect_equal(local_r2$points$value, fit$local_r2)
  expect_equal(residuals$points$value, fit$residuals)
  expect_equal(fitted$points$value, fit$fitted.values)
  expect_equal(response$points$value, fit$fitted.values + fit$residuals)
  expect_true(all(intercept$points$selected))
  expect_true(any(!named$points$selected))
})

test_that("spatial value errors list valid choices", {
  fit <- make_heatmap_fit()$fit
  expect_error(
    gwrs:::extract_gwrs_spatial_values(fit, "not_a_variable"),
    "Valid choices"
  )
  expect_error(
    gwrs:::extract_gwrs_spatial_values(list(), "income"),
    "inherit"
  )
})

test_that("coordinate validation rejects missing and insufficient locations", {
  fit <- make_heatmap_fit()$fit
  missing <- fit
  missing$neighbors$query_coords[1L, 1L] <- NA_real_
  expect_error(
    gwrs:::extract_gwrs_spatial_values(missing, "income"),
    "finite numeric matrix"
  )

  two <- fit
  two$neighbors$query_coords <- matrix(
    rep(c(0, 0, 1, 1), length.out = 2L * nrow(two$coefficients)),
    ncol = 2L,
    byrow = TRUE
  )
  expect_error(
    gwrs:::extract_gwrs_spatial_values(two, "income"),
    "three unique"
  )
})

test_that("duplicate locations are aggregated deterministically", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  fit$neighbors$query_coords[2L, ] <- fit$neighbors$query_coords[1L, ]
  fit$neighbors$coords[2L, ] <- fit$neighbors$coords[1L, ]

  result <- plot_heatmap(
    fit, "income", resolution = 30L, return_data = TRUE
  )
  expect_equal(nrow(result$points), nrow(fit$coefficients))
  expect_equal(nrow(result$interpolation_points), nrow(fit$coefficients) - 1L)
  expect_equal(result$settings$duplicate_locations, 1L)
  expect_equal(
    result$interpolation_points$value_1[1L],
    mean(fit$coefficients[1:2, "income"])
  )
})

test_that("selection is distinct from coefficient value", {
  fit <- make_heatmap_fit(lambda = 0.35)$fit
  sparse <- gwrs:::extract_gwrs_spatial_values(fit, "income")
  expect_equal(
    sparse$points$selected,
    abs(fit$coefficients[, "income"]) > 1e-10
  )

  ridge <- make_heatmap_fit(lambda = 0.35, alpha = 0)$fit
  ridge_values <- gwrs:::extract_gwrs_spatial_values(ridge, "income")
  expect_true(all(ridge_values$points$selected))
})
