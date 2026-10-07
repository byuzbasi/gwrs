test_that("warm path agrees with independent scalar fits", {
  example <- make_spatial_example(n = 65, p = 4, seed = 41)
  neighbors <- gwr_neighbors(example$coords, k = 28)
  control <- serial_control()
  anchor <- gwr_sl_lambda_max(
    example$x, example$y, neighbors = neighbors, alpha = 0.7,
    control = control
  )
  lambda <- anchor * c(1, 0.3, 0.08)
  path <- gwr_sl_path(
    example$x, example$y, neighbors = neighbors, lambda = lambda,
    alpha = 0.7, d = 0.4, keep_coefficients = TRUE, control = control
  )
  scalar <- gwr_sl_fit(
    example$x, example$y, neighbors = neighbors, lambda = lambda[2L],
    alpha = 0.7, d = 0.4, control = control
  )
  expect_equal(path$fitted.values[, 2L], scalar$fitted.values,
               tolerance = 2e-7)
  expect_equal(unname(path$coefficients[, , 2L]),
               unname(scalar$coefficients),
               tolerance = 2e-7)
  expect_lt(max(path$path_summary$max_kkt_violation), 1e-6)
})

test_that("screening is corrected to the unscreened KKT solution", {
  example <- make_spatial_example(n = 52, p = 5, seed = 43)
  neighbors <- gwr_neighbors(example$coords, k = 27)
  lambda <- c(0.5, 0.2, 0.08, 0.03)
  screened <- gwr_sl_path(
    example$x, example$y, neighbors = neighbors, lambda = lambda,
    alpha = 0.75, d = 0.3, screening = TRUE, control = serial_control()
  )
  full <- gwr_sl_path(
    example$x, example$y, neighbors = neighbors, lambda = lambda,
    alpha = 0.75, d = 0.3, screening = FALSE, control = serial_control()
  )
  expect_equal(screened$fitted.values, full$fitted.values,
               tolerance = 2e-7)
})

test_that("spatial CV evaluates a common grid and refits the selected model", {
  example <- make_spatial_example(n = 60, p = 3, seed = 47)
  fold <- spatial_folds(example$coords, n_folds = 3, seed = 4)
  result <- cv_gwr_sl(
    example$x, example$y, example$coords, k = 15, fold = fold,
    lambda = c(0.3, 0.1, 0.03), alpha = 0.7, d = c(0, 0.5),
    control = serial_control(tolerance = 1e-7)
  )
  expect_s3_class(result, "gwrs_cv")
  expect_equal(nrow(result$cv), 6)
  expect_true(all(c("rmse", "mse", "mae") %in% names(result$cv)))
  expect_s3_class(result$fit, "gwrs_fit")
  expect_equal(result$fit$lambda, result$best$lambda)
  expect_equal(result$fit$d, result$best$d)
})

test_that("new-location prediction returns matching local coefficients", {
  example <- make_spatial_example(n = 70, p = 3, seed = 53)
  neighbors <- gwr_neighbors(example$coords, k = 24)
  fit <- gwr_sl_fit(
    example$x, example$y, neighbors = neighbors, lambda = 0.08,
    alpha = 0.7, d = 0.4, control = serial_control()
  )
  prediction <- predict(
    fit,
    newdata = example$x[1:6, , drop = FALSE],
    newcoords = example$coords[1:6, , drop = FALSE],
    k = 24,
    type = "both"
  )
  expect_length(prediction$fit, 6)
  expect_identical(dim(prediction$coefficients), c(6L, 4L))
  expect_equal(
    prediction$fit,
    rowSums(cbind(1, example$x[1:6, ]) * prediction$coefficients),
    tolerance = 1e-10
  )
})
