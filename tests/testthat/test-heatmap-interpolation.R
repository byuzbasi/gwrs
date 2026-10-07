test_that("kernel interpolation equals a hand calculation", {
  source <- rbind(c(0, 0), c(1, 0), c(0, 1))
  query <- matrix(c(0.25, 0.25), nrow = 1L)
  values <- matrix(c(1, 2, 4), ncol = 1L)
  colnames(values) <- "value_1"
  bandwidth <- 0.7
  distance <- sqrt(rowSums((source - query[1L, ])^2))
  weights <- exp(-0.5 * (distance / bandwidth)^2)
  expected <- sum(weights * values[, 1L]) / sum(weights)

  result <- gwrs:::gwrs_interpolate_value_matrix(
    source, query, values,
    interpolation = "kernel",
    distance = "euclidean",
    kernel = "gaussian",
    bandwidth = bandwidth,
    idw_power = 2,
    max_neighbors = Inf,
    chunk_size = 10L
  )
  expect_equal(
    unname(result$values[1L, 1L]), expected, tolerance = 1e-12
  )
})

test_that("IDW returns the observation at an exact location", {
  source <- rbind(c(0, 0), c(1, 0), c(0, 1))
  values <- matrix(c(1, 2, 4), ncol = 1L)
  result <- gwrs:::gwrs_interpolate_value_matrix(
    source, source[2L, , drop = FALSE], values,
    interpolation = "idw",
    distance = "euclidean",
    kernel = "gaussian",
    bandwidth = NULL,
    idw_power = 2,
    max_neighbors = Inf,
    chunk_size = 10L
  )
  expect_equal(result$values[1L, 1L], 2, tolerance = 0)
})

test_that("all supported kernels return finite deterministic results", {
  source <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  query <- rbind(c(0.2, 0.3), c(0.8, 0.7))
  values <- matrix(c(1, 2, 4, 8), ncol = 1L)
  for (kernel in c("gaussian", "bisquare", "exponential")) {
    first <- gwrs:::gwrs_interpolate_value_matrix(
      source, query, values, "kernel", "euclidean", kernel,
      bandwidth = 2, idw_power = 2, max_neighbors = Inf,
      chunk_size = 10L
    )
    second <- gwrs:::gwrs_interpolate_value_matrix(
      source, query, values, "kernel", "euclidean", kernel,
      bandwidth = 2, idw_power = 2, max_neighbors = Inf,
      chunk_size = 10L
    )
    expect_true(all(is.finite(first$values)))
    expect_identical(first, second)
  }
})

test_that("automatic and user bandwidths are reported", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  automatic <- plot_heatmap(
    fit, "income", resolution = 30L, return_data = TRUE
  )
  supplied <- suppressWarnings(plot_heatmap(
    fit, "income", bandwidth = 0.75, resolution = 30L,
    return_data = TRUE
  ))

  expect_true(is.finite(automatic$bandwidth))
  expect_gt(automatic$bandwidth, 0)
  expect_equal(supplied$bandwidth, 0.75)
  expect_equal(attr(automatic$plot, "gwrs_bandwidth"), automatic$bandwidth)
  expect_error(
    plot_heatmap(fit, "income", bandwidth = 0, resolution = 30L),
    "positive finite"
  )
})

test_that("public IDW and kernel surfaces are reproducible", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  first <- plot_heatmap(
    fit, "income", resolution = 32L, return_data = TRUE
  )
  second <- plot_heatmap(
    fit, "income", resolution = 32L, return_data = TRUE
  )
  idw <- plot_heatmap(
    fit, "income", interpolation = "idw", resolution = 32L,
    return_data = TRUE
  )

  expect_identical(first$grid$value, second$grid$value)
  expect_true(all(is.finite(idw$grid$value[idw$grid$interpolation_valid])))
  expect_null(idw$bandwidth)
})

test_that("geographic distances use earth-surface metres", {
  nearest <- gwrs:::gwrs_knn_query(
    source = rbind(c(0, 0), c(2, 0)),
    query = matrix(c(1, 0), nrow = 1L),
    k = 1L,
    distance = "geographic"
  )
  expect_equal(nearest$distance[1L, 1L], 111195, tolerance = 250)

  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  result <- plot_heatmap(
    fit, "income", distance = "geographic", crs = 4326,
    resolution = 30L, return_data = TRUE
  )
  expect_identical(result$settings$distance_units, "metres")
  expect_gt(result$bandwidth, 1000)
})

test_that("unsupported hyperbolic interpolation fails explicitly", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  expect_error(
    plot_heatmap(
      fit, "income", distance = "hyperbolic", resolution = 30L
    ),
    "does not provide a valid mapping"
  )
})
