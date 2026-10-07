test_that("numeric boundaries mask every returned grid cell", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  boundary <- rbind(
    c(0.10, 0.10), c(0.90, 0.10), c(0.80, 0.85),
    c(0.20, 0.90)
  )
  expect_warning(
    result <- plot_heatmap(
      fit, "income", boundary = boundary, resolution = 35L,
      return_data = TRUE
    ),
    "no CRS"
  )
  grid_sf <- sf::st_as_sf(
    result$grid, coords = c("x", "y"), crs = sf::st_crs(result$boundary)
  )
  expect_true(all(lengths(sf::st_intersects(grid_sf, result$boundary)) > 0L))
  expect_lt(result$settings$inside_cells, result$settings$grid_cells)
})

test_that("mask false retains the complete boundary bounding box", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  boundary <- rbind(
    c(0.10, 0.10), c(0.90, 0.10), c(0.55, 0.90)
  )
  expect_warning(
    result <- plot_heatmap(
      fit, "income", boundary = boundary, mask = FALSE,
      resolution = 30L, return_data = TRUE
    ),
    "no CRS"
  )
  expect_equal(result$settings$inside_cells, result$settings$grid_cells)
  expect_equal(nrow(result$grid), result$settings$grid_cells)
  expect_true(all(result$grid$inside_boundary))
})

test_that("boundary CRS is transformed to the coordinate CRS", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  ring <- rbind(
    c(-0.05, -0.05), c(1.05, -0.05), c(1.05, 1.05),
    c(-0.05, 1.05), c(-0.05, -0.05)
  )
  boundary_4326 <- sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326)
  boundary_3857 <- sf::st_transform(boundary_4326, 3857)
  result <- plot_heatmap(
    fit, "income", boundary = boundary_3857, crs = 4326,
    distance = "geographic", resolution = 30L, return_data = TRUE
  )
  expect_true(isTRUE(sf::st_crs(result$boundary) == sf::st_crs(4326)))
})

test_that("multiple polygon outlines are retained for overlays", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  left <- rbind(
    c(0, 0), c(0.5, 0), c(0.5, 1), c(0, 1), c(0, 0)
  )
  right <- rbind(
    c(0.5, 0), c(1, 0), c(1, 1), c(0.5, 1), c(0.5, 0)
  )
  boundary <- sf::st_sfc(
    sf::st_polygon(list(left)), sf::st_polygon(list(right)), crs = 3857
  )
  result <- plot_heatmap(
    fit, "income", boundary = boundary, crs = 3857,
    resolution = 30L, return_data = TRUE
  )
  expect_length(result$boundary, 2L)
  expect_silent(ggplot2::ggplot_build(result$plot))
})

test_that("CRS omissions and incompatible distance choices are explicit", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  ring <- rbind(
    c(-0.05, -0.05), c(1.05, -0.05), c(1.05, 1.05),
    c(-0.05, 1.05), c(-0.05, -0.05)
  )
  known_boundary <- sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326)
  expect_error(
    plot_heatmap(
      fit, "income", boundary = known_boundary, resolution = 30L
    ),
    "model coordinates do not"
  )

  missing_boundary <- sf::st_sfc(sf::st_polygon(list(ring)))
  expect_warning(
    plot_heatmap(
      fit, "income", boundary = missing_boundary, crs = 4326,
      distance = "geographic", resolution = 30L
    ),
    "being assigned"
  )
  expect_warning(
    plot_heatmap(
      fit, "income", crs = 4326, distance = "euclidean",
      resolution = 30L
    ),
    "longitude/latitude"
  )
})

test_that("resolution and spatial controls are validated", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  expect_error(
    plot_heatmap(fit, "income", resolution = 24L),
    "between 25 and 1000"
  )
  expect_error(
    plot_heatmap(fit, "income", resolution = 30.5),
    "integer"
  )
  expect_warning(
    gwrs:::validate_gwrs_resolution(501L),
    "more than 250,000"
  )
  expect_error(
    plot_heatmap(fit, "income", resolution = 30L, max_neighbors = -Inf),
    "positive integer"
  )
})
