test_that("neighbor structures use the intended target-contiguous layout", {
  example <- make_spatial_example(n = 30)
  neighbors <- gwr_neighbors(example$coords, k = 8)
  expect_s3_class(neighbors, "gwrs_neighbors")
  expect_identical(dim(neighbors$index), c(8L, 30L))
  expect_identical(dim(neighbors$distance), c(8L, 30L))
  expect_identical(neighbors$index[1L, ], seq_len(30L))
  expect_equal(neighbors$distance[1L, ], rep(0, 30))

  cross <- gwr_neighbors(
    example$coords[1:20, ],
    query_coords = example$coords[21:30, ],
    k = 5
  )
  expect_identical(dim(cross$index), c(5L, 10L))
  expect_identical(cross$n_train, 20L)
  expect_identical(cross$n_target, 10L)
  expect_true(all(cross$index >= 1L & cross$index <= 20L))
})

test_that("neighbor structures retain optional CRS metadata", {
  skip_if_not_installed("sf")
  example <- make_spatial_example(n = 20)
  coordinates <- example$coords
  attr(coordinates, "crs") <- sf::st_crs(3857)
  neighbors <- gwr_neighbors(coordinates, k = 8)
  expect_true(isTRUE(neighbors$crs == sf::st_crs(3857)))
  expect_true(isTRUE(neighbors$query_crs == sf::st_crs(3857)))

  points <- sf::st_as_sf(
    data.frame(x = example$coords[, 1L], y = example$coords[, 2L]),
    coords = c("x", "y"), crs = 4326
  )
  sf_neighbors <- gwr_neighbors(points, k = 8)
  expect_true(isTRUE(sf_neighbors$crs == sf::st_crs(4326)))
})

test_that("spatial folds are reproducible without changing caller RNG state", {
  example <- make_spatial_example(n = 50)
  set.seed(901)
  state <- .Random.seed
  first <- spatial_folds(example$coords, n_folds = 4, seed = 11)
  expect_identical(.Random.seed, state)
  second <- spatial_folds(example$coords, n_folds = 4, seed = 11)
  expect_identical(first, second)
  expect_identical(sort(unique(first)), 1:4)
})
