test_that("kernel graph Laplacian is sparse, symmetric, and conservative", {
  example <- make_spatial_example(n = 40, seed = 61)
  neighbors <- gwr_neighbors(example$coords, k = 12)
  laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
  expect_s4_class(laplacian, "dgCMatrix")
  expect_true(Matrix::isSymmetric(laplacian, tol = 1e-12))
  expect_lt(max(abs(Matrix::rowSums(laplacian))), 1e-10)
  expect_true(all(Matrix::diag(laplacian) >= 0))
})

test_that("graph Laplacian supports dimensions beyond a 32-bit element count", {
  skip_if(.Machine$sizeof.pointer < 8L, "requires a 64-bit platform")
  expect_identical(gwrs:::cpp_arma_uword_bits(), 64L)
  n <- 65536L
  neighbors <- structure(
    list(
      index = matrix(seq_len(n), nrow = 1L),
      distance = matrix(0, nrow = 1L, ncol = n),
      n = n,
      n_train = n,
      n_target = n
    ),
    class = "gwrs_neighbors"
  )

  laplacian <- gwr_graph_laplacian(neighbors, kernel = "boxcar")
  expect_s4_class(laplacian, "dgCMatrix")
  expect_identical(dim(laplacian), c(n, n))
  expect_identical(Matrix::nnzero(laplacian), 0L)
})

test_that("GWR-KTDD satisfies the joint first-order equations", {
  example <- make_spatial_example(n = 24, p = 2, seed = 67)
  neighbors <- gwr_neighbors(example$coords, k = 12)
  laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
  source <- sin(seq_len(24) / 4)
  source <- source - mean(source)
  lambda_diffusion <- 2
  gamma <- 0.08
  fit <- gwr_ktdd_fit(
    example$x, example$y, neighbors = neighbors, kernel = "gaussian",
    diffusion_operator = laplacian, source = source,
    coefficient_laplacian = laplacian,
    lambda_diffusion = lambda_diffusion, gamma = gamma,
    solver = "matrix_free",
    coefficient_tolerance = 1e-10,
    coefficient_max_iterations = 4000L,
    standardize = FALSE, linear_tolerance = 1e-10,
    linear_max_iterations = 4000,
    control = serial_control(tolerance = 1e-8, max_iterations = 2000)
  )
  expect_true(fit$diagnostics$converged)
  expect_true(fit$diagnostics$coefficient_linear_converged)
  expect_true(fit$diagnostics$linear_converged)
  expect_equal(mean(fit$process), 0, tolerance = 1e-10)

  design <- cbind(1, example$x)
  coefficients <- fit$coefficients
  process <- fit$process
  gradient_coefficients <- matrix(0, nrow = 24, ncol = 3)
  gradient_process <- numeric(24)
  local_loss <- 0
  for (target in seq_len(24)) {
    index <- neighbors$index[, target]
    weight <- adaptive_weights(neighbors, target, "gaussian")
    residual <- as.double(
      design[index, , drop = FALSE] %*% coefficients[target, ] +
        process[index] - example$y[index]
    )
    gradient_coefficients[target, ] <- crossprod(
      design[index, , drop = FALSE], weight * residual
    )[, 1L]
    gradient_process[index] <- gradient_process[index] + weight * residual
    local_loss <- local_loss + 0.5 * sum(weight * residual^2)
  }
  gradient_coefficients <- gradient_coefficients +
    gamma * as.matrix(laplacian %*% coefficients)
  physical_residual <- as.double(laplacian %*% process - source)
  gradient_process <- gradient_process + lambda_diffusion *
    as.double(Matrix::t(laplacian) %*% physical_residual)
  gradient_process <- gradient_process - mean(gradient_process)

  expect_lt(max(abs(gradient_coefficients)), 2e-5)
  expect_lt(max(abs(gradient_process)), 2e-5)
  diffusion_loss <- 0.5 * lambda_diffusion * sum(physical_residual^2)
  smoothness_loss <- 0.5 * gamma *
    sum(coefficients * as.matrix(laplacian %*% coefficients))
  expect_equal(
    fit$objective,
    local_loss + diffusion_loss + smoothness_loss,
    tolerance = 1e-8
  )
})

test_that("KTDD serial and parallel reverse-neighbor reductions agree", {
  example <- make_spatial_example(n = 36, p = 2, seed = 71)
  neighbors <- gwr_neighbors(example$coords, k = 16)
  laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
  process <- sin(2 * pi * example$coords[, 1L])
  process <- process - mean(process)
  rectangular_operator <- laplacian[-nrow(laplacian), , drop = FALSE]
  source <- as.double(rectangular_operator %*% process)
  y <- example$y + process
  arguments <- list(
    x = example$x,
    y = y,
    neighbors = neighbors,
    kernel = "gaussian",
    diffusion_operator = rectangular_operator,
    source = source,
    coefficient_laplacian = laplacian,
    lambda_diffusion = 1.5,
    gamma = 0.06,
    solver = "matrix_free",
    coefficient_tolerance = 1e-10,
    coefficient_max_iterations = 3000L,
    standardize = FALSE,
    linear_tolerance = 1e-10,
    linear_max_iterations = 3000L
  )
  serial <- do.call(
    gwr_ktdd_fit,
    c(arguments, list(
      control = gwrs_control(
        tolerance = 1e-8, max_iterations = 2000L,
        n_threads = 1L, grain_size = 2L
      )
    ))
  )
  parallel <- do.call(
    gwr_ktdd_fit,
    c(arguments, list(
      control = gwrs_control(
        tolerance = 1e-8, max_iterations = 2000L,
        n_threads = 2L, grain_size = 2L
      )
    ))
  )

  expect_true(serial$diagnostics$converged)
  expect_true(parallel$diagnostics$converged)
  expect_true(serial$diagnostics$linear_converged)
  expect_true(parallel$diagnostics$linear_converged)
  expect_equal(parallel$coefficients, serial$coefficients, tolerance = 1e-7)
  expect_equal(parallel$process, serial$process, tolerance = 1e-7)
  expect_equal(parallel$objective, serial$objective, tolerance = 1e-8)
})

test_that("KTDD factorized and matrix-free coefficient solvers agree", {
  example <- make_spatial_example(n = 30, p = 2, seed = 72)
  neighbors <- gwr_neighbors(example$coords, k = 14)
  laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
  process <- sin(2 * pi * example$coords[, 1L])
  process <- process - mean(process)
  arguments <- list(
    x = example$x,
    y = example$y + process,
    neighbors = neighbors,
    kernel = "gaussian",
    diffusion_operator = laplacian,
    source = as.double(laplacian %*% process),
    coefficient_laplacian = laplacian,
    lambda_diffusion = 2,
    gamma = 0.08,
    standardize = FALSE,
    linear_tolerance = 1e-10,
    linear_max_iterations = 4000L,
    control = serial_control(tolerance = 1e-8, max_iterations = 2000L)
  )
  factorized <- do.call(
    gwr_ktdd_fit,
    c(arguments, list(
      solver = "factorized",
      factor_cache_limit_mib = Inf
    ))
  )
  matrix_free <- do.call(
    gwr_ktdd_fit,
    c(arguments, list(
      solver = "matrix_free",
      coefficient_tolerance = 1e-10,
      coefficient_max_iterations = 4000L
    ))
  )

  expect_true(factorized$diagnostics$converged)
  expect_true(matrix_free$diagnostics$converged)
  expect_true(matrix_free$diagnostics$coefficient_linear_converged)
  expect_equal(matrix_free$coefficients, factorized$coefficients,
               tolerance = 5e-7)
  expect_equal(matrix_free$process, factorized$process, tolerance = 5e-7)
  expect_equal(matrix_free$objective, factorized$objective, tolerance = 1e-9)
})

test_that("KTDD automatic solver protects the projected factor cache", {
  small <- gwrs:::resolve_ktdd_solver(
    "auto", factor_cache_limit_mib = 512,
    n = 100, p = 5, k = 20
  )
  large <- gwrs:::resolve_ktdd_solver(
    "auto", factor_cache_limit_mib = 512,
    n = 250000, p = 400, k = 100
  )

  expect_identical(small$resolved, "factorized")
  expect_identical(large$resolved, "matrix_free")
  expect_identical(large$armadillo_index_bits, 64L)
  expect_gt(large$factor_cache_mib, 512)
  expect_lt(large$approximate_working_mib, large$factor_cache_mib)
  expect_error(
    gwrs:::resolve_ktdd_solver(
      "factorized", factor_cache_limit_mib = 512,
      n = 250000, p = 400, k = 100
    ),
    "above.*factor_cache_limit_mib"
  )
})

test_that("KTDD rejects a signed matrix that is not a graph Laplacian", {
  example <- make_spatial_example(n = 12, p = 2, seed = 73)
  neighbors <- gwr_neighbors(example$coords, k = 8)
  laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
  signed_block <- Matrix::Matrix(
    matrix(c(
      1, 1, -2,
      1, 1, -2,
      -2, -2, 4
    ), nrow = 3, byrow = TRUE),
    sparse = TRUE
  )
  invalid <- Matrix::bdiag(
    signed_block,
    Matrix::Diagonal(nrow(example$x) - 3L, x = 0)
  )

  expect_error(
    gwr_ktdd_fit(
      example$x, example$y,
      neighbors = neighbors,
      kernel = "gaussian",
      diffusion_operator = laplacian,
      source = rep(0, nrow(example$x)),
      coefficient_laplacian = invalid,
      lambda_diffusion = 1,
      gamma = 0.1,
      control = serial_control()
    ),
    "nonpositive off-diagonal"
  )
})

test_that("KTDD centers only when the PDE operator has a constant null space", {
  example <- make_spatial_example(n = 20, p = 2, seed = 79)
  neighbors <- gwr_neighbors(example$coords, k = 12)
  coefficient_laplacian <- gwr_graph_laplacian(
    neighbors, kernel = "gaussian"
  )
  diffusion_operator <- Matrix::Diagonal(nrow(example$x), x = 1)
  process <- 0.4 + 0.2 * sin(2 * pi * example$coords[, 1L])
  fit <- gwr_ktdd_fit(
    example$x, example$y + process,
    neighbors = neighbors,
    kernel = "gaussian",
    diffusion_operator = diffusion_operator,
    source = process,
    coefficient_laplacian = coefficient_laplacian,
    lambda_diffusion = 10,
    gamma = 0.05,
    standardize = FALSE,
    linear_tolerance = 1e-10,
    linear_max_iterations = 3000L,
    control = serial_control(tolerance = 1e-8, max_iterations = 2000L)
  )

  expect_true(fit$diagnostics$converged)
  expect_true(fit$diagnostics$linear_converged)
  expect_false(fit$diagnostics$process_centered)
  expect_gt(abs(mean(fit$process)), 0.1)
})
