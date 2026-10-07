ridge_native_path <- function(x, y, target_x, neighbors, lambda,
                              kernel = "bisquare", threads = 1L,
                              diagnostics = FALSE, adaptive = TRUE, bandwidth = 1) {
  gwrs:::cpp_gwr_sl_path_predict(x, y, target_x, neighbors$index, neighbors$distance,
    gwrs:::kernel_code(kernel)$code, adaptive, bandwidth, lambda, 0, 0,
    1e-7, 1L, TRUE, TRUE, diagnostics, threads, 2L)
}

ridge_reference <- function(x, y, target_x, neighbors, lambda, kernel = "bisquare",
                            dual = FALSE, adaptive = TRUE, bandwidth = 1) {
  n <- nrow(target_x); p <- ncol(x)
  b <- matrix(NA_real_, n, p + 1L)
  for (i in seq_len(n)) {
    ids <- neighbors$index[, i]
    d <- neighbors$distance[, i]
    h <- if (adaptive) max(d) * (1 + 1e-10) else bandwidth
    r <- d / h
    w <- switch(kernel, bisquare = ifelse(r < 1, (1 - r^2)^2, 0),
                gaussian = exp(-r^2 / 2), boxcar = as.numeric(r <= 1))
    w <- w / sum(w)
    xm <- colSums(x[ids, , drop = FALSE] * w); ym <- sum(y[ids] * w)
    z <- sweep(x[ids, , drop = FALSE], 2, xm) * sqrt(w)
    t <- (y[ids] - ym) * sqrt(w)
    beta <- if (dual) crossprod(z, solve(tcrossprod(z) + diag(lambda, nrow(z)), t)) else
      qr.solve(rbind(z, diag(sqrt(lambda), p)), c(t, rep(0, p)))
    b[i, ] <- c(ym - sum(xm * beta), beta)
  }
  list(coefficients = b, predictions = rowSums(cbind(1, target_x) * b))
}

test_that("weak-penalty HD Ridge matches an independent dual solution", {
  set.seed(20260904)
  x <- matrix(rnorm(48 * 64), 48, 64)
  coords <- cbind(runif(48), runif(48))
  y <- 2 + 3 * x[, 1] - 2 * x[, 2] + rnorm(48, sd = .2)
  prep <- gwrs:::standardize_design(x[1:40, ], TRUE)
  target <- gwrs:::apply_standardization(x[41:48, ], prep$center, prep$scale)
  nn <- gwr_neighbors(coords[1:40, ], query_coords = coords[41:48, ], k = 40, include_self = FALSE)
  lambda <- exp(seq(log(1e4), log(1e-4), length.out = 4))
  path <- ridge_native_path(prep$x, y[1:40], target, nn, lambda, "boxcar")
  for (j in seq_along(lambda)) {
    reference <- ridge_reference(prep$x, y[1:40], target, nn, lambda[j], "boxcar", dual = TRUE)
    expect_equal(path$predictions[, j], reference$predictions, tolerance = 1e-8)
    expect_equal(unname(t(path$coefficients[, , j])), reference$coefficients, tolerance = 1e-8)
  }
  expect_true(all(path$converged == 1L))
  expect_true(all(path$iterations == 1L))
  expect_true(all(path$status == 0L))
  expect_lt(max(path$kkt), 1e-9)
})

test_that("primal and dual weighted Ridge agree with augmented QR", {
  for (p in c(3L, 28L)) {
    e <- make_spatial_example(n = 30, p = p, seed = 400 + p)
    nn <- gwr_neighbors(e$coords, query_coords = e$coords[1:4, ], k = 16, include_self = FALSE)
    for (kernel in c("bisquare", "gaussian", "boxcar")) {
      path <- ridge_native_path(e$x, e$y, e$x[1:4, ], nn, c(.8, .0001), kernel)
      for (j in 1:2) {
        reference <- ridge_reference(e$x, e$y, e$x[1:4, ], nn, c(.8, .0001)[j], kernel)
        expect_equal(unname(t(path$coefficients[, , j])), reference$coefficients, tolerance = 1e-8)
      }
    }
  }
})

test_that("single-fit, path and prediction preserve original coefficient scale", {
  e <- make_spatial_example(n = 18, p = 22, seed = 503)
  e$x <- sweep(e$x, 2, seq(.5, 3, length.out = 22), "*") + 4
  nn <- gwr_neighbors(e$coords, k = 12)
  control <- gwrs_control(n_threads = 1, max_iterations = 1, diagnostics = "none")
  fit <- gwr_ridge_fit(e$x, e$y, neighbors = nn, lambda = .0001, control = control)
  path <- gwr_sl_path(e$x, e$y, neighbors = nn, alpha = 0, d = 0,
    lambda = c(.2, .0001), keep_coefficients = TRUE, control = control)
  expect_equal(unname(fit$coefficients), unname(path$coefficients[, , 2]), tolerance = 1e-8)
  pred <- predict(fit, newdata = e$x[1:3, ], newcoords = e$coords[1:3, ], k = 12, type = "both")
  expect_equal(pred$fit, fit$fitted.values[1:3], tolerance = 1e-8)
  expect_equal(rowSums(cbind(1, e$x[1:3, ]) * pred$coefficients), pred$fit, tolerance = 1e-10)
  expect_true(all(fit$diagnostics$converged))
  expect_equal(dim(fit$gwr_center), c(18L, 22L))
})

test_that("duplicated and locally constant columns need no Ridge pseudoinverse", {
  e <- make_spatial_example(n = 20, p = 3, seed = 607)
  e$x <- cbind(e$x, e$x[, 1], 1, e$x[, 2] + 1e-11 * e$x[, 3])
  nn <- gwr_neighbors(e$coords, query_coords = e$coords[1:3, ], k = 5, include_self = FALSE)
  path <- ridge_native_path(e$x, e$y, e$x[1:3, ], nn, c(1, .0001))
  ref <- ridge_reference(e$x, e$y, e$x[1:3, ], nn, .0001)
  expect_equal(unname(t(path$coefficients[, , 2])), ref$coefficients, tolerance = 1e-8)
  expect_true(all(path$status == 0))
  expect_equal(path$coefficients[2, , ], path$coefficients[5, , ], tolerance = 1e-8)
  expect_equal(as.double(path$coefficients[6, , ]), rep(0, 6), tolerance = 1e-12)
})

test_that("the dual path supports 1101 columns in a small smoke fixture", {
  e <- make_spatial_example(n = 18, p = 1101, seed = 709)
  nn <- gwr_neighbors(e$coords, query_coords = e$coords[1:2, ], k = 14, include_self = FALSE)
  path <- ridge_native_path(e$x, e$y, e$x[1:2, ], nn, c(1, .0001))
  ref <- ridge_reference(e$x, e$y, e$x[1:2, ], nn, .0001, dual = TRUE)
  expect_equal(path$predictions[, 2], ref$predictions, tolerance = 1e-8)
  expect_identical(dim(path$coefficients), c(1102L, 2L, 2L))
  expect_true(all(path$converged == 1L))
})

test_that("serial/parallel and prefix/full direct Ridge paths agree", {
  e <- make_spatial_example(n = 26, p = 31, seed = 811)
  nn <- gwr_neighbors(e$coords, k = 16)
  a <- ridge_native_path(e$x, e$y, e$x, nn, c(1, .03, .0001))
  b <- ridge_native_path(e$x, e$y, e$x, nn, c(1, .03, .0001), threads = 2L)
  prefix <- ridge_native_path(e$x, e$y, e$x, nn, c(1, .03))
  expect_equal(a$coefficients, b$coefficients, tolerance = 1e-12)
  expect_equal(a$predictions[, 1:2], prefix$predictions, tolerance = 1e-12)
})

test_that("fixed bandwidth and zero-weight failure retain their contract", {
  e <- make_spatial_example(n = 20, p = 4, seed = 907)
  nn <- gwr_neighbors(e$coords, query_coords = e$coords[1:2, ] + 10, k = 12, include_self = FALSE)
  good <- ridge_native_path(e$x, e$y, e$x[1:2, ], nn, .1, "gaussian", adaptive = FALSE, bandwidth = 10)
  ref <- ridge_reference(e$x, e$y, e$x[1:2, ], nn, .1, "gaussian", adaptive = FALSE, bandwidth = 10)
  expect_equal(as.double(good$predictions), ref$predictions, tolerance = 1e-8)
  bad <- ridge_native_path(e$x, e$y, e$x[1:2, ], nn, .1, adaptive = FALSE, bandwidth = .001)
  expect_true(all(!is.finite(bad$predictions)))
  expect_true(all(bad$converged == 0L))
  expect_true(all(bad$status == 3L))
})

test_that("unpenalized endpoint and optional path diagnostics remain available", {
  e <- make_spatial_example(n = 22, p = 3, seed = 1009)
  nn <- gwr_neighbors(e$coords, k = 15)
  a <- ridge_native_path(e$x, e$y, e$x, nn, c(.2, 0), diagnostics = TRUE)
  ordinary <- gwr_fit(e$x, e$y, neighbors = nn, standardize = FALSE,
                      control = gwrs_control(n_threads = 1, diagnostics = "none"))
  expect_equal(a$predictions[, 2], ordinary$fitted.values, tolerance = 1e-8)
  expect_true(all(is.finite(a$hat_diag_path)))
  expect_true(all(is.finite(a$post_fitted_path)))
  expect_equal(a$predictions[, 1], ridge_native_path(e$x, e$y, e$x, nn, .2)$predictions[, 1], tolerance = 1e-9)
})
