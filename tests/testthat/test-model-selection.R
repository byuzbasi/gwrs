test_that("joint search agrees with independent dense GWR and criterion algebra", {
  d <- selection_example()
  n <- nrow(d)
  result <- gwr_select(y ~ x1 + x2, d, c("east", "north"), k = c(18, n),
    f_tests = FALSE, control = gwrs_control(n_threads = 1), verbose = FALSE)
  expect_equal(nrow(result$ranking), 34L)
  expect_true(result$search_status$exhaustive_optimum)
  expect_identical(result$settings$criterion, "aicc")
  expect_true(all(grepl("weights_not_supplied", result$ranking$moran_status)))
  expect_true(all(result$ranking$eligible))
  criteria <- c("aic", "aicc", "bic", "gcv", "loo_cv")
  for (i in seq_len(nrow(result$ranking))) {
    row <- result$ranking[i, ]
    x <- model.matrix(as.formula(paste("y ~", row$terms)), d)
    if (row$family == "OLS") {
      s <- tcrossprod(qr.Q(qr(x)))
    } else {
      nb <- gwr_neighbors(as.matrix(d[, c("east", "north")]), k = row$k)
      s <- matrix(0, n, n)
      for (j in seq_len(n)) {
        idx <- nb$index[, j]
        u <- nb$distance[, j] / (max(nb$distance[, j]) * (1 + 1e-10))
        w <- switch(row$kernel, bisquare = (1 - u^2)^2, gaussian = exp(-u^2 / 2),
          exponential = exp(-u), tricube = (1 - u^3)^3, boxcar = rep(1, length(u)))
        local <- x[idx, , drop = FALSE]
        s[j, idx] <- x[j, , drop = FALSE] %*% solve(crossprod(local, local * w), t(local * w))
      }
    }
    residual <- d$y - drop(s %*% d$y)
    rss <- sum(residual^2)
    tr <- sum(diag(s))
    dev <- n * (log(rss / n) + log(2 * pi) + 1)
    expected <- c(aic = dev + 2 * (tr + 1), aicc = dev + 2 * n * (tr + 1) / (n - tr - 2),
      bic = dev + log(n) * (tr + 1), gcv = (rss / n) / (1 - tr / n)^2,
      loo_cv = sum((residual / (1 - diag(s)))^2))
    expect_equal(unlist(row[criteria], use.names = FALSE), unname(expected), tolerance = 1e-8)
    expect_equal(row$trStS, sum(s^2), tolerance = 1e-9)
  }
  gwr <- subset(result$ranking, family == "GWR")
  expect_equal(result$selected$aicc, min(gwr$aicc))
  expect_equal(result$selected_overall$aicc, min(result$ranking$aicc))
  expect_equal(nrow(result$best_by_criterion), 10L)
  expect_equal(coef(result), coef(result$best_model))
  expect_output(print(result), "best GWR")
  expect_s3_class(summary(result), "summary.gwrs_selection")
})

test_that("fixed kernels use full support and both bandwidth types are distinct", {
  d <- selection_example()
  result <- gwr_select(y ~ x1, d, c("east", "north"),
    kernels = c("gaussian", "exponential"), bandwidth = c("fixed", "adaptive"),
    bw = c(.4, 1.5), k = 20, f_tests = FALSE, verbose = FALSE,
    control = gwrs_control(n_threads = 1))
  fixed <- subset(result$ranking, bandwidth_type == "fixed")
  expect_true(all(fixed$k == nrow(d)))
  expect_true(all(fixed$support == "all_observations"))
  for (i in seq_len(nrow(fixed))) {
    ref <- gwr_fit(as.matrix(d["x1"]), d$y, coords = as.matrix(d[c("east", "north")]),
      k = nrow(d), kernel = fixed$kernel[i], bandwidth = fixed$bw[i],
      standardize = FALSE, control = gwrs_control(n_threads = 1))
    expect_equal(fixed$aicc[i], gwr_diagnostics(ref)$global$aicc, tolerance = 1e-10)
  }
  plan <- gwr_select(y ~ x1, d, c("east", "north"), kernels = "gaussian",
    bandwidth = "fixed", dry_run = TRUE)
  expect_equal(length(plan$settings$bw), 30L)
  expect_equal(plan$expected_candidates, 32L)
  expect_error(small_selection(max_models = 3), "max_models")
  expect_error(small_selection(memory_limit_mb = 1e-6), "workspace")
})

test_that("factor groups, hierarchy, required terms and predictions are preserved", {
  d <- selection_example(n = 36)
  d$group <- factor(rep(c("a", "b", "c"), 12))
  plan <- gwr_select(y ~ x1 * group, d, c("east", "north"), kernels = "gaussian",
    k = 36, required = "x1:group", dry_run = TRUE)
  expect_equal(plan$expected_candidates, 2L)
  expect_true(all(plan$grids$n_columns == 5))
  expect_false(any(plan$grids$terms == "1"))
  out <- gwr_select(y ~ x1 + group, d, c("east", "north"), kernels = "gaussian",
    k = 36, required = "group", f_tests = FALSE, verbose = FALSE)
  expect_true(all(out$ranking$n_columns %in% c(2, 3)))
  expect_true(all(grepl("group", out$ranking$terms)))
  pred <- predict(out, newdata = d[1:3, ], newcoords = as.matrix(d[1:3, c("east", "north")]))
  direct <- predict(out$best_model, newdata = d[1:3, ], newcoords = as.matrix(d[1:3, c("east", "north")]))
  expect_equal(pred, direct)
  expect_error(gwr_select(y ~ x1:group, d, c("east", "north"), dry_run = TRUE), "main effects")
  expect_error(gwr_select(y ~ x1 + offset(x2), d, c("east", "north")), "Offsets")
  d$x2[1] <- NA
  expect_error(small_selection(d), "missing")
})

test_that("invalid fits remain visible and unfinished selection is provisional", {
  d <- selection_example()
  d$x2 <- d$x1
  result <- small_selection(d)
  bad <- subset(result$ranking, n_columns == 2)
  expect_true(all(!bad$eligible))
  expect_true(all(nzchar(bad$reason)))
  expect_equal(result$selected$n_columns, 1L)
  none <- gwr_select(y ~ x1 + x2, d, c("east", "north"), required = c("x1", "x2"),
    kernels = "gaussian", k = 20, f_tests = FALSE, verbose = FALSE)
  expect_null(none$best_model)
  expect_identical(none$search_status$status, "no_valid_gwr")
  expect_error(coef(none), "No eligible")
  empty <- small_selection(max_evaluations = 0)
  expect_false(empty$search_status$complete)
  expect_equal(empty$search_status$evaluated, 0)
  expect_null(empty$best_model)
  limited <- small_selection(max_evaluations = 3)
  expect_equal(limited$search_status$evaluated, 3)
  expect_identical(limited$search_status$stop_reason, "invocation_budget")
})

test_that("diagnostics use the selected ordinary GWR and common user weights", {
  d <- selection_example()
  n <- nrow(d)
  w <- matrix(0, n, n)
  w[cbind(seq_len(n), c(seq.int(2, n), 1))] <- 1
  result <- gwr_select(y ~ x1, d, c("east", "north"), kernels = "gaussian", k = 20,
    moran_weights = w, moran_permutations = 19, moran_seed = 83,
    control = gwrs_control(n_threads = 1), verbose = FALSE)
  fit <- result$best_model
  tests <- gwr_diagnostic_tests(fit)
  expect_equal(result$selected$f1_p, tests$global_tests$p_value[1])
  expect_equal(result$selected$f2_p, tests$global_tests$p_value[2])
  expect_equal(nrow(result$coefficient_tests), 2L)
  m <- gwr_moran(fit, weights = w, permutations = 19, seed = 83)
  z <- residuals(fit) - mean(residuals(fit))
  expect_equal(m$I, n / sum(w) * sum(z * drop(w %*% z)) / sum(z^2), tolerance = 1e-12)
  expect_equal(result$selected$moran_I, m$I)
  expect_equal(result$selected$moran_p, m$p_value)
  fit$control$n_threads <- 2L
  expect_identical(m, gwr_moran(fit, weights = w, permutations = 19, seed = 83))
  both <- gwr_moran_diagnostics(fit, weights = Matrix::Matrix(w, sparse = TRUE), permutations = 19, seed = 83)
  expect_equal(both$global, m)
  expect_equal(both$local$local_I, drop(z * (w %*% z)) / mean(z^2))
  w[1, ] <- 0
  expect_error(gwr_moran(fit, weights = w), "islands")
  island <- gwr_local_moran(fit, weights = w, zero_policy = "keep", permutations = 0)
  expect_true(is.na(island$local_I[1]))
  expect_error(gwr_moran(fit, weights = diag(n)), "zero diagonal")
  w[2, 3] <- -1
  expect_error(gwr_moran(fit, weights = w), "nonnegative")
  dimnames(w) <- list(rev(rownames(d)), rev(rownames(d)))
  expect_error(small_selection(moran_weights = w), "row order")
})

test_that("verified resume exactly reproduces exhaustive and stepwise searches", {
  skip_if_not_installed("jsonlite")
  clean <- function(x) {
    x <- x$ranking[order(x$ranking$candidate_id), ]
    x$elapsed_seconds <- NULL
    rownames(x) <- NULL
    x
  }
  for (method in c("exhaustive", "stepwise")) {
    path <- tempfile("selection-checkpoint-")
    full <- small_selection(search = method)
    partial <- small_selection(search = method, checkpoint = path, checkpoint_every = 2,
      max_evaluations = 4)
    expect_false(partial$search_status$complete)
    expect_equal(gwr_select_verify(path)$durable_candidates, 4L)
    resumed <- small_selection(search = method, checkpoint = path, resume = TRUE)
    expect_equal(clean(resumed), clean(full))
    expect_true(resumed$search_status$complete)
    expect_equal(resumed$search_status$new_evaluations, nrow(full$ranking) - 4L)
    again <- small_selection(search = method, checkpoint = path, resume = TRUE)
    expect_equal(again$search_status$new_evaluations, 0)
    expect_equal(clean(again), clean(full))
    expect_true(gwr_select_verify(path)$valid)
    changed <- selection_example()
    changed$y[1] <- changed$y[1] + .1
    expect_error(small_selection(changed, search = method, checkpoint = path, resume = TRUE), "signature differs")
    expect_error(small_selection(search = method, checkpoint = path, resume = TRUE, criterion = "bic"), "signature differs")
    expect_error(small_selection(search = method, checkpoint = path), "already exists")
    pack <- list.files(path, pattern = "^pack-", full.names = TRUE)[1]
    con <- file(pack, "ab")
    writeBin(as.raw(0), con)
    close(con)
    expect_error(gwr_select_verify(path), "checksum mismatch")
  }
  d <- selection_example(p = 3)
  out <- gwr_select(y ~ x1 + x2 + x3, d, c("east", "north"), kernels = "gaussian",
    k = 22, search = "stepwise", f_tests = FALSE, verbose = FALSE)
  expect_true(out$search_status$complete)
  expect_false(out$search_status$exhaustive_optimum)
  expect_true(any(out$ranking$n_terms == 2))
})
