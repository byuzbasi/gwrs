# Small deterministic usage example, not a simulation result or model comparison.
library(gwrs)

# Unit-square coordinates are dimensionless; all distances use this same scale.
locations <- expand.grid(east = seq(0, 1, length.out = 8),
                         north = seq(0, 1, length.out = 8))
index <- seq_len(nrow(locations))
locations$x1 <- sin(index / 5)
locations$x2 <- cos(index / 7)
locations$y <- 1 + (0.8 + locations$east) * locations$x1 -
  0.7 * locations$x2 + 0.05 * sin(index)
train <- locations[1:48, ]
test <- locations[49:64, ]
coords <- as.matrix(train[, c("east", "north")])
x <- as.matrix(train[, c("x1", "x2")])
control <- gwrs_control(n_threads = 1, diagnostics = "none")

# These illustrative settings are shared across methods, not the paper's grids.
fits <- list(
  ridge = gwr_ridge(y ~ x1 + x2, train, c("east", "north"),
    k = 18, lambda = 0.08, control = control),
  lasso = gwr_lasso(y ~ x1 + x2, train, c("east", "north"),
    k = 18, lambda = 0.08, control = control),
  elastic_net = gwr_en(y ~ x1 + x2, train, c("east", "north"),
    k = 18, lambda = 0.08, alpha = 0.5, control = control),
  scad = gwr_scad(y ~ x1 + x2, train, c("east", "north"),
    k = 18, lambda = 0.08, gamma = 3.7, control = control),
  mcp = gwr_mcp(y ~ x1 + x2, train, c("east", "north"),
    k = 18, lambda = 0.08, gamma = 3, control = control)
)
stopifnot(all(vapply(fits, function(fit) all(is.finite(coef(fit))), logical(1))))
head(coef(fits$scad))

# Every method uses the same observations, folds, neighborhood and lambda grid.
# Predictor standardization is fitted separately within each training fold.
fold <- spatial_folds(coords, n_folds = 3, seed = 42)
lambda <- c(0.2, 0.08, 0.02)
cv <- lapply(names(fits), function(method) {
  cv_gwr_penalized(x, train$y, coords, k = 18, fold = fold,
    lambda = lambda, penalty = method, alpha = 0.5, refit = TRUE,
    control = control)
})
names(cv) <- names(fits)
stopifnot(all(vapply(cv, function(z) nrow(z$best) == 1L && !is.null(z$fit), logical(1))))
lapply(cv, function(z) z$best)

# New-location predictions use training observations and training neighbors only.
prediction <- lapply(cv, function(z) predict(z$fit,
  newdata = as.matrix(test[, c("x1", "x2")]),
  newcoords = as.matrix(test[, c("east", "north")]), k = 18))
stopifnot(all(vapply(prediction, function(z) length(z) == 16L && all(is.finite(z)), logical(1))))
head(as.data.frame(prediction))
cat("FIVE_METHOD_USAGE_PASS: fits, shared-fold CV, coefficients and new-location predictions\n")
