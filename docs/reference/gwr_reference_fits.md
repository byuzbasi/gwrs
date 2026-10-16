# Fit geographically weighted lasso, ridge, or elastic net

The three wrappers use the same normalized local least-squares loss,
neighbor structure and preprocessing convention. Ridge applies a
quadratic penalty, lasso an absolute-value penalty, and elastic net a
mixture controlled by alpha. Positive ridge uses a direct Cholesky
solve; an unsuccessful solve is reported as a failed fit. Lasso and
elastic net use coordinate updates. Coefficients are returned on the
original predictor scale.

## Usage

``` r
gwr_lasso_fit(
  x,
  y,
  coords = NULL,
  neighbors = NULL,
  k = NULL,
  kernel = c("bisquare", "gaussian", "exponential", "tricube", "boxcar"),
  bandwidth = NULL,
  lambda,
  standardize = TRUE,
  control = gwrs_control()
)

gwr_ridge_fit(
  x,
  y,
  coords = NULL,
  neighbors = NULL,
  k = NULL,
  kernel = c("bisquare", "gaussian", "exponential", "tricube", "boxcar"),
  bandwidth = NULL,
  lambda,
  standardize = TRUE,
  control = gwrs_control()
)

gwr_en_fit(
  x,
  y,
  coords = NULL,
  neighbors = NULL,
  k = NULL,
  kernel = c("bisquare", "gaussian", "exponential", "tricube", "boxcar"),
  bandwidth = NULL,
  lambda,
  alpha = 0.5,
  standardize = TRUE,
  control = gwrs_control()
)
```

## Arguments

- x:

  Numeric design matrix without an intercept.

- y:

  Numeric response vector.

- coords:

  Optional coordinate matrix used to construct neighbors.

- neighbors:

  Optional reusable object from \[gwr_neighbors()\].

- k:

  Number of neighbors when \`neighbors\` is not supplied.

- kernel:

  Spatial kernel.

- bandwidth:

  Fixed distance bandwidth. The default \`NULL\` uses the farthest
  retained neighbor as an adaptive bandwidth at each target.

- lambda:

  Nonnegative regularization strength.

- standardize:

  Whether predictors are standardized globally before the local fits.
  Coefficients are returned on the original scale.

- control:

  Numerical controls from \[gwrs_control()\].

- alpha:

  Mixing parameter in \`\[0, 1\]\`.

## Value

An object of class \`gwr_sl_fit\` and \`gwrs_fit\`.

## Examples

``` r
quake <- datasets::quakes[1:30, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
neighbors <- gwr_neighbors(coords, k = 18)
control <- gwrs_control(n_threads = 1)

lasso <- gwr_lasso_fit(
  x, quake$mag, neighbors = neighbors, lambda = 0.05, control = control
)
ridge <- gwr_ridge_fit(
  x, quake$mag, neighbors = neighbors, lambda = 0.05, control = control
)
elastic_net <- gwr_en_fit(
  x, quake$mag, neighbors = neighbors, lambda = 0.05,
  alpha = 0.6, control = control
)
c(lasso = lasso$alpha, ridge = ridge$alpha, elastic_net = elastic_net$alpha)
```
