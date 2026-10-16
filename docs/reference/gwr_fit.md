# Fit an unpenalized geographically weighted regression

Fit an unpenalized geographically weighted regression

## Usage

``` r
gwr_fit(
  x,
  y,
  coords = NULL,
  neighbors = NULL,
  k = NULL,
  kernel = c("bisquare", "gaussian", "exponential", "tricube", "boxcar"),
  bandwidth = NULL,
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

- standardize:

  Whether predictors are standardized globally before the local fits.
  Coefficients are returned on the original scale.

- control:

  Numerical controls from \[gwrs_control()\].

## Value

An object of class \`gwrs_fit\`.

## Examples

``` r
quake <- datasets::quakes[1:32, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
fit <- gwr_fit(
  x, quake$mag, coords = coords, k = 20,
  control = gwrs_control(n_threads = 1)
)
summary(fit)
```
