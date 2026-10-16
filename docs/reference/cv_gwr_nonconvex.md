# Spatial cross-validation for GWR-SCAD and GWR-MCP

Evaluates a common decreasing lambda path using held-out spatial blocks.
Standardization and nearest-neighbor construction are repeated inside
each training fold, preventing transformation and neighborhood leakage.

## Usage

``` r
cv_gwr_nonconvex(
  x,
  y,
  coords,
  k,
  fold,
  kernel = c("bisquare", "gaussian", "exponential", "tricube", "boxcar"),
  bandwidth = NULL,
  lambda = NULL,
  n_lambda = 50L,
  lambda_min_ratio = 0.001,
  penalty = c("scad", "mcp"),
  gamma = NULL,
  standardize = TRUE,
  screening = TRUE,
  metric = c("rmse", "mse", "mae"),
  refit = TRUE,
  control = gwrs_control(),
  cv_cache = NULL
)

cv_gwr_scad(x, y, coords, k, fold, ..., gamma = 3.7)

cv_gwr_mcp(x, y, coords, k, fold, ..., gamma = 3)
```

## Arguments

- x:

  Numeric design matrix without an intercept.

- y:

  Numeric response vector.

- coords:

  Optional coordinate matrix used to construct neighbors.

- k:

  Number of neighbors when \`neighbors\` is not supplied.

- fold:

  User-supplied spatial fold labels, one per observation. Create a
  reproducible assignment with \[spatial_folds()\] or supply
  application- specific blocks/buffers.

- kernel:

  Spatial kernel.

- bandwidth:

  Fixed distance bandwidth. The default \`NULL\` uses the farthest
  retained neighbor as an adaptive bandwidth at each target.

- lambda:

  Decreasing vector of nonnegative regularization strengths. When
  \`NULL\`, a logarithmic sequence is generated from
  \[gwr_sl_lambda_max()\].

- n_lambda:

  Number of automatically generated values.

- lambda_min_ratio:

  Ratio of the smallest to largest generated value.

- penalty:

  Nonconvex penalty, either \`"scad"\` or \`"mcp"\`.

- gamma:

  Concavity parameter. The default is \`3.7\` for SCAD and \`3\` for
  MCP. SCAD requires \`gamma \> 2\`; MCP requires \`gamma \> 1\`.

- standardize:

  Whether predictors are standardized globally before the local fits.
  Coefficients are returned on the original scale.

- screening:

  Whether to use sequential screening. A full KKT scan is always
  performed, so screened variables that violate optimality are restored
  before a solution is accepted.

- metric:

  Criterion used for selecting the final pair: root mean squared error
  (\`"rmse"\`), mean squared error (\`"mse"\`), or mean absolute error
  (\`"mae"\`). All three are returned.

- refit:

  Whether to fit the selected model on all observations.

- control:

  Numerical controls from \[gwrs_control()\].

- cv_cache:

  Optional read-only preparation from
  [`gwr_nonconvex_cv_cache`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_cv_cache.md).
  Inputs must match exactly; responses and penalty paths are not cached.

- ...:

  Additional arguments passed by \[cv_gwr_scad()\] or \[cv_gwr_mcp()\]
  to \`cv_gwr_nonconvex()\`.

## Value

An object inheriting from \`gwrs_nonconvex_cv\` and \`gwrs_cv\`.

## Details

With the guarded-block solver, a lambda is eligible only if all held-out
observations have finite, converged predictions. Otherwise all its error
scores are \`NA\`; no observation is removed. \`cv\$observation_count\`
retains the planned full denominator and \`cv\$valid\` reports
eligibility. If none are valid, the object has zero rows in \`best\`, a
\`NULL\` fit and a warning. Refit uses the existing single-lambda
initialization, not the CV path's warm start, and is independently
checked for convergence.

## Examples

``` r
quake <- datasets::quakes[1:30, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
fold <- spatial_folds(coords, n_folds = 3, seed = 42)
cv <- cv_gwr_mcp(
  x, quake$mag, coords, k = 8, fold = fold,
  lambda = c(0.12, 0.04), refit = FALSE,
  control = gwrs_control(n_threads = 1, diagnostics = "none")
)
cv$best
```
