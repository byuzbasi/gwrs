# Warm-started GWR-SCAD and GWR-MCP paths

Fits a decreasing nonconvex regularization path. Local neighborhoods and
centered design columns are reused. Warm starts, persistent active sets,
sequential candidate screening, and a complete inactive-coordinate
objective-improvement scan reduce work while retaining a coordinatewise
minimum check at every accepted solution.

## Usage

``` r
gwr_nonconvex_path(
  x,
  y,
  coords = NULL,
  neighbors = NULL,
  k = NULL,
  kernel = c("bisquare", "gaussian", "exponential", "tricube", "boxcar"),
  bandwidth = NULL,
  lambda = NULL,
  n_lambda = 50L,
  lambda_min_ratio = 0.001,
  penalty = c("scad", "mcp"),
  gamma = NULL,
  standardize = TRUE,
  screening = TRUE,
  keep_coefficients = FALSE,
  return_diagnostics = FALSE,
  ebic_gamma = 0.5,
  control = gwrs_control()
)

gwr_scad_path(x, y, ..., gamma = 3.7)

gwr_mcp_path(x, y, ..., gamma = 3)
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

- keep_coefficients:

  Whether to retain the full locations-by-coefficients-by-lambda array.
  The default avoids a potentially large allocation.

- return_diagnostics:

  Whether to compute conditional smoother traces, information criteria,
  and post-selection refit diagnostics along the path. These use \`O(n
  \* length(lambda))\` memory and do not require retaining the
  coefficient cube.

- ebic_gamma:

  EBIC model-space multiplier in \`\[0, 1\]\`.

- control:

  Numerical controls from \[gwrs_control()\].

- ...:

  Additional arguments passed by \[gwr_scad_path()\] or
  \[gwr_mcp_path()\] to \`gwr_nonconvex_path()\`.

## Value

An object inheriting from \`gwrs_nonconvex_path\` and \`gwrs_path\`.

## Details

In guarded-block mode, a failed target stops its remaining lambda path.
\`path_state\` labels converged, failed and unattempted target/lambda
entries; \`solver_diagnostics\` retains objectives, stopping diagnostics
and compact pair/block counters. Failed coefficients, if requested, are
diagnostic iterates only; invalid fitted values and incomplete-path
statistical summaries are unavailable. A warning identifies partial
paths. The coordinate mode retains the original behavior.

## Examples

``` r
quake <- datasets::quakes[1:30, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
path <- gwr_mcp_path(
  x, quake$mag, coords = coords, k = 18,
  lambda = c(0.2, 0.08, 0.02),
  control = gwrs_control(n_threads = 1, diagnostics = "none")
)
path
path$path_summary
```
