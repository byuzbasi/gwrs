# Fast nonconvex geographically weighted regression

Fits Gaussian geographically weighted regression with either the
smoothly clipped absolute deviation (SCAD) or minimax concave penalty
(MCP). At target location \`i\`, the profiled local objective is
\$\$\frac{1}{2c_i}(y-X\beta_i)^T W_i(y-X\beta_i)+ \sum\_{m=1}^p
p\_{\lambda,\gamma}(\|\beta\_{im}\|),\$\$ where \`c_i\` is the sum of
local kernel weights. The intercept is locally profiled and is never
penalized. Predictor coefficients are optimized in a parallel
RcppArmadillo core using exact one-coordinate SCAD/MCP minimizers.

## Usage

``` r
gwr_nonconvex_fit(
  x,
  y,
  coords = NULL,
  neighbors = NULL,
  k = NULL,
  kernel = c("bisquare", "gaussian", "exponential", "tricube", "boxcar"),
  bandwidth = NULL,
  lambda,
  penalty = c("scad", "mcp"),
  gamma = NULL,
  standardize = TRUE,
  control = gwrs_control()
)

gwr_scad_fit(x, y, ..., gamma = 3.7)

gwr_mcp_fit(x, y, ..., gamma = 3)
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

- penalty:

  Nonconvex penalty, either \`"scad"\` or \`"mcp"\`.

- gamma:

  Concavity parameter. The default is \`3.7\` for SCAD and \`3\` for
  MCP. SCAD requires \`gamma \> 2\`; MCP requires \`gamma \> 1\`.

- standardize:

  Whether predictors are standardized globally before the local fits.
  Coefficients are returned on the original scale.

- control:

  Numerical controls from \[gwrs_control()\].

- ...:

  Additional arguments passed by \[gwr_scad_fit()\] or \[gwr_mcp_fit()\]
  to \`gwr_nonconvex_fit()\`.

## Value

An object inheriting from \`gwr_nonconvex_fit\` and \`gwrs_fit\`.

## Details

With \`control = gwrs_control(nonconvex_solver = "guarded_block")\`, the
same objective also admits guarded pair and small flat-tail QR
proposals. This is an optimizer choice, not another penalty. Solutions
can differ because the objective is nonconvex. A failed guarded fit
raises a \`gwrs_nonconvex_convergence_error\` condition containing
\`solver_result\`; it is not used to construct statistical diagnostics
or predictions.

## References

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. \*Journal of the American
Statistical Association\*, 96, 1348–1360.
[doi:10.1198/016214501753382273](https://doi.org/10.1198/016214501753382273)
.

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. \*The Annals of Statistics\*, 38, 894–942.
[doi:10.1214/09-AOS729](https://doi.org/10.1214/09-AOS729) .

## Examples

``` r
quake <- datasets::quakes[1:32, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
fit <- gwr_mcp_fit(
  x, quake$mag, coords = coords, k = 20, lambda = 0.05,
  control = gwrs_control(n_threads = 1, diagnostics = "none")
)
fit
head(coef(fit))
```
