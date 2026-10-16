# Fast local multicollinearity diagnostics

Computes local condition numbers, variance-inflation factors, maximum
absolute local correlations, effective sample sizes, and positive kernel
counts directly from the reusable nearest-neighbor structure.

## Usage

``` r
gwr_local_collinearity(
  object,
  threshold_cn = 30,
  threshold_vif = 10,
  center = TRUE,
  scale = TRUE,
  return_pairwise = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted \`gwrs_fit\` retaining its design matrix.

- threshold_cn:

  Condition-number flag threshold.

- threshold_vif:

  VIF flag threshold.

- center:

  Whether predictors are locally centered.

- scale:

  Whether predictors are locally scaled.

- return_pairwise:

  Whether all pairwise local correlations are returned.

- ...:

  Unused.

## Value

A \`gwrs_local_collinearity\` object.

## Details

With the default centering, CN is the square root of the ratio of
largest to smallest eigenvalues of the weighted predictor correlation
matrix. The intercept is excluded. For comparable global diagnostics use
\[gwr_global_collinearity()\]. The rank tolerance is the largest
eigenvalue (or one, if larger) times \`1e-10\`. Singular designs return
infinite CN/VIF. Constant local predictors have undefined correlations
and infinite CN/VIF. Even with \`scale = FALSE\`, the covariance is
converted to a correlation matrix before computing CN/VIF. These
predictor diagnostics differ from the local design condition number in
\[gwr_diagnostics()\], which includes the intercept and uses the model's
working scale.

## See also

\[summary.gwrs_local_collinearity()\],
\[plot.gwrs_local_collinearity()\]

## Examples

``` r
quake <- datasets::quakes[1:36, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
fit <- gwr_fit(
  x, quake$mag, coords = coords, k = 24,
  control = gwrs_control(n_threads = 1)
)
collinearity <- gwr_local_collinearity(fit, return_pairwise = TRUE)
head(collinearity$condition_number)
head(collinearity$vif_long)
```
