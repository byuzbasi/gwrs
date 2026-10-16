# Ridge, Lasso, Elastic Net, SCAD and MCP with gwrs

## Scope and interfaces

This guide covers geographically weighted ridge, lasso, elastic-net,
SCAD and MCP regression. All five estimators fit a local intercept and
local slopes. Ridge shrinks coefficients; Lasso, Elastic Net, SCAD and
MCP can also set local slopes to zero. SCAD and MCP use nonconvex
penalties, so convergence does not by itself establish a global optimum
or recovery of the true active set.

| Estimator | Formula     | Matrix          |
|:----------|:------------|:----------------|
| GWR-Ridge | gwr_ridge() | gwr_ridge_fit() |
| GWR-Lasso | gwr_lasso() | gwr_lasso_fit() |
| GWR-ENET  | gwr_en()    | gwr_en_fit()    |
| GWR-SCAD  | gwr_scad()  | gwr_scad_fit()  |
| GWR-MCP   | gwr_mcp()   | gwr_mcp_fit()   |

## Run a small example

The installed example constructs a deterministic grid on the unit
square, uses the first 48 locations for training and the last 16 for
prediction, and fits all five estimators. Coordinates and distances are
dimensionless. The example makes no statistical comparison or
performance claim.

The demonstration fixes the adaptive bisquare neighborhood at 18
training neighbors. It uses penalty 0.08 for the initial fits, Elastic
Net mixing 0.5, and the default SCAD/MCP concavity constants 3.7 and 3.
Spatial CV then uses three common folds (seed 42) and the explicit
penalty grid 0.2, 0.08, 0.02. These are example settings, not the
settings used in the manuscript.

``` r

example_file <- system.file("examples", "five-methods.R", package = "gwrs")
source(example_file, local = knitr::knit_global())
#> FIVE_METHOD_USAGE_PASS: fits, shared-fold CV, coefficients and new-location predictions
```

Inspect the selected penalty for each method:

``` r

lapply(cv, function(z) z$best)
#> $ridge
#>   lambda d      rmse       mse       mae convergence_rate
#> 3   0.02 0 0.3094922 0.0957854 0.2395486                1
#> 
#> $lasso
#>   lambda d      rmse       mse       mae convergence_rate
#> 3   0.02 0 0.3676823 0.1351902 0.2905299                1
#> 
#> $elastic_net
#>   lambda d     rmse       mse       mae convergence_rate
#> 3   0.02 0 0.358257 0.1283481 0.2776355                1
#> 
#> $scad
#>   lambda      rmse       mse       mae convergence_rate
#> 3   0.02 0.3362073 0.1130353 0.2535229                1
#> 
#> $mcp
#>   lambda      rmse       mse       mae convergence_rate
#> 3   0.02 0.3362073 0.1130353 0.2535229                1
```

Inspect SCAD coefficients and new-location predictions:

``` r

head(coef(cv$scad$fit))
#>      (Intercept)        x1         x2
#> [1,]   1.0116238 0.8978659 -0.6488478
#> [2,]   1.0082939 0.9520375 -0.6512968
#> [3,]   1.0044745 1.0496238 -0.6647299
#> [4,]   0.9990337 1.2171170 -0.7095883
#> [5,]   0.9680645 1.4490780 -0.7876932
#> [6,]   0.9215754 1.7015540 -0.8839481
head(as.data.frame(prediction))
#>         ridge      lasso elastic_net        scad         mcp
#> 1  0.43369460  0.7426695  0.56788760  0.22530960  0.22530958
#> 2  0.36477466  0.6552055  0.49266148  0.16179572  0.16179572
#> 3  0.20863631  0.3947105  0.29357838  0.04253704  0.04253704
#> 4  0.01394306  0.1098925  0.05861234 -0.13117765 -0.13117765
#> 5 -0.31935765 -0.3125693 -0.31622899 -0.42455670 -0.42455670
#> 6 -0.55231892 -0.5637020 -0.55797957 -0.62960456 -0.62960456
```

## Paths, tuning and interpretation

Use
[`gwr_ridge_path()`](../articles/gwrs-introduction.html#paths-tuning-and-interpretation),
[`gwr_lasso_path()`](../articles/gwrs-introduction.html#paths-tuning-and-interpretation),
[`gwr_en_path()`](../articles/gwrs-introduction.html#paths-tuning-and-interpretation),
[`gwr_scad_path()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_path.md)
or
[`gwr_mcp_path()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_path.md)
when the complete penalty path is needed. The predictors passed to
matrix interfaces exclude the intercept; the package handles it
separately. Predictor standardization is enabled by default and
coefficients are transformed back to the original predictor units.

``` r

path <- gwr_scad_path(x, train$y, coords = coords, k = 18,
  lambda = c(0.2, 0.08, 0.02), control = control)
path$lambda
#> [1] 0.20 0.08 0.02
```

[`cv_gwr_penalized()`](../articles/gwrs-introduction.html#paths-tuning-and-interpretation)
accepts an explicit `penalty` and returns CV diagnostics, the selected
tuning row and, with `refit = TRUE`, a fitted model. Dedicated
[`cv_gwr_scad()`](https://byuzbasi.github.io/gwrs/reference/cv_gwr_nonconvex.md)
and
[`cv_gwr_mcp()`](https://byuzbasi.github.io/gwrs/reference/cv_gwr_nonconvex.md)
interfaces are also available. Training-fold standardization and
training-only neighbors prevent preprocessing and neighborhood leakage.
Hold folds fixed for a paired method comparison.

The fixed penalty grid in this example is specified before CV. An
adaptive bandwidth or mixing-parameter search requires an explicitly
designed outer search; this guide does not perform such a search.
Examine convergence, candidate eligibility and numerical diagnostics
before interpreting results. The default SCAD/MCP solver remains
`coordinate`.

## Joint ordinary-GWR selection

Use
[`gwr_select()`](https://byuzbasi.github.io/gwrs/reference/gwr_select.md)
to compare unpenalized GWR predictor subsets, kernels and bandwidths on
a common sample. The [joint-selection
guide](https://byuzbasi.github.io/gwrs/articles/gwrs-selection.md)
provides an executable example, diagnostic limits and verified resume.
This search selects a predictor subset for each candidate model; local
penalized fits can instead have different active sets at different
locations.

F reference probabilities are exploratory and have no established
calibration, including for unpenalized fits. F1 and F2 share one
residual-sum-of-squares ratio; they do not provide independent
confirmations. These diagnostics do not rank candidate models or
validate a selected coefficient map.

## Paper-specific source versions

This guide describes the current package. Reproducing the paper requires
its recorded simulation (0.4.0) and county-analysis (0.4.0.9006) source
versions, settings, inputs and software records. Use the accompanying
reproduction archive for those studies. This demonstration is
intentionally small and is not a substitute for either numerical
experiment.
