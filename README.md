# gwrs: Fast Geographically Weighted Regression and Selection

Geographically weighted regression with **Ridge, Lasso, Elastic Net, SCAD and MCP**.

`gwrs` provides local regression coefficients, regularization paths, spatial
cross-validation and prediction at new locations. Its compiled kernels use
RcppArmadillo and RcppParallel, with reusable nearest-neighbor structures.
The package contains the GWR-SCAD and GWR-MCP estimators studied in
*Nonconvex Geographically Weighted Regression for Local Variable Selection*.
SCAD and MCP are established penalties; the paper concerns their geographically
weighted estimators and computational treatment.

| Estimator | Formula interface | Matrix interface | Purpose |
|:--|:--|:--|:--|
| GWR-Ridge | `gwr_ridge()` | `gwr_ridge_fit()` | Quadratic shrinkage |
| GWR-Lasso | `gwr_lasso()` | `gwr_lasso_fit()` | Shrinkage and local selection |
| GWR-Elastic Net | `gwr_en()` | `gwr_en_fit()` | Mixed quadratic and absolute-value penalties |
| GWR-SCAD | `gwr_scad()` | `gwr_scad_fit()` | Nonconvex local variable selection |
| GWR-MCP | `gwr_mcp()` | `gwr_mcp_fit()` | Nonconvex local variable selection |

## Installation

This repository is private. Cloning and repository-based installation require
access granted by the owner. From an authorized local checkout:

```sh
R CMD INSTALL .
```

R >= 4.1, a C++17 compiler and the dependencies in `DESCRIPTION` are required.
To install a supplied source archive instead:

```r
install.packages("gwrs_0.4.0.9013.tar.gz", repos = NULL, type = "source")
```

Install the declared dependencies before the source archive if they are absent.
The source-archive command does not download missing dependencies automatically.

## Fit the five methods

The example uses deterministic, dimensionless unit-square coordinates.
The small neighborhood and penalty choices illustrate the API and are not the
paper's simulation or application settings.

```r
library(gwrs)
data <- expand.grid(east = seq(0, 1, length.out = 8),
                    north = seq(0, 1, length.out = 8))
i <- seq_len(nrow(data))
data$x1 <- sin(i / 5)
data$x2 <- cos(i / 7)
data$y <- 1 + (0.8 + data$east) * data$x1 - 0.7 * data$x2 + 0.05 * sin(i)
train <- data[1:48, ]
test <- data[49:64, ]
control <- gwrs_control(n_threads = 1, diagnostics = "none")

ridge <- gwr_ridge(y ~ x1 + x2, train, c("east", "north"),
                   k = 18, lambda = 0.08, control = control)
lasso <- gwr_lasso(y ~ x1 + x2, train, c("east", "north"),
                   k = 18, lambda = 0.08, control = control)
enet <- gwr_en(y ~ x1 + x2, train, c("east", "north"),
               k = 18, lambda = 0.08, alpha = 0.5, control = control)
scad <- gwr_scad(y ~ x1 + x2, train, c("east", "north"),
                 k = 18, lambda = 0.08, gamma = 3.7, control = control)
mcp <- gwr_mcp(y ~ x1 + x2, train, c("east", "north"),
               k = 18, lambda = 0.08, gamma = 3, control = control)
head(coef(scad))
```

Here `k` is the number of retained neighbors, and `lambda` is the penalty
strength. The default kernel is adaptive bisquare. Predictor standardization
is enabled by default; coefficients are returned on the original predictor
scale. The intercept is not penalized. For Elastic Net, `alpha` controls the
mix of penalties. For SCAD/MCP, `gamma` controls penalty concavity; it is not a
bandwidth parameter. Nonconvex fits are local numerical solutions and do not
carry a universal global-optimality guarantee.

## Spatial cross-validation and prediction

Supply the same folds and observations when comparing methods. The following
example tunes MCP over a prespecified grid; use `penalty = "ridge"`, `"lasso"`,
`"elastic_net"` or `"scad"` for the other methods. With Elastic Net, supply the
chosen `alpha`; this example does not jointly tune bandwidth or alpha.

```r
coords <- as.matrix(train[, c("east", "north")])
x <- as.matrix(train[, c("x1", "x2")])
fold <- spatial_folds(coords, n_folds = 3, seed = 42)
cv <- cv_gwr_penalized(x, train$y, coords, k = 18, fold = fold,
  lambda = c(0.2, 0.08, 0.02), penalty = "mcp", refit = TRUE,
  control = control)
cv$best

pred <- predict(cv$fit,
  newdata = as.matrix(test[, c("x1", "x2")]),
  newcoords = as.matrix(test[, c("east", "north")]), k = 18)
head(pred)
```

Standardization and neighborhoods are constructed within each training fold.
New-location predictions use training observations only. Check convergence and
the returned CV diagnostics before interpreting fitted coefficients.

The [complete five-method example](inst/examples/five-methods.R) runs all five
fits, uses common CV folds and predicts at held-out locations. The
[user guide source](vignettes/gwrs-introduction.Rmd) is included with the package;
after installing a build containing vignettes, open it with:

```r
vignette("gwrs-introduction", package = "gwrs")
citation("gwrs")
```

## County dataset

The package includes the documented **2020–2024 ACS five-year county dataset**:
3,107 counties and 91 columns, including the study response, 26 predictors,
90% margins of error, controlled-MOE flags, geographic identifiers, projected
coordinates and spatial folds.

```r
data(acs2024_counties, package = "gwrs")
head(acs2024_counties[, c("fips", "county_name", "median_household_income_2024")])
help("acs2024_counties", package = "gwrs")
citation("gwrs")
```

The documentation describes sample exclusions, units and transformations.
The income MOE is in dollars even though its column name starts with the
log-response name. Please acknowledge the U.S. Census Bureau, the accompanying
paper and `gwrs` when using this prepared dataset. The complete variable
dictionary is installed in `system.file("extdata", package = "gwrs")`.

## Reproducing the paper

Current package version: **0.4.0.9013**. The paper's simulation used the recorded
**0.4.0** implementation and the county analysis used **0.4.0.9006**. Their
versioned source archives and analysis scripts are preserved separately; the
current package version must not be substituted silently in a historical run.
See [reproducibility records](reproducibility/README.md). A private repository
link does not grant access to reviewers; the separate code archive can be
supplied with the manuscript.

Author and maintainer: Bahadir Yuzbasi. License: GPL (>= 3).
