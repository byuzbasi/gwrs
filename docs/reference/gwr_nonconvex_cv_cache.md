# Prepare reusable training-fold inputs for SCAD/MCP cross-validation

Prepares full-data and per-training-fold scaling and neighbor structures
for reuse by SCAD and MCP. Pass it as `cv_cache` to
[`cv_gwr_nonconvex`](https://byuzbasi.github.io/gwrs/reference/cv_gwr_nonconvex.md),
[`cv_gwr_scad`](https://byuzbasi.github.io/gwrs/reference/cv_gwr_nonconvex.md)
or
[`cv_gwr_mcp`](https://byuzbasi.github.io/gwrs/reference/cv_gwr_nonconvex.md).

## Usage

``` r
gwr_nonconvex_cv_cache(x, coords, k, fold, standardize = TRUE)
```

## Arguments

- x:

  Numeric predictor matrix with observations in rows.

- coords:

  Numeric coordinate matrix or an sf point object.

- k:

  Retained neighbor count, no larger than any training-fold sample.

- fold:

  One non-missing spatial-fold label per observation.

- standardize:

  Whether to center and scale predictors within each training fold.

## Value

A read-only `gwrs_nonconvex_cv_cache` environment with locked bindings.

## Details

Predictors, coordinates, neighbor count, fold labels and standardization
must match exactly, including row order. Each fold uses its own
training-only transformation. Responses, penalties, coefficients,
warm-start state and lambda paths are never cached. Automatic lambda
anchors retain their existing method-specific full-data definition.
Kernel weights are not cached, allowing preparation to be reused across
kernels and bandwidths. Holding all prepared folds costs additional
memory; `cv_cache = NULL` retains the one-fold-at-a-time route. No
global cache is created.
