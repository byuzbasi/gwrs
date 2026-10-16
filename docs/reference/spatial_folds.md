# Construct reproducible spatially contiguous folds

Uses k-means clustering of standardized coordinates to form spatial
blocks. The spatial CV routines never substitute ordinary random folds
when a spatial fold assignment is absent.

## Usage

``` r
spatial_folds(coords, n_folds = 5L, seed = 1L, nstart = 20L)
```

## Arguments

- coords:

  Numeric coordinate matrix.

- n_folds:

  Number of spatial folds.

- seed:

  Integer seed used by k-means initialization. The caller's random
  number state is restored before returning.

- nstart:

  Number of k-means starts.

## Value

An integer fold vector with one value per coordinate row.

## Examples

``` r
coords <- as.matrix(datasets::quakes[1:30, c("long", "lat")])
fold <- spatial_folds(coords, n_folds = 3, seed = 42)
table(fold)
```
