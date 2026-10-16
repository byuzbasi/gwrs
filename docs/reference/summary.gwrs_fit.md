# Summarize a fitted gwrs model

Summarize a fitted gwrs model

## Usage

``` r
# S3 method for class 'gwrs_fit'
summary(object, ...)
```

## Arguments

- object:

  A fitted \`gwrs_fit\` object.

- ...:

  Unused.

## Value

An object of class \`summary.gwrs_fit\`.

## Examples

``` r
quake <- datasets::quakes[1:32, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
fit <- gwr_lasso_fit(
  x, quake$mag, coords = coords, k = 20, lambda = 0.05,
  control = gwrs_control(n_threads = 1)
)
summary(fit)
```
