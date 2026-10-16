# Formula interfaces for reference penalized GWR models

These are concise formula wrappers around \[gwr_penalized()\].

## Usage

``` r
gwr_lasso(formula, data, coords, ...)

gwr_ridge(formula, data, coords, ...)

gwr_en(formula, data, coords, ..., alpha = 0.5)
```

## Arguments

- formula:

  Model formula with an intercept.

- data:

  Data frame containing model variables.

- coords:

  Optional coordinate matrix used to construct neighbors.

- ...:

  Additional arguments passed to \[gwr_penalized()\].

- alpha:

  Elastic-net mixing parameter in \`\[0, 1\]\`.

## Value

An object inheriting from \`gwrs_fit\`.

## Examples

``` r
quake <- datasets::quakes[1:30, ]
control <- gwrs_control(n_threads = 1)
lasso <- gwr_lasso(
  mag ~ depth + stations, quake, c("long", "lat"),
  k = 18, lambda = 0.05, control = control
)
ridge <- gwr_ridge(
  mag ~ depth + stations, quake, c("long", "lat"),
  k = 18, lambda = 0.05, control = control
)
elastic_net <- gwr_en(
  mag ~ depth + stations, quake, c("long", "lat"),
  k = 18, lambda = 0.05, alpha = 0.6, control = control
)
c(lasso$method, ridge$method, elastic_net$method)
```
