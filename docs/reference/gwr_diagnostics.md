# Extract statistical diagnostics from a gwrs fit

Statistical diagnostics are computed without constructing a dense
spatial weight or smoother matrix. For sparse penalties,
effective-degrees-of- freedom quantities are conditional on the fitted
active sets.

## Usage

``` r
gwr_diagnostics(object, ...)
```

## Arguments

- object:

  A fitted \`gwrs_fit\` object.

- ...:

  Unused.

## Value

A \`gwrs_diagnostics\` object.

## Examples

``` r
quake <- datasets::quakes[1:32, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
fit <- gwr_fit(
  x, quake$mag, coords = coords, k = 20,
  control = gwrs_control(n_threads = 1)
)
diagnostics <- gwr_diagnostics(fit)
diagnostics$global[c("rss", "r_squared", "aicc", "gcv")]
head(diagnostics$local)
```
