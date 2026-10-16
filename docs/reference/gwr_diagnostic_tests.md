# GWR-style diagnostic F-tests

Computes F1 (local model versus global OLS), F2 (joint coefficient
stationarity), and coefficient-specific F3 spatial-variability
diagnostics. Exact first and second operator moments are computed in
neighbor blocks. The F reference distributions are moment-matched
approximations for a fixed, identified Gaussian GWR smoother, not exact
finite-sample distributions. All reference p-values are exploratory and
have no established calibration. Penalized fits additionally require
\`allow_approximate = TRUE\`.

## Usage

``` r
gwr_diagnostic_tests(
  object,
  allow_approximate = FALSE,
  ...,
  block_size = 64L,
  memory_limit_mb = 512
)
```

## Arguments

- object:

  A least-squares \`gwrs_fit\` retaining \`x\` and \`y\`.

- allow_approximate:

  Permit conditional active-set tests for a penalized fit.

- ...:

  Unused.

- block_size:

  Positive integer number of operator columns per block. Smaller blocks
  reduce workspace without changing the target moments.

- memory_limit_mb:

  Positive finite workspace budget in MiB. A conservative estimate is
  checked before operator allocation. This is not a hard process RSS
  limit; the R session, allocator and numerical libraries use other
  memory. No stochastic or diagonal-only fallback is used when the
  budget is exceeded.

## Value

A \`gwrs_diagnostic_tests\` object. \`diagnostics\` includes full
operator moments, coefficient moments and workspace settings.
Unidentified global OLS references or degenerate variability produce
unavailable (\`NA\`) tests. A singular/non-positive local Hessian or a
nonconvex penalty knot raises an informative error instead of silently
substituting a pseudoinverse.

## Details

F3 tests spatial variation of an entire coefficient surface; it is not a
pointwise test of a zero coefficient. Exact moments do not correct
selection bias, response-dependent tuning or residual spatial
dependence. Selection maps should label gray as "not selected", with
inferential decisions displayed separately. Exact all-coefficient second
moments can be expensive even though no dense spatial smoother is
stored.

For fixed weights, let nu be the OLS residual degrees of freedom, delta1
the local residual trace, and eta1 = nu - delta1. When the denominators
are positive, F2 = nu/eta1 - (delta1/eta1) \* F1. Thus F1 and F2 are
decreasing affine transforms of one RSS ratio. Separate moment-matched F
references can give different p-values; these are not independent
evidence. Numerator/denominator quadratic forms are generally dependent,
so exact operator moments alone do not establish an F law, correct test
size, or valid inference after model selection. F3 calibration with all
coefficients constant would not by itself justify testing one constant
coefficient with spatially varying nuisance effects.

## References

Leung, Y., Mei, C.-L. and Zhang, W.-X. (2000). Statistical Tests for
Spatial Nonstationarity Based on the Geographically Weighted Regression
Model. Environment and Planning A, 32, 9-32.
[doi:10.1068/a3162](https://doi.org/10.1068/a3162) .

## Examples

``` r
quake <- datasets::quakes[1:36, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
fit <- gwr_fit(
  x, quake$mag, coords = coords, k = 24,
  control = gwrs_control(n_threads = 1)
)
tests <- gwr_diagnostic_tests(fit)
tests$global_tests
tests$coefficient_tests
```
