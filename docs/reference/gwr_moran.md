# Sparse residual Moran diagnostics

Computes global Moran's I and local Moran (LISA) statistics directly
from a fitted model's nearest-neighbor graph. Self-neighbor weights are
removed and no dense spatial matrix is formed. Permutations use
deterministic independent C++ random streams and can be evaluated in
parallel.

## Usage

``` r
gwr_moran(
  object,
  permutations = 999L,
  row_standardize = TRUE,
  seed = 1L,
  ...,
  weights = NULL,
  zero_policy = c("error", "keep"),
  memory_limit_mb = 512
)

gwr_local_moran(
  object,
  permutations = 999L,
  row_standardize = TRUE,
  seed = 1L,
  alpha = 0.05,
  ...,
  weights = NULL,
  zero_policy = c("error", "keep"),
  memory_limit_mb = 512
)

gwr_moran_diagnostics(
  object,
  permutations = 999L,
  row_standardize = TRUE,
  seed = 1L,
  alpha = 0.05,
  ...,
  weights = NULL,
  zero_policy = c("error", "keep"),
  memory_limit_mb = 512
)
```

## Arguments

- object:

  A fitted \`gwrs_fit\`.

- permutations:

  Number of random permutations; zero skips inference.

- row_standardize:

  Whether non-self spatial weights are row-standardized.

- seed:

  Nonnegative deterministic C++ permutation seed.

- ...:

  Unused.

- weights:

  Optional nonnegative numeric or sparse spatial weight matrix, in
  observation order, with a zero diagonal. When supplied it replaces the
  fitted kernel graph. Directed weights are allowed; no symmetrization
  occurs.

- zero_policy:

  For explicit weights, \`"error"\` rejects islands; \`"keep"\` retains
  them in global centering and permutations and reports local NA.

- memory_limit_mb:

  Workspace guard for the explicit-weight calculation.

- alpha:

  Significance level used for LISA cluster labels.

## Value

\`gwr_moran()\` returns a \`gwrs_moran\` object. \`gwr_local_moran()\`
returns a \`gwrs_local_moran\` data frame. \`gwr_moran_diagnostics()\`
returns both as \`global\` and \`local\`, reusing one C++ permutation
calculation.

## Details

The existing permutation null randomly relabels the entire centered
residual vector (total randomization); it is not conditional LISA with
the focal residual fixed. Two-sided global exceedances are centered at
\`-1/(n-1)\`; local exceedances are centered at each location's
simulated mean. Both use the plus-one permutation correction. LISA
cluster labels use raw, unadjusted p-values. Regression residuals need
not be exchangeable, so these permutation diagnostics are exploratory,
not exact residual-model tests. No fitted model or tuning parameters are
re-estimated during permutations.

## References

Anselin, L. (1995). Local Indicators of Spatial Association–LISA.
Geographical Analysis, 27, 93-115.
[doi:10.1111/j.1538-4632.1995.tb00338.x](https://doi.org/10.1111/j.1538-4632.1995.tb00338.x)
.

## Examples

``` r
quake <- datasets::quakes[1:36, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
fit <- gwr_fit(
  x, quake$mag, coords = coords, k = 24,
  control = gwrs_control(n_threads = 1)
)
global <- gwr_moran(fit, permutations = 19, seed = 42)
local <- gwr_local_moran(fit, permutations = 19, seed = 42)
global
head(local)
plot(local, significant_only = TRUE)
```
