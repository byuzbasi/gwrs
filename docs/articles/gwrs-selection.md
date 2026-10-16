# Joint model selection for ordinary GWR

[`gwr_select()`](https://byuzbasi.github.io/gwrs/reference/gwr_select.md)
accepts a response and candidate predictor scope, then compares ordinary
Gaussian GWR fits across predictor subsets, kernels and bandwidths. Each
candidate is unpenalized. The default ranking criterion is AICc; AIC,
BIC, GCV and conditional leave-one-out squared error are reported
alongside it. The selected GWR is the best eligible candidate under the
chosen criterion within the evaluated search space.

## A small reproducible demonstration

The following simulated observations have unitless planar coordinates.
They illustrate the interface and are not empirical evidence about a
real region. Real applications should supply suitable projected
coordinates; the package uses Euclidean distance and does not transform
coordinates.

``` r

set.seed(1401)
dat <- data.frame(east = runif(28), north = runif(28),
                  x1 = rnorm(28), x2 = rnorm(28))
dat$y <- 1 + (1 + dat$east) * dat$x1 + rnorm(28, sd = 0.5)
selection <- gwr_select(y ~ x1 + x2, dat, c("east", "north"),
  kernels = c("bisquare", "gaussian"), k = c(20, 28),
  control = gwrs_control(n_threads = 1), verbose = FALSE)
selection
#> Joint unpenalized GWR selection
#>   search:exhaustive / complete
#>   evaluated: 16   criterion: aicc 
#>   best GWR: x1 
#>   kernel: gaussian   bandwidth: adaptive 28 
#>   criterion: 39.99189 
#>   Moran: not_calculated: weights_not_supplied 
#>   local coefficient inference: calculated 
#>   flags: grid_boundary 
#>   best including OLS: OLS - x1
knitr::kable(selection$best_by_criterion[, c("criterion", "comparison",
  "terms", "kernel", "k", "criterion_value")], digits = 3)
```

|     | criterion | comparison | terms | kernel   |   k | criterion_value |
|:----|:----------|:-----------|:------|:---------|----:|----------------:|
| 6   | aic       | GWR        | x1    | gaussian |  28 |          38.780 |
| 61  | aic       | overall    | x1    | gaussian |  28 |          38.780 |
| 62  | aicc      | GWR        | x1    | gaussian |  28 |          39.992 |
| 2   | aicc      | overall    | x1    | NA       |  NA |          39.860 |
| 63  | bic       | GWR        | x1    | gaussian |  28 |          43.204 |
| 21  | bic       | overall    | x1    | NA       |  NA |          42.857 |
| 64  | gcv       | GWR        | x1    | gaussian |  28 |           0.219 |
| 65  | gcv       | overall    | x1    | gaussian |  28 |           0.219 |
| 66  | loo_cv    | GWR        | x1    | gaussian |  28 |           6.192 |
| 22  | loo_cv    | overall    | x1    | NA       |  NA |           6.186 |

`best_model` is a fitted GWR object supporting
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html); `best_overall` also
allows OLS to win. The corresponding metadata are `selected` and
`selected_overall`. `ranking` retains every evaluated candidate,
including failed or invalid fits and their reasons. `summary(selection)`
shows leading candidates and the selected model’s conditional F1/F2/F3
diagnostics.

## Local coefficient p-values and BH adjustment

By default, only the selected GWR is refitted with full diagnostics.
Candidate search fits retain standard diagnostics, so local inference is
not repeated for every subset/kernel/bandwidth combination. The table
below contains two-sided coefficient-zero reference p-values and
BH-adjusted p-values. These differ from F3, which concerns spatial
constancy of a whole surface.

``` r

local <- selection$local_inference
local$status
#> [1] "calculated"
display <- head(local$local[, c("location", "coefficient", "estimate",
  "standard_error", "statistic", "p_value", "p_adjusted")])
display$p_value <- format.pval(display$p_value, digits = 3)
display$p_adjusted <- format.pval(display$p_adjusted, digits = 3)
knitr::kable(display, digits = 4)
```

| location | coefficient | estimate | standard_error | statistic | p_value  | p_adjusted |
|---------:|:------------|---------:|---------------:|----------:|:---------|:-----------|
|        1 | (Intercept) |   1.0735 |         0.0861 |   12.4712 | 2.50e-12 | 3.28e-12   |
|        2 | (Intercept) |   1.0875 |         0.0861 |   12.6253 | 1.91e-12 | 3.28e-12   |
|        3 | (Intercept) |   1.0905 |         0.0860 |   12.6745 | 1.76e-12 | 3.28e-12   |
|        4 | (Intercept) |   1.0750 |         0.0859 |   12.5079 | 2.35e-12 | 3.28e-12   |
|        5 | (Intercept) |   1.0587 |         0.0859 |   12.3274 | 3.23e-12 | 3.28e-12   |
|        6 | (Intercept) |   1.0590 |         0.0859 |   12.3344 | 3.19e-12 | 3.28e-12   |

``` r

local$settings[c("adjust", "adjust_scope", "exploratory",
                 "reference_calibrated", "selection_conditioned")]
#> $adjust
#> [1] "BH"
#> 
#> $adjust_scope
#> [1] "coefficient"
#> 
#> $exploratory
#> [1] TRUE
#> 
#> $reference_calibrated
#> [1] FALSE
#> 
#> $selection_conditioned
#> [1] TRUE
```

P-values are formatted for display; `local$local` retains the numeric
values.

The default BH family is the locations for each design-matrix
coefficient separately, including the intercept. Factor contrasts are
separate coefficients. To adjust over every local coefficient test in
the selected model instead, call the same inference function with an
explicit scope:

``` r

all_local_tests <- gwr_local_inference(selection, adjust = "BH",
                                       adjust_scope = "global")
all_local_tests$summary
#>   coefficient selected_percent significant_percent median_abs_statistic
#> 1 (Intercept)              100                 100             12.40802
#> 2          x1              100                 100             15.30688
#>   min_adjusted_p available_tests unavailable_tests
#> 1   3.170602e-12              28                 0
#> 2   6.681090e-14              28                 0
```

All local p-values are exploratory. Neither BH scope accounts for the
response-based selection of predictors, kernel and bandwidth; BH does
not establish calibration of the underlying reference p-values or remove
smoothing bias and dependence between overlapping neighborhoods. No
post-selection or calibrated FDR guarantee is claimed. Local p-values
never affect model ranking. The selected model and whether the search
completed are retained in `local$settings$selection`.

Set `local_inference = FALSE` in
[`gwr_select()`](https://byuzbasi.github.io/gwrs/reference/gwr_select.md)
to skip full diagnostics and local inference; the result records
`local_inference_status = "disabled"`. No eligible GWR gives NULL
inference with an explicit status. Incomplete searches identify their
inference as provisional. Unavailable reference degrees of freedom or
coefficient statistics produce NA p-values and significance flags; the
function does not silently switch to a normal approximation.
`df = "normal"` is an explicit option for standalone inference.

## Define and inspect the search

Factors and multi-column formula terms are selected as whole groups.
Interactions require their main effects in the full scope and retain all
declared lower-order parents. `required` accepts term labels or a
one-sided formula. Required terms also constrain the OLS comparison;
intercept-only OLS is compared only when no terms are required.
Intercept-only GWR is not fitted.

All candidates use the same full-scope observations and design coding.
Missing or nonfinite values cause an error, with no silent sample
deletion or imputation. Standardization defaults to `FALSE`; setting it
to `TRUE` is an explicit preprocessing choice passed to the existing
fitting function.

``` r

design <- gwr_select(y ~ x1 + x2, dat, c("east", "north"), dry_run = TRUE)
design
#> GWR selection design (no fits)
#>   n: 28   scope terms: 2 
#>   search: exhaustive   candidates: 389 
#>   neighbor workspace estimate (MiB): 0.01196289 
#>   Moran: not_calculated: weights_not_supplied
```

The default exhaustive search considers every nonempty hierarchy-valid
subset and all five kernels. With adaptive bandwidth and no `k`, each
subset searches all integers from its expanded predictor count plus two
through the sample size. `k` includes the target observation. Adaptive
Gaussian and exponential fits retain only these `k` neighbors, matching
the existing package semantics.

For fixed-distance selection, set `bandwidth = "fixed"` and provide `bw`
in coordinate units. All observations are retained in the neighbor
representation; compact kernels assign zero weight outside their radius.
Without `bw`, a 30-value logarithmic grid spans the smallest positive
pairwise distance to twice the largest. A finite grid gives no
continuous-bandwidth guarantee. Both modes can be requested with
`bandwidth = c("adaptive", "fixed")`.

`max_models` rejects oversized exhaustive designs before fitting. For a
larger scope, explicitly choose `search = "stepwise"`: full-model and
required-plus-single-term seeds are followed by best-improvement
additions and deletions. Each visited subset searches the complete
requested kernel/bandwidth grid. This deterministic search has no
global-optimum guarantee.

## Diagnostics and scientific limits

Eligibility requires full local predictor rank, direct converged solves,
finite fitted quantities and positive residual effective degrees of
freedom. Each ranking criterion separately requires a finite score. Low
effective sample size, high leverage, collinearity and grid-boundary
solutions are reported as descriptive flags. Diagnostic p-values do not
affect ranking.

F3 tests spatial variability of coefficients, not whether they are zero.
F reference distributions are moment-matched approximations. All F and
Moran p-values are raw, conditional and exploratory; selection and
multiplicity have not been adjusted for. Conditional LOO holds the
fitted neighborhoods and weights fixed. It does not replace external or
nested spatial validation that repeats model selection inside training
folds.

For a fixed smoother, F1 and F2 are decreasing affine transformations of
the same residual-sum-of-squares ratio. Their separate moment-matched
reference distributions may give different p-values. These are not two
independent pieces of evidence, and exact operator moments do not
establish calibration of either reference distribution. F3 likewise has
no established calibration under response-based selection or spatially
varying nuisance coefficients.

Moran is reported as `not_calculated: weights_not_supplied` by default.
To compute it, supply a scientifically chosen common matrix as
`moran_weights`. Rows and columns must follow data order; supplied names
must match exactly. Weights must be finite, nonnegative and have zero
diagonal. Row standardization is controlled explicitly. Islands fail by
default; `moran_zero_policy = "keep"` retains them in global
centering/permutations and marks their local statistics unavailable. No
neighbor graph is chosen for you.

## Bounded execution and verified resume

A work unit is one candidate fit and its requested diagnostics. The
following template uses your already prepared `dat` and a new versioned
output directory. It is not executed when building this guide.

``` r

run_dir <- "results/gwr-selection-v01"
partial <- gwr_select(y ~ x1 + x2, dat, c("east", "north"),
  k = c(20, 28), checkpoint = run_dir, max_evaluations = 10,
  control = gwrs_control(n_threads = 1))
gwr_select_verify(run_dir)
finished <- gwr_select(y ~ x1 + x2, dat, c("east", "north"),
  k = c(20, 28), checkpoint = run_dir, resume = TRUE,
  control = gwrs_control(n_threads = 1))
gwr_select_verify(run_dir)
```

Invocation limits produce explicitly provisional results. Resume
verifies data, settings (including Moran seeds/weights and thread
controls), software, file sizes and SHA-256 hashes before reusing
immutable candidate packs. Uncommitted orphan packs are retained but
ignored. A single owned lock prevents concurrent writers. Stale locks
are never silently removed: confirm the owner is no longer running, then
archive its `run.lock` directory before resuming.

`progress.json` and append-only `progress.tsv` distinguish calculated
and durable counts. Search ETA uses measured throughput for new
candidates in the current invocation and excludes final refits/output
verification; stepwise total work and ETA remain unknown. Progress
updates occur between candidates, at most every 30 seconds, so a long
individual fit may delay the next heartbeat. Final result and ranking
snapshots are numbered without overwriting earlier snapshots. Use
[`gwr_select_verify()`](https://byuzbasi.github.io/gwrs/reference/gwr_select_verify.md)
for integrity and completeness checks; progress alone is not validation.
`memory_limit_mb` guards estimated neighbor and diagnostic workspaces,
not total process memory.
