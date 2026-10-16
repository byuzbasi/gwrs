# Exploratory local coefficient inference

Computes local standard errors, test statistics, raw p-values, adjusted
p-values, and significance flags from a fit created with
\`gwrs_control(diagnostics = "full")\`. Sparse fits default to an
unpenalized refit on each selected local active set; all reported
p-values remain exploratory because selection and overlapping
neighborhoods induce dependence.

## Usage

``` r
gwr_local_inference(
  object,
  alpha = 0.05,
  df = c("edf", "residual", "normal"),
  adjust = c("BH", "none", "holm", "bonferroni", "BY", "fdr"),
  adjust_scope = c("coefficient", "global"),
  estimator = c("auto", "post_selection", "penalized"),
  ...
)
```

## Arguments

- object:

  A Gaussian least-squares \`gwrs_fit\` with full diagnostics, or a
  \`gwrs_selection\` whose selected GWR retains full diagnostics.

- alpha:

  Significance level.

- df:

  Reference degrees-of-freedom rule.

- adjust:

  Multiplicity adjustment passed to \[stats::p.adjust()\].

- adjust_scope:

  Adjust separately by coefficient or over all local tests.

- estimator:

  Use the penalized active-set derivative, the local post-selection
  refit, or choose automatically.

- ...:

  Unused.

## Value

A \`gwrs_local_inference\` object with local estimates, standard errors,
statistics, raw \`p_value\`, adjusted \`p_adjusted\`, \`significant\`,
coefficient summaries, settings, notes and a calculation status.

## Details

The two-sided reference calculation concerns whether each local
coefficient is zero. It is distinct from the F3 spatial-constancy
diagnostic. By default, BH adjustment is applied separately to the
locations of each design-matrix coefficient, including the intercept.
Factor contrasts are separate coefficients. \`adjust_scope = "global"\`
instead adjusts over all local coefficient tests in the fitted model.
Neither scope includes tests from other candidate models.

The default \`df = "edf"\` uses effective residual degrees of freedom;
\`"residual"\` uses n minus the smoother trace. If the requested degrees
of freedom are unavailable or nonpositive, p-values are NA with a
warning. There is no automatic normal fallback; \`df = "normal"\`
requests it explicitly. Unavailable tests also have NA significance
flags and NA summary values, rather than apparent nonsignificance or
infinite minima.

All results are exploratory. BH does not establish calibration of the
underlying reference p-values or account for response-based predictor,
kernel, bandwidth or penalty selection, smoothing bias, or dependence
between overlapping neighborhoods. Selection context is retained when
available; no calibrated post-selection or FDR guarantee is claimed.

## Examples

``` r
quake <- datasets::quakes[1:36, ]
coords <- as.matrix(quake[, c("long", "lat")])
x <- as.matrix(quake[, c("depth", "stations")])
fit <- gwr_lasso_fit(
  x, quake$mag, coords = coords, k = 24, lambda = 0.05,
  control = gwrs_control(n_threads = 1, diagnostics = "full")
)
inference <- gwr_local_inference(fit, adjust = "BH")
inference$summary
head(inference$local)
```
