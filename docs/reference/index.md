# Package index

## Select an ordinary GWR model

Joint subset, kernel and bandwidth selection on a common sample.

- [`gwr_fit()`](https://byuzbasi.github.io/gwrs/reference/gwr_fit.md) :
  Fit an unpenalized geographically weighted regression
- [`gwr_select()`](https://byuzbasi.github.io/gwrs/reference/gwr_select.md)
  : Joint subset, kernel and bandwidth selection for unpenalized GWR
- [`gwr_select_verify()`](https://byuzbasi.github.io/gwrs/reference/gwr_select_verify.md)
  : Verify a joint GWR selection checkpoint
- [`summary(`*`<gwrs_selection>`*`)`](https://byuzbasi.github.io/gwrs/reference/summary.gwrs_selection.md)
  : Summarize joint GWR selection

## Fit a model

Formula and matrix interfaces for the five documented estimators.

- [`gwr_lasso()`](https://byuzbasi.github.io/gwrs/reference/gwr_reference_formulas.md)
  [`gwr_ridge()`](https://byuzbasi.github.io/gwrs/reference/gwr_reference_formulas.md)
  [`gwr_en()`](https://byuzbasi.github.io/gwrs/reference/gwr_reference_formulas.md)
  : Formula interfaces for reference penalized GWR models
- [`gwr_lasso_fit()`](https://byuzbasi.github.io/gwrs/reference/gwr_reference_fits.md)
  [`gwr_ridge_fit()`](https://byuzbasi.github.io/gwrs/reference/gwr_reference_fits.md)
  [`gwr_en_fit()`](https://byuzbasi.github.io/gwrs/reference/gwr_reference_fits.md)
  : Fit geographically weighted lasso, ridge, or elastic net
- [`gwr_nonconvex()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex.md)
  [`gwr_scad()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex.md)
  [`gwr_mcp()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex.md)
  : Formula interfaces for nonconvex GWR
- [`gwr_nonconvex_fit()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_fit.md)
  [`gwr_scad_fit()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_fit.md)
  [`gwr_mcp_fit()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_fit.md)
  : Fast nonconvex geographically weighted regression

## Tune SCAD and MCP

- [`gwr_nonconvex_path()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_path.md)
  [`gwr_scad_path()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_path.md)
  [`gwr_mcp_path()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_path.md)
  : Warm-started GWR-SCAD and GWR-MCP paths
- [`cv_gwr_nonconvex()`](https://byuzbasi.github.io/gwrs/reference/cv_gwr_nonconvex.md)
  [`cv_gwr_scad()`](https://byuzbasi.github.io/gwrs/reference/cv_gwr_nonconvex.md)
  [`cv_gwr_mcp()`](https://byuzbasi.github.io/gwrs/reference/cv_gwr_nonconvex.md)
  : Spatial cross-validation for GWR-SCAD and GWR-MCP
- [`gwr_nonconvex_cv_cache()`](https://byuzbasi.github.io/gwrs/reference/gwr_nonconvex_cv_cache.md)
  : Prepare reusable training-fold inputs for SCAD/MCP cross-validation

## Predict and inspect

- [`predict(`*`<gwrs_fit>`*`)`](https://byuzbasi.github.io/gwrs/reference/predict.gwrs_fit.md)
  : Predict from a fitted geographically weighted model
- [`summary(`*`<gwrs_fit>`*`)`](https://byuzbasi.github.io/gwrs/reference/summary.gwrs_fit.md)
  : Summarize a fitted gwrs model

## Diagnostics and interpretation

Descriptive checks and explicitly exploratory reference calculations.

- [`gwr_diagnostics()`](https://byuzbasi.github.io/gwrs/reference/gwr_diagnostics.md)
  : Extract statistical diagnostics from a gwrs fit
- [`gwr_local_inference()`](https://byuzbasi.github.io/gwrs/reference/gwr_local_inference.md)
  : Exploratory local coefficient inference
- [`gwr_diagnostic_tests()`](https://byuzbasi.github.io/gwrs/reference/gwr_diagnostic_tests.md)
  : GWR-style diagnostic F-tests
- [`gwr_moran()`](https://byuzbasi.github.io/gwrs/reference/gwr_moran.md)
  [`gwr_local_moran()`](https://byuzbasi.github.io/gwrs/reference/gwr_moran.md)
  [`gwr_moran_diagnostics()`](https://byuzbasi.github.io/gwrs/reference/gwr_moran.md)
  : Sparse residual Moran diagnostics
- [`gwr_local_collinearity()`](https://byuzbasi.github.io/gwrs/reference/gwr_local_collinearity.md)
  : Fast local multicollinearity diagnostics
- [`gwr_global_collinearity()`](https://byuzbasi.github.io/gwrs/reference/gwr_global_collinearity.md)
  : Global predictor collinearity diagnostics

## Spatial inputs and numerical controls

- [`gwr_neighbors()`](https://byuzbasi.github.io/gwrs/reference/gwr_neighbors.md)
  : Build a reusable nearest-neighbor structure
- [`spatial_folds()`](https://byuzbasi.github.io/gwrs/reference/spatial_folds.md)
  : Construct reproducible spatially contiguous folds
- [`gwrs_control()`](https://byuzbasi.github.io/gwrs/reference/gwrs_control.md)
  : Control numerical fitting in \`gwrs\`

## County data

- [`acs2020_2024_counties`](https://byuzbasi.github.io/gwrs/reference/acs2020_2024_counties.md)
  : County analysis data from the 2020–2024 ACS five-year Data Profiles
- [`acs2024_counties`](https://byuzbasi.github.io/gwrs/reference/acs2024_counties.md)
  : Compatibility name for the 2020–2024 ACS county data
