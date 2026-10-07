# Protocol v1: large-n, moderate-p USA county application

## Research question

Which social, demographic, labour-market, occupational, industrial,
transportation, and housing characteristics have location-specific predictive
associations with county median household income, and do local SCAD and MCP
provide more accurate and more stable spatially varying supports than local
Lasso, Elastic Net, Ridge, and unpenalized GWR?

This is a predictive and descriptive question. In a real data set the true
coefficient surfaces are not observed. Consequently, the application will not
claim that heterogeneity is known in advance or that fitted associations are
causal. Evidence for heterogeneity must come from held-out spatial prediction,
coefficient/support maps, stability diagnostics, and valid nonstationarity
tests.

## Population and spatial domain

The source population is all county and county-equivalent records in the 2024
ACS 5-year Data Profiles. The intended analysis domain is the contiguous 48
states and the District of Columbia. Alaska, Hawaii, Puerto Rico, and the other
territories are excluded by the state FIPS codes fixed in
`config/study-v1.json`. This matches the earlier county display domain while
keeping the new raw snapshot and all results separate.

County geometry will come from the U.S. Census Bureau's 2024 1:5,000,000
cartographic boundary files. Geometry processing, coordinate representation,
projection, and spatial folds are deliberately deferred until the acquisition
and sample audit is accepted.

## Outcome and predictors

The response is the natural logarithm of 2024 median household income in
inflation-adjusted U.S. dollars (`DP03_0062E`). The logarithm is prespecified to
place proportional rather than absolute income errors on a comparable scale
and to reduce the influence of the long upper tail. Counties with a missing,
sentinel-coded, zero, or negative response will be ineligible; they will not be
imputed.

Exactly 26 candidate predictors were chosen from ACS Data Profiles before
examining any fitted result. Percentages are used in their published percentage
units; mean commute time remains in minutes. No standardization is performed
during acquisition. If modelling proceeds, centring and scaling will be learned
from each training fold only.

The catalogue intentionally omits one or more categories from occupational and
other compositional families. This prevents a set of shares that sums exactly
to 100 from becoming linearly dependent with the intercept. Correlation and
local conditioning will still be diagnosed rather than assumed away.

## Leakage policy

Variables are excluded before modelling when they:

1. are another direct summary of the household-income distribution;
2. are constructed from household income, including poverty thresholds and
   housing-cost-to-income ratios;
3. represent mean values of income components; or
4. complete an exact percentage composition already represented by the kept
   categories.

These exclusions are scientific design rules, not data-driven screening. The
full decision trail is retained in `config/leakage-audit-v1.csv`.

## ACS estimates and margins of error

For every response or predictor estimate, the corresponding ACS 90% margin of
error (MOE) is downloaded. Estimate availability determines complete-case
eligibility. Missing MOE values are reported separately and do not silently
remove a county. No reliability threshold is used in v1.

ACS numeric special values are interpreted according to the Census Bureau
policy frozen in `config/acs-special-values-v1.csv`. In particular,
`-555555555` in an MOE field denotes an estimate controlled to an independent
population or housing estimate. Consistent with Census guidance, audit v2
stores its MOE as zero and also writes a row-level `*_moe90_controlled` flag so
that a controlled MOE cannot be confused with an ordinary reported zero. All
other documented negative special values and JSON nulls remain missing. This
correction changes neither estimate availability nor sample eligibility.

The audit reports absolute MOE distributions and, for nonzero estimates, the
ratio of MOE to the absolute estimate. That ratio is descriptive only: it is
not a model weight, exclusion rule, or transformation. Any later use of MOEs
in estimation would change the method and would require a separately approved
protocol.

## Planned comparison after the audit gate

The primary methods will be local SCAD and local MCP. Planned comparators are
local Lasso, Elastic Net, Ridge, unpenalized GWR, and their appropriate global
counterparts. All methods must use the same eligible counties, frozen outer
spatial folds, and training-only preprocessing. Hyperparameters must be chosen
inside the corresponding outer-training samples.

F1, F2, and coefficient-wise F3 diagnostics will be computed only for the
unpenalized fixed GWR fit when their assumptions and residual degrees of
freedom are valid. Moran and LISA diagnostics will use a prespecified common
graph and out-of-fold residuals. On selection maps, grey will mean an exactly
stored zero coefficient; it will not mean an insignificant coefficient.

No bandwidth grid, folds, penalty grid, uncertainty threshold, neighbourhood
graph, or production resource request is fixed in this acquisition version.
Those decisions require the completed sample and geometry audit.

## Gates

1. Metadata gate: all 54 estimate/MOE codes (one response and 26 predictors)
   must match the official 2024 metadata, and the predictor count must be 26.
2. Acquisition gate: all raw files must be nonempty, parseable, checksummed,
   and recorded without exposing the API key.
3. Sample gate: county keys must be unique; all exclusions and missing values
   must be enumerated; controlled MOEs must remain identifiable; no imputation
   is allowed.
4. Review gate: modelling remains locked until the sample and MOE report is
   reviewed and the geometry/fold protocol is explicitly approved.
