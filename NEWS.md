# gwrs 0.4.0.9015

- Adds selected-model local coefficient inference to `gwr_select()` by
  default, with raw two-sided reference p-values and BH adjustment separately
  across locations for each coefficient. Set `local_inference = FALSE` to
  skip it. Search candidates retain standard diagnostics; only the selected
  GWR final refit uses full diagnostics. Model ranking is unchanged.
- `gwr_local_inference()` also accepts a selection result, retains selection
  context, and explicitly labels reference calibration as unestablished.
  BH does not account for response-based model/tuning selection.
- Invalid reference degrees of freedom now return NA p-values with a warning
  instead of silently switching to normal references. An explicit
  `df = "normal"` remains available. Missing tests have NA significance flags
  and finite-or-NA summaries, with available/unavailable counts.
- Clarifies that local coefficient-zero tests differ from F3 spatial
  constancy, and documents per-coefficient versus global adjustment families.

# gwrs 0.4.0.9014

- Adds `gwr_select()` for joint predictor-subset, kernel and bandwidth selection
  in ordinary Gaussian GWR. AICc is the default; AIC, BIC, GCV and conditional
  leave-one-out error are also reported on the same observations.
- Supports exhaustive or deterministic stepwise search, common user-supplied
  Moran weights, design inspection, bounded execution and verified resume.
- Improves candidate-grid construction and contiguous C++ F-moment access.
  Numerical moment definitions are preserved. Small benchmarks measured about
  2.1--2.4 times faster F diagnostics at n=320; this is not a whole-search claim.
- Labels all F references exploratory and uncalibrated. F1 and F2 transform
  the same RSS ratio and do not provide independent confirmations.
- Adds an executable joint-selection guide and corresponding website reference
  pages, retaining the ACS period-name alias and historical reproduction engines.

# gwrs 0.4.0.9013

- Adds `acs2020_2024_counties` as the primary name for the unchanged 2020–2024 ACS county table, with period-estimate documentation; `acs2024_counties` remains available for compatibility.

- Focused documentation and examples for GWR-Ridge, GWR-Lasso, GWR-ENET, GWR-SCAD and GWR-MCP.
- Retains the validated coordinate cache, regularization paths and spatial prediction implementation.
- Documentation packaging does not change numerical functions, defaults, tuning or stopping rules.
