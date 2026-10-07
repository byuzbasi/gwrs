# Spatial protocol v1

This stage fixes the geometry, distance coordinates, common diagnostic graph,
and spatial fold assignment for the 2024 ACS county application. It does not
fit or tune any statistical model.

## Geometry and distance coordinates

The analysis domain remains the contiguous 48 states and the District of
Columbia. County records are linked to the official 2024 Census cartographic
boundaries by the five-character county FIPS key. The linkage must be one to
one and exhaustive for all 3,109 domain counties; only the 3,107 counties that
passed `audit-v2` enter the analysis geometry.

All geometries are made valid and transformed from NAD83 geographic
coordinates to NAD83 / Conus Albers (EPSG:5070). Distances are therefore
Euclidean metres in a continental equal-area projection. One representative
point per eligible county is generated with `sf::st_point_on_surface()` after
projection. This keeps each modelling point inside its county polygon and
avoids offshore centroids for irregular coastal counties.

## Frozen spatial folds

Five spatial folds are generated once with `gwrs::spatial_folds()` from the
standardized EPSG:5070 representative-point coordinates. The frozen seed is
20260907 and `nstart` is 100. These folds are regional blocks, not independent
random folds. Later model comparisons must reuse the saved assignments and
must learn every transformation and tuning parameter from training data only.
The current stage only checks determinism, coverage, minimum fold size, and
the prespecified maximum fold-size ratio.

## Common Moran and LISA graph

The base graph is binary queen contiguity among eligible county polygons with
a one-metre numerical snap tolerance. Filtering to the analysis sample leaves
Nantucket County, Massachusetts, and San Juan County, Washington, without a
polygon neighbour. Each isolate is connected symmetrically to its nearest
non-isolated eligible county using the frozen projected representative-point
distance. Every repair is written to a separate table.

The saved directed edge table contains both binary weights and row-standardized
weights. Later global Moran and local Moran (LISA) calculations must use this
same graph for the out-of-fold residuals of every compared method. This avoids
confounding residual diagnostics with method-specific GWR bandwidth graphs.
Local p-values and any multiplicity adjustment will be fixed at the modelling
protocol gate; no inferential calculation occurs here.

## Diagnostic boundary

F1, F2, and coefficient-wise F3 tests are reserved for the unpenalized fixed
Gaussian GWR fit. They will not be reported as classical post-selection tests
for SCAD, MCP, Lasso, Elastic Net, or Ridge. Selection maps will later use grey
for an exactly stored zero coefficient, separately from inferential maps.

## Reproducibility and safety

The complete stage is assembled under a temporary sibling directory and moved
atomically to `spatial-v1` only after all validation gates pass. Inputs,
outputs, software versions, external spatial-library versions, file sizes, and
SHA-256 checksums are recorded. An existing valid release is verified and
resumed read-only; an invalid existing release is never overwritten. Raw
downloads, `audit-v1`, `audit-v2`, and all earlier study outputs remain
unchanged. Model fitting stays locked after this stage.

The model-data CSV and RDS retain the complete estimate, MOE, and controlled-
MOE fields. The portable county GeoPackage intentionally contains identifiers,
coordinates, folds, the response, and the 26 predictor estimates only. Later
diagnostic and selection results will be joined by FIPS rather than expanding
the GeoPackage to a driver-sensitive number of fields.
