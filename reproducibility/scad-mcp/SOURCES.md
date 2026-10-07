# Source and data provenance

The simulation uses `packages/gwrs_0.4.0.tar.gz`; the county application and
input reconstruction use `packages/gwrs_0.4.0.9006.tar.gz`. These are minimal
source distributions of the versions used in the analyses. Their R code,
native code and NAMESPACE are preserved byte-for-byte; manifests are in
`packages/numerical-source-*.json`. The source archives were repackaged to
exclude author installation logs, machine metadata and unrelated documentation.
The new archive checksums therefore differ from earlier source-archive checksums.

The current package in the repository root includes later computational work
and the documented ACS data object. It is separate from the paper's reference
engines. Do not substitute it when aiming to reproduce the reported analyses.
The portable wrappers change paths, environment isolation, progress reporting
and integrity checks. They do not change estimator functions, data-generating
processes, model grids, seeds, selection thresholds or statistical summaries.
Simulation thread count and ACS process count are explicit resource settings.

The study seeds and full simulation task plan are recorded under
`simulation/reference/`. ACS settings are under
`real-data/study-v3/config/`. Each execution records its actual software stack;
there is no claim of a universally identical dependency stack on all platforms.
The historical simulation used R 4.3.0 on Linux; local portability validation
used R 4.6.0 on macOS. See the validation report for the tested environment.

## Census inputs

The data are U.S. Census Bureau **2020--2024 ACS five-year Data Profiles**,
downloaded on 7 September 2026, with the 2024 county/state cartographic
boundaries. Original public files and download manifests are under
`data-preparation/data/raw/2026-09-07/`. Variable codes, official labels and
units are recorded in `data-preparation/config/variable-catalog-v1.csv`.

- [ACS API Data Profiles](https://api.census.gov/data/2024/acs/acs5/profile.html)
- [ACS program](https://www.census.gov/programs-surveys/acs)
- [Cartographic boundary files](https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html)

The sample contains 3,107 counties in the contiguous 48 states plus DC.
De Baca County, New Mexico, is excluded for missing income; Loving County,
Texas, for missing mean commute time. Estimates are not imputed. The income
response is log-transformed, and the 26 predictor estimates retain their units.
Controlled MOE sentinel values are represented by zero with explicit flags.
The income MOE column retains dollar units despite its log-response prefix.
MOEs are not used as model weights. No new exclusions or transformations were
introduced during packaging.

Please acknowledge the Census Bureau as the original data source, Bahadir
Yuzbasi's accompanying manuscript for the prepared study dataset, and the
`gwrs` package via `citation("gwrs")`. The paper has no assigned DOI in this
distribution; no journal acceptance or publication is implied.
