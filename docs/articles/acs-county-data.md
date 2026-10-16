# Explore the ACS county dataset

The package includes `acs2020_2024_counties`, the documented **2020–2024
ACS five-year county analysis table** used in the accompanying SCAD/MCP
paper. It contains **3,107 counties and 91 columns**. The table
preserves the study values, column names and row order.

This is a single cross-sectional table with one row per county. The 2024
ACS five-year product contains Census Bureau period estimates based on
survey data collected during 2020–2024. The analysis uses these
published estimates directly: separate annual observations are neither
stacked nor arithmetically averaged. The 2024 suffix in the retained
income column names identifies the ACS product and the
inflation-adjustment year, not a single-year observation. The legacy
dataset name `acs2024_counties` remains available for compatibility.

## Load and inspect

``` r

library(gwrs)
data(acs2020_2024_counties)
dim(acs2020_2024_counties)
head(acs2020_2024_counties[, c("fips", "county_name",
                         "median_household_income_2024")])
help("acs2020_2024_counties", package = "gwrs")
```

The 91 fields contain geographic identifiers, the original and log
income response, 26 predictors, Census 90% margins of error,
controlled-MOE flags, projected coordinates and the saved five-fold
spatial assignment. The [complete variable
reference](https://byuzbasi.github.io/gwrs/reference/acs2020_2024_counties.md)
describes each column, unit and transformation.

## Coordinates and uncertainty

`easting_m` and `northing_m` are projected coordinates in metres under
NAD83 / Conus Albers (EPSG:5070). They are not longitude and latitude.
Retain the same coordinate units when fitting and predicting.

The response `log_median_household_income_2024` is the natural logarithm
of the income estimate. Its companion income margin of error remains in
**dollars**, even though the column name begins with the log-response
name. Margins of error for percentage predictors are in percentage
points. Controlled Census MOE sentinels are represented by zero and
their companion flags; these flags do not mean that the source estimate
has no uncertainty.

The documented analysis sample uses the study exclusions. Do not treat
it as an unfiltered table of every US county or silently replace its
observations when reproducing the application.

## Data dictionary and attribution

``` r

system.file("extdata", package = "gwrs")
citation("gwrs")
```

Please acknowledge the **U.S. Census Bureau**, the accompanying paper
and `gwrs` when using this prepared table. The [source and
reconstruction
guide](https://github.com/byuzbasi/gwrs/blob/main/reproducibility/scad-mcp/SOURCES.md)
records the Census sources, exclusions and processing steps. The
downloadable reproduction bundle includes the public input snapshot and
reconstruction code.
