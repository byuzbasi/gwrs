# Single cross-sectional Census period-estimate table for 2020-2024.
# No stacking of annual observations or averaging of annual estimates.
# Run from the repository root. This converts the reviewed CSV without changes.
source_path <- "reproducibility/scad-mcp/real-data/study-v1/frozen/input/spatial/data/county-model-data-v1.csv"
provenance <- jsonlite::read_json("inst/extdata/acs2020_2024_counties_provenance.json")
stopifnot(digest::digest(file = source_path, algo = "sha256") == provenance$input_sha256)
acs2020_2024_counties <- read.csv(source_path, stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = c(fips = "character", state_fips = "character", county_fips = "character"))
stopifnot(identical(dim(acs2020_2024_counties), c(3107L, 91L)),
          !anyNA(acs2020_2024_counties), !anyDuplicated(acs2020_2024_counties$fips))
for (name in names(acs2020_2024_counties)[vapply(acs2020_2024_counties, is.character, logical(1))]) {
  Encoding(acs2020_2024_counties[[name]]) <- "UTF-8"
}
dir.create("data", showWarnings = FALSE)
path <- "data/acs2020_2024_counties.rda"
if (file.exists(path)) {
  previous <- new.env(parent = emptyenv())
  load(path, envir = previous)
  stopifnot(identical(previous$acs2020_2024_counties, acs2020_2024_counties))
}
save(acs2020_2024_counties, file = path, compress = "xz", version = 2)
restored <- new.env(parent = emptyenv())
load(path, envir = restored)
stopifnot(identical(restored$acs2020_2024_counties, acs2020_2024_counties))
cat("ACS_DATASET_EXACT_EQUALITY_PASS: 3107 x 91; no data values changed\n")
