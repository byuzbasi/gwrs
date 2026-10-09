# Compatibility entry point; existing legacy data are preserved.
source("data-raw/acs2020_2024_counties.R")
legacy <- new.env(parent = emptyenv())
load("data/acs2024_counties.rda", envir = legacy)
stopifnot(identical(legacy$acs2024_counties, acs2020_2024_counties))
