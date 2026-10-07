args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(Sys.getenv("GWRS_REPRO_ROOT"), mustWork = TRUE)
prepared <- normalizePath(args[[1L]], mustWork = TRUE)
reference <- file.path(root, "real-data/study-v1/frozen/input/spatial")
model_path <- "data/county-model-data-v1.csv"
new <- read.csv(file.path(prepared, "spatial-v1", model_path), check.names = FALSE,
                colClasses = c(fips = "character", state_fips = "character", county_fips = "character"))
old <- read.csv(file.path(reference, model_path), check.names = FALSE,
                colClasses = c(fips = "character", state_fips = "character", county_fips = "character"))
stopifnot(identical(names(new), names(old)), identical(new$fips, old$fips), nrow(new) == 3107L)
for (name in names(old)) {
  if (is.numeric(old[[name]])) {
    stopifnot(isTRUE(all.equal(new[[name]], old[[name]], tolerance = 1e-12)))
  } else stopifnot(identical(new[[name]], old[[name]]))
}
for (name in c("spatial-folds-v1.csv", "common-neighbor-graph-directed-v1.csv")) {
  new <- read.csv(file.path(prepared, "spatial-v1/tables", name), colClasses = "character")
  old <- read.csv(file.path(reference, "tables", name), colClasses = "character")
  stopifnot(identical(new, old))
}
cat("ACS_INPUT_RECONSTRUCTION_PASS: 3107 rows, all fields, frozen folds and common graph\n")
