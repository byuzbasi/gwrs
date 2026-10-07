args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
app <- dirname(dirname(script))
source(file.path(app, "code", "common.R"), local = TRUE)
source(file.path(app, "code", "spatial-common.R"), local = TRUE)
require_spatial_packages()

# Spatial-v1 is immutable. Its original observer compared a decimal CSV
# round-trip with an attributed RDS matrix bit for bit, then recomputed k-means
# from the rounded CSV coordinates. On this release, the maximum round-trip
# error is about 5e-9 metres. The partition is unchanged, but k-means assigns
# different integer labels to the same five clusters. This wrapper verifies
# those facts explicitly and evaluates the original observer with only its two
# observer-side comparisons corrected in memory.

release_root <- file.path(app, "spatial-v1")
release <- validate_spatial_release(app, release_root)
if (nrow(release$outputs) != 12L || nrow(release$inputs) != 19L) {
  abort("Spatial-v1 has an unexpected input/output manifest size.")
}

model_data <- utils::read.csv(
  file.path(release_root, "data", "county-model-data-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = c(
    fips = "character", state_fips = "character",
    county_fips = "character"
  )
)
bundle <- readRDS(file.path(
  release_root, "data", "spatial-design-v1.rds"
))
folds <- utils::read.csv(
  file.path(release_root, "tables", "spatial-folds-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = c(fips = "character", state_fips = "character")
)
config <- read_json_file(file.path(app, "config", "spatial-design-v1.json"))

csv_coordinates <- as.matrix(
  model_data[, c("easting_m", "northing_m"), drop = FALSE]
)
storage.mode(csv_coordinates) <- "double"
rds_coordinates <- bundle$coordinates
if (!identical(dim(csv_coordinates), dim(rds_coordinates)) ||
    any(!is.finite(csv_coordinates)) || any(!is.finite(rds_coordinates))) {
  abort("CSV and RDS coordinate matrices have incompatible dimensions.")
}
coordinate_tolerance_m <- 1e-6
maximum_coordinate_delta_m <- max(abs(csv_coordinates - rds_coordinates))
if (!is.finite(maximum_coordinate_delta_m) ||
    maximum_coordinate_delta_m > coordinate_tolerance_m) {
  abort("CSV coordinate round-trip exceeds the one-micrometre gate.")
}

n_folds <- as.integer(config$spatial_folds$n_folds)
saved_fold <- as.integer(folds$fold)
rds_fold <- gwrs::spatial_folds(
  rds_coordinates,
  n_folds = n_folds,
  seed = as.integer(config$spatial_folds$seed),
  nstart = as.integer(config$spatial_folds$nstart)
)
if (!identical(saved_fold, rds_fold) ||
    !identical(saved_fold, as.integer(bundle$folds))) {
  abort("Authoritative RDS coordinates do not reproduce saved fold labels.")
}

csv_fold <- gwrs::spatial_folds(
  csv_coordinates,
  n_folds = n_folds,
  seed = as.integer(config$spatial_folds$seed),
  nstart = as.integer(config$spatial_folds$nstart)
)
label_map <- vapply(seq_len(n_folds), function(index) {
  labels <- unique(csv_fold[saved_fold == index])
  if (length(labels) != 1L) {
    abort("CSV round-trip changes spatial cluster membership.")
  }
  as.integer(labels)
}, integer(1L))
if (anyDuplicated(label_map) ||
    !identical(as.integer(label_map[saved_fold]), csv_fold)) {
  abort("CSV and RDS fold partitions are not equivalent up to label names.")
}

validator_path <- file.path(app, "code", "04_validate_spatial_design.R")
validator_lines <- readLines(validator_path, warn = FALSE)
if (length(validator_lines) < 6L) abort("Original observer is incomplete.")
validator_text <- paste(validator_lines[-seq_len(5L)], collapse = "\n")

replace_once <- function(text, old, new, label) {
  positions <- gregexpr(old, text, fixed = TRUE)[[1L]]
  if (length(positions) != 1L || positions[[1L]] <= 0L) {
    abort("Original observer no longer has the expected ", label, " block.")
  }
  sub(old, new, text, fixed = TRUE)
}

old_coordinate_check <- paste0(
  "!isTRUE(all.equal(\n",
  "      unname(coordinates), unname(bundle$coordinates), tolerance = 0\n",
  "    ))"
)
new_coordinate_check <- paste0(
  "!(identical(dim(coordinates), dim(bundle$coordinates)) &&\n",
  "      max(abs(coordinates - bundle$coordinates)) <= 1e-6)"
)
old_fold_source <- paste0(
  "recomputed_fold <- gwrs::spatial_folds(\n",
  "  coordinates,\n"
)
new_fold_source <- paste0(
  "recomputed_fold <- gwrs::spatial_folds(\n",
  "  bundle$coordinates,\n"
)
validator_text <- replace_once(
  validator_text, old_coordinate_check, new_coordinate_check,
  "coordinate-comparison"
)
validator_text <- replace_once(
  validator_text, old_fold_source, new_fold_source,
  "fold-coordinate-source"
)

validation_environment <- new.env(parent = globalenv())
validation_environment$app <- app
eval(parse(text = validator_text), envir = validation_environment)

message(
  "Observer-v2 diagnostics: maximum CSV/RDS coordinate delta=",
  format(maximum_coordinate_delta_m, scientific = TRUE, digits = 8L),
  " m; CSV label map=",
  paste(seq_len(n_folds), label_map, sep = "->", collapse = ","), "."
)
message(
  "Spatial-v1 observer-v2 verification: OK. Frozen release unchanged; ",
  "no model fitted."
)
