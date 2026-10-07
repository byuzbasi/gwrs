args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
app <- dirname(dirname(script))
source(file.path(app, "code", "common.R"), local = TRUE)
source(file.path(app, "code", "spatial-common.R"), local = TRUE)
require_spatial_packages()

config <- read_json_file(file.path(app, "config", "spatial-design-v1.json"))
release_root <- file.path(app, config$output_root)
release <- validate_spatial_release(app, release_root)
completion <- release$completion
if (nrow(release$outputs) != 12L || nrow(release$inputs) != 19L) {
  abort("Spatial release has an unexpected input/output manifest size.")
}

data_dir <- file.path(release_root, "data")
tables_dir <- file.path(release_root, "tables")
model_data <- utils::read.csv(
  file.path(data_dir, "county-model-data-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = c(
    fips = "character", state_fips = "character",
    county_fips = "character"
  )
)
bundle <- readRDS(file.path(data_dir, "spatial-design-v1.rds"))
counties <- sf::st_read(
  file.path(data_dir, "analysis-counties-v1.gpkg"), quiet = TRUE
)
states <- sf::st_read(
  file.path(data_dir, "conus-states-v1.gpkg"), quiet = TRUE
)
folds <- utils::read.csv(
  file.path(tables_dir, "spatial-folds-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = c(fips = "character", state_fips = "character")
)
fold_summary <- utils::read.csv(
  file.path(tables_dir, "fold-summary-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
edges <- utils::read.csv(
  file.path(tables_dir, "common-neighbor-graph-directed-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = c(from_fips = "character", to_fips = "character")
)
repairs <- utils::read.csv(
  file.path(tables_dir, "isolate-repairs-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = c(isolate_fips = "character", attached_fips = "character")
)
graph_summary <- utils::read.csv(
  file.path(tables_dir, "graph-summary-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
geometry_audit <- utils::read.csv(
  file.path(tables_dir, "geometry-audit-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
catalog <- utils::read.csv(
  file.path(app, "config", "variable-catalog-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

expected_n <- as.integer(config$expected_eligible_counties)
epsg <- as.integer(config$projection$epsg)
allowed_bundle_names <- c(
  "schema", "model_data", "coordinates", "folds", "analysis_counties",
  "states", "common_graph", "metadata"
)
if (bundle$schema != "gwrs-usa-counties-acs2024-spatial-bundle-v1" ||
    !setequal(names(bundle), allowed_bundle_names) ||
    !isTRUE(bundle$metadata$model_fitting_locked) ||
    completion$model_fits_created != 0L ||
    !isTRUE(completion$model_fitting_locked)) {
  abort("Spatial bundle schema or model-fitting lock is invalid.")
}

if (nrow(model_data) != expected_n || nrow(counties) != expected_n ||
    nrow(folds) != expected_n || nrow(bundle$model_data) != expected_n ||
    nrow(states) != 49L || anyDuplicated(model_data$fips) ||
    !identical(model_data$fips, counties$fips) ||
    !identical(model_data$fips, folds$fips) ||
    !identical(model_data$fips, bundle$model_data$fips)) {
  abort("Spatial release dimensions, keys, or row order are inconsistent.")
}
if (sf::st_crs(counties)$epsg != epsg || sf::st_crs(states)$epsg != epsg ||
    any(!sf::st_is_valid(counties)) || any(sf::st_is_empty(counties)) ||
    any(!sf::st_is_valid(states)) || any(sf::st_is_empty(states))) {
  abort("Saved GeoPackage geometry failed CRS, validity, or emptiness checks.")
}

coordinates <- as.matrix(model_data[, c("easting_m", "northing_m")])
storage.mode(coordinates) <- "double"
if (any(!is.finite(coordinates)) || anyDuplicated(as.data.frame(coordinates)) ||
    !isTRUE(all.equal(
      unname(coordinates), unname(bundle$coordinates), tolerance = 0
    ))) {
  abort("Saved representative-point coordinates are inconsistent.")
}
points <- sf::st_as_sf(
  model_data, coords = c("easting_m", "northing_m"), crs = epsg,
  remove = FALSE
)
covered <- sf::st_covered_by(points, counties)
if (any(!vapply(seq_along(covered), function(index) {
      index %in% covered[[index]]
    }, logical(1L)))) {
  abort("A saved representative point is outside its county.")
}

fold <- as.integer(folds$fold)
n_folds <- as.integer(config$spatial_folds$n_folds)
fold_sizes <- tabulate(fold, nbins = n_folds)
recomputed_fold <- gwrs::spatial_folds(
  coordinates,
  n_folds = n_folds,
  seed = as.integer(config$spatial_folds$seed),
  nstart = as.integer(config$spatial_folds$nstart)
)
if (!identical(fold, as.integer(bundle$folds)) ||
    !identical(fold, recomputed_fold) ||
    !identical(fold_sizes, as.integer(fold_summary$n)) ||
    !identical(fold_sizes, as.integer(unlist(
      completion$spatial_fold_sizes, use.names = FALSE
    ))) ||
    min(fold_sizes) < as.integer(config$spatial_folds$minimum_fold_size) ||
    max(fold_sizes) / min(fold_sizes) >
      as.double(config$spatial_folds$maximum_size_ratio)) {
  abort("Frozen fold assignment failed deterministic or size validation.")
}

if (nrow(edges) != completion$common_graph_directed_edges ||
    nrow(edges) != graph_summary$directed_edges ||
    graph_summary$vertices != expected_n ||
    graph_summary$initial_isolates !=
      as.integer(config$common_diagnostic_graph$expected_initial_isolates) ||
    graph_summary$final_isolates != 0L ||
    graph_summary$final_components !=
      as.integer(config$common_diagnostic_graph$expected_final_components) ||
    anyDuplicated(edges[c("from_index", "to_index")]) ||
    any(edges$from_index == edges$to_index) ||
    any(!edges$edge_type %in% c("queen", "isolate_nearest")) ||
    any(edges$binary_weight != 1)) {
  abort("Common diagnostic graph structure is invalid.")
}
edge_key <- paste(edges$from_index, edges$to_index, sep = ":")
reverse_key <- paste(edges$to_index, edges$from_index, sep = ":")
if (!setequal(edge_key, reverse_key)) {
  abort("Common diagnostic graph is not symmetric.")
}
weight_sum <- tapply(
  edges$row_standard_weight, edges$from_index, sum
)
if (length(weight_sum) != expected_n ||
    any(abs(weight_sum - 1) > 1e-12)) {
  abort("Common diagnostic graph weights are not row-standardized.")
}

adjacency <- split(edges$to_index, edges$from_index)
visited <- rep(FALSE, expected_n)
queue <- 1L
visited[[1L]] <- TRUE
while (length(queue)) {
  current <- queue[[1L]]
  queue <- queue[-1L]
  next_vertices <- adjacency[[as.character(current)]]
  next_vertices <- next_vertices[!visited[next_vertices]]
  if (length(next_vertices)) {
    visited[next_vertices] <- TRUE
    queue <- c(queue, next_vertices)
  }
}
if (!all(visited)) abort("Common diagnostic graph is disconnected.")

expected_isolates <- sort(unlist(
  config$common_diagnostic_graph$expected_isolate_fips,
  use.names = FALSE
))
if (nrow(repairs) != length(expected_isolates) ||
    !identical(sort(repairs$isolate_fips), expected_isolates) ||
    sum(edges$edge_type == "isolate_nearest") != 2L * nrow(repairs) ||
    any(!is.finite(repairs$distance_m)) || any(repairs$distance_m <= 0)) {
  abort("Isolate-repair records are incomplete or inconsistent.")
}

response <- catalog$alias[catalog$role == "response"]
predictors <- catalog$alias[catalog$role == "predictor"]
if (length(response) != 1L || length(predictors) != 26L ||
    any(!c(response, predictors) %in% names(model_data)) ||
    any(!is.finite(model_data[[response]])) ||
    any(!is.finite(as.matrix(model_data[, predictors, drop = FALSE])))) {
  abort("Model-ready response or predictor matrix is invalid.")
}
if (geometry_audit$fips_matched !=
      as.integer(config$expected_domain_counties) ||
    geometry_audit$eligible_analysis_counties != expected_n ||
    geometry_audit$analysis_epsg != epsg ||
    geometry_audit$invalid_geometries_after_repair != 0L ||
    geometry_audit$empty_geometries != 0L ||
    geometry_audit$representative_points_outside_own_county != 0L ||
    geometry_audit$duplicate_representative_points != 0L) {
  abort("Geometry audit contains a failed gate.")
}

message(
  "Spatial release verified read-only: n=", expected_n,
  ", p=26, folds=", paste(fold_sizes, collapse = ","),
  ", graph edges=", nrow(edges), "."
)
message("No model fitted; geometry, folds, graph, and manifests are valid.")
