args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
app <- dirname(dirname(script))
source(file.path(app, "code", "common.R"), local = TRUE)
source(file.path(app, "code", "spatial-common.R"), local = TRUE)
require_spatial_packages()

main <- function() {
project_root <- find_gwrs_project_root(app)

config_path <- file.path(app, "config", "spatial-design-v1.json")
config <- read_json_file(config_path)
if (!isTRUE(config$model_fitting_locked) ||
    config$schema != "gwrs-usa-counties-acs2024-spatial-design-v1" ||
    config$output_root != "spatial-v1") {
  abort("Spatial design configuration is invalid or unlocked.")
}
release_name <- config$output_root
release_root <- file.path(app, release_name)
if (file.exists(release_root)) {
  verified <- validate_spatial_release(app, release_root)
  message(
    "RESUME spatial-v1: ", verified$completion$eligible_counties,
    " counties; all signatures and outputs are valid. No model fitted."
  )
  quit(save = "no", status = 0L)
}

stage_root <- file.path(
  app, paste0(".", release_name, "-build-", Sys.getpid())
)
if (file.exists(stage_root)) abort("Staging path already exists: ", stage_root)
dir.create(file.path(stage_root, "data"), recursive = TRUE)
dir.create(file.path(stage_root, "tables"), recursive = TRUE)
dir.create(file.path(stage_root, "validation"), recursive = TRUE)
on.exit({
  if (file.exists(stage_root)) unlink(stage_root, recursive = TRUE)
}, add = TRUE)

study_path <- file.path(app, "config", "study-v1.json")
catalog_path <- file.path(app, "config", "variable-catalog-v1.csv")
study <- read_json_file(study_path)
catalog <- utils::read.csv(
  catalog_path, stringsAsFactors = FALSE, check.names = FALSE
)
design_gate_path <- file.path(app, "validation", "DESIGN_VALIDATED_v3.json")
design_manifest_path <- file.path(
  app, "validation", "design-manifest-sha256-v3.csv"
)
audit_gate_path <- file.path(
  app, "validation", "audit-v2", "AUDIT_COMPLETED.json"
)
audit_manifest_path <- file.path(
  app, "validation", "audit-v2", "output-manifest-sha256-v2.csv"
)
design_gate <- read_json_file(design_gate_path)
audit_gate <- read_json_file(audit_gate_path)
if (design_gate$schema !=
      "gwrs-usa-counties-acs2024-design-validation-v3" ||
    !isTRUE(design_gate$model_fitting_locked) ||
    !isTRUE(design_gate$geometry_and_folds_locked) ||
    design_gate$design_manifest_sha256 != sha256_file(design_manifest_path)) {
  abort("Design-v3 input gate is invalid.")
}
invisible(validate_record_manifest(design_manifest_path, app))
if (audit_gate$schema !=
      "gwrs-usa-counties-acs2024-sample-audit-v2" ||
    !isTRUE(audit_gate$model_fitting_locked) ||
    audit_gate$output_manifest_sha256 != sha256_file(audit_manifest_path)) {
  abort("Audit-v2 input gate is invalid.")
}
invisible(validate_record_manifest(audit_manifest_path, app))

expected_domain <- as.integer(config$expected_domain_counties)
expected_eligible <- as.integer(config$expected_eligible_counties)
if (audit_gate$domain_counties != expected_domain ||
    audit_gate$eligible_counties != expected_eligible ||
    nrow(catalog) != 27L || sum(catalog$role == "predictor") != 26L) {
  abort("Spatial configuration and audit-v2 dimensions disagree.")
}

description_version <- unname(read.dcf(file.path(project_root, "DESCRIPTION"))[
  1L, "Version"
])
installed_version <- as.character(utils::packageVersion("gwrs"))
if (!identical(installed_version, description_version)) {
  abort(
    "Installed gwrs version differs from the project source: ",
    installed_version, " versus ", description_version
  )
}
source_environment <- new.env(parent = globalenv())
sys.source(file.path(project_root, "R", "cv.R"), envir = source_environment)
if (!identical(
      formals(source_environment$spatial_folds),
      formals(gwrs::spatial_folds)
    ) ||
    !identical(
      body(source_environment$spatial_folds), body(gwrs::spatial_folds)
    )) {
  abort("Installed gwrs::spatial_folds differs from R/cv.R.")
}
message("SPATIAL-V1 phase=preflight status=complete")

raw_dir <- file.path(app, "data", "raw", study$raw_snapshot)
county_zip <- file.path(raw_dir, "cb_2024_us_county_5m.zip")
state_zip <- file.path(raw_dir, "cb_2024_us_state_5m.zip")
validate_raw_files(
  app, raw_dir,
  c(
    "acs5-profile-2024-counties-batch-01.json",
    "acs5-profile-2024-counties-batch-02.json",
    basename(county_zip), basename(state_zip)
  )
)
audit_path <- file.path(
  app, "data", "processed", "audit-v2", "county-profile-audit-v2.csv"
)
audit <- utils::read.csv(
  audit_path, stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = c(
    fips = "character", state_fips = "character",
    county_fips = "character"
  )
)
if (nrow(audit) != expected_domain || anyDuplicated(audit$fips) ||
    sum(audit$eligible_complete_case) != expected_eligible) {
  abort("Audit-v2 sample keys or dimensions are invalid.")
}
message("SPATIAL-V1 phase=input status=complete")

county_uri <- paste0(
  "/vsizip/", normalizePath(county_zip, mustWork = TRUE),
  "/cb_2024_us_county_5m.shp"
)
state_uri <- paste0(
  "/vsizip/", normalizePath(state_zip, mustWork = TRUE),
  "/cb_2024_us_state_5m.shp"
)
counties <- sf::st_read(county_uri, quiet = TRUE)
states <- sf::st_read(state_uri, quiet = TRUE)
excluded_state_fips <- unlist(study$excluded_state_fips, use.names = FALSE)
counties <- counties[
  !counties$STATEFP %in% excluded_state_fips, , drop = FALSE
]
states <- states[!states$STATEFP %in% excluded_state_fips, , drop = FALSE]
if (nrow(counties) != expected_domain || nrow(states) != 49L ||
    anyDuplicated(counties$GEOID) || any(sf::st_is_empty(counties))) {
  abort("Unexpected or invalid CONUS+DC boundary inventory.")
}
boundary_index <- match(audit$fips, counties$GEOID)
if (anyNA(boundary_index) ||
    length(setdiff(counties$GEOID, audit$fips)) != 0L) {
  abort("County FIPS linkage is not exhaustive and one-to-one.")
}
counties <- counties[boundary_index, , drop = FALSE]
if (!identical(counties$GEOID, audit$fips)) {
  abort("County boundary order does not match audit-v2.")
}
for (column in names(audit)) counties[[column]] <- audit[[column]]

epsg <- as.integer(config$projection$epsg)
counties <- sf::st_transform(sf::st_make_valid(counties), epsg)
states <- sf::st_transform(sf::st_make_valid(states), epsg)
if (any(!sf::st_is_valid(counties)) || any(sf::st_is_empty(counties)) ||
    sf::st_crs(counties)$epsg != epsg || sf::st_crs(states)$epsg != epsg) {
  abort("Projected geometry failed validity or CRS checks.")
}
message("SPATIAL-V1 phase=geometry status=complete")
analysis_counties <- counties[counties$eligible_complete_case, , drop = FALSE]
if (nrow(analysis_counties) != expected_eligible) {
  abort("Eligible analysis geometry has the wrong row count.")
}
representative_points <- suppressWarnings(
  sf::st_point_on_surface(analysis_counties)
)
coordinates <- sf::st_coordinates(representative_points)[, c("X", "Y"),
                                                          drop = FALSE]
colnames(coordinates) <- c("easting_m", "northing_m")
storage.mode(coordinates) <- "double"
attr(coordinates, "crs") <- sf::st_crs(epsg)
covered <- sf::st_covered_by(representative_points, analysis_counties)
if (any(!is.finite(coordinates)) || anyDuplicated(as.data.frame(coordinates)) ||
    any(!vapply(seq_along(covered), function(index) {
      index %in% covered[[index]]
    }, logical(1L)))) {
  abort("Representative points are non-unique, non-finite, or outside counties.")
}

n_folds <- as.integer(config$spatial_folds$n_folds)
fold_seed <- as.integer(config$spatial_folds$seed)
fold_nstart <- as.integer(config$spatial_folds$nstart)
fold <- gwrs::spatial_folds(
  coordinates, n_folds = n_folds, seed = fold_seed, nstart = fold_nstart
)
fold_repeat <- gwrs::spatial_folds(
  coordinates, n_folds = n_folds, seed = fold_seed, nstart = fold_nstart
)
fold_sizes <- tabulate(fold, nbins = n_folds)
if (!identical(fold, fold_repeat) ||
    !identical(sort(unique(fold)), seq_len(n_folds)) ||
    min(fold_sizes) < as.integer(config$spatial_folds$minimum_fold_size) ||
    max(fold_sizes) / min(fold_sizes) >
      as.double(config$spatial_folds$maximum_size_ratio)) {
  abort("Frozen spatial folds failed determinism, coverage, or size gates.")
}

graph <- suppressWarnings(common_queen_graph(
  analysis_counties, coordinates, analysis_counties$fips,
  analysis_counties$county_name,
  snap_metres = as.double(config$common_diagnostic_graph$snap_metres)
))
expected_isolates <- sort(unlist(
  config$common_diagnostic_graph$expected_isolate_fips,
  use.names = FALSE
))
if (graph$summary$initial_isolates !=
      as.integer(config$common_diagnostic_graph$expected_initial_isolates) ||
    graph$summary$final_components !=
      as.integer(config$common_diagnostic_graph$expected_final_components) ||
    !identical(sort(graph$repairs$isolate_fips), expected_isolates)) {
  abort("Common graph differs from the prespecified isolate/component audit.")
}
message("SPATIAL-V1 phase=folds_graph status=complete")

state_lookup <- sf::st_drop_geometry(states)[, c("STATEFP", "STUSPS", "NAME")]
if (anyDuplicated(state_lookup$STATEFP)) abort("State lookup is not unique.")
state_index <- match(analysis_counties$state_fips, state_lookup$STATEFP)
if (anyNA(state_index)) abort("An eligible county has no state lookup.")

analysis_counties$state_abbreviation <- state_lookup$STUSPS[state_index]
analysis_counties$state_name <- state_lookup$NAME[state_index]
analysis_counties$easting_m <- coordinates[, "easting_m"]
analysis_counties$northing_m <- coordinates[, "northing_m"]
analysis_counties$spatial_fold <- fold

scientific_columns <- c(
  "median_household_income_2024", catalog$alias,
  paste0(catalog$alias, "__moe90"),
  paste0(catalog$alias, "__moe90_controlled")
)
if (!all(scientific_columns %in% names(analysis_counties))) {
  abort("Analysis geometry omits a catalogue estimate or MOE field.")
}
model_data <- data.frame(
  fips = analysis_counties$fips,
  county_name = analysis_counties$county_name,
  state_fips = analysis_counties$state_fips,
  state_abbreviation = analysis_counties$state_abbreviation,
  state_name = analysis_counties$state_name,
  county_fips = analysis_counties$county_fips,
  easting_m = coordinates[, "easting_m"],
  northing_m = coordinates[, "northing_m"],
  spatial_fold = fold,
  sf::st_drop_geometry(analysis_counties)[, scientific_columns, drop = FALSE],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
response_alias <- catalog$alias[catalog$role == "response"]
predictor_aliases <- catalog$alias[catalog$role == "predictor"]
if (anyDuplicated(model_data$fips) || nrow(model_data) != expected_eligible ||
    any(!is.finite(model_data[[response_alias]])) ||
    any(!is.finite(as.matrix(model_data[, predictor_aliases, drop = FALSE])))) {
  abort("Model-ready estimates failed key, dimension, or finite-value checks.")
}
message("SPATIAL-V1 phase=model_data status=complete")

map_columns <- unique(c(
  "fips", "county_name", "state_fips", "state_abbreviation", "state_name",
  "county_fips", "median_household_income_2024", response_alias,
  predictor_aliases, "easting_m", "northing_m", "spatial_fold"
))
map_counties <- analysis_counties[, map_columns, drop = FALSE]
if (ncol(sf::st_drop_geometry(map_counties)) != length(map_columns) ||
    anyDuplicated(names(map_counties))) {
  abort("Portable map geometry has missing or duplicate attributes.")
}

fold_table <- model_data[, c(
  "fips", "county_name", "state_fips", "state_abbreviation",
  "easting_m", "northing_m", "spatial_fold"
)]
names(fold_table)[names(fold_table) == "spatial_fold"] <- "fold"
fold_summary <- do.call(rbind, lapply(seq_len(n_folds), function(index) {
  selected <- fold == index
  data.frame(
    fold = index,
    n = sum(selected),
    state_count = length(unique(model_data$state_fips[selected])),
    easting_min_m = min(coordinates[selected, "easting_m"]),
    easting_max_m = max(coordinates[selected, "easting_m"]),
    northing_min_m = min(coordinates[selected, "northing_m"]),
    northing_max_m = max(coordinates[selected, "northing_m"]),
    stringsAsFactors = FALSE
  )
}))
geometry_audit <- data.frame(
  boundary_counties_all_products = nrow(sf::st_read(county_uri, quiet = TRUE)),
  boundary_counties_domain = nrow(counties),
  audit_counties_domain = nrow(audit),
  fips_matched = sum(counties$GEOID == audit$fips),
  eligible_analysis_counties = nrow(analysis_counties),
  conus_dc_states = nrow(states),
  source_epsg = 4269L,
  analysis_epsg = epsg,
  invalid_geometries_after_repair = sum(!sf::st_is_valid(counties)),
  empty_geometries = sum(sf::st_is_empty(counties)),
  representative_points_outside_own_county = 0L,
  duplicate_representative_points = anyDuplicated(as.data.frame(coordinates)),
  stringsAsFactors = FALSE
)

spatial_bundle <- list(
  schema = "gwrs-usa-counties-acs2024-spatial-bundle-v1",
  model_data = model_data,
  coordinates = coordinates,
  folds = fold,
  analysis_counties = analysis_counties,
  states = states,
  common_graph = graph,
  metadata = list(
    configuration = config,
    data_vintage = study$data_vintage,
    response = response_alias,
    predictors = predictor_aliases,
    gwrs_version = installed_version,
    projection = sf::st_crs(epsg),
    model_fitting_locked = TRUE
  )
)

input_paths <- c(
  file.path(project_root, c("DESCRIPTION", "R/cv.R")),
  file.path(app, c(
  "SPATIAL_PROTOCOL_v1.md",
  "config/study-v1.json", "config/variable-catalog-v1.csv",
  "config/spatial-design-v1.json", "code/common.R",
  "code/spatial-common.R", "code/03_prepare_spatial_design.R",
  "code/04_validate_spatial_design.R", "tests/test-spatial-design.R",
  "validation/DESIGN_VALIDATED_v3.json",
  "validation/design-manifest-sha256-v3.csv",
  "validation/audit-v2/AUDIT_COMPLETED.json",
  "validation/audit-v2/output-manifest-sha256-v2.csv",
  "data/processed/audit-v2/county-profile-audit-v2.csv",
  paste0("data/raw/", study$raw_snapshot, "/source_manifest_v4.csv"),
  paste0("data/raw/", study$raw_snapshot, "/cb_2024_us_county_5m.zip"),
  paste0("data/raw/", study$raw_snapshot, "/cb_2024_us_state_5m.zip")
)))
if (any(!file.exists(input_paths))) {
  abort("A spatial-release input is absent: ",
        input_paths[!file.exists(input_paths)][1L])
}
input_manifest <- file_records(input_paths, project_root)
input_manifest_path <- file.path(
  stage_root, "validation", "input-manifest-sha256-v1.csv"
)
write_csv_atomic_once(input_manifest, input_manifest_path)
message("SPATIAL-V1 phase=input_manifest status=complete")

data_paths <- file.path(stage_root, "data", c(
  "county-model-data-v1.csv", "spatial-design-v1.rds",
  "analysis-counties-v1.gpkg", "conus-states-v1.gpkg"
))
table_paths <- file.path(stage_root, "tables", c(
  "spatial-folds-v1.csv", "fold-summary-v1.csv",
  "common-neighbor-graph-directed-v1.csv", "isolate-repairs-v1.csv",
  "graph-summary-v1.csv", "geometry-audit-v1.csv"
))
write_csv_atomic_once(model_data, data_paths[[1L]])
message("SPATIAL-V1 phase=write_csv status=complete")
write_rds_atomic_new(spatial_bundle, data_paths[[2L]])
message("SPATIAL-V1 phase=write_rds status=complete")
write_gpkg_atomic_new(
  map_counties, data_paths[[3L]], layer = "analysis_counties"
)
message("SPATIAL-V1 phase=write_county_gpkg status=complete")
write_gpkg_atomic_new(states, data_paths[[4L]], layer = "conus_states")
message("SPATIAL-V1 phase=write_state_gpkg status=complete")
write_csv_atomic_once(fold_table, table_paths[[1L]])
write_csv_atomic_once(fold_summary, table_paths[[2L]])
write_csv_atomic_once(graph$directed_edges, table_paths[[3L]])
write_csv_atomic_once(graph$repairs, table_paths[[4L]])
write_csv_atomic_once(graph$summary, table_paths[[5L]])
write_csv_atomic_once(geometry_audit, table_paths[[6L]])

session_path <- file.path(
  stage_root, "validation", "session-info-v1.txt"
)
session_text <- c(
  capture.output(sessionInfo()), "", "Library paths:", .libPaths(), "",
  "Package locations:",
  paste("gwrs", find.package("gwrs")),
  paste("sf", find.package("sf")),
  paste("spdep", find.package("spdep")), "",
  "sf external software:", capture.output(sf::sf_extSoftVersion())
)
write_text_atomic_once(session_text, session_path)

manifested_stage_paths <- c(
  input_manifest_path, data_paths, table_paths, session_path
)
output_manifest <- stage_file_records(
  manifested_stage_paths, stage_root, release_name
)
output_manifest_path <- file.path(
  stage_root, "validation", "output-manifest-sha256-v1.csv"
)
write_csv_atomic_once(output_manifest, output_manifest_path)
independent_stage_records <- stage_file_records(
  manifested_stage_paths, stage_root, release_name
)
if (!identical(output_manifest, independent_stage_records)) {
  abort("Staged output manifest is not reproducible.")
}

completion <- list(
  schema = "gwrs-usa-counties-acs2024-spatial-release-v1",
  domain_counties = nrow(audit),
  eligible_counties = nrow(model_data),
  predictor_count = length(predictor_aliases),
  analysis_epsg = epsg,
  coordinate_unit = config$projection$coordinate_unit,
  representative_point_method = config$representative_points$method,
  spatial_fold_method = config$spatial_folds$method,
  spatial_fold_seed = fold_seed,
  spatial_fold_sizes = as.list(as.integer(fold_sizes)),
  common_graph_vertices = graph$summary$vertices,
  common_graph_directed_edges = graph$summary$directed_edges,
  common_graph_initial_isolates = graph$summary$initial_isolates,
  common_graph_final_components = graph$summary$final_components,
  isolate_fips = as.list(graph$repairs$isolate_fips),
  gwrs_version = installed_version,
  input_manifest_sha256 = sha256_file(input_manifest_path),
  output_manifest_sha256 = sha256_file(output_manifest_path),
  model_fits_created = 0L,
  model_fitting_locked = TRUE
)
write_json_atomic_once(
  completion,
  file.path(stage_root, "validation", "SPATIAL_DESIGN_COMPLETED.json")
)

if (!file.rename(stage_root, release_root)) {
  abort("Atomic spatial-release directory move failed: ", release_root)
}
verified <- validate_spatial_release(app, release_root)
message(
  "Spatial design v1 complete: ", verified$completion$eligible_counties,
  " counties, folds=", paste(fold_sizes, collapse = ","),
  ", directed graph edges=", graph$summary$directed_edges, "."
)
message("No model fitted; modelling protocol remains locked.")
}

main()
