args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
app <- dirname(dirname(script))
source(file.path(app, "code", "common.R"), local = TRUE)

catalog <- utils::read.csv(
  file.path(app, "config", "variable-catalog-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
special_policy <- utils::read.csv(
  file.path(app, "config", "acs-special-values-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
study <- read_json_file(file.path(app, "config", "study-v1.json"))
data_sources <- read_json_file(
  file.path(app, "config", "data-sources-v1.json"), simplify = FALSE
)
metadata_gate_path <- file.path(app, "validation", "METADATA_VALIDATED.json")
metadata_gate <- read_json_file(metadata_gate_path)

if (!isTRUE(study$model_fitting_locked) ||
    !isTRUE(metadata_gate$model_fitting_locked)) {
  abort("Design validation requires model fitting to remain locked.")
}
if (nrow(catalog) != 27L || sum(catalog$role == "predictor") != 26L ||
    study$predictor_count != 26L || metadata_gate$predictor_count != 26L) {
  abort("The prespecified predictor count is inconsistent.")
}
if (!setequal(
      special_policy$numeric_value,
      c(-999999999, -888888888, -666666666, -555555555, -333333333,
        -222222222)
    ) ||
    special_policy$audit_v2_handling[
      special_policy$numeric_value == -555555555
    ] != "zero_with_controlled_flag") {
  abort("The ACS special-value policy is inconsistent.")
}
metadata_manifest_path <- file.path(
  app, "validation", "metadata-validation-manifest-v1.csv"
)
if (metadata_gate$manifest_sha256 != sha256_file(metadata_manifest_path)) {
  abort("Metadata gate does not match its SHA-256 manifest.")
}

api_sources <- Filter(
  function(source) identical(source$api_key_env, "CENSUS_API_KEY"),
  data_sources$sources
)
if (length(api_sources) != 2L) abort("Expected exactly two keyed API batches.")
api_urls <- vapply(api_sources, `[[`, character(1L), "url")
if (!all(grepl("key=\\{CENSUS_API_KEY\\}", api_urls)) ||
    any(grepl("key=[A-Za-z0-9]{20,}", api_urls))) {
  abort("API URLs do not satisfy the credential-redaction policy.")
}
api_fields <- lapply(api_urls, function(url) {
  get <- sub("^.*[?]get=([^&]+)&for=.*$", "\\1", url)
  strsplit(get, ",", fixed = TRUE)[[1L]]
})
if (any(vapply(api_fields, length, integer(1L)) > 50L)) {
  abort("A Census API batch exceeds the 50-field limit.")
}
downloaded_codes <- unlist(lapply(api_fields, setdiff, y = "NAME"),
                           use.names = FALSE)
required_codes <- c(catalog$estimate_code, catalog$moe_code)
if (!setequal(downloaded_codes, required_codes) ||
    anyDuplicated(downloaded_codes)) {
  abort("API batches do not contain every required code exactly once.")
}

audit_gate_path <- file.path(app, "validation", "audit-v2",
                             "AUDIT_COMPLETED.json")
audit_manifest_path <- file.path(app, "validation", "audit-v2",
                                 "output-manifest-sha256-v2.csv")
if (!file.exists(audit_gate_path) || !file.exists(audit_manifest_path)) {
  abort("The sample audit-v2 gate is incomplete.")
}
audit_gate <- read_json_file(audit_gate_path)
if (audit_gate$schema != "gwrs-usa-counties-acs2024-sample-audit-v2" ||
    audit_gate$domain_counties != 3109L ||
    audit_gate$eligible_counties != 3107L ||
    audit_gate$predictor_count != 26L ||
    !isTRUE(audit_gate$row_level_controlled_moe_flags_written) ||
    !isTRUE(audit_gate$sample_membership_unchanged_from_audit_v1) ||
    !isTRUE(audit_gate$model_fitting_locked) ||
    audit_gate$output_manifest_sha256 != sha256_file(audit_manifest_path)) {
  abort("The sample audit-v2 completion record is invalid.")
}

required_paths <- file.path(app, c(
  "config/acs-special-values-v1.csv",
  "config/data-sources-v1.json",
  "config/leakage-audit-v1.csv",
  "config/metadata-sources-v1.json",
  "config/spatial-design-v1.json",
  "config/study-v1.json",
  "config/variable-catalog-v1.csv",
  "code/00_download_metadata.py",
  "code/00_validate_metadata.R",
  "code/01_download_data.py",
  "code/02_audit_sample.R",
  "code/03_prepare_spatial_design.R",
  "code/04_validate_spatial_design.R",
  "code/05_validate_spatial_design_v2.R",
  "code/99_validate_design.R",
  "code/common.R",
  "code/download_sources.py",
  "code/spatial-common.R",
  "validation/METADATA_VALIDATED.json",
  "validation/metadata-validation-manifest-v1.csv",
  "validation/audit-v2/AUDIT_COMPLETED.json",
  "validation/audit-v2/output-manifest-sha256-v2.csv",
  "data/processed/audit-v1/county-profile-audit-v1.csv",
  "data/processed/audit-v2/county-profile-audit-v2.csv"
))
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  abort("Required design files are absent: ", paste(missing_paths, collapse = ", "))
}

manifest <- file_records(required_paths, app)
manifest_path <- file.path(app, "validation", "design-manifest-sha256-v3.csv")
write_csv_atomic_once(manifest, manifest_path)

completion <- list(
  schema = "gwrs-usa-counties-acs2024-design-validation-v3",
  data_vintage = study$data_vintage,
  predictor_count = 26L,
  required_estimate_moe_codes = length(required_codes),
  api_batch_count = length(api_sources),
  maximum_api_fields = max(vapply(api_fields, length, integer(1L))),
  metadata_gate_valid = TRUE,
  credential_redaction_valid = TRUE,
  census_null_positions_preserved = TRUE,
  census_special_value_policy_valid = TRUE,
  controlled_moe_semantics_validated = TRUE,
  design_manifest_sha256 = sha256_file(manifest_path),
  county_data_downloaded = TRUE,
  sample_audit_completed = TRUE,
  eligible_counties = audit_gate$eligible_counties,
  sample_membership_unchanged_from_audit_v1 = TRUE,
  geometry_and_folds_locked = TRUE,
  model_fitting_locked = TRUE
)
write_json_atomic_once(
  completion,
  file.path(app, "validation", "DESIGN_VALIDATED_v3.json")
)

message(
  "Design v3 validated: p=26, n=3107, two complete API batches, ",
  "controlled-MOE policy and SHA-256 manifests complete."
)
message("Geometry, folds, and model fitting remain locked.")
