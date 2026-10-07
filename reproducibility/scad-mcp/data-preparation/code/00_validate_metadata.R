args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
app <- dirname(dirname(script))
source(file.path(app, "code", "common.R"), local = TRUE)

raw_dir <- file.path(app, "data", "raw", "2026-09-07")
validation_dir <- file.path(app, "validation")
catalog_path <- file.path(app, "config", "variable-catalog-v1.csv")
catalog <- utils::read.csv(catalog_path, stringsAsFactors = FALSE,
                           check.names = FALSE)

required_columns <- c(
  "position", "role", "domain", "alias", "estimate_code", "moe_code",
  "unit", "transformation", "official_group", "official_label",
  "selection_basis"
)
if (!identical(names(catalog), required_columns)) {
  abort("Unexpected variable-catalog schema.")
}
if (nrow(catalog) != 27L || sum(catalog$role == "response") != 1L ||
    sum(catalog$role == "predictor") != 26L) {
  abort("The catalogue must contain one response and exactly 26 predictors.")
}
if (anyDuplicated(catalog$alias) || anyDuplicated(catalog$estimate_code) ||
    anyDuplicated(catalog$moe_code)) {
  abort("Aliases and ACS estimate/MOE codes must be unique.")
}
if (!identical(sort(catalog$position), 0:26)) {
  abort("Variable positions must be exactly 0 through 26.")
}

metadata_files <- setNames(
  file.path(
    raw_dir,
    paste0("acs5-profile-2024-", c("DP02", "DP03", "DP04", "DP05"),
           "-metadata.json")
  ),
  c("DP02", "DP03", "DP04", "DP05")
)
raw_receipt <- validate_raw_files(app, raw_dir, basename(metadata_files))

metadata <- lapply(metadata_files, function(path) {
  value <- read_json_file(path, simplify = FALSE)
  if (is.null(value$variables) || !length(value$variables)) {
    abort("Official metadata has no variable collection: ", path)
  }
  value$variables
})

audit_rows <- lapply(seq_len(nrow(catalog)), function(index) {
  row <- catalog[index, , drop = FALSE]
  group_variables <- metadata[[row$official_group]]
  estimate <- group_variables[[row$estimate_code]]
  moe <- group_variables[[row$moe_code]]
  if (is.null(estimate) || is.null(moe)) {
    abort("Variable or MOE code is absent from official metadata: ",
          row$estimate_code, " / ", row$moe_code)
  }
  expected_moe_label <- if (startsWith(row$official_label, "Percent!!")) {
    sub("^Percent!!", "Percent Margin of Error!!", row$official_label)
  } else if (startsWith(row$official_label, "Estimate!!")) {
    sub("^Estimate!!", "Margin of Error!!", row$official_label)
  } else {
    abort("Unsupported official estimate label: ", row$estimate_code)
  }
  data.frame(
    position = row$position,
    role = row$role,
    alias = row$alias,
    group = row$official_group,
    estimate_code = row$estimate_code,
    estimate_predicate_type = estimate$predicateType,
    estimate_label_match = identical(estimate$label, row$official_label),
    moe_code = row$moe_code,
    moe_predicate_type = moe$predicateType,
    moe_label_match = identical(moe$label, expected_moe_label),
    stringsAsFactors = FALSE
  )
})
audit <- do.call(rbind, audit_rows)

if (!all(audit$estimate_label_match) || !all(audit$moe_label_match)) {
  bad <- audit[!audit$estimate_label_match | !audit$moe_label_match, ,
               drop = FALSE]
  print(bad)
  abort("At least one catalogue label differs from official ACS metadata.")
}
if (!all(audit$estimate_predicate_type %in% c("int", "float")) ||
    !all(audit$moe_predicate_type %in% c("int", "float"))) {
  abort("All selected ACS estimates and MOEs must be numeric.")
}

audit_path <- file.path(validation_dir, "metadata-audit-v1.csv")
write_csv_atomic_once(audit, audit_path)

manifest_inputs <- c(
  metadata_files,
  catalog_path,
  file.path(app, "config", "leakage-audit-v1.csv"),
  file.path(app, "config", "study-v1.json"),
  audit_path
)
manifest <- file_records(manifest_inputs, app)
manifest_path <- file.path(validation_dir,
                           "metadata-validation-manifest-v1.csv")
write_csv_atomic_once(manifest, manifest_path)

completion <- list(
  schema = "gwrs-usa-counties-acs2024-metadata-validation-v1",
  data_vintage = "2020-2024 ACS 5-year Data Profiles",
  catalog_rows = nrow(catalog),
  response_count = sum(catalog$role == "response"),
  predictor_count = sum(catalog$role == "predictor"),
  estimate_moe_code_count = 2L * nrow(catalog),
  estimate_labels_match = all(audit$estimate_label_match),
  moe_labels_match = all(audit$moe_label_match),
  raw_source_count = nrow(raw_receipt),
  manifest_sha256 = sha256_file(manifest_path),
  model_fitting_locked = TRUE
)
write_json_atomic_once(
  completion,
  file.path(validation_dir, "METADATA_VALIDATED.json")
)

message(
  "Official metadata verified: one response, 26 predictors, 54 estimate/MOE codes."
)
message("Model fitting remains locked.")

