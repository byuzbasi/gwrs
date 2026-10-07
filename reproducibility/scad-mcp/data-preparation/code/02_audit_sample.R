args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
app <- dirname(dirname(script))
source(file.path(app, "code", "common.R"), local = TRUE)

metadata_gate <- file.path(app, "validation", "METADATA_VALIDATED.json")
if (!file.exists(metadata_gate)) {
  abort("Metadata gate is absent; run 00_validate_metadata.R first.")
}
gate <- read_json_file(metadata_gate)
if (!isTRUE(gate$estimate_labels_match) || !isTRUE(gate$moe_labels_match) ||
    gate$predictor_count != 26L) {
  abort("Metadata gate is invalid.")
}

study <- read_json_file(file.path(app, "config", "study-v1.json"))
catalog <- utils::read.csv(
  file.path(app, "config", "variable-catalog-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
special_policy <- utils::read.csv(
  file.path(app, "config", "acs-special-values-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
expected_special_values <- c(
  -999999999, -888888888, -666666666, -555555555, -333333333,
  -222222222
)
if (!setequal(special_policy$numeric_value, expected_special_values) ||
    anyDuplicated(special_policy$numeric_value) ||
    anyDuplicated(special_policy$status) ||
    special_policy$audit_v2_handling[
      special_policy$numeric_value == -555555555
    ] != "zero_with_controlled_flag") {
  abort("The ACS special-value policy is incomplete or inconsistent.")
}
raw_dir <- file.path(app, "data", "raw", study$raw_snapshot)
batch_files <- file.path(
  raw_dir,
  c(
    "acs5-profile-2024-counties-batch-01.json",
    "acs5-profile-2024-counties-batch-02.json"
  )
)
boundary_files <- file.path(
  raw_dir,
  c("cb_2024_us_county_5m.zip", "cb_2024_us_state_5m.zip")
)
source_receipt <- validate_raw_files(
  app, raw_dir, basename(c(batch_files, boundary_files))
)

batch_1 <- parse_census_json(batch_files[[1L]])
batch_2 <- parse_census_json(batch_files[[2L]])
county_all <- merge_census_batches(batch_1, batch_2)

required_codes <- c(catalog$estimate_code, catalog$moe_code)
if (!all(required_codes %in% names(county_all))) {
  abort(
    "Downloaded batches omit catalogue codes: ",
    paste(setdiff(required_codes, names(county_all)), collapse = ", ")
  )
}
county_all$fips <- paste0(county_all$state, county_all$county)
if (anyDuplicated(county_all$fips) || any(nchar(county_all$fips) != 5L)) {
  abort("ACS county FIPS keys are not unique five-character identifiers.")
}

excluded_state_fips <- unlist(study$excluded_state_fips, use.names = FALSE)
in_domain <- !county_all$state %in% excluded_state_fips
county <- county_all[in_domain, , drop = FALSE]
county <- county[order(county$fips), , drop = FALSE]
if (nrow(county) < 3000L || nrow(county) > 3200L) {
  abort("Unexpected contiguous-USA+DC county count: ", nrow(county))
}

estimate_status <- lapply(catalog$estimate_code, function(code) {
  acs_numeric_status(county[[code]])
})
estimate_values <- lapply(catalog$estimate_code, function(code) {
  acs_estimate_numeric(county[[code]])
})
names(estimate_values) <- catalog$alias
names(estimate_status) <- catalog$alias
moe_status <- lapply(catalog$moe_code, function(code) {
  acs_numeric_status(county[[code]])
})
moe_values <- lapply(catalog$moe_code, function(code) {
  acs_moe_numeric(county[[code]])
})
names(moe_values) <- paste0(catalog$alias, "__moe90")
names(moe_status) <- paste0(catalog$alias, "__moe90")

response_index <- which(catalog$role == "response")
response_raw <- estimate_values[[response_index]]
response_alias <- catalog$alias[[response_index]]
predictor_catalog <- catalog[catalog$role == "predictor", , drop = FALSE]

missing_matrix <- vapply(
  seq_len(nrow(catalog)),
  function(index) !is.finite(estimate_values[[index]]),
  logical(nrow(county))
)
colnames(missing_matrix) <- catalog$alias
nonpositive_response <- is.finite(response_raw) & response_raw <= 0
eligible <- rowSums(missing_matrix) == 0L & !nonpositive_response

exclusion_reason <- vapply(seq_len(nrow(county)), function(index) {
  reasons <- paste0("missing_estimate:",
                    colnames(missing_matrix)[missing_matrix[index, ]])
  if (nonpositive_response[[index]]) {
    reasons <- c(reasons, "nonpositive_response")
  }
  paste(reasons, collapse = ";")
}, character(1L))

processed <- data.frame(
  fips = county$fips,
  county_name = county$NAME,
  state_fips = county$state,
  county_fips = county$county,
  eligible_complete_case = eligible,
  exclusion_reason = exclusion_reason,
  median_household_income_2024 = response_raw,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (index in seq_len(nrow(catalog))) {
  alias <- catalog$alias[[index]]
  if (index == response_index) {
    processed[[alias]] <- ifelse(
      is.finite(response_raw) & response_raw > 0,
      log(response_raw),
      NA_real_
    )
  } else {
    processed[[alias]] <- estimate_values[[index]]
  }
  processed[[paste0(alias, "__moe90")]] <- moe_values[[index]]
  processed[[paste0(alias, "__moe90_controlled")]] <-
    moe_status[[index]] == "sentinel_controlled_moe"
}

quantile_or_na <- function(value, probability) {
  value <- value[is.finite(value)]
  if (!length(value)) return(NA_real_)
  unname(stats::quantile(value, probability, names = FALSE, type = 7L))
}

variable_audit <- do.call(rbind, lapply(seq_len(nrow(catalog)), function(index) {
  estimate <- estimate_values[[index]]
  moe <- moe_values[[index]]
  estimate_state <- estimate_status[[index]]
  moe_state <- moe_status[[index]]
  relative <- moe / abs(estimate)
  relative[!is.finite(relative) | abs(estimate) <= 0] <- NA_real_
  data.frame(
    role = catalog$role[[index]],
    alias = catalog$alias[[index]],
    estimate_code = catalog$estimate_code[[index]],
    moe_code = catalog$moe_code[[index]],
    n_domain = nrow(county),
    estimate_missing_n = sum(!is.finite(estimate)),
    estimate_json_null_n = sum(estimate_state == "json_null"),
    estimate_special_value_n = sum(estimate_state != "observed"),
    estimate_zero_n = sum(is.finite(estimate) & estimate == 0),
    estimate_median = quantile_or_na(estimate, 0.50),
    estimate_q05 = quantile_or_na(estimate, 0.05),
    estimate_q95 = quantile_or_na(estimate, 0.95),
    moe_missing_n = sum(!is.finite(moe)),
    moe_json_null_n = sum(moe_state == "json_null"),
    moe_controlled_n = sum(moe_state == "sentinel_controlled_moe"),
    moe_unavailable_n = sum(
      !moe_state %in% c("observed", "sentinel_controlled_moe")
    ),
    moe_zero_n = sum(is.finite(moe) & moe == 0),
    moe_median = quantile_or_na(moe, 0.50),
    moe_q90 = quantile_or_na(moe, 0.90),
    relative_moe_median = quantile_or_na(relative, 0.50),
    relative_moe_q90 = quantile_or_na(relative, 0.90),
    moe_ge_abs_estimate_n = sum(
      is.finite(moe) & is.finite(estimate) & moe >= abs(estimate)
    ),
    stringsAsFactors = FALSE
  )
}))

status_levels <- c(
  "observed", "json_null", "blank", "nonnumeric",
  "sentinel_no_data_small_sample",
  "sentinel_not_applicable_or_available",
  "sentinel_insufficient_sample", "sentinel_controlled_moe",
  "sentinel_open_ended_median_moe", "sentinel_moe_not_computable",
  "other_negative"
)
sentinel_rows <- function(index, field_role) {
  state <- if (field_role == "estimate") {
    estimate_status[[index]]
  } else {
    moe_status[[index]]
  }
  code <- if (field_role == "estimate") {
    catalog$estimate_code[[index]]
  } else {
    catalog$moe_code[[index]]
  }
  data.frame(
    role = catalog$role[[index]],
    alias = catalog$alias[[index]],
    field_role = field_role,
    code = code,
    status = status_levels,
    n = tabulate(match(state, status_levels), nbins = length(status_levels)),
    audit_v2_handling = ifelse(
      status_levels == "observed", "numeric_value",
      ifelse(
        status_levels == "sentinel_controlled_moe",
        "zero_with_controlled_flag", "missing"
      )
    ),
    stringsAsFactors = FALSE
  )
}
sentinel_audit <- do.call(
  rbind,
  unlist(lapply(seq_len(nrow(catalog)), function(index) {
    list(sentinel_rows(index, "estimate"), sentinel_rows(index, "moe"))
  }), recursive = FALSE)
)
if (any(xtabs(n ~ alias + field_role, sentinel_audit) != nrow(county))) {
  abort("ACS special-value counts do not cover every county and field.")
}

sample_flow <- data.frame(
  stage = c(
    "All ACS county and county-equivalent rows",
    "Excluded state/territory FIPS rows",
    "Contiguous 48 states plus District of Columbia",
    "Excluded for missing/nonpositive required estimates",
    "Complete-case candidate analysis sample"
  ),
  n = c(
    nrow(county_all),
    sum(!in_domain),
    nrow(county),
    sum(!eligible),
    sum(eligible)
  ),
  stringsAsFactors = FALSE
)

processed_dir <- file.path(app, "data", "processed", "audit-v2")
tables_dir <- file.path(app, "tables", "audit-v2")
validation_dir <- file.path(app, "validation", "audit-v2")
processed_path <- file.path(processed_dir, "county-profile-audit-v2.csv")
exclusions_path <- file.path(tables_dir, "county-exclusions-v2.csv")
variable_audit_path <- file.path(tables_dir, "variable-missingness-moe-v2.csv")
sentinel_audit_path <- file.path(tables_dir, "acs-special-values-v2.csv")
sample_flow_path <- file.path(tables_dir, "sample-flow-v2.csv")
source_receipt_path <- file.path(validation_dir, "source-receipt-v2.csv")

write_csv_atomic_once(processed, processed_path)
write_csv_atomic_once(
  processed[!eligible, c("fips", "county_name", "state_fips",
                         "exclusion_reason"), drop = FALSE],
  exclusions_path
)
write_csv_atomic_once(variable_audit, variable_audit_path)
write_csv_atomic_once(sentinel_audit, sentinel_audit_path)
write_csv_atomic_once(sample_flow, sample_flow_path)
write_csv_atomic_once(source_receipt, source_receipt_path)

if (anyDuplicated(processed$fips) || nrow(processed) != nrow(county)) {
  abort("Processed audit output failed key/dimension validation.")
}
if (sum(processed$eligible_complete_case) != sample_flow$n[[5L]]) {
  abort("Eligible sample count is internally inconsistent.")
}
if (nrow(predictor_catalog) != 26L) {
  abort("Predictor count changed during the audit.")
}

prior_path <- file.path(
  app, "data", "processed", "audit-v1", "county-profile-audit-v1.csv"
)
if (!file.exists(prior_path)) {
  abort("The preserved audit-v1 sample is required for membership parity.")
}
prior <- utils::read.csv(
  prior_path, stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = "character"
)
if (!identical(prior$fips, processed$fips) ||
    !identical(prior$eligible_complete_case == "TRUE", eligible) ||
    !identical(prior$exclusion_reason, processed$exclusion_reason)) {
  abort("Audit-v2 unexpectedly changed sample membership or exclusions.")
}

output_paths <- c(
  processed_path, exclusions_path, variable_audit_path, sentinel_audit_path,
  sample_flow_path, source_receipt_path
)
output_manifest <- file_records(output_paths, app)
output_manifest_path <- file.path(validation_dir,
                                  "output-manifest-sha256-v2.csv")
write_csv_atomic_once(output_manifest, output_manifest_path)

completion <- list(
  schema = "gwrs-usa-counties-acs2024-sample-audit-v2",
  all_acs_counties = nrow(county_all),
  domain_counties = nrow(county),
  eligible_counties = sum(eligible),
  excluded_for_required_estimates = sum(!eligible),
  response_count = 1L,
  predictor_count = nrow(predictor_catalog),
  imputation_performed = FALSE,
  moe_based_exclusion_performed = FALSE,
  controlled_moe_values_set_to_zero = sum(vapply(
    moe_status,
    function(state) sum(state == "sentinel_controlled_moe"),
    integer(1L)
  )),
  unavailable_moe_values_set_to_missing = sum(vapply(
    moe_status,
    function(state) sum(
      !state %in% c("observed", "sentinel_controlled_moe")
    ),
    integer(1L)
  )),
  row_level_controlled_moe_flags_written = TRUE,
  sample_membership_unchanged_from_audit_v1 = TRUE,
  special_value_policy = "config/acs-special-values-v1.csv",
  output_manifest_sha256 = sha256_file(output_manifest_path),
  model_fitting_locked = TRUE
)
write_json_atomic_once(
  completion,
  file.path(validation_dir, "AUDIT_COMPLETED.json")
)

message(
  "Sample audit v2 complete: ", sum(eligible), "/", nrow(county),
  " domain counties eligible; controlled MOEs are zero with row flags."
)
message("Review tables/audit-v2 before defining geometry, folds, or models.")
