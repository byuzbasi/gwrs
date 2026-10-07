args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
app <- dirname(dirname(script))
source(file.path(app, "code", "common.R"), local = TRUE)
require_package("jsonlite")

catalog <- utils::read.csv(
  file.path(app, "config", "variable-catalog-v1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
stopifnot(
  nrow(catalog) == 27L,
  sum(catalog$role == "response") == 1L,
  sum(catalog$role == "predictor") == 26L,
  !anyDuplicated(catalog$alias),
  !anyDuplicated(catalog$estimate_code),
  !anyDuplicated(catalog$moe_code)
)

data_config <- read_json_file(
  file.path(app, "config", "data-sources-v1.json"), simplify = FALSE
)
sources <- data_config$sources
api_sources <- Filter(function(source) !is.null(source$api_key_env), sources)
stopifnot(length(api_sources) == 2L)
api_urls <- vapply(api_sources, `[[`, character(1L), "url")
stopifnot(
  all(grepl("key=\\{CENSUS_API_KEY\\}", api_urls)),
  !any(grepl("key=[A-Za-z0-9]{20,}", api_urls))
)

api_fields <- lapply(api_urls, function(url) {
  get <- sub("^.*[?]get=([^&]+)&for=.*$", "\\1", url)
  strsplit(get, ",", fixed = TRUE)[[1L]]
})
stopifnot(all(vapply(api_fields, length, integer(1L)) <= 50L))
downloaded_codes <- unlist(lapply(api_fields, function(fields) {
  setdiff(fields, "NAME")
}), use.names = FALSE)
required_codes <- c(catalog$estimate_code, catalog$moe_code)
stopifnot(
  setequal(downloaded_codes, required_codes),
  !anyDuplicated(downloaded_codes)
)

fixture_dir <- tempfile("gwrs-acs-fixture-")
dir.create(fixture_dir)
on.exit(unlink(fixture_dir, recursive = TRUE), add = TRUE)
batch_1_path <- file.path(fixture_dir, "batch-1.json")
batch_2_path <- file.path(fixture_dir, "batch-2.json")
writeLines(
  jsonlite::toJSON(
    list(
      as.list(c("NAME", "A", "AM", "state", "county")),
      as.list(c("Alpha County", "10", "1", "01", "001")),
      list("Beta County", NULL, "2", "01", "003")
    ),
    auto_unbox = TRUE, null = "null"
  ),
  batch_1_path
)
writeLines(
  jsonlite::toJSON(
    list(
      as.list(c("NAME", "B", "BM", "state", "county")),
      as.list(c("Alpha County", "20", "2", "01", "001")),
      as.list(c("Beta County", "30", "3", "01", "003"))
    ),
    auto_unbox = TRUE, null = "null"
  ),
  batch_2_path
)
fixture <- merge_census_batches(
  parse_census_json(batch_1_path), parse_census_json(batch_2_path)
)
stopifnot(
  nrow(fixture) == 2L,
  identical(fixture$NAME, c("Alpha County", "Beta County")),
  identical(acs_estimate_numeric(fixture$A), c(10, NA_real_)),
  identical(acs_estimate_numeric(fixture$B), c(20, 30))
)

special_fixture <- c(
  "1.5", "-999999999", "-888888888", "-666666666",
  "-555555555.0", "-333333333", "-222222222", NA_character_
)
stopifnot(
  identical(
    acs_numeric_status(special_fixture),
    c(
      "observed", "sentinel_no_data_small_sample",
      "sentinel_not_applicable_or_available",
      "sentinel_insufficient_sample", "sentinel_controlled_moe",
      "sentinel_open_ended_median_moe", "sentinel_moe_not_computable",
      "json_null"
    )
  ),
  identical(
    acs_estimate_numeric(special_fixture),
    c(1.5, rep(NA_real_, 7L))
  ),
  identical(
    acs_moe_numeric(special_fixture),
    c(1.5, NA_real_, NA_real_, NA_real_, 0, NA_real_, NA_real_,
      NA_real_)
  )
)

atomic_path <- file.path(fixture_dir, "atomic.txt")
write_text_atomic_once("same", atomic_path)
write_text_atomic_once("same", atomic_path)
overwrite_error <- try(write_text_atomic_once("different", atomic_path),
                       silent = TRUE)
stopifnot(inherits(overwrite_error, "try-error"))

message(
  "Acquisition smoke tests passed: catalogue, API batching, redaction, ",
  "Census null-position parsing, sentinel handling, and overwrite guard."
)
