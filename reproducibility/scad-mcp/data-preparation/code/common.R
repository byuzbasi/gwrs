abort <- function(...) stop(..., call. = FALSE)

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    abort("Required package is unavailable: ", package)
  }
}

sha256_file <- function(path) {
  require_package("digest")
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

read_json_file <- function(path, simplify = TRUE) {
  require_package("jsonlite")
  jsonlite::fromJSON(path, simplifyVector = simplify)
}

write_text_atomic_once <- function(text, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), "."),
                        tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(text, temporary, useBytes = TRUE)
  if (file.exists(path)) {
    if (!identical(sha256_file(path), sha256_file(temporary))) {
      abort("Existing output differs; refusing to overwrite: ", path)
    }
    return(invisible(path))
  }
  if (!file.rename(temporary, path)) {
    abort("Atomic move failed: ", path)
  }
  invisible(path)
}

write_csv_atomic_once <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), "."),
                        tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(data, temporary, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8", quote = TRUE)
  if (file.exists(path)) {
    if (!identical(sha256_file(path), sha256_file(temporary))) {
      abort("Existing output differs; refusing to overwrite: ", path)
    }
    return(invisible(path))
  }
  if (!file.rename(temporary, path)) {
    abort("Atomic move failed: ", path)
  }
  invisible(path)
}

write_json_atomic_once <- function(object, path) {
  require_package("jsonlite")
  text <- jsonlite::toJSON(
    object, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"
  )
  write_text_atomic_once(c(text, ""), path)
}

parse_census_json <- function(path) {
  payload <- read_json_file(path, simplify = FALSE)
  if (!is.list(payload) || length(payload) < 2L) {
    abort("Census response has no data rows: ", path)
  }
  preserve_cells <- function(row) {
    if (!is.list(row)) row <- as.list(row)
    vapply(row, function(cell) {
      if (is.null(cell) || !length(cell)) return(NA_character_)
      if (length(cell) != 1L || is.list(cell)) {
        abort("Nested or repeated cell in Census response: ", path)
      }
      as.character(cell[[1L]])
    }, character(1L), USE.NAMES = FALSE)
  }
  header <- preserve_cells(payload[[1L]])
  if (anyNA(header) || any(!nzchar(header))) {
    abort("Census response has a missing header field: ", path)
  }
  rows <- lapply(payload[-1L], preserve_cells)
  widths <- vapply(rows, length, integer(1L))
  if (!length(header) || any(widths != length(header)) || anyDuplicated(header)) {
    abort("Invalid Census response shape: ", path)
  }
  matrix_data <- do.call(rbind, rows)
  result <- as.data.frame(matrix_data, stringsAsFactors = FALSE,
                          check.names = FALSE)
  names(result) <- header
  result
}

merge_census_batches <- function(left, right) {
  keys <- c("state", "county")
  if (!all(c(keys, "NAME") %in% names(left)) ||
      !all(c(keys, "NAME") %in% names(right))) {
    abort("Every Census batch must contain NAME, state, and county.")
  }
  if (anyDuplicated(left[keys]) || anyDuplicated(right[keys])) {
    abort("Duplicate state/county key in Census batch.")
  }
  overlap <- intersect(setdiff(names(left), c(keys, "NAME")),
                       setdiff(names(right), c(keys, "NAME")))
  if (length(overlap)) {
    abort("Census batches repeat data columns: ", paste(overlap, collapse = ", "))
  }
  result <- merge(left, right, by = keys, all = FALSE,
                  suffixes = c(".left", ".right"), sort = TRUE)
  if (nrow(result) != nrow(left) || nrow(result) != nrow(right)) {
    abort("Census batches do not have identical county coverage.")
  }
  if (!identical(result$NAME.left, result$NAME.right)) {
    abort("County names differ between Census batches.")
  }
  result$NAME <- result$NAME.left
  result$NAME.left <- NULL
  result$NAME.right <- NULL
  result
}

acs_numeric_status <- function(value) {
  raw <- as.character(value)
  trimmed <- trimws(raw)
  number <- suppressWarnings(as.numeric(trimmed))
  status <- rep("observed", length(trimmed))
  status[is.na(raw)] <- "json_null"
  status[!is.na(raw) & !nzchar(trimmed)] <- "blank"
  status[!is.na(raw) & nzchar(trimmed) & !is.finite(number)] <- "nonnumeric"

  sentinel_status <- c(
    "-999999999" = "sentinel_no_data_small_sample",
    "-888888888" = "sentinel_not_applicable_or_available",
    "-666666666" = "sentinel_insufficient_sample",
    "-555555555" = "sentinel_controlled_moe",
    "-333333333" = "sentinel_open_ended_median_moe",
    "-222222222" = "sentinel_moe_not_computable"
  )
  for (sentinel in names(sentinel_status)) {
    status[is.finite(number) & number == as.numeric(sentinel)] <-
      unname(sentinel_status[[sentinel]])
  }
  recognized <- status != "observed"
  status[!recognized & is.finite(number) & number < 0] <- "other_negative"
  status
}

acs_estimate_numeric <- function(value) {
  number <- suppressWarnings(as.numeric(trimws(as.character(value))))
  status <- acs_numeric_status(value)
  number[status != "observed" | !is.finite(number)] <- NA_real_
  number
}

acs_moe_numeric <- function(value) {
  number <- suppressWarnings(as.numeric(trimws(as.character(value))))
  status <- acs_numeric_status(value)
  number[status == "sentinel_controlled_moe"] <- 0
  number[!status %in% c("observed", "sentinel_controlled_moe") |
           !is.finite(number)] <- NA_real_
  number
}

# Retained for compatibility with the acquisition-v1 smoke test and artifacts.
# New analysis code should call the role-specific conversion functions above.
acs_numeric <- acs_estimate_numeric

latest_manifest_rows <- function(raw_dir) {
  manifests <- list.files(raw_dir, pattern = "^source_manifest_v[0-9]+[.]csv$",
                          full.names = TRUE)
  if (!length(manifests)) abort("No source manifest found in ", raw_dir)
  version <- as.integer(sub(".*source_manifest_v([0-9]+)[.]csv$", "\\1",
                            manifests))
  manifests <- manifests[order(version)]
  combined <- do.call(rbind, lapply(seq_along(manifests), function(index) {
    value <- utils::read.csv(manifests[[index]], stringsAsFactors = FALSE,
                             check.names = FALSE)
    value$.manifest_order <- index
    value
  }))
  combined
}

validate_raw_files <- function(app, raw_dir, filenames) {
  rows <- latest_manifest_rows(raw_dir)
  result <- lapply(filenames, function(filename) {
    candidates <- rows[rows$filename == filename, , drop = FALSE]
    if (!nrow(candidates)) abort("Unmanifested raw file: ", filename)
    row <- candidates[which.max(candidates$.manifest_order), , drop = FALSE]
    if (!row$status %in% c("downloaded", "preserved_existing")) {
      abort("Raw source is incomplete: ", filename, " [", row$status, "]")
    }
    path <- file.path(raw_dir, filename)
    if (!file.exists(path)) abort("Manifested raw file is missing: ", path)
    if (as.double(file.info(path)$size) != as.double(row$bytes) ||
        sha256_file(path) != row$sha256) {
      abort("Raw source checksum/size mismatch: ", path)
    }
    data.frame(
      filename = filename,
      bytes = as.double(row$bytes),
      sha256 = row$sha256,
      relative_path = substring(path, nchar(app) + 2L),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

file_records <- function(paths, root) {
  paths <- sort(normalizePath(paths, mustWork = TRUE))
  data.frame(
    path = substring(paths, nchar(normalizePath(root, mustWork = TRUE)) + 2L),
    bytes = as.double(file.info(paths)$size),
    sha256 = vapply(paths, sha256_file, character(1L)),
    stringsAsFactors = FALSE
  )
}
