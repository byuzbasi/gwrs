# Descriptive validation-horizon analysis from completed OOF predictions.
# No model is fitted or retuned by this script.

arguments <- commandArgs(TRUE)
if (length(arguments) != 3L || !arguments[[1L]] %in% c("create", "verify")) {
  stop("usage: Rscript --vanilla distance-profile.R create|verify FULL_RUN OUTPUT")
}
action <- arguments[[1L]]
full_run <- normalizePath(arguments[[2L]], mustWork = TRUE)
output_requested <- arguments[[3L]]
if (!grepl("^/", output_requested)) stop("Output path must be absolute")

file_option <- commandArgs(FALSE)
file_option <- sub("^--file=", "", file_option[grepl("^--file=", file_option)])
stopifnot(length(file_option) == 1L)
app <- normalizePath(file.path(dirname(file_option), ".."), mustWork = TRUE)
upstream <- normalizePath(file.path(app, "../study-v1"), mustWork = TRUE)

required_packages <- c("digest", "ggplot2", "jsonlite", "knitr")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1L), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing distance-profile package(s): ", paste(missing_packages, collapse = ", "))
}

hash_file <- function(path) {
  unname(digest::digest(path, algo = "sha256", file = TRUE))
}

file_records <- function(paths, base) {
  paths <- sort(normalizePath(paths, mustWork = TRUE))
  base <- normalizePath(base, mustWork = TRUE)
  data.frame(
    file = substring(paths, nchar(base) + 2L),
    bytes = unname(file.info(paths)$size),
    sha256 = vapply(paths, hash_file, ""),
    stringsAsFactors = FALSE
  )
}

atomic_text <- function(value, path) {
  stopifnot(!file.exists(path))
  temporary <- tempfile(".publish-", dirname(path))
  writeLines(enc2utf8(as.character(value)), temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("Atomic publication failed: ", path)
  invisible(path)
}

atomic_csv <- function(value, path) {
  stopifnot(!file.exists(path))
  temporary <- tempfile(".publish-", dirname(path))
  utils::write.csv(value, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Atomic publication failed: ", path)
  invisible(path)
}

verify_full <- function() {
  required <- c(
    file.path(full_run, "COMPLETED"),
    file.path(full_run, "manifest-sha256.csv"),
    file.path(full_run, "numeric/oof-predictions.rds")
  )
  stopifnot(all(file.exists(required)))
  command <- c(
    file.path(app, "code/run.py"), "verify", "--mode", "full",
    "--output", full_run
  )
  status <- system2("python3", shQuote(command), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
    stop("Full-run verification failed:\n", paste(status, collapse = "\n"))
  }
  invisible(required)
}

required_output <- c(
  "REPORT.md", "session-info.txt",
  "tables/oof-nearest-training-distance.csv",
  "tables/oof-performance-by-distance-band.csv",
  "tables/oof-performance-by-distance-decile.csv",
  "tables/oof-paired-loss-vs-global-ols-by-distance.csv",
  "figures/oof-distance-distribution.png",
  "figures/oof-rmse-by-distance.png"
)

verify_output <- function(output) {
  verify_full()
  output <- normalizePath(output, mustWork = TRUE)
  manifest_path <- file.path(output, "manifest-sha256.csv")
  completion_path <- file.path(output, "COMPLETED")
  stopifnot(
    file.exists(manifest_path), file.exists(completion_path),
    all(file.exists(file.path(output, required_output)))
  )
  recorded <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  recorded$bytes <- as.numeric(recorded$bytes)
  files <- list.files(
    output, recursive = TRUE, full.names = TRUE,
    include.dirs = FALSE, all.files = TRUE, no.. = TRUE
  )
  files <- files[!basename(files) %in% c("manifest-sha256.csv", "COMPLETED")]
  stopifnot(identical(recorded, file_records(files, output)))
  completion <- readLines(completion_path, warn = FALSE)
  stopifnot(
    "Schema: gwrs-usa-counties-acs2024-study-v3-distance-v1" %in% completion,
    paste0(
      "Source-Manifest-SHA256: ",
      hash_file(file.path(full_run, "manifest-sha256.csv"))
    ) %in% completion,
    paste0("Manifest-SHA256: ", hash_file(manifest_path)) %in% completion
  )
  cat("Study-v3 OOF nearest-training-distance diagnostics: OK\n")
}

if (action == "verify") {
  verify_output(output_requested)
  quit(save = "no", status = 0L)
}

verify_full()
if (file.exists(output_requested)) stop("Preserving existing output: ", output_requested)
dir.create(dirname(output_requested), recursive = TRUE, showWarnings = FALSE)
building <- tempfile(paste0(".", basename(output_requested), "-building-"), dirname(output_requested))
dir.create(building)
dir.create(file.path(building, "tables"))
dir.create(file.path(building, "figures"))

cfg <- jsonlite::fromJSON(file.path(app, "config/numerics.json"))
model <- utils::read.csv(
  file.path(upstream, "frozen/input/spatial/data/county-model-data-v1.csv"),
  colClasses = c(fips = "character"), check.names = FALSE
)
predictions <- readRDS(file.path(full_run, "numeric/oof-predictions.rds"))
stopifnot(
  nrow(model) == cfg$n,
  nrow(predictions) == 12L * cfg$n,
  setequal(unique(predictions$sample_id), model$fips),
  setequal(unique(predictions$outer_fold), 1:5)
)

coords <- as.matrix(model[, c("easting_m", "northing_m")])
storage.mode(coords) <- "double"
fold <- as.integer(model$spatial_fold)

nearest_distance <- rep(NA_real_, nrow(model))
for (outer in 1:5) {
  test <- which(fold == outer)
  train <- which(fold != outer)
  train_sq <- rowSums(coords[train, , drop = FALSE]^2)
  test_sq <- rowSums(coords[test, , drop = FALSE]^2)
  squared <- outer(test_sq, train_sq, "+") -
    2 * tcrossprod(coords[test, , drop = FALSE], coords[train, , drop = FALSE])
  squared[squared < 0 & squared > -1e-6] <- 0
  stopifnot(all(is.finite(squared)), all(squared >= 0))
  nearest_distance[test] <- sqrt(apply(squared, 1L, min))
}
stopifnot(all(is.finite(nearest_distance)), all(nearest_distance > 0))

distance <- data.frame(
  sample_id = model$fips,
  outer_fold = fold,
  nearest_training_km = nearest_distance / 1000,
  stringsAsFactors = FALSE
)
breaks_km <- c(as.double(cfg$postanalysis$distance_bands_km), Inf)
labels <- c("[0,50)", "[50,100)", "[100,200)", "[200,400)", "[400,Inf)")
stopifnot(length(labels) == length(breaks_km) - 1L)
distance$distance_band_km <- cut(
  distance$nearest_training_km, breaks = breaks_km, labels = labels,
  include.lowest = TRUE, right = FALSE
)
distance$distance_decile <- pmin(
  10L,
  pmax(1L, ceiling(rank(distance$nearest_training_km, ties.method = "average") /
    (nrow(distance) / 10)))
)
atomic_csv(distance, file.path(building, "tables/oof-nearest-training-distance.csv"))

predictions$distance_index <- match(predictions$sample_id, distance$sample_id)
stopifnot(!anyNA(predictions$distance_index))
predictions$nearest_training_km <- distance$nearest_training_km[predictions$distance_index]
predictions$distance_band_km <- distance$distance_band_km[predictions$distance_index]
predictions$distance_decile <- distance$distance_decile[predictions$distance_index]
predictions$error <- predictions$prediction - predictions$observed

metric_table <- function(value, group_name) {
  groups <- split(value, interaction(value$method, value[[group_name]], drop = TRUE))
  result <- do.call(rbind, lapply(groups, function(part) {
    error <- part$error
    sse <- sum(error^2)
    data.frame(
      method = part$method[[1L]],
      group = as.character(part[[group_name]][[1L]]),
      n = nrow(part),
      distance_min_km = min(part$nearest_training_km),
      distance_median_km = stats::median(part$nearest_training_km),
      distance_max_km = max(part$nearest_training_km),
      rmse = sqrt(mean(error^2)),
      mae = mean(abs(error)),
      bias = mean(error),
      predictive_R2 = 1 - sse / sum((part$observed - part$baseline)^2),
      stringsAsFactors = FALSE
    )
  }))
  names(result)[names(result) == "group"] <- group_name
  rownames(result) <- NULL
  result
}

band_metrics <- metric_table(predictions, "distance_band_km")
decile_metrics <- metric_table(predictions, "distance_decile")
atomic_csv(
  band_metrics,
  file.path(building, "tables/oof-performance-by-distance-band.csv")
)
atomic_csv(
  decile_metrics,
  file.path(building, "tables/oof-performance-by-distance-decile.csv")
)

ols <- predictions[predictions$method == "global-ols", c(
  "sample_id", "error", "nearest_training_km", "distance_band_km"
)]
names(ols)[names(ols) == "error"] <- "ols_error"
comparison <- merge(
  predictions[predictions$method != "global-ols", c(
    "sample_id", "method", "error", "nearest_training_km", "distance_band_km"
  )],
  ols[, c("sample_id", "ols_error")], by = "sample_id", all.x = TRUE,
  sort = FALSE
)
stopifnot(!anyNA(comparison$ols_error))
comparison$squared_loss_difference <- comparison$error^2 - comparison$ols_error^2
comparison$absolute_loss_difference <- abs(comparison$error) - abs(comparison$ols_error)
paired_groups <- split(
  comparison,
  interaction(comparison$method, comparison$distance_band_km, drop = TRUE)
)
paired <- do.call(rbind, lapply(paired_groups, function(part) {
  data.frame(
    method = part$method[[1L]],
    distance_band_km = as.character(part$distance_band_km[[1L]]),
    n = nrow(part),
    mean_squared_loss_difference_vs_ols = mean(part$squared_loss_difference),
    mean_absolute_loss_difference_vs_ols = mean(part$absolute_loss_difference),
    fraction_lower_squared_loss_than_ols = mean(part$squared_loss_difference < 0),
    stringsAsFactors = FALSE
  )
}))
rownames(paired) <- NULL
atomic_csv(
  paired,
  file.path(building, "tables/oof-paired-loss-vs-global-ols-by-distance.csv")
)

distance_plot <- ggplot2::ggplot(
  distance, ggplot2::aes(x = nearest_training_km, fill = factor(outer_fold))
) +
  ggplot2::geom_histogram(bins = 45, alpha = 0.75, position = "identity") +
  ggplot2::facet_wrap(~outer_fold, scales = "free_y") +
  ggplot2::labs(
    title = "OOF validation horizon by regional fold",
    subtitle = "Distance from each held-out county centroid to its nearest training county centroid",
    x = "Nearest training distance (km)", y = "Count", fill = "Outer fold"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(
  file.path(building, "figures/oof-distance-distribution.png"), distance_plot,
  width = 10, height = 6.5, dpi = 320, bg = "white"
)
ggplot2::ggsave(
  file.path(building, "figures/oof-distance-distribution.pdf"), distance_plot,
  width = 10, height = 6.5, device = grDevices::cairo_pdf
)

decile_metrics$distance_decile <- as.integer(decile_metrics$distance_decile)
rmse_plot <- ggplot2::ggplot(
  decile_metrics,
  ggplot2::aes(x = distance_median_km, y = rmse, color = method)
) +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::geom_point(size = 1.4) +
  ggplot2::facet_wrap(~ifelse(grepl("^local-", method), "Local", "Global")) +
  ggplot2::labs(
    title = "OOF RMSE across the regional-transfer distance horizon",
    subtitle = "Distance deciles are descriptive; models are not refitted within distance groups",
    x = "Median nearest-training distance (km)", y = "OOF RMSE", color = "Method"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
  file.path(building, "figures/oof-rmse-by-distance.png"), rmse_plot,
  width = 11, height = 7, dpi = 320, bg = "white"
)
ggplot2::ggsave(
  file.path(building, "figures/oof-rmse-by-distance.pdf"), rmse_plot,
  width = 11, height = 7, device = grDevices::cairo_pdf
)

report <- c(
  "# OOF validation-horizon diagnostics", "",
  "These diagnostics use the already completed nested regional-block OOF predictions.",
  "No model is fitted, tuned, or selected by this analysis.", "",
  paste0("Minimum nearest-training distance: ", signif(min(distance$nearest_training_km), 5), " km."),
  paste0("Median nearest-training distance: ", signif(stats::median(distance$nearest_training_km), 5), " km."),
  paste0("Maximum nearest-training distance: ", signif(max(distance$nearest_training_km), 5), " km."), "",
  "Distance-stratified RMSE describes how the existing regional-transfer errors vary with the validation horizon.",
  "It is not a substitute for the separately planned location/infill-CV estimand.", "",
  paste0("Source run: `", full_run, "`"),
  paste0("Source manifest SHA-256: `", hash_file(file.path(full_run, "manifest-sha256.csv")), "`")
)
atomic_text(report, file.path(building, "REPORT.md"))
atomic_text(c(
  paste0("created_utc: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  paste0("source_run: ", full_run),
  capture.output(sessionInfo()), capture.output(.libPaths())
), file.path(building, "session-info.txt"))

files <- list.files(
  building, recursive = TRUE, full.names = TRUE,
  include.dirs = FALSE, all.files = TRUE, no.. = TRUE
)
manifest <- file_records(files, building)
atomic_csv(manifest, file.path(building, "manifest-sha256.csv"))
atomic_text(c(
  "Schema: gwrs-usa-counties-acs2024-study-v3-distance-v1",
  paste0("Source-Manifest-SHA256: ", hash_file(file.path(full_run, "manifest-sha256.csv"))),
  paste0("Manifest-SHA256: ", hash_file(file.path(building, "manifest-sha256.csv")))
), file.path(building, "COMPLETED"))
if (!file.rename(building, output_requested)) {
  stop("Could not atomically publish the distance-profile directory")
}
verify_output(output_requested)
cat("Study-v3 distance-profile rendering complete: ", output_requested, "\n", sep = "")
