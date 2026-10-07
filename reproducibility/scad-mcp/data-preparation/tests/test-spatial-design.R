args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
app <- dirname(dirname(script))
source(file.path(app, "code", "common.R"), local = TRUE)
source(file.path(app, "code", "spatial-common.R"), local = TRUE)
require_spatial_packages()

square <- function(xmin, ymin, xmax, ymax) {
  sf::st_polygon(list(rbind(
    c(xmin, ymin), c(xmax, ymin), c(xmax, ymax), c(xmin, ymax),
    c(xmin, ymin)
  )))
}
fixture <- sf::st_sf(
  fips = c("00001", "00002", "00003"),
  name = c("Left", "Middle", "Island"),
  geometry = sf::st_sfc(
    square(0, 0, 1, 1), square(1, 0, 2, 1), square(4, 0, 5, 1),
    crs = 5070
  )
)
coordinates <- rbind(c(0.5, 0.5), c(1.5, 0.5), c(4.5, 0.5))
graph <- suppressWarnings(common_queen_graph(
  fixture, coordinates, fixture$fips, fixture$name, snap_metres = 0.001
))
stopifnot(
  graph$summary$vertices == 3L,
  graph$summary$initial_undirected_edges == 1,
  graph$summary$final_undirected_edges == 2,
  graph$summary$initial_isolates == 1L,
  graph$summary$final_isolates == 0L,
  graph$summary$initial_components == 2L,
  graph$summary$final_components == 1L,
  nrow(graph$repairs) == 1L,
  graph$repairs$isolate_fips == "00003",
  graph$repairs$attached_fips == "00002",
  nrow(graph$directed_edges) == 4L,
  sum(graph$directed_edges$edge_type == "isolate_nearest") == 2L
)
weight_sum <- tapply(
  graph$directed_edges$row_standard_weight,
  graph$directed_edges$from_index,
  sum
)
stopifnot(all(abs(weight_sum - 1) < 1e-12))

fold_coordinates <- as.matrix(expand.grid(x = 1:6, y = 1:5))
fold_1 <- gwrs::spatial_folds(
  fold_coordinates, n_folds = 5L, seed = 20260907L, nstart = 20L
)
fold_2 <- gwrs::spatial_folds(
  fold_coordinates, n_folds = 5L, seed = 20260907L, nstart = 20L
)
stopifnot(
  identical(fold_1, fold_2),
  identical(sort(unique(fold_1)), 1:5),
  all(tabulate(fold_1, nbins = 5L) > 0L)
)

fixture_root <- tempfile("gwrs-spatial-manifest-")
dir.create(fixture_root)
on.exit(unlink(fixture_root, recursive = TRUE), add = TRUE)
fixture_file <- file.path(fixture_root, "value.txt")
writeLines("fixed", fixture_file)
manifest <- file_records(fixture_file, fixture_root)
manifest_path <- file.path(fixture_root, "manifest.csv")
write_csv_atomic_once(manifest, manifest_path)
validated <- validate_record_manifest(manifest_path, fixture_root)
stopifnot(nrow(validated) == 1L, validated$path == "value.txt")

message(
  "Spatial design unit tests passed: isolate repair, graph weights, ",
  "deterministic folds, and SHA-256 manifest validation."
)
