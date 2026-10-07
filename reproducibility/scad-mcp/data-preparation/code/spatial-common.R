if (!exists("abort", mode = "function") ||
    !exists("sha256_file", mode = "function")) {
  stop("Source code/common.R before code/spatial-common.R.", call. = FALSE)
}

require_spatial_packages <- function() {
  for (package in c("sf", "spdep", "gwrs", "digest", "jsonlite")) {
    require_package(package)
  }
  invisible(TRUE)
}

find_gwrs_project_root <- function(start) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        file.exists(file.path(current, "R", "cv.R"))) {
      description <- tryCatch(
        read.dcf(file.path(current, "DESCRIPTION")), error = identity
      )
      if (!inherits(description, "error") &&
          identical(unname(description[1L, "Package"]), "gwrs")) {
        return(current)
      }
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  abort("Could not locate the gwrs project root above: ", start)
}

validate_record_manifest <- function(path, root) {
  if (!file.exists(path)) abort("Manifest is absent: ", path)
  records <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c("path", "bytes", "sha256")
  if (!identical(names(records), required) || !nrow(records) ||
      anyDuplicated(records$path) || anyNA(records)) {
    abort("Invalid file manifest structure: ", path)
  }
  if (any(grepl("^/", records$path)) ||
      any(vapply(strsplit(records$path, "/", fixed = TRUE),
                 function(parts) ".." %in% parts, logical(1L)))) {
    abort("Manifest contains a path outside its root: ", path)
  }
  root <- normalizePath(root, mustWork = TRUE)
  full <- file.path(root, records$path)
  if (any(!file.exists(full))) {
    abort("Manifested file is absent: ", records$path[!file.exists(full)][1L])
  }
  resolved <- normalizePath(full, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (any(!startsWith(resolved, prefix))) {
    abort("Manifested file resolves outside its root: ", path)
  }
  actual_bytes <- as.double(file.info(full)$size)
  actual_sha256 <- unname(vapply(full, sha256_file, character(1L)))
  if (any(actual_bytes != as.double(records$bytes)) ||
      any(actual_sha256 != records$sha256)) {
    abort("File size or SHA-256 mismatch in manifest: ", path)
  }
  records
}

write_rds_atomic_new <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) abort("Refusing to overwrite RDS: ", path)
  temporary <- tempfile(
    pattern = paste0(basename(path), "."), tmpdir = dirname(path),
    fileext = ".rds"
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3L, compress = "xz")
  if (!file.rename(temporary, path)) abort("Atomic RDS move failed: ", path)
  invisible(path)
}

write_gpkg_atomic_new <- function(object, path, layer) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) abort("Refusing to overwrite GeoPackage: ", path)
  temporary <- tempfile(
    pattern = paste0(basename(path), "."), tmpdir = dirname(path),
    fileext = ".gpkg"
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  sf::st_write(
    object, temporary, layer = layer, quiet = TRUE, append = FALSE,
    delete_dsn = FALSE
  )
  if (!file.exists(temporary) || file.info(temporary)$size <= 0) {
    abort("GeoPackage writer produced no data: ", path)
  }
  if (!file.rename(temporary, path)) {
    abort("Atomic GeoPackage move failed: ", path)
  }
  invisible(path)
}

stage_file_records <- function(paths, stage_root, release_name) {
  stage_root <- normalizePath(stage_root, mustWork = TRUE)
  paths <- sort(normalizePath(paths, mustWork = TRUE))
  prefix <- paste0(stage_root, .Platform$file.sep)
  if (any(!startsWith(paths, prefix))) {
    abort("A staged output resolves outside the staging root.")
  }
  relative <- substring(paths, nchar(stage_root) + 2L)
  data.frame(
    path = file.path(release_name, relative),
    bytes = as.double(file.info(paths)$size),
    sha256 = unname(vapply(paths, sha256_file, character(1L))),
    stringsAsFactors = FALSE
  )
}

common_queen_graph <- function(geometry,
                               coordinates,
                               ids,
                               labels,
                               snap_metres = 1) {
  require_spatial_packages()
  coordinates <- as.matrix(coordinates)
  storage.mode(coordinates) <- "double"
  ids <- as.character(ids)
  labels <- as.character(labels)
  n <- length(ids)
  if (!inherits(geometry, "sf") || nrow(geometry) != n ||
      nrow(coordinates) != n || ncol(coordinates) != 2L ||
      length(labels) != n || any(!is.finite(coordinates)) ||
      anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    abort("Geometry, coordinates, identifiers, and labels are inconsistent.")
  }
  if (is.na(sf::st_crs(geometry)) || sf::st_is_longlat(geometry)) {
    abort("Queen graph requires a projected polygon geometry with a CRS.")
  }
  if (!is.numeric(snap_metres) || length(snap_metres) != 1L ||
      !is.finite(snap_metres) || snap_metres < 0) {
    abort("`snap_metres` must be one nonnegative finite number.")
  }

  neighbours <- spdep::poly2nb(
    geometry, queen = TRUE, row.names = ids, snap = snap_metres
  )
  initial_neighbours <- neighbours
  initial_degree <- spdep::card(initial_neighbours)
  isolates <- which(initial_degree == 0L)
  nonisolated <- which(initial_degree > 0L)
  if (length(isolates) && !length(nonisolated)) {
    abort("Queen graph has no non-isolated vertex for isolate repair.")
  }
  for (index in isolates) {
    # `spdep` represents a no-neighbour entry by the sentinel integer `0L`.
    # Replace it by a genuinely empty vector before adding a valid index.
    neighbours[[index]] <- integer()
  }

  repairs <- lapply(isolates, function(index) {
    difference <- sweep(
      coordinates[nonisolated, , drop = FALSE], 2L,
      coordinates[index, ], FUN = "-"
    )
    distance_squared <- rowSums(difference^2)
    target <- nonisolated[[which.min(distance_squared)]]
    neighbours[[index]] <<- sort(unique(c(neighbours[[index]], target)))
    neighbours[[target]] <<- sort(unique(c(neighbours[[target]], index)))
    data.frame(
      isolate_index = index,
      isolate_fips = ids[[index]],
      isolate_name = labels[[index]],
      attached_index = target,
      attached_fips = ids[[target]],
      attached_name = labels[[target]],
      distance_m = sqrt(min(distance_squared)),
      rule = "nearest_nonisolated_representative_point_symmetric",
      stringsAsFactors = FALSE
    )
  })
  repairs <- if (length(repairs)) {
    do.call(rbind, repairs)
  } else {
    data.frame(
      isolate_index = integer(), isolate_fips = character(),
      isolate_name = character(), attached_index = integer(),
      attached_fips = character(), attached_name = character(),
      distance_m = numeric(), rule = character(),
      stringsAsFactors = FALSE
    )
  }

  degree <- vapply(neighbours, length, integer(1L))
  if (any(degree == 0L)) abort("Isolate repair left a zero-degree vertex.")
  symmetric <- all(vapply(seq_len(n), function(index) {
    all(vapply(neighbours[[index]], function(other) {
      index %in% neighbours[[other]]
    }, logical(1L)))
  }, logical(1L)))
  if (!symmetric) abort("The repaired neighbour graph is not symmetric.")

  initial_components <- spdep::n.comp.nb(initial_neighbours)$nc
  final_components <- spdep::n.comp.nb(neighbours)$nc
  if (final_components != 1L) {
    abort("The repaired neighbour graph is not connected.")
  }

  from <- rep(seq_len(n), degree)
  to <- unlist(neighbours, use.names = FALSE)
  pair_key <- paste(pmin(from, to), pmax(from, to), sep = ":")
  repair_key <- if (nrow(repairs)) {
    paste(
      pmin(repairs$isolate_index, repairs$attached_index),
      pmax(repairs$isolate_index, repairs$attached_index), sep = ":"
    )
  } else {
    character()
  }
  difference <- coordinates[from, , drop = FALSE] -
    coordinates[to, , drop = FALSE]
  directed <- data.frame(
    from_index = from,
    to_index = to,
    from_fips = ids[from],
    to_fips = ids[to],
    edge_type = ifelse(
      pair_key %in% repair_key, "isolate_nearest", "queen"
    ),
    distance_m = sqrt(rowSums(difference^2)),
    binary_weight = 1,
    row_standard_weight = 1 / degree[from],
    stringsAsFactors = FALSE
  )
  directed <- directed[order(directed$from_index, directed$to_index), ,
                       drop = FALSE]
  row.names(directed) <- NULL
  if (anyDuplicated(directed[c("from_index", "to_index")]) ||
      any(directed$from_index == directed$to_index)) {
    abort("Directed graph contains a duplicate or self edge.")
  }
  row_sums <- tapply(
    directed$row_standard_weight, directed$from_index, sum
  )
  if (length(row_sums) != n || any(abs(row_sums - 1) > 1e-12)) {
    abort("Row-standardized graph weights do not sum to one.")
  }

  summary <- data.frame(
    vertices = n,
    queen = TRUE,
    snap_metres = snap_metres,
    initial_undirected_edges = sum(initial_degree) / 2,
    final_undirected_edges = sum(degree) / 2,
    directed_edges = nrow(directed),
    initial_isolates = length(isolates),
    final_isolates = sum(degree == 0L),
    initial_components = initial_components,
    final_components = final_components,
    minimum_degree = min(degree),
    median_degree = stats::median(degree),
    maximum_degree = max(degree),
    binary_edges = TRUE,
    row_standardized_weights = TRUE,
    stringsAsFactors = FALSE
  )
  list(
    neighbours = neighbours,
    directed_edges = directed,
    repairs = repairs,
    summary = summary
  )
}

validate_spatial_release <- function(app, release_root) {
  require_spatial_packages()
  app <- normalizePath(app, mustWork = TRUE)
  project_root <- find_gwrs_project_root(app)
  release_root <- normalizePath(release_root, mustWork = TRUE)
  completion_path <- file.path(
    release_root, "validation", "SPATIAL_DESIGN_COMPLETED.json"
  )
  input_manifest_path <- file.path(
    release_root, "validation", "input-manifest-sha256-v1.csv"
  )
  output_manifest_path <- file.path(
    release_root, "validation", "output-manifest-sha256-v1.csv"
  )
  if (!file.exists(completion_path)) {
    abort("Spatial release has no completion record: ", release_root)
  }
  completion <- read_json_file(completion_path)
  if (completion$schema !=
        "gwrs-usa-counties-acs2024-spatial-release-v1" ||
      !isTRUE(completion$model_fitting_locked) ||
      completion$input_manifest_sha256 != sha256_file(input_manifest_path) ||
      completion$output_manifest_sha256 != sha256_file(output_manifest_path)) {
    abort("Spatial release completion record is invalid.")
  }
  inputs <- validate_record_manifest(input_manifest_path, project_root)
  outputs <- validate_record_manifest(output_manifest_path, app)
  list(completion = completion, inputs = inputs, outputs = outputs)
}
