# Reproducible interpolation benchmark; run after installing gwrs.
if (!"package:gwrs" %in% search()) library(gwrs)

set.seed(20260818)
n_source <- 1000L
source <- cbind(runif(n_source), runif(n_source))
values <- cbind(
  value_1 = sin(4 * source[, 1L]) + cos(5 * source[, 2L]),
  selected_1 = as.double(source[, 1L] > 0.25)
)
bandwidth <- gwrs:::gwrs_auto_bandwidth(source, "euclidean")$bandwidth
cell_counts <- c(1000L, 10000L, 100000L)
max_neighbors <- 200L
chunk_size <- 2000L

make_query <- function(n_cell) {
  nx <- ceiling(sqrt(n_cell))
  grid <- expand.grid(
    x = seq(0, 1, length.out = nx),
    y = seq(0, 1, length.out = nx),
    KEEP.OUT.ATTRS = FALSE
  )
  as.matrix(grid[seq_len(n_cell), , drop = FALSE])
}

rows <- list()
results <- list()
row_index <- 0L
for (n_cell in cell_counts) {
  query <- make_query(n_cell)
  for (method in c("kernel", "idw")) {
    gc(FALSE)
    timing <- system.time({
      result <- gwrs:::gwrs_interpolate_value_matrix(
        source = source,
        query = query,
        values = values,
        interpolation = method,
        distance = "euclidean",
        kernel = "gaussian",
        bandwidth = if (method == "kernel") bandwidth else NULL,
        idw_power = 2,
        max_neighbors = max_neighbors,
        chunk_size = chunk_size
      )
    })
    row_index <- row_index + 1L
    rows[[row_index]] <- data.frame(
      method = method,
      grid_cells = n_cell,
      elapsed_seconds = unname(timing[["elapsed"]]),
      result_mib = as.numeric(object.size(result)) / 1024^2,
      finite_fraction = mean(is.finite(result$values)),
      stringsAsFactors = FALSE
    )
    results[[paste(method, n_cell, sep = "_")]] <- result
  }
}

summary <- do.call(rbind, rows)
rownames(summary) <- NULL
block_working_mib <- min(max(cell_counts), chunk_size) * max_neighbors *
  (8 + 4 + 8) / 1024^2

benchmark <- list(
  dimensions = c(
    source_locations = n_source,
    interpolated_columns = ncol(values),
    max_neighbors = max_neighbors,
    chunk_size = chunk_size
  ),
  automatic_bandwidth = bandwidth,
  timings = summary,
  approximate_block_distance_index_weight_mib = block_working_mib,
  full_distance_matrix_gib_avoided_at_100k =
    8 * max(cell_counts) * n_source / 1024^3,
  session = c(
    R = R.version.string,
    platform = R.version$platform,
    gwrs = as.character(utils::packageVersion("gwrs"))
  )
)
print(benchmark)
