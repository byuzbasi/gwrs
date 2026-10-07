# Reuse the complete study-v2 estimator and correct only the post-analysis
# common-graph root.

study_v2_app <- Sys.getenv("GWRS_STUDY_V2_APP")
stopifnot(nzchar(study_v2_app), dir.exists(study_v2_app))
sys.source(file.path(study_v2_app, "code/core.R"), envir = environment())

study_v3_base_common_graph <- usa_common_graph
usa_common_graph <- function(ctx, ids) {
  input_root <- normalizePath(
    Sys.getenv("GWRS_STUDY_INPUT_ROOT"), mustWork = TRUE
  )
  expected_suffix <- file.path("frozen", "input")
  stopifnot(endsWith(input_root, expected_suffix))
  upstream_app <- normalizePath(
    file.path(input_root, "../.."), mustWork = TRUE
  )
  corrected <- ctx
  corrected$app <- upstream_app
  graph <- study_v3_base_common_graph(corrected, ids)
  stopifnot(
    normalizePath(graph$source, mustWork = TRUE) == normalizePath(
      file.path(
        input_root, "spatial/tables",
        "common-neighbor-graph-directed-v1.csv"
      ),
      mustWork = TRUE
    )
  )
  graph
}
