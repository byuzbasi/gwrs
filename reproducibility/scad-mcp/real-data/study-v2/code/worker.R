# Execute the frozen study-v1 worker with study-v2 paths and adapters.
arguments <- commandArgs(TRUE)
stopifnot(length(arguments) == 7L)
study_v2_app <- normalizePath(arguments[[2L]], mustWork = TRUE)
study_v1_worker <- normalizePath(
  file.path(study_v2_app, "../study-v1/code/worker.R"), mustWork = TRUE
)
sys.source(study_v1_worker, envir = globalenv())
