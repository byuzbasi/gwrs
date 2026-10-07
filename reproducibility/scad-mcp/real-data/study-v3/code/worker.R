# Execute the frozen study-v2/study-v1 worker stack with study-v3 paths.
arguments <- commandArgs(TRUE)
stopifnot(length(arguments) == 7L)
study_v3_app <- normalizePath(arguments[[2L]], mustWork = TRUE)
study_v2_worker <- normalizePath(
  file.path(study_v3_app, "../study-v2/code/worker.R"), mustWork = TRUE
)
sys.source(study_v2_worker, envir = globalenv())
