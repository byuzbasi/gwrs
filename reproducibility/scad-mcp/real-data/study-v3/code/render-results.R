# Reuse the validated study-v2 publication renderer. Because this file is the
# Rscript entry point, the inherited renderer resolves `app` to study-v3 and
# validates the study-v3 run before creating any derived artifact.

file_option <- commandArgs(FALSE)
file_option <- sub("^--file=", "", file_option[grepl("^--file=", file_option)])
stopifnot(length(file_option) == 1L)
study_v3_app <- normalizePath(
  file.path(dirname(file_option), ".."), mustWork = TRUE
)
study_v2_renderer <- normalizePath(
  file.path(study_v3_app, "../study-v2/code/render-results.R"),
  mustWork = TRUE
)
sys.source(study_v2_renderer, envir = globalenv())
