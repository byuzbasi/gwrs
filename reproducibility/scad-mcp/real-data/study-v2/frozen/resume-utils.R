# Reuse the validated resume utilities, changing only RDS compression.
study_v1_app <- Sys.getenv("GWRS_STUDY_V1_APP")
stopifnot(nzchar(study_v1_app), dir.exists(study_v1_app))
sys.source(
  file.path(study_v1_app, "frozen/resume-utils.R"),
  envir = environment()
)

gwrs_sim_atomic_rds <- function(object, path, replace = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !replace) {
    stop("Refusing to overwrite `", path, "`.", call. = FALSE)
  }
  temporary <- tempfile(".partial-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3, compress = "gzip")
  if (!file.rename(temporary, path)) {
    if (!replace || !file.exists(path)) {
      stop("Could not atomically create `", path, "`.", call. = FALSE)
    }
    unlink(path)
    if (!file.rename(temporary, path)) {
      stop("Could not replace derived state `", path, "`.", call. = FALSE)
    }
  }
  invisible(path)
}
