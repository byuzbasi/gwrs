# study-v2 adapter over the hash-bound study-v1 execution engine.

study_v2_app <- normalizePath(app, mustWork = TRUE)
study_v1_app <- normalizePath(file.path(study_v2_app, "../study-v1"), mustWork = TRUE)
Sys.setenv(
  GWRS_STUDY_V1_APP = study_v1_app,
  GWRS_STUDY_INPUT_ROOT = file.path(study_v1_app, "frozen/input")
)
source(file.path(study_v1_app, "code/bootstrap.R"))

study_v2_base_load <- usa_study_v1_load
usa_study_v1_load <- function(app, mode) {
  Sys.setenv(
    GWRS_STUDY_V1_APP = normalizePath(file.path(app, "../study-v1"), mustWork = TRUE),
    GWRS_STUDY_INPUT_ROOT = normalizePath(
      file.path(app, "../study-v1/frozen/input"), mustWork = TRUE
    )
  )
  value <- study_v2_base_load(app, mode)
  if (mode == "smoke") value$ctx$cfg$q <- as.double(value$policy$smoke$q)
  # study-v1's read-only preflight expects a fixed-k compatibility field.
  # study-v2 does not use this field to construct its grid or fit a model;
  # every candidate is resolved from q inside its own training split. The
  # stability floor is supplied here solely so the inherited preflight checks
  # the attainable lower bound rather than the obsolete full-sample k value.
  value$ctx$cfg$k <- as.integer(min(
    value$ctx$cfg$n, 2L * (value$ctx$cfg$p + 1L)
  ))
  namespace <- asNamespace("gwrs")
  value$ctx$g$compute_gwr_model_diagnostics <- get(
    "compute_gwr_model_diagnostics", envir = namespace
  )
  value$ctx$g$compute_gwr_nonconvex_model_diagnostics <- get(
    "compute_gwr_nonconvex_model_diagnostics", envir = namespace
  )
  value
}
