# study-v3 is a path-correction adapter over the hash-bound study-v2 design.

study_v3_app <- normalizePath(app, mustWork = TRUE)
study_v2_app <- normalizePath(file.path(study_v3_app, "../study-v2"), mustWork = TRUE)
study_v1_app <- normalizePath(file.path(study_v3_app, "../study-v1"), mustWork = TRUE)
Sys.setenv(
  GWRS_STUDY_V2_APP = study_v2_app,
  GWRS_STUDY_V1_APP = study_v1_app,
  GWRS_STUDY_INPUT_ROOT = file.path(study_v1_app, "frozen/input")
)
source(file.path(study_v2_app, "code/bootstrap.R"))

# The inherited loader already applies the study-v2 q-fraction design. Reset
# the bound paths on every call so subprocess environments cannot redirect
# scientific inputs.
study_v3_base_load <- usa_study_v1_load
usa_study_v1_load <- function(app, mode) {
  v3 <- normalizePath(app, mustWork = TRUE)
  v2 <- normalizePath(file.path(v3, "../study-v2"), mustWork = TRUE)
  v1 <- normalizePath(file.path(v3, "../study-v1"), mustWork = TRUE)
  Sys.setenv(
    GWRS_STUDY_V2_APP = v2,
    GWRS_STUDY_V1_APP = v1,
    GWRS_STUDY_INPUT_ROOT = file.path(v1, "frozen/input")
  )
  study_v3_base_load(v3, mode)
}
