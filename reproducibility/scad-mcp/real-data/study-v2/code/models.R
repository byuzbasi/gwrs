# Preserve study-v1 numerical validation and append selected-fit diagnostics.

study_v1_app <- Sys.getenv("GWRS_STUDY_V1_APP")
sys.source(file.path(study_v1_app, "code/models.R"), envir = environment())
study_v2_base_install_models <- usa_study_v1_install_models

usa_study_v1_install_models <- function(e) {
  study_v2_base_install_models(e)
  validated_evaluate <- e$usa_evaluate
  e$usa_evaluate <- function(ctx, train, data, split, candidate,
                             until = NULL, keep = FALSE,
                             full_diagnostics = FALSE) {
    ctx$native_capture$selected_model_diagnostics <- NULL
    value <- validated_evaluate(
      ctx, train, data, split, candidate, until, keep, full_diagnostics
    )
    value$selected_model_diagnostics <-
      ctx$native_capture$selected_model_diagnostics
    value
  }
}
