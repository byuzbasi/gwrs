study_v2_app <- Sys.getenv("GWRS_STUDY_V2_APP")
stopifnot(nzchar(study_v2_app), dir.exists(study_v2_app))
sys.source(file.path(study_v2_app, "code/diagnostics.R"), envir = environment())
