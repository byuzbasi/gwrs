# The approved spatial diagnostic definitions are unchanged from study-v1.
study_v1_app <- Sys.getenv("GWRS_STUDY_V1_APP")
sys.source(file.path(study_v1_app, "code/diagnostics.R"), envir = environment())
