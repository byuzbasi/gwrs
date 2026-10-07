args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(Sys.getenv("GWRS_REPRO_ROOT"), mustWork = TRUE)
target <- args[[2L]]
version <- if (target == "simulation") "0.4.0" else "0.4.0.9006"
required <- c("Rcpp", "RcppArmadillo", "RcppParallel", "FNN", "Matrix",
              "digest", "jsonlite", "ggplot2", "knitr")
required <- c(required, if (target == "simulation") c("ragg", "svglite") else "sf")
if (target == "data") required <- c(required, "spdep")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (args[[1L]] == "install" && length(missing)) {
  install.packages(missing, lib = file.path(root, "runtime/dependencies"),
                   repos = "https://cloud.r-project.org", dependencies = NA)
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
}
if (length(missing)) stop("Missing declared dependencies: ", paste(missing, collapse = ", "),
                         ". Run setup with --install-dependencies, or install them explicitly.")
if (args[[1L]] == "verify") {
  lib <- file.path(root, "runtime", paste0("gwrs-", version), "library")
  stopifnot(normalizePath(find.package("gwrs")) == normalizePath(file.path(lib, "gwrs")),
            as.character(packageVersion("gwrs")) == version)
  if (version == "0.4.0.9006") stopifnot(gwrs:::cpp_arma_uword_bits() == 64L)
  info <- list(R = R.version.string, platform = R.version$platform,
    packages = as.list(vapply(c("gwrs", required), function(p) as.character(packageVersion(p)), "")))
  jsonlite::write_json(info, file.path(dirname(lib), "software.json"), pretty = TRUE, auto_unbox = TRUE)
}
cat("DEPENDENCIES_PASS", target, "R", as.character(getRversion()), "\n")
