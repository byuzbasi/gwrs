# Run from any directory; dependencies and gwrs must already be installed.
# This builds documentation locally and does not install or publish anything.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
stopifnot(length(script) == 1L)
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
stopifnot(requireNamespace("pkgdown", quietly = TRUE))
stopifnot(requireNamespace("gwrs", quietly = TRUE))
expected <- read.dcf(file.path(root, "DESCRIPTION"))[1, "Version"]
if (as.character(utils::packageVersion("gwrs")) != expected) {
  stop("Install the current gwrs source before building its executable guide.")
}
pkg <- pkgdown::as_pkgdown(root)
# Curate website reference coverage without changing the package help or API.
topics <- c(
  "gwr_fit", "gwr_select", "gwr_select_verify", "summary.gwrs_selection",
  "gwr_diagnostics", "gwr_diagnostic_tests", "gwr_moran", "gwr_local_inference",
  "gwr_local_collinearity", "gwr_global_collinearity",
  "gwr_reference_formulas", "gwr_reference_fits", "gwr_nonconvex",
  "gwr_nonconvex_fit", "gwr_nonconvex_path", "cv_gwr_nonconvex",
  "gwr_nonconvex_cv_cache",
  "predict.gwrs_fit", "summary.gwrs_fit", "gwr_neighbors", "spatial_folds",
  "gwrs_control", "acs2020_2024_counties", "acs2024_counties", "gwrs-package"
)
stopifnot(all(topics %in% pkg$topics$name))
pkg$topics <- pkg$topics[pkg$topics$name %in% topics, ]
# Reference examples remain displayed; only the two small guides are executed.
pkgdown::build_site(pkg, examples = FALSE, new_process = FALSE,
                    install = FALSE, devel = FALSE, preview = FALSE)
# Version assets so local previews cannot reuse another package's cached logo/CSS.
asset_version <- gsub(".", "-", expected, fixed = TRUE)
logo_asset <- paste0("gwrs-logo-", asset_version, ".svg")
style_asset <- paste0("gwrs-style-", asset_version, ".css")
script_asset <- paste0("gwrs-site-", asset_version, ".js")
stopifnot(file.copy(file.path(root, "man/figures/logo.svg"),
                   file.path(pkg$dst_path, logo_asset), overwrite = TRUE))
stopifnot(file.copy(file.path(root, "pkgdown/extra.css"),
                   file.path(pkg$dst_path, style_asset), overwrite = TRUE))
# Search indexes use canonical URLs. Keep local preview navigation on localhost.
script <- readLines(file.path(pkg$dst_path, "pkgdown.js"), warn = FALSE)
search_navigation <- 'window.location.href = s.path + "?q=" + q + "#" + s.id;'
stopifnot(any(grepl(search_navigation, script, fixed = TRUE)))
local_navigation <- paste(
  'var destination = s.path;',
  'if (["localhost", "127.0.0.1"].includes(window.location.hostname)) {',
  '  destination = new URL(s.path, window.location.href).pathname.replace(/^\\/gwrs\\//, "/");',
  '}',
  'window.location.href = destination + "?q=" + encodeURIComponent(q) + "#" + (s.id || "");',
  sep = "\n")
script <- gsub(search_navigation, local_navigation, script, fixed = TRUE)
writeLines(script, file.path(pkg$dst_path, script_asset))
# Link shared path/CV wrappers to the executable guide. Those wrappers share
# broader package help topics; the website retains its curated model scope.
pages <- list.files(pkg$dst_path, recursive = TRUE, pattern = "\\.(html|md)$")
for (page in pages) {
  filename <- file.path(pkg$dst_path, page)
  text <- readLines(filename, warn = FALSE, encoding = "UTF-8")
  depth <- length(strsplit(page, "/", fixed = TRUE)[[1]]) - 1L
  prefix <- paste(rep("../", depth), collapse = "")
  target <- paste0(prefix,
    "articles/gwrs-introduction.html#paths-tuning-and-interpretation")
  revised <- gsub(
    "https://rdrr.io/pkg/gwrs/man/(cv_gwr_penalized|gwr_penalized_path)\\.html",
    target, text
  )
  revised <- gsub('href="vignettes/gwrs-selection.Rmd"',
    'href="articles/gwrs-selection.html"', revised, fixed = TRUE)
  revised <- gsub("reproducibility/gwrs_0.4.0.9015.tar.gz",
    "https://github.com/byuzbasi/gwrs/blob/main/reproducibility/gwrs_0.4.0.9015.tar.gz",
    revised, fixed = TRUE)
  revised <- gsub("logo.svg", logo_asset, revised, fixed = TRUE)
  revised <- gsub("extra.css", style_asset, revised, fixed = TRUE)
  revised <- gsub("pkgdown.js", script_asset, revised, fixed = TRUE)
  if (!identical(text, revised)) writeLines(revised, filename, useBytes = TRUE)
}
cat("GWRS_LOCAL_DOCUMENTATION_BUILD_PASS\n")
