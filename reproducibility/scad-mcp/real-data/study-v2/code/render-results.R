# Read-only publication rendering for a completed ACS county study-v2 run.
# No estimator is called here. All outputs are derived from validated shards.

arguments <- commandArgs(TRUE)
if (length(arguments) != 3L ||
    !arguments[[1L]] %in% c("create", "verify")) {
  stop(
    "usage: Rscript --vanilla code/render-results.R create|verify ",
    "FULL_RUN_ABSOLUTE OUTPUT_ABSOLUTE"
  )
}
action <- arguments[[1L]]
full_run <- normalizePath(arguments[[2L]], mustWork = TRUE)
output_requested <- arguments[[3L]]
if (!grepl("^/", output_requested)) stop("Output path must be absolute")

file_option <- commandArgs(FALSE)
file_option <- sub("^--file=", "", file_option[grepl("^--file=", file_option)])
stopifnot(length(file_option) == 1L)
app <- normalizePath(file.path(dirname(file_option), ".."), mustWork = TRUE)
upstream <- normalizePath(file.path(app, "../study-v1"), mustWork = TRUE)

required_packages <- c("digest", "ggplot2", "jsonlite", "knitr", "sf")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1L), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing publication package(s): ", paste(missing_packages, collapse = ", "))
}

hash_file <- function(path) {
  unname(digest::digest(path, algo = "sha256", file = TRUE))
}

file_records <- function(paths, base) {
  paths <- sort(normalizePath(paths, mustWork = TRUE))
  base <- normalizePath(base, mustWork = TRUE)
  data.frame(
    file = substring(paths, nchar(base) + 2L),
    bytes = unname(file.info(paths)$size),
    sha256 = vapply(paths, hash_file, ""),
    stringsAsFactors = FALSE
  )
}

atomic_text <- function(value, path) {
  stopifnot(!file.exists(path))
  temporary <- tempfile(".publish-", dirname(path))
  writeLines(enc2utf8(as.character(value)), temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("Atomic publication failed: ", path)
  invisible(path)
}

atomic_csv <- function(value, path) {
  stopifnot(!file.exists(path))
  temporary <- tempfile(".publish-", dirname(path))
  utils::write.csv(value, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Atomic publication failed: ", path)
  invisible(path)
}

write_latex <- function(value, caption, path, digits = 4L) {
  table <- knitr::kable(
    value, format = "latex", booktabs = TRUE, caption = caption,
    digits = digits, row.names = FALSE
  )
  if (requireNamespace("kableExtra", quietly = TRUE)) {
    table <- kableExtra::kable_styling(
      table, latex_options = c("hold_position", "scale_down")
    )
  }
  atomic_text(as.character(table), path)
}

verify_full <- function() {
  required <- c(
    file.path(full_run, "COMPLETED"),
    file.path(full_run, "manifest-sha256.csv"),
    file.path(full_run, "numeric/oof-metrics.rds"),
    file.path(full_run, "numeric/oof-predictions.rds"),
    file.path(full_run, "numeric/oof-lisa-local.csv"),
    file.path(full_run, "numeric/oof-moran-global.csv"),
    file.path(full_run, "numeric/selection-full.rds"),
    file.path(full_run, "numeric/full-selection-summary.rds"),
    file.path(full_run, "numeric/gwr-full-local-inference.rds"),
    file.path(full_run, "numeric/gwr-full-f1-f2.csv"),
    file.path(full_run, "numeric/gwr-full-f3.csv")
  )
  stopifnot(all(file.exists(required)))
  completed <- readLines(required[[1L]], warn = FALSE)
  stopifnot(
    "Mode: full" %in% completed,
    any(grepl(
      "^Schema: gwrs-usa-counties-acs2024-study-v1-run-v1$",
      completed
    ))
  )
  command <- c(
    file.path(app, "code/run.py"), "verify", "--mode", "full",
    "--output", full_run
  )
  status <- system2("python3", shQuote(command), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
    stop("Full-run verification failed:\n", paste(status, collapse = "\n"))
  }
  invisible(required)
}

required_render_files <- c(
  "REPORT.md", "session-info.txt",
  "tables/oof-performance.csv", "tables/oof-performance-by-fold.csv",
  "tables/full-data-diagnostics.csv", "tables/full-tuning-curve.csv",
  "tables/full-bandwidth-selection.csv", "tables/outer-bandwidth-selection.csv",
  "tables/local-selection-by-variable.csv", "tables/local-selection-stability.csv",
  "tables/local-selection-oof-agreement.csv", "tables/global-selection.csv",
  "tables/oof-moran-global.csv", "tables/gwr-f1-f2.csv", "tables/gwr-f3.csv",
  "figures/oof-performance.png", "figures/full-tuning-curves.png",
  "figures/selected-bandwidths.png", "figures/oof-lisa-all-methods.png",
  "figures/gwr-local-inference-main.png", "figures/local-active-counts.png"
)

verify_render <- function(output) {
  verify_full()
  output <- normalizePath(output, mustWork = TRUE)
  manifest_path <- file.path(output, "manifest-sha256.csv")
  completion_path <- file.path(output, "COMPLETED")
  stopifnot(
    file.exists(manifest_path), file.exists(completion_path),
    all(file.exists(file.path(output, required_render_files)))
  )
  recorded <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  recorded$bytes <- as.numeric(recorded$bytes)
  files <- list.files(
    output, recursive = TRUE, full.names = TRUE,
    include.dirs = FALSE, all.files = TRUE, no.. = TRUE
  )
  files <- files[!basename(files) %in% c("manifest-sha256.csv", "COMPLETED")]
  stopifnot(identical(recorded, file_records(files, output)))
  completion <- readLines(completion_path, warn = FALSE)
  stopifnot(
    "Schema: gwrs-usa-counties-acs2024-study-v2-render-v1" %in% completion,
    paste0(
      "Source-Manifest-SHA256: ",
      hash_file(file.path(full_run, "manifest-sha256.csv"))
    ) %in% completion,
    paste0("Manifest-SHA256: ", hash_file(manifest_path)) %in% completion
  )
  cat("Study-v2 publication tables, diagnostics, maps and checksums: OK\n")
}

if (action == "verify") {
  verify_render(output_requested)
  quit(save = "no", status = 0L)
}

verify_full()
if (file.exists(output_requested)) stop("Preserving existing output: ", output_requested)
parent <- dirname(output_requested)
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
building <- tempfile(paste0(".", basename(output_requested), "-building-"), parent)
dir.create(building)
dir.create(file.path(building, "figures"))
dir.create(file.path(building, "tables"))

numeric <- file.path(full_run, "numeric")
cfg <- jsonlite::fromJSON(file.path(app, "config/numerics.json"))
policy <- jsonlite::fromJSON(file.path(app, "config/policy.json"))
model <- utils::read.csv(
  file.path(upstream, "frozen/input/spatial/data/county-model-data-v1.csv"),
  colClasses = c(fips = "character"), check.names = FALSE
)
counties <- sf::st_read(
  file.path(upstream, "frozen/input/spatial/data/analysis-counties-v1.gpkg"),
  quiet = TRUE
)
stopifnot(
  nrow(model) == cfg$n, nrow(counties) == cfg$n,
  identical(as.character(counties$fips), model$fips)
)

method_order <- c(
  "local-gwr", "local-ridge", "local-lasso", "local-en",
  "local-scad", "local-mcp", "global-ols", "global-ridge",
  "global-lasso", "global-en", "global-scad", "global-mcp"
)
local_methods <- method_order[grepl("^local-", method_order)]
local_selectors <- c("local-lasso", "local-en", "local-scad", "local-mcp")
global_selectors <- c("global-lasso", "global-en", "global-scad", "global-mcp")
main_variables <- cfg$postanalysis$main_map_predictors

read_payload <- function(task) {
  envelope <- readRDS(file.path(numeric, "tasks", paste0(task, ".rds")))
  stopifnot(identical(envelope$id, task), is.list(envelope$payload))
  envelope$payload
}

named_stat <- function(value, name) {
  answer <- unname(value[[name]])
  if (!length(answer)) NA_real_ else as.double(answer)
}

# Primary predictive assessment: pooled and fold-specific nested spatial OOF.
predictions <- readRDS(file.path(numeric, "oof-predictions.rds"))
metrics <- readRDS(file.path(numeric, "oof-metrics.rds"))
stopifnot(
  nrow(metrics) == length(method_order),
  setequal(metrics$method, method_order),
  nrow(predictions) == length(method_order) * cfg$n
)
metrics <- metrics[match(method_order, metrics$method), , drop = FALSE]
ols_rmse <- metrics$rmse[metrics$method == "global-ols"]
gwr_rmse <- metrics$rmse[metrics$method == "local-gwr"]
metrics$rmse_improvement_vs_global_ols_pct <-
  100 * (ols_rmse - metrics$rmse) / ols_rmse
metrics$rmse_improvement_vs_local_gwr_pct <-
  100 * (gwr_rmse - metrics$rmse) / gwr_rmse
metrics$rank_rmse <- rank(metrics$rmse, ties.method = "min")

fold_groups <- split(
  predictions,
  interaction(predictions$method, predictions$outer_fold, drop = TRUE)
)
fold_metrics <- do.call(rbind, lapply(fold_groups, function(value) {
  error <- value$prediction - value$observed
  sse <- sum(error^2)
  data.frame(
    method = value$method[[1L]], outer_fold = value$outer_fold[[1L]],
    n = nrow(value), rmse = sqrt(mean(error^2)),
    mae = mean(abs(error)), bias = mean(error),
    predictive_R2 = 1 - sse / sum((value$observed - value$baseline)^2),
    stringsAsFactors = FALSE
  )
}))
rownames(fold_metrics) <- NULL
fold_metrics <- fold_metrics[order(
  match(fold_metrics$method, method_order), fold_metrics$outer_fold
), , drop = FALSE]

atomic_csv(metrics, file.path(building, "tables/oof-performance.csv"))
atomic_csv(fold_metrics, file.path(building, "tables/oof-performance-by-fold.csv"))
write_latex(
  metrics, "Nested spatial out-of-fold predictive performance",
  file.path(building, "tables/oof-performance.tex"), 4L
)

performance_plot <- ggplot2::ggplot(
  metrics,
  ggplot2::aes(
    x = reorder(method, rmse), y = rmse,
    fill = ifelse(grepl("^local-", method), "Local", "Global")
  )
) +
  ggplot2::geom_col(width = 0.78) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c(Local = "#2166AC", Global = "#B2182B")) +
  ggplot2::labs(
    title = "Nested spatial out-of-fold RMSE",
    subtitle = "Lower is better; all methods use the same 3,107 held-out counties",
    x = NULL, y = "RMSE (log-income units)", fill = "Model scope"
  ) +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(
  file.path(building, "figures/oof-performance.png"), performance_plot,
  width = 9, height = 6.5, dpi = 320, bg = "white"
)
ggplot2::ggsave(
  file.path(building, "figures/oof-performance.pdf"), performance_plot,
  width = 9, height = 6.5, device = grDevices::cairo_pdf
)

# Full-data tuning curves and the selected q/lambda/alpha combinations.
selection_object <- readRDS(file.path(numeric, "selection-full.rds"))
full_selection <- readRDS(file.path(numeric, "full-selection-summary.rds"))
stopifnot(is.data.frame(selection_object$candidates), nrow(full_selection) == 12L)
candidate_groups <- split(selection_object$candidates, selection_object$candidates$id)
tuning <- do.call(rbind, lapply(candidate_groups, function(value) {
  eligible <- value$valid & is.finite(value$mse)
  stopifnot(any(eligible))
  selected <- which(eligible)[which.min(value$mse[eligible])]
  value[selected, , drop = FALSE]
}))
rownames(tuning) <- NULL
tuning$family <- ifelse(
  tuning$method == "en",
  paste0(tuning$where, "-en (alpha=", tuning$alpha, ")"),
  paste(tuning$where, tuning$method, sep = "-")
)
tuning$selected <- mapply(
  function(id, index) any(full_selection$id == id & full_selection$index == index),
  tuning$id, tuning$index
)
tuning$q_boundary <- with(
  tuning,
  where == "local" & (abs(q - min(cfg$q)) < 1e-14 | abs(q - 1) < 1e-14)
)
atomic_csv(tuning, file.path(building, "tables/full-tuning-curve.csv"))
atomic_csv(full_selection, file.path(building, "tables/full-selection.csv"))
write_latex(
  full_selection, "Full-data spatial-CV tuning selections",
  file.path(building, "tables/full-selection.tex"), 5L
)

local_tuning <- tuning[tuning$where == "local", , drop = FALSE]
tuning_plot <- ggplot2::ggplot(
  local_tuning,
  ggplot2::aes(x = q, y = mse, color = family, group = family)
) +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::geom_point(size = 1.4) +
  ggplot2::geom_point(
    data = local_tuning[local_tuning$selected, , drop = FALSE],
    shape = 21, size = 3.3, stroke = 0.8, fill = "white"
  ) +
  ggplot2::scale_x_log10(breaks = cfg$q) +
  ggplot2::facet_wrap(~family, scales = "free_y", ncol = 3) +
  ggplot2::labs(
    title = "Full-data spatial-CV bandwidth profiles",
    subtitle = "Circled points are selected; q=1 retains bisquare distance weights",
    x = "Requested training-neighborhood fraction q (log scale)",
    y = "Pooled validation MSE", color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(
  file.path(building, "figures/full-tuning-curves.png"), tuning_plot,
  width = 12, height = 9, dpi = 320, bg = "white"
)
ggplot2::ggsave(
  file.path(building, "figures/full-tuning-curves.pdf"), tuning_plot,
  width = 12, height = 9, device = grDevices::cairo_pdf
)

geometry_row <- function(value, scope, outer = NA_integer_) {
  candidate <- value$candidate[1L, , drop = FALSE]
  geometry <- value$geometry
  stopifnot(!is.null(geometry), identical(candidate$where, "local"))
  data.frame(
    scope = scope, outer_fold = outer,
    method = paste(candidate$where, candidate$method, sep = "-"),
    q = as.double(geometry$requested_q),
    training_n = as.integer(geometry$training_n),
    k = as.integer(geometry$realized_k),
    k_fraction = as.double(geometry$realized_fraction),
    bandwidth_min_km = named_stat(geometry$bandwidth_m, "minimum") / 1000,
    bandwidth_median_km = named_stat(geometry$bandwidth_m, "median") / 1000,
    bandwidth_max_km = named_stat(geometry$bandwidth_m, "maximum") / 1000,
    effective_n_min = named_stat(geometry$effective_weight_count, "minimum"),
    effective_n_median = named_stat(geometry$effective_weight_count, "median"),
    effective_n_max = named_stat(geometry$effective_weight_count, "maximum"),
    q_lower_boundary = abs(geometry$requested_q - min(cfg$q)) < 1e-14,
    q_upper_boundary = abs(geometry$requested_q - 1) < 1e-14,
    stringsAsFactors = FALSE
  )
}

full_bandwidth <- do.call(rbind, lapply(local_methods, function(method) {
  geometry_row(read_payload(paste0("full-refit-", method)), "full")
}))
outer_bandwidth <- do.call(rbind, lapply(seq_len(5L), function(outer) {
  do.call(rbind, lapply(local_methods, function(method) {
    geometry_row(
      read_payload(paste0("o", outer, "-refit-", method)),
      "outer_refit", outer
    )
  }))
}))
atomic_csv(full_bandwidth, file.path(building, "tables/full-bandwidth-selection.csv"))
atomic_csv(outer_bandwidth, file.path(building, "tables/outer-bandwidth-selection.csv"))
write_latex(
  full_bandwidth, "Selected full-data local bandwidth geometry",
  file.path(building, "tables/full-bandwidth-selection.tex"), 3L
)

bandwidth_plot <- ggplot2::ggplot(
  full_bandwidth,
  ggplot2::aes(x = reorder(method, q), y = q, fill = method)
) +
  ggplot2::geom_col(show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0("k=", k)), hjust = -0.08, size = 3.4
  ) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  ggplot2::labs(
    title = "Selected local bandwidths",
    subtitle = "q is converted separately within every training split",
    x = NULL, y = "Selected neighborhood fraction q"
  ) +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(
  file.path(building, "figures/selected-bandwidths.png"), bandwidth_plot,
  width = 9, height = 5.5, dpi = 320, bg = "white"
)
ggplot2::ggsave(
  file.path(building, "figures/selected-bandwidths.pdf"), bandwidth_plot,
  width = 9, height = 5.5, device = grDevices::cairo_pdf
)

# Local and global variable-selection summaries. True-support metrics are not
# reported because this is observational data and the support is unknown.
coefficient_rows <- list()
selection_rows <- list()
active_rows <- list()
full_selector_state <- list()
for (method in local_methods) {
  value <- read_payload(paste0("full-refit-", method))
  stopifnot(identical(value$test_ids, model$fips), nrow(value$coefficients) == cfg$n)
  slopes <- value$coefficients[, -1L, drop = FALSE]
  selected <- if (is.null(value$selected)) {
    matrix(TRUE, nrow(slopes), ncol(slopes), dimnames = dimnames(slopes))
  } else {
    value$selected
  }
  stopifnot(identical(dim(selected), dim(slopes)))
  if (method %in% local_selectors) full_selector_state[[method]] <- selected
  for (variable_index in seq_len(ncol(slopes))) {
    variable <- colnames(slopes)[[variable_index]]
    chosen <- selected[, variable_index]
    coefficient <- slopes[, variable_index]
    if (method %in% local_selectors) {
      active_coefficient <- coefficient[chosen]
      positive <- if (length(active_coefficient)) mean(active_coefficient > 0) else NA_real_
      negative <- if (length(active_coefficient)) mean(active_coefficient < 0) else NA_real_
      selection_rows[[length(selection_rows) + 1L]] <- data.frame(
        method = method, variable = variable,
        selected_n = sum(chosen), total_n = length(chosen),
        selected_fraction = mean(chosen), zero_fraction = mean(!chosen),
        positive_fraction_among_selected = positive,
        negative_fraction_among_selected = negative,
        sign_stability_among_selected = if (all(is.na(c(positive, negative)))) {
          NA_real_
        } else max(positive, negative, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    if (variable %in% main_variables) {
      coefficient_rows[[length(coefficient_rows) + 1L]] <- data.frame(
        fips = model$fips, method = method, variable = variable,
        coefficient = coefficient, selected = chosen,
        stringsAsFactors = FALSE
      )
    }
  }
  if (method %in% local_selectors) {
    active <- rowSums(selected)
    active_rows[[length(active_rows) + 1L]] <- data.frame(
      method = method, fips = model$fips, active_count = active,
      stringsAsFactors = FALSE
    )
  }
}
coefficient_map <- do.call(rbind, coefficient_rows)
local_selection <- do.call(rbind, selection_rows)
active_counts <- do.call(rbind, active_rows)
active_summary <- do.call(rbind, lapply(split(active_counts, active_counts$method), function(value) {
  quantile_value <- stats::quantile(value$active_count, c(0, 0.25, 0.5, 0.75, 1))
  data.frame(
    method = value$method[[1L]], n = nrow(value),
    minimum = unname(quantile_value[[1L]]), q1 = unname(quantile_value[[2L]]),
    median = unname(quantile_value[[3L]]), mean = mean(value$active_count),
    q3 = unname(quantile_value[[4L]]), maximum = unname(quantile_value[[5L]]),
    stringsAsFactors = FALSE
  )
}))
rownames(active_summary) <- NULL

fold_selection_rows <- list()
agreement_counts <- list()
for (outer in seq_len(5L)) {
  for (method in local_selectors) {
    value <- read_payload(paste0("o", outer, "-refit-", method))
    selected <- value$selected
    stopifnot(!is.null(selected), nrow(selected) == length(value$test_ids))
    full_selected <- full_selector_state[[method]][
      match(value$test_ids, model$fips), , drop = FALSE
    ]
    for (variable_index in seq_len(ncol(selected))) {
      variable <- colnames(value$coefficients)[variable_index + 1L]
      current <- selected[, variable_index]
      reference <- full_selected[, variable_index]
      fold_selection_rows[[length(fold_selection_rows) + 1L]] <- data.frame(
        method = method, outer_fold = outer, variable = variable,
        n = length(current), selected_fraction = mean(current),
        stringsAsFactors = FALSE
      )
      key <- paste(method, variable, sep = "::")
      old <- agreement_counts[[key]]
      if (is.null(old)) old <- c(n = 0, agree = 0, both = 0, union = 0)
      agreement_counts[[key]] <- old + c(
        n = length(current), agree = sum(current == reference),
        both = sum(current & reference), union = sum(current | reference)
      )
    }
  }
}
fold_selection <- do.call(rbind, fold_selection_rows)
stability <- do.call(rbind, lapply(
  split(fold_selection, interaction(
    fold_selection$method, fold_selection$variable, drop = TRUE
  )),
  function(value) data.frame(
    method = value$method[[1L]], variable = value$variable[[1L]],
    mean_selected_fraction = mean(value$selected_fraction),
    sd_selected_fraction = stats::sd(value$selected_fraction),
    min_selected_fraction = min(value$selected_fraction),
    max_selected_fraction = max(value$selected_fraction),
    range_selected_fraction = diff(range(value$selected_fraction)),
    stringsAsFactors = FALSE
  )
))
rownames(stability) <- NULL
agreement <- do.call(rbind, lapply(names(agreement_counts), function(key) {
  parts <- strsplit(key, "::", fixed = TRUE)[[1L]]
  value <- agreement_counts[[key]]
  data.frame(
    method = parts[[1L]], variable = parts[[2L]], n = value[["n"]],
    agreement_fraction = value[["agree"]] / value[["n"]],
    selected_jaccard = if (value[["union"]] == 0) 1 else
      value[["both"]] / value[["union"]],
    stringsAsFactors = FALSE
  )
}))
rownames(agreement) <- NULL

global_selection_rows <- list()
for (method in global_selectors) {
  value <- read_payload(paste0("full-refit-", method))
  selected <- as.logical(value$selected[1L, ])
  coefficient <- as.double(value$coefficients[1L, -1L])
  global_selection_rows[[length(global_selection_rows) + 1L]] <- data.frame(
    method = method, variable = colnames(value$coefficients)[-1L],
    selected = selected, coefficient = coefficient,
    stringsAsFactors = FALSE
  )
}
global_selection <- do.call(rbind, global_selection_rows)

atomic_csv(local_selection, file.path(building, "tables/local-selection-by-variable.csv"))
atomic_csv(active_summary, file.path(building, "tables/local-active-count-summary.csv"))
atomic_csv(fold_selection, file.path(building, "tables/local-selection-by-fold.csv"))
atomic_csv(stability, file.path(building, "tables/local-selection-stability.csv"))
atomic_csv(agreement, file.path(building, "tables/local-selection-oof-agreement.csv"))
atomic_csv(global_selection, file.path(building, "tables/global-selection.csv"))
atomic_csv(coefficient_map, file.path(building, "tables/main-coefficient-surfaces.csv"))
atomic_csv(
  local_selection[local_selection$variable %in% main_variables, , drop = FALSE],
  file.path(building, "tables/main-selection-rates.csv")
)
write_latex(
  active_summary, "Distribution of active local predictors by selector",
  file.path(building, "tables/local-active-count-summary.tex"), 3L
)

coefficient_map$method <- factor(coefficient_map$method, levels = local_methods)
for (variable in main_variables) {
  values <- coefficient_map[coefficient_map$variable == variable, , drop = FALSE]
  geometry <- counties[match(values$fips, counties$fips), "geom"]
  geometry$method <- values$method
  geometry$coefficient <- values$coefficient
  geometry$selected <- values$selected
  shown <- geometry[geometry$selected, , drop = FALSE]
  plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = geometry, fill = "grey80", color = NA) +
    ggplot2::geom_sf(data = shown, ggplot2::aes(fill = coefficient), color = NA) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, name = "Coefficient"
    ) +
    ggplot2::facet_wrap(~method, ncol = 3) +
    ggplot2::coord_sf(datum = NA) +
    ggplot2::labs(
      title = paste("Local coefficient surface:", variable),
      subtitle = paste(
        "Gray is an exact zero for selector methods;",
        "GWR and ridge are shown everywhere"
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  stem <- paste0("coefficient-", gsub("_", "-", variable))
  ggplot2::ggsave(
    file.path(building, "figures", paste0(stem, ".png")), plot,
    width = 12, height = 7.5, dpi = 320, bg = "white"
  )
  ggplot2::ggsave(
    file.path(building, "figures", paste0(stem, ".pdf")), plot,
    width = 12, height = 7.5, device = grDevices::cairo_pdf
  )
}

active_geometry <- counties[match(active_counts$fips, counties$fips), "geom"]
active_geometry$method <- factor(active_counts$method, levels = local_selectors)
active_geometry$active_count <- active_counts$active_count
active_plot <- ggplot2::ggplot(active_geometry) +
  ggplot2::geom_sf(ggplot2::aes(fill = active_count), color = NA) +
  ggplot2::scale_fill_viridis_c(option = "C", name = "Active") +
  ggplot2::facet_wrap(~method, ncol = 2) +
  ggplot2::coord_sf(datum = NA) +
  ggplot2::labs(
    title = "Number of locally selected predictors",
    subtitle = "Counts range from zero to 26 at each county",
    x = NULL, y = NULL
  ) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(
  file.path(building, "figures/local-active-counts.png"), active_plot,
  width = 11, height = 7.5, dpi = 320, bg = "white"
)
ggplot2::ggsave(
  file.path(building, "figures/local-active-counts.pdf"), active_plot,
  width = 11, height = 7.5, device = grDevices::cairo_pdf
)

# Conditional full-data fit diagnostics. Local rows come from the gwrs
# active-set derivative. Global rows use the corresponding centered-design
# derivative, with nonconvex information criteria withheld if local convexity
# fails at the selected stationary point.
diagnostic_fields <- c(
  "n", "rss", "tss", "rmse", "mae", "bias", "r_squared",
  "adjusted_r_squared", "trS", "trStS", "enp", "edf", "sigma2",
  "deviance", "logLik", "model_df", "aic", "aicc", "bic", "gcv", "loo_cv"
)

global_active_diagnostics <- function(value, method) {
  y <- as.double(value$observed)
  prediction <- as.double(value$predictions)
  x <- as.matrix(model[, colnames(value$coefficients)[-1L], drop = FALSE])
  center <- colMeans(x)
  z <- sweep(x, 2L, center, "-")
  scale <- sqrt(colMeans(z^2))
  stopifnot(all(is.finite(scale)), all(scale > 0))
  z <- sweep(z, 2L, scale, "/")
  beta_working <- as.double(value$coefficients[1L, -1L]) * scale
  selected <- if (is.null(value$selected)) {
    rep(TRUE, ncol(z))
  } else as.logical(value$selected[1L, ])
  active <- which(selected)
  lambda <- as.double(tail(value$score$lambda, 1L))
  alpha <- as.double(value$candidate$alpha[[1L]])
  curvature <- rep(0, length(active))
  if (length(active) && grepl("scad$", method)) {
    absolute <- abs(beta_working[active])
    curved <- absolute > lambda & absolute < cfg$scad_gamma * lambda
    curvature[curved] <- -1 / (cfg$scad_gamma - 1)
  } else if (length(active) && grepl("mcp$", method)) {
    absolute <- abs(beta_working[active])
    curved <- absolute > 0 & absolute < cfg$mcp_gamma * lambda
    curvature[curved] <- -1 / cfg$mcp_gamma
  } else if (length(active) && grepl("ridge$|en$", method)) {
    curvature[] <- lambda * (1 - alpha)
  }
  if (length(active)) {
    gram <- crossprod(z[, active, drop = FALSE]) / nrow(z)
    hessian <- gram + diag(curvature, nrow = length(active))
    minimum_eigenvalue <- min(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values)
    valid_derivative <- is.finite(minimum_eigenvalue) && minimum_eigenvalue > 1e-10
  } else {
    gram <- matrix(numeric(), 0L, 0L)
    hessian <- gram
    minimum_eigenvalue <- Inf
    valid_derivative <- TRUE
  }
  tr_s <- tr_st_s <- NA_real_
  hat_diag <- rep(NA_real_, nrow(z))
  if (valid_derivative) {
    if (length(active)) {
      inverse <- solve(hessian)
      influence <- inverse %*% gram
      tr_s <- 1 + sum(diag(influence))
      tr_st_s <- 1 + sum(diag(inverse %*% gram %*% inverse %*% gram))
      hat_diag <- 1 / nrow(z) + rowSums(
        (z[, active, drop = FALSE] %*% inverse) * z[, active, drop = FALSE]
      ) / nrow(z)
    } else {
      tr_s <- tr_st_s <- 1
      hat_diag[] <- 1 / nrow(z)
    }
  }
  residual <- y - prediction
  n <- length(y)
  rss <- sum(residual^2)
  tss <- sum((y - mean(y))^2)
  enp <- if (valid_derivative) 2 * tr_s - tr_st_s else NA_real_
  edf <- n - enp
  model_df <- tr_s + 1
  deviance <- n * (log(2 * pi) + 1 + log(rss / n))
  ic_valid <- valid_derivative && is.finite(edf) && edf > 0 &&
    is.finite(model_df) && n - model_df - 1 > 0 && all(hat_diag < 1)
  answer <- list(
    n = n, rss = rss, tss = tss, rmse = sqrt(rss / n),
    mae = mean(abs(residual)), bias = mean(prediction - y),
    r_squared = 1 - rss / tss,
    adjusted_r_squared = if (ic_valid) 1 - (rss / edf) / (tss / (n - 1)) else NA_real_,
    trS = tr_s, trStS = tr_st_s, enp = enp, edf = edf,
    sigma2 = if (ic_valid) rss / edf else NA_real_,
    deviance = deviance, logLik = -deviance / 2,
    model_df = model_df,
    aic = if (ic_valid) deviance + 2 * model_df else NA_real_,
    aicc = if (ic_valid) deviance + 2 * model_df +
      2 * model_df * (model_df + 1) / (n - model_df - 1) else NA_real_,
    bic = if (ic_valid) deviance + log(n) * model_df else NA_real_,
    gcv = if (ic_valid) n * rss / (n - tr_s)^2 else NA_real_,
    loo_cv = if (ic_valid) sum((residual / (1 - hat_diag))^2) else NA_real_
  )
  sparse <- method %in% global_selectors
  ebic <- if (ic_valid && sparse) {
    answer$bic + 2 * 0.5 * lchoose(ncol(z), length(active))
  } else NA_real_
  list(
    values = answer, information_criteria_valid = ic_valid,
    minimum_active_hessian_eigenvalue = minimum_eigenvalue,
    locally_convex_fraction = if (grepl("scad$|mcp$", method)) {
      as.numeric(valid_derivative)
    } else NA_real_,
    active_count = length(active), ebic_0.5 = ebic
  )
}

diagnostic_rows <- list()
for (method in method_order) {
  value <- read_payload(paste0("full-refit-", method))
  candidate <- value$candidate[1L, , drop = FALSE]
  if (grepl("^local-", method)) {
    diagnostic <- if (method == "local-gwr") {
      value$full_gwr$model_diagnostics
    } else value$selected_model_diagnostics
    stopifnot(is.list(diagnostic), is.list(diagnostic$global))
    global <- diagnostic$global
    selected_count <- if (is.null(value$selected)) cfg$p else mean(rowSums(value$selected))
    minimum_eigenvalue <- if (
      "minimum_active_hessian_eigenvalue" %in% names(diagnostic$local)
    ) min(diagnostic$local$minimum_active_hessian_eigenvalue) else NA_real_
    convex_fraction <- diagnostic$validity$locally_convex_fraction
    if (is.null(convex_fraction)) convex_fraction <- NA_real_
    selected_k <- diagnostic$selected_k
    if (is.null(selected_k)) selected_k <- candidate$k
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
      method = method, scope = "local", q = as.double(candidate$q),
      k = as.integer(selected_k),
      alpha = as.double(candidate$alpha),
      selected_lambda = as.double(tail(value$score$lambda, 1L)),
      active_count = selected_count,
      information_criteria_valid = isTRUE(diagnostic$validity$information_criteria),
      exact_linear = isTRUE(diagnostic$validity$exact_linear),
      conditional = isTRUE(diagnostic$validity$conditional),
      minimum_active_hessian_eigenvalue = minimum_eigenvalue,
      locally_convex_fraction = as.double(convex_fraction),
      ebic_0.5 = NA_real_,
      diagnostic_basis = diagnostic$validity$basis,
      as.data.frame(as.list(global[diagnostic_fields]), check.names = FALSE),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  } else {
    derived <- global_active_diagnostics(value, method)
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
      method = method, scope = "global", q = NA_real_, k = 0L,
      alpha = as.double(candidate$alpha),
      selected_lambda = as.double(tail(value$score$lambda, 1L)),
      active_count = derived$active_count,
      information_criteria_valid = derived$information_criteria_valid,
      exact_linear = method %in% c("global-ols", "global-ridge"),
      conditional = method != "global-ols",
      minimum_active_hessian_eigenvalue =
        derived$minimum_active_hessian_eigenvalue,
      locally_convex_fraction = derived$locally_convex_fraction,
      ebic_0.5 = derived$ebic_0.5,
      diagnostic_basis = "centered full-data active-set derivative",
      as.data.frame(derived$values, check.names = FALSE),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
}
full_diagnostics <- do.call(rbind, diagnostic_rows)
rownames(full_diagnostics) <- NULL
atomic_csv(full_diagnostics, file.path(building, "tables/full-data-diagnostics.csv"))
write_latex(
  full_diagnostics[, c(
    "method", "q", "k", "active_count", "rmse", "r_squared",
    "adjusted_r_squared", "enp", "aic", "aicc", "bic", "gcv"
  )],
  "Conditional descriptive full-data fit diagnostics",
  file.path(building, "tables/full-data-diagnostics.tex"), 3L
)

# Global and local residual spatial diagnostics.
moran <- utils::read.csv(
  file.path(numeric, "oof-moran-global.csv"), stringsAsFactors = FALSE
)
lisa <- utils::read.csv(
  file.path(numeric, "oof-lisa-local.csv"),
  stringsAsFactors = FALSE, colClasses = c(sample_id = "character")
)
stopifnot(nrow(moran) == 12L, nrow(lisa) == 12L * cfg$n)
atomic_csv(moran, file.path(building, "tables/oof-moran-global.csv"))
write_latex(
  moran, "Out-of-fold residual Moran diagnostics",
  file.path(building, "tables/oof-moran-global.tex"), 4L
)

lisa$method <- factor(lisa$method, levels = method_order)
lisa_geometry <- counties[match(lisa$sample_id, counties$fips), "geom"]
lisa_geometry$method <- lisa$method
lisa_geometry$cluster <- factor(
  lisa$cluster_BH_005,
  levels = c("Not Significant", "High-High", "Low-Low", "High-Low", "Low-High")
)
lisa_plot <- ggplot2::ggplot(lisa_geometry) +
  ggplot2::geom_sf(ggplot2::aes(fill = cluster), color = NA) +
  ggplot2::scale_fill_manual(values = c(
    "Not Significant" = "grey85", "High-High" = "#B2182B",
    "Low-Low" = "#2166AC", "High-Low" = "#EF8A62",
    "Low-High" = "#67A9CF"
  ), drop = FALSE, name = "BH LISA cluster") +
  ggplot2::facet_wrap(~method, ncol = 4) +
  ggplot2::coord_sf(datum = NA) +
  ggplot2::labs(
    title = "OOF residual LISA on the common queen graph",
    subtitle = "9,999 common-seed permutations; BH adjustment within method",
    x = NULL, y = NULL
  ) +
  ggplot2::theme_void(base_size = 9) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
  file.path(building, "figures/oof-lisa-all-methods.png"), lisa_plot,
  width = 14, height = 9, dpi = 320, bg = "white"
)
ggplot2::ggsave(
  file.path(building, "figures/oof-lisa-all-methods.pdf"), lisa_plot,
  width = 14, height = 9, device = grDevices::cairo_pdf
)

# GWR F1/F2/F3 and local inference. F3 multiplicity is across all 27 surfaces.
f1_f2 <- utils::read.csv(
  file.path(numeric, "gwr-full-f1-f2.csv"), stringsAsFactors = FALSE
)
f3 <- utils::read.csv(
  file.path(numeric, "gwr-full-f3.csv"), stringsAsFactors = FALSE
)
f1_f2$reject_raw_005 <- f1_f2$p_value < 0.05
f3$p_value_BH <- stats::p.adjust(f3$p_value, method = "BH")
f3$reject_raw_005 <- f3$p_value < 0.05
f3$reject_BH_005 <- f3$p_value_BH < 0.05
atomic_csv(f1_f2, file.path(building, "tables/gwr-f1-f2.csv"))
atomic_csv(f3, file.path(building, "tables/gwr-f3.csv"))
write_latex(
  f1_f2, "GWR global F1 and F2 diagnostics",
  file.path(building, "tables/gwr-f1-f2.tex"), 4L
)
write_latex(
  f3, "Coefficient-specific GWR F3 diagnostics with BH adjustment",
  file.path(building, "tables/gwr-f3.tex"), 4L
)

inference <- readRDS(file.path(numeric, "gwr-full-local-inference.rds"))
gwr_local <- inference$local
gwr_local$sample_id <- model$fips[gwr_local$location]
gwr_main <- gwr_local[gwr_local$coefficient %in% main_variables, , drop = FALSE]
gwr_geometry <- counties[match(gwr_main$sample_id, counties$fips), "geom"]
gwr_geometry$coefficient <- factor(gwr_main$coefficient, levels = main_variables)
gwr_geometry$estimate <- gwr_main$estimate
gwr_geometry$significant <- gwr_main$significant
gwr_significant <- gwr_geometry[gwr_geometry$significant, , drop = FALSE]
gwr_plot <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = gwr_geometry, fill = "grey94", color = NA) +
  ggplot2::geom_sf(
    data = gwr_significant, ggplot2::aes(fill = estimate), color = NA
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, name = "GWR coefficient"
  ) +
  ggplot2::facet_wrap(~coefficient, ncol = 2) +
  ggplot2::coord_sf(datum = NA) +
  ggplot2::labs(
    title = "Unpenalized GWR local inference",
    subtitle = "Pale gray is BH-adjusted nonsignificance within a coefficient surface",
    x = NULL, y = NULL
  ) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(legend.position = "right")
ggplot2::ggsave(
  file.path(building, "figures/gwr-local-inference-main.png"), gwr_plot,
  width = 10, height = 7.5, dpi = 320, bg = "white"
)
ggplot2::ggsave(
  file.path(building, "figures/gwr-local-inference-main.pdf"), gwr_plot,
  width = 10, height = 7.5, device = grDevices::cairo_pdf
)

# A concise, self-contained interpretation layer. It deliberately keeps OOF
# prediction separate from conditional full-data description.
best <- metrics[which.min(metrics$rmse), , drop = FALSE]
ols <- metrics[metrics$method == "global-ols", , drop = FALSE]
gwr <- metrics[metrics$method == "local-gwr", , drop = FALSE]
f1 <- f1_f2[f1_f2$test == "F1", , drop = FALSE]
f2 <- f1_f2[f1_f2$test == "F2", , drop = FALSE]
selected_bandwidth_text <- paste(
  paste0(full_bandwidth$method, ": q=", signif(full_bandwidth$q, 4),
         " (k=", full_bandwidth$k, ")"),
  collapse = "; "
)
boundary_methods <- full_bandwidth$method[
  full_bandwidth$q_lower_boundary | full_bandwidth$q_upper_boundary
]
report <- c(
  "# ACS 2024 county study-v2: full-range local-bandwidth results", "",
  "## Primary predictive evidence", "",
  paste0(
    "All twelve methods are compared on the same 3,107 nested spatial ",
    "out-of-fold predictions. The lowest RMSE is **", best$method,
    "** (", signif(best$rmse, 5), "), with predictive R-squared ",
    signif(best$predictive_R2, 5), "."
  ),
  paste0(
    "Global OLS has RMSE ", signif(ols$rmse, 5),
    " and local GWR has RMSE ", signif(gwr$rmse, 5),
    ". Positive values in the two improvement columns mean lower RMSE than ",
    "the named reference."
  ), "",
  "## Bandwidth audit", "",
  selected_bandwidth_text,
  if (length(boundary_methods)) paste0(
    "A q-grid boundary was selected for: ", paste(boundary_methods, collapse = ", "),
    ". This is reported as a boundary result rather than silently extending the grid."
  ) else "No local method selected either endpoint of the approved q grid.",
  "The q=1 candidate uses all training observations but retains bisquare distance weights; it is not global OLS.",
  "", "## Heterogeneity and spatial residual evidence", "",
  paste0(
    "GWR F1 p-value: ", signif(f1$p_value, 4),
    "; GWR F2 p-value: ", signif(f2$p_value, 4),
    ". After BH correction across all 27 F3 tests, ", sum(f3$reject_BH_005),
    " coefficient surfaces reject stationarity at 0.05."
  ),
  "F1/F2/F3 address conditional in-sample GWR heterogeneity; they do not establish superior out-of-fold prediction. Moran and LISA use only nested spatial OOF residuals.",
  "", "## Variable selection", "",
  "Selection rates, fold-wise rates, sign stability, and OOF-versus-full selection agreement are reported. Gray map regions are exact zeros under the frozen threshold, not nonsignificant estimates.",
  "MCC, F1 score, sensitivity, and specificity are not reported for this observational application because the true active set is unknown.",
  "", "## Full-data fit criteria", "",
  "RSS, RMSE, MAE, R-squared, adjusted R-squared, ENP, EDF, AIC, AICc, BIC, GCV, and leave-one-out criteria are descriptive full-data diagnostics. Penalized rows condition on the spatial-CV-selected hyperparameters and active sets. Nonconvex information criteria are withheld whenever the selected active Hessian is not locally convex.",
  "These conditional criteria are not used to replace the primary nested spatial OOF ranking. EBIC(0.5) is supplied only for global sparse fits where a single active-set size has a conventional interpretation.",
  "", "## Reproducibility", "",
  paste0("Source run: `", full_run, "`"),
  paste0(
    "Source manifest SHA-256: `",
    hash_file(file.path(full_run, "manifest-sha256.csv")), "`"
  )
)
atomic_text(report, file.path(building, "REPORT.md"))

session <- c(
  paste0("created_utc: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  paste0("source_run: ", full_run),
  paste0(
    "source_manifest_sha256: ",
    hash_file(file.path(full_run, "manifest-sha256.csv"))
  ),
  capture.output(sessionInfo()), capture.output(.libPaths())
)
atomic_text(session, file.path(building, "session-info.txt"))
files <- list.files(
  building, recursive = TRUE, full.names = TRUE,
  include.dirs = FALSE, all.files = TRUE, no.. = TRUE
)
manifest <- file_records(files, building)
atomic_csv(manifest, file.path(building, "manifest-sha256.csv"))
atomic_text(c(
  "Schema: gwrs-usa-counties-acs2024-study-v2-render-v1",
  paste0(
    "Source-Manifest-SHA256: ",
    hash_file(file.path(full_run, "manifest-sha256.csv"))
  ),
  paste0("Manifest-SHA256: ", hash_file(file.path(building, "manifest-sha256.csv")))
), file.path(building, "COMPLETED"))
if (!file.rename(building, output_requested)) {
  stop("Could not atomically publish the rendering directory")
}
verify_render(output_requested)
cat("Study-v2 publication rendering complete: ", output_requested, "\n", sep = "")
