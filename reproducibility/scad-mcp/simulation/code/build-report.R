# Render a validated completed simulation; no model fitting.
# Usage: Rscript --vanilla code/build-report.R ANALYSIS_ROOT
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1L]], mustWork = TRUE)
source_root <- file.path(root, "study")
summary_root <- file.path(root, "summary-v1")
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
for (package in c("ggplot2", "knitr", "ragg", "svglite")) {
  stopifnot(requireNamespace(package, quietly = TRUE))
}
library(ggplot2)
source(file.path(dirname(script_path), "resume-utils.R"))
source(file.path(dirname(script_path), "summarize-penalized-selection-study.R"))

stopifnot(gwrs_selection_summary_verify_output(summary_root))
completion <- read.dcf(file.path(summary_root, "COMPLETED"))
stopifnot(identical(unname(completion[1, "Manifest-SHA256"]),
                    gwrs_sim_sha256(file.path(summary_root, "manifest-sha256.csv"))))
read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE,
                                  check.names = FALSE, na.strings = c("", "NA"))
snapshot <- function(directory) {
  files <- sort(list.files(directory, recursive = TRUE, all.files = TRUE,
                          no.. = TRUE, include.dirs = FALSE))
  paths <- file.path(directory, files)
  data.frame(file = files, bytes = as.double(file.info(paths)$size),
             sha256 = unname(gwrs_sim_sha256(paths)), stringsAsFactors = FALSE)
}
source_before <- snapshot(source_root)
summary_before <- snapshot(summary_root)
metrics <- read_csv(file.path(summary_root, "metric-summary.csv"))
pairs <- read_csv(file.path(summary_root, "paired-differences.csv"))
convergence <- read_csv(file.path(summary_root, "convergence-summary.csv"))
design <- read_csv(file.path(summary_root, "design-summary.csv"))
bandwidth <- read_csv(file.path(summary_root, "bandwidth-summary.csv"))
raw <- read_csv(file.path(root, "study/task-summary.csv"))
stopifnot(nrow(raw) == 15200L, nrow(design) == 39L,
          sum(design$completed_tasks) == 1900L,
          all(design$planned_tasks == design$completed_tasks))
sections <- c("main", "reference", "stress")
methods <- gwrs_selection_summary_methods
method_labels <- c(gwr = "GWR", ridge = "Ridge", lasso = "Lasso",
                   elastic_net = "Elastic Net", scad = "SCAD", mcp = "MCP",
                   oracle_union_gwr = "Oracle union", oracle_local_gwr = "Oracle local")
primary <- c("new_location_signal_rmse", "active_coefficient_rmse", "mcc")
metric_labels <- c(new_location_signal_rmse = "New-location signal RMSE\nLower is better",
                   active_coefficient_rmse = "Active-coefficient RMSE\nLower is better",
                   mcc = "Selection MCC\nHigher is better")
design <- design[order(match(design$section, sections), design$scenario_id,
                       design$n, design$p, design$predictor_rho,
                       design$spatial_rho, design$target_snr), ]
design$cell_code <- unlist(lapply(sections, function(section) {
  paste0(c(main = "M", reference = "R", stress = "S")[[section]],
         sprintf("%02d", seq_len(sum(design$section == section))))
}))
code <- setNames(design$cell_code, design$cell_id)
add_code <- function(data) {
  data$cell_code <- unname(code[data$cell_id])
  stopifnot(!anyNA(data$cell_code))
  data
}
metrics <- add_code(metrics)
pairs <- add_code(pairs)
convergence <- add_code(convergence)
bandwidth <- add_code(bandwidth)
raw <- add_code(raw)
primary_metrics <- metrics[metrics$metric %in% primary, ]
direct <- pairs[pairs$metric %in% primary &
                  pairs$target_method %in% c("scad", "mcp") &
                  pairs$reference_method == "lasso", ]
stopifnot(nrow(primary_metrics) == 39L * 8L * 3L * 2L,
          nrow(direct) == 39L * 2L * 3L * 2L)

# Reconcile every primary summary and direct paired row against task-level data.
# This checks existing calculations; it introduces no new estimand or exclusion.
groups <- split(raw, paste(raw$cell_id, raw$method, sep = "|"))
equal_numbers <- function(left, right) {
  isTRUE(all.equal(as.double(left), as.double(right), tolerance = 1e-11))
}
for (i in seq_len(nrow(primary_metrics))) {
  row <- primary_metrics[i, ]
  observed <- groups[[paste(row$cell_id, row$method, sep = "|")]]
  if (row$analysis_set == "all_converged") observed <- observed[observed$all_converged, ]
  check <- gwrs_selection_summary_stats(observed[[row$metric]])
  stopifnot(all(vapply(names(check), function(name) equal_numbers(row[[name]], check[[name]]),
                      logical(1L))))
}
for (i in seq_len(nrow(direct))) {
  row <- direct[i, ]
  target <- groups[[paste(row$cell_id, row$target_method, sep = "|")]]
  reference <- groups[[paste(row$cell_id, row$reference_method, sep = "|")]]
  reference <- reference[match(target$task_id, reference$task_id), ]
  stopifnot(identical(target$task_id, reference$task_id))
  keep <- if (row$analysis_set == "all_converged") {
    target$all_converged & reference$all_converged
  } else rep(TRUE, nrow(target))
  check <- gwrs_selection_summary_pair_stats(
    target[[row$metric]][keep] - reference[[row$metric]][keep])
  stopifnot(all(vapply(names(check), function(name) equal_numbers(row[[name]], check[[name]]),
                      logical(1L))))
}
cat("Primary rows reconciled:", nrow(primary_metrics), "Paired rows:", nrow(direct), "\n")

output <- file.path(root, "report", "v1")
output_link <- Sys.readlink(output)
stopifnot(!file.exists(output), !dir.exists(output),
          is.na(output_link) || !nzchar(output_link))
dir.create(dirname(output), showWarnings = FALSE)
staging <- tempfile(".report-v1-", tmpdir = dirname(output))
stopifnot(dir.create(staging))
for (directory in c("tables", "figures", "validation", "code")) {
  stopifnot(dir.create(file.path(staging, directory)))
}
# On any failure leave the clearly named staging directory for inspection.
write_csv <- function(data, name, directory = "tables") {
  write.csv(data, file.path(staging, directory, name), row.names = FALSE, na = "")
}
fmt <- function(x, digits = 4L) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
write_table <- function(data, filename, caption) {
  table <- knitr::kable(data, format = "latex", booktabs = TRUE, longtable = TRUE,
                       row.names = FALSE, caption = caption, escape = TRUE)
  writeLines(as.character(table), file.path(staging, "tables", paste0(filename, ".tex")))
}
write_csv(design, "design-key.csv")
write_csv(primary_metrics, "primary-metrics.csv")
write_csv(direct, "scad-mcp-vs-lasso.csv")
write_csv(pairs[pairs$metric %in% primary, ], "all-prespecified-primary-comparisons.csv")
write_csv(convergence, "convergence-and-timing-by-cell.csv")
write_csv(bandwidth, "common-bandwidth.csv")

decisions <- do.call(rbind, lapply(split(direct, interaction(
    direct$section, direct$analysis_set, direct$target_method, direct$metric, drop = TRUE)),
  function(data) data.frame(
    section = data$section[1], analysis_set = data$analysis_set[1],
    target_method = data$target_method[1], reference_method = "lasso", metric = data$metric[1],
    cells = nrow(data), target_better = sum(data$favors == "target"),
    lasso_better = sum(data$favors == "reference"),
    inconclusive = sum(data$favors == "inconclusive"),
    not_estimable = sum(data$favors == "not_estimable"))))
stopifnot(all(rowSums(decisions[c("target_better", "lasso_better", "inconclusive", "not_estimable")]) ==
                decisions$cells))
write_csv(decisions, "paired-decision-counts.csv")
warning_rows <- raw[!raw$all_converged, c("task_id", "cell_code", "cell_id", "section", "method",
                                        "replication", "convergence_rate", "cv_convergence_rate",
                                        "maximum_stationarity_violation", "maximum_coordinate_gap")]
write_csv(warning_rows, "convergence-warnings.csv")
method_convergence <- do.call(rbind, lapply(methods, function(method) {
  data <- raw[raw$method == method, ]
  data.frame(method = method, fits = nrow(data), all_converged = sum(data$all_converged),
             with_warning = sum(!data$all_converged),
             mean_location_convergence_rate = mean(data$convergence_rate))
}))
write_csv(method_convergence, "method-convergence.csv")
write_table(method_convergence, "method-convergence", "Final-fit convergence; warnings are retained in primary analyses.")
key <- c("cell_id", "target_method", "reference_method", "metric")
main_pairs <- direct[direct$analysis_set == "all_completed", ]
sens_pairs <- direct[direct$analysis_set == "all_converged", ]
sensitivity <- merge(main_pairs, sens_pairs, by = key, suffixes = c("_primary", "_sensitivity"))
sensitivity$decision_changed <- sensitivity$favors_primary != sensitivity$favors_sensitivity
sensitivity$paired_fits_excluded <- sensitivity$n_total_primary - sensitivity$n_total_sensitivity
stopifnot(all(sensitivity$paired_fits_excluded >= 0))
write_csv(sensitivity, "convergence-sensitivity-comparison.csv")
write_csv(main_pairs[main_pairs$favors == "reference", ], "lasso-favored-primary-cells.csv")
null <- metrics[metrics$scenario_id == "null_slopes" & metrics$analysis_set == "all_completed" &
                  metrics$metric %in% c("false_positive_rate", "specificity", "false_discovery_rate",
                                        "mean_selected_predictors", primary), ]
write_csv(null, "null-slope-reference.csv")

for (section in sections) {
  d <- design[design$section == section, ]
  write_table(d[c("cell_code", "scenario_id", "n", "p", "predictor_rho", "spatial_rho", "target_snr",
                  "planned_tasks")], paste0("design-", section), paste("Design key:", section))
  base <- primary_metrics[primary_metrics$section == section &
                            primary_metrics$analysis_set == "all_completed", ]
  base <- base[order(match(base$cell_code, d$cell_code), match(base$method, methods)), ]
  ids <- unique(base[c("cell_code", "method")])
  printable <- data.frame(Cell = ids$cell_code, Method = unname(method_labels[ids$method]))
  counts <- character(nrow(ids))
  for (metric in primary) {
    one <- base[base$metric == metric, ]
    stopifnot(identical(one$cell_code, ids$cell_code), identical(one$method, ids$method))
    printable[[c(new_location_signal_rmse = "Signal RMSE (MCSE)",
                  active_coefficient_rmse = "Active RMSE (MCSE)", mcc = "MCC (MCSE)")[[metric]]]] <-
      ifelse(is.finite(one$mean), paste0(fmt(one$mean), " (", fmt(one$mcse), ")"), "NA")
    counts <- if (metric == primary[1]) as.character(one$n_finite) else paste0(counts, "/", one$n_finite)
  }
  printable[["Finite r: signal/active/MCC"]] <- counts
  write_table(printable, paste0("primary-", section),
              paste("Primary all-completed metrics:", section, "; MCSE in parentheses; no cross-cell pooling."))
  decision_table <- decisions[decisions$section == section & decisions$analysis_set == "all_completed",
                              c("target_method", "metric", "target_better", "lasso_better", "inconclusive", "not_estimable")]
  write_table(decision_table, paste0("paired-counts-", section),
              paste("Descriptive cell counts from unadjusted 95 percent paired Monte Carlo intervals:", section))
}

theme_set(theme_minimal(base_size = 11, base_family = "sans") +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
                  plot.caption = element_text(hjust = 0, size = 9),
                  legend.position = "bottom", strip.text = element_text(face = "bold")))
figure_index <- list()
save_figure <- function(plot, name, width, height, source_file, note) {
  path <- file.path(staging, "figures", name)
  ggsave(paste0(path, ".png"), plot, width = width, height = height, dpi = 180,
         device = ragg::agg_png, bg = "white", limitsize = FALSE)
  ggsave(paste0(path, ".svg"), plot, width = width, height = height,
         device = svglite::svglite, bg = "white", limitsize = FALSE)
  figure_index[[length(figure_index) + 1L]] <<- data.frame(
    figure = name, png = paste0("figures/", name, ".png"), svg = paste0("figures/", name, ".svg"),
    data_source = source_file, note = note)
}
for (section in sections) {
  for (analysis_set in c("all_completed", "all_converged")) {
    data <- direct[direct$section == section & direct$analysis_set == analysis_set, ]
    cell_order <- design$cell_code[design$section == section]
    data$cell <- factor(data$cell_code, levels = rev(cell_order))
    data$endpoint <- factor(data$metric, levels = primary, labels = unname(metric_labels[primary]))
    data$method <- factor(data$target_method, levels = c("scad", "mcp"), labels = c("SCAD", "MCP"))
    finite <- data[is.finite(data$mean_difference), ]
    missing <- unique(data[!is.finite(data$mean_difference), c("cell", "endpoint")])
    p <- ggplot(finite, aes(x = mean_difference, y = cell, color = method)) +
      geom_vline(xintercept = 0, color = "#677383", linetype = "dashed", linewidth = .4) +
      geom_errorbar(aes(xmin = ci_lower_95, xmax = ci_upper_95), orientation = "y",
                    width = .18, position = position_dodge(width = .45), linewidth = .4) +
      geom_point(aes(shape = method), position = position_dodge(width = .45), size = 1.8) +
      facet_wrap(~endpoint, nrow = 1, scales = "free_x", drop = FALSE) +
      scale_y_discrete(drop = FALSE) + scale_color_manual(values = c("#0072B2", "#D55E00")) +
      labs(title = paste("SCAD and MCP versus Lasso |", section),
           subtitle = paste("Analysis:", analysis_set, "| Each row is a separate design cell"),
           x = "Paired mean difference: method minus Lasso", y = "Design cell", color = NULL, shape = NULL,
           caption = paste("Bars: unadjusted 95% paired Monte Carlo t intervals. RMSE: left favors SCAD/MCP; MCC: right favors SCAD/MCP.",
                           "\nNA denotes a structurally undefined metric. Cell definitions: tables/design-key.csv."))
    if (nrow(missing)) p <- p + geom_text(data = missing, aes(x = 0, y = cell, label = "NA"),
                                         inherit.aes = FALSE, color = "#697386", size = 3)
    save_figure(p, paste0("paired-", section, "-", analysis_set), 13,
                max(4, 2.8 + .24 * length(cell_order)), "tables/scad-mcp-vs-lasso.csv",
                "No replication pooling; both methods must converge for paired sensitivity.")
  }
}

gain <- pairs[pairs$analysis_set == "all_completed" & pairs$reference_method == "gwr" &
                pairs$metric == "new_location_signal_rmse", ]
gain$cell <- factor(gain$cell_code, levels = rev(design$cell_code))
gain$method <- factor(gain$target_method, levels = methods[-1], labels = method_labels[methods[-1]])
p <- ggplot(gain, aes(method, cell, fill = mean_difference)) +
  geom_tile(color = "white", linewidth = .25) +
  geom_point(data = gain[gain$favors == "inconclusive", ], shape = 1, size = 1.5, color = "#222222") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       name = "RMSE difference\nmethod - GWR") +
  labs(title = "New-location prediction versus unpenalized GWR", x = NULL, y = "Design cell",
       subtitle = "All 39 cells; primary all-completed analysis",
       caption = "Blue: lower signal RMSE than GWR. Circle: interval includes zero. No cross-cell pooling.\nOracles use true supports; they are simulation benchmarks. M/R/S denote main/reference/stress cells.") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")
save_figure(p, "prediction-vs-gwr", 10, 12.5,
            "tables/all-prespecified-primary-comparisons.csv", "Fill is the existing paired mean difference; not a percentage gain.")

cv <- convergence
cv$cell <- factor(cv$cell_code, levels = rev(design$cell_code))
cv$method_label <- factor(cv$method, levels = methods, labels = method_labels[methods])
cv$warnings <- cv$n_total - cv$n_all_converged
p <- ggplot(cv, aes(method_label, cell, fill = warnings)) +
  geom_tile(color = "#DDDDDD", linewidth = .2) +
  geom_text(aes(label = ifelse(warnings == 0, "", warnings)), size = 3) +
  scale_fill_gradient(low = "#FAFAFA", high = "#D55E00", name = "Fits with\na warning") +
  labs(title = "Final-fit convergence audit", subtitle = "Nonzero labels count fits with at least one nonconverged local estimate",
       x = NULL, y = "Design cell", caption = "Unlabelled cells have zero warnings, not missing data. All valid fits remain in the primary analysis.") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")
save_figure(p, "convergence-audit", 10, 12.5, "tables/convergence-and-timing-by-cell.csv",
            "Final convergence only; CV convergence is retained separately in the tables.")

time_fields <- c("tuning_elapsed_seconds_mean", "fit_elapsed_seconds_mean", "prediction_elapsed_seconds_mean")
times <- do.call(rbind, lapply(time_fields, function(field) data.frame(
  cell_code = convergence$cell_code, section = convergence$section, method = convergence$method,
  component = sub("_elapsed_seconds_mean", "", field), seconds = convergence[[field]])))
write_csv(times, "method-time-components.csv")
positive <- times[is.finite(times$seconds) & times$seconds > 0, ]
positive$x <- match(positive$method, methods) +
  .5 * (match(positive$cell_code, design$cell_code) - 20) / 39
positive$section <- factor(positive$section, levels = sections)
p <- ggplot(positive, aes(x, seconds, color = section)) + geom_point(size = 1.6, alpha = .75) +
  facet_grid(section ~ component) + scale_y_log10() +
  scale_x_continuous(breaks = seq_along(methods), labels = unname(method_labels), limits = c(.5, 8.5)) +
  labs(title = "Method timing components", subtitle = "Each dot is one design-cell mean; no cross-cell averaging",
       x = NULL, y = "Seconds (log scale)",
       caption = "Zero/nonfinite components are omitted from this log-scale plot only and retained in the CSV.\nCommon bandwidth and neighbor-search costs are separate in common-bandwidth.csv and convergence-and-timing-by-cell.csv.") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none")
save_figure(p, "timing-components", 14, 11, "tables/method-time-components.csv",
            "Log scale is display only. No shared preprocessing cost is attributed repeatedly to methods.")
write_csv(do.call(rbind, figure_index), "figure-index.csv", directory = ".")

notes <- c(
  "# Penalized GWR production results: the validated simulation",
  "", "## Verification and scope", "",
  "The imported archive, source files, scientific signature and completed task shards passed local validation.",
  "There are 39 design cells, 1,900 datasets and 15,200 fits. No model was refitted.",
  "The frozen analysis uses the same data, spatial folds and GWR-selected bandwidth across methods.",
  "Main (32 x 50), reference (2 x 100) and stress (5 x 20) cells are never pooled into a single performance average.",
  "", "## Main comparisons against Lasso", "",
  as.character(knitr::kable(decisions[decisions$section == "main" & decisions$analysis_set == "all_completed",
    c("target_method", "metric", "target_better", "lasso_better", "inconclusive")], format = "pipe", row.names = FALSE)),
  "", "Counts describe separate cells using the existing unadjusted 95% paired Monte Carlo t intervals.",
  "They are not a pooled test, a multiplicity-adjusted significance claim, or a probability of winning on a new dataset.",
  "SCAD/MCP improve support selection relative to Lasso across this main design, but do not uniformly dominate it in prediction or active-coefficient recovery.",
  "See tables/lasso-favored-primary-cells.csv for every primary cell where the interval favors Lasso.",
  "A direct SCAD-versus-MCP inferential ranking was not prespecified and is not added here.",
  "", "## Convergence and sensitivity", "",
  as.character(knitr::kable(method_convergence, format = "pipe", row.names = FALSE, digits = 7)),
  "", paste("Among the 234 direct primary comparison rows,", sum(sensitivity$decision_changed),
             "change interval classification in the all-converged sensitivity analysis."),
  "The primary analysis retains valid completed fits with warnings. Paired sensitivity requires both compared fits to converge.",
  "Final-fit convergence flags are distinct from CV/target convergence indicators, which remain in the source and diagnostic tables.",
  "", "## Null reference and measurement rules", "",
  "MCC and active-coefficient RMSE are undefined for null-slope truth and remain NA, not zero.",
  "GWR, Ridge and the oracles are not evaluated as variable selectors; their selection metrics remain NA.",
  "See tables/null-slope-reference.csv for false positives and other defined null-scenario metrics.",
  "All numeric tables preserve n_total and n_finite. Display tables round to four decimals; CSV values retain numeric precision.",
  "", "## Timing and limits", "",
  "Runtime components are reported separately; shared bandwidth/search cost is not multiplied by the number of methods.",
  "The timing figure uses a log scale only for display; zero/nonfinite values remain in its source table.",
  "These results apply to the frozen Gaussian design and common GWR bandwidth rule; no external-package speed comparison or general optimality claim follows.",
  "Oracle methods use known true supports and are simulation reference benchmarks.",
  "Full coefficient surfaces were not retained. No spatial coefficient maps or new fitting were generated for this report.",
  "The original ARF runtime (R 4.3.0) is preserved in source; this reporting session records local R independently.",
  "", "## Reproducibility and file index", "",
  "- figure-index.csv maps every PNG/SVG figure to its exact numeric source table.",
  "- tables/design-key.csv maps compact M/R/S labels to the unchanged design cell identifiers.",
  "- tables/primary-*.tex and paired-counts-*.tex are knitr::kable longtables; load longtable and booktabs in LaTeX.",
  "- tables/primary-metrics.csv retains both analysis sets; all-prespecified-primary-comparisons.csv retains all frozen primary contrasts.",
  "- validation/ contains numerical reconciliation, input snapshots and session information.",
  "- code/build-report.R is the exact reporting script used; manifest-sha256.csv records all output sizes and SHA256 hashes.")
writeLines(notes, file.path(staging, "README.md"))
tex <- c("% Generated from validated the validated simulation; no additional models or tests.",
         "\\paragraph{Simulation results.}",
         "The study comprised 39 design cells, 1,900 datasets and 15,200 geographically weighted model fits.",
         "All completed valid fits were retained in the primary analysis, with a separate all-converged sensitivity analysis.",
         "The main, reference and stress designs were summarized separately.",
         "The paired differences used common data, spatial folds and bandwidths.",
         "Cell-level conclusions are based on unadjusted 95\\% Monte Carlo intervals and do not imply a pooled or multiplicity-adjusted test.")
for (method in c("scad", "mcp")) {
  d <- decisions[decisions$section == "main" & decisions$analysis_set == "all_completed" &
                   decisions$target_method == method, ]
  get_count <- function(metric, column) d[d$metric == metric, column]
  tex <- c(tex, sprintf(
    "%s was favored over Lasso in %d of 32 main cells for new-location signal RMSE, %d for active-coefficient RMSE, and %d for MCC.",
    toupper(method), get_count(primary[1], "target_better"), get_count(primary[2], "target_better"), get_count(primary[3], "target_better")))
}
tex <- c(tex, "Neither nonconvex estimator uniformly dominated Lasso. Undefined null-support MCC values were retained as missing, not replaced by zero.")
writeLines(tex, file.path(staging, "results-summary.tex"))
stopifnot(file.copy(script_path, file.path(staging, "code/build-report.R")))
write_csv(source_before, "source-before.csv", "validation")
write_csv(summary_before, "summary-before.csv", "validation")
stopifnot(gwrs_sim_verify_completion(source_root))
source_after <- snapshot(source_root)
summary_after <- snapshot(summary_root)
stopifnot(identical(source_before, source_after), identical(summary_before, summary_after))
write_csv(source_after, "source-after.csv", "validation")
write_csv(summary_after, "summary-after.csv", "validation")
write_csv(data.frame(check = c("primary_rows", "direct_paired_rows", "source_preserved", "summary_preserved"),
                      value = c(nrow(primary_metrics), nrow(direct), 1, 1), passed = TRUE),
          "checks.csv", "validation")
writeLines(capture.output(sessionInfo()), file.path(staging, "validation/session-info.txt"))
manifest <- snapshot(staging)
write.csv(manifest, file.path(staging, "manifest-sha256.csv"), row.names = FALSE)
stopifnot(identical(manifest$sha256, unname(gwrs_sim_sha256(file.path(staging, manifest$file)))))
writeLines(c("Schema: gwrs-paper-report-v1", "Study: penalized-selection-v2",
             paste0("Completed-UTC: ", gwrs_sim_utc()),
             paste0("Manifest-SHA256: ", gwrs_sim_sha256(file.path(staging, "manifest-sha256.csv")))),
           file.path(staging, "COMPLETED"))
stopifnot(!dir.exists(output), file.rename(staging, output))
cat("Report:", output, "\n")
cat("Figures:", length(figure_index), "(PNG + SVG). Models fitted: 0.\n")
