# Auditable variable-selection metrics for gwrs simulation studies.

gwrs_selection_safe_ratio <- function(numerator, denominator) {
  if (length(denominator) != 1L || !is.finite(denominator) ||
      denominator == 0) {
    return(NA_real_)
  }
  as.double(numerator) / as.double(denominator)
}

gwrs_selection_metrics <- function(estimated,
                                   truth,
                                   selection_tolerance = 1e-8,
                                   truth_tolerance = 1e-12,
                                   truth_support = NULL,
                                   evaluate = TRUE) {
  estimated <- as.matrix(estimated)
  truth <- as.matrix(truth)
  if (!identical(dim(estimated), dim(truth)) ||
      any(!is.finite(estimated)) || any(!is.finite(truth))) {
    stop("`estimated` and `truth` must be finite matrices of equal size.",
         call. = FALSE)
  }
  if (length(selection_tolerance) != 1L ||
      !is.finite(selection_tolerance) || selection_tolerance < 0 ||
      length(truth_tolerance) != 1L || !is.finite(truth_tolerance) ||
      truth_tolerance < 0) {
    stop("Selection and truth tolerances must be nonnegative scalars.",
         call. = FALSE)
  }

  metric_names <- c(
    "true_positive", "false_positive", "false_negative", "true_negative",
    "true_positive_rate", "false_positive_rate", "specificity",
    "precision", "false_discovery_rate", "f1", "mcc", "support_iou",
    "selection_accuracy",
    "mean_selected_predictors"
  )
  if (!isTRUE(evaluate)) {
    values <- as.list(rep.int(NA_real_, length(metric_names)))
    names(values) <- metric_names
    return(c(list(selection_evaluable = FALSE), values))
  }

  selected <- abs(estimated) > selection_tolerance
  if (is.null(truth_support)) {
    active_truth <- abs(truth) > truth_tolerance
  } else {
    active_truth <- as.matrix(truth_support)
    if (!identical(dim(active_truth), dim(truth)) || anyNA(active_truth)) {
      stop("`truth_support` must be a nonmissing logical matrix matching `truth`.",
           call. = FALSE)
    }
    storage.mode(active_truth) <- "logical"
  }
  true_positive <- sum(selected & active_truth)
  false_positive <- sum(selected & !active_truth)
  false_negative <- sum(!selected & active_truth)
  true_negative <- sum(!selected & !active_truth)
  precision <- gwrs_selection_safe_ratio(
    true_positive, true_positive + false_positive
  )
  recall <- gwrs_selection_safe_ratio(
    true_positive, true_positive + false_negative
  )
  f1 <- gwrs_selection_safe_ratio(
    2 * true_positive,
    2 * true_positive + false_positive + false_negative
  )

  predicted_positive <- true_positive + false_positive
  truth_positive <- true_positive + false_negative
  truth_negative <- true_negative + false_positive
  predicted_negative <- true_negative + false_negative
  mcc_denominator <- sqrt(
    as.double(predicted_positive) * as.double(truth_positive) *
      as.double(truth_negative) * as.double(predicted_negative)
  )
  mcc <- if (truth_positive > 0 && truth_negative > 0 &&
             (predicted_positive == 0 || predicted_negative == 0)) {
    # A constant prediction is uninformative when the truth contains both
    # classes. Zero retains the failed replicate in paired summaries.
    0
  } else if (!is.finite(mcc_denominator) || mcc_denominator == 0) {
    # With a one-class truth (for example all-null slopes), correlation is
    # unidentified and FPR remains the directly auditable measure.
    NA_real_
  } else {
    (as.double(true_positive) * as.double(true_negative) -
       as.double(false_positive) * as.double(false_negative)) /
      mcc_denominator
  }

  list(
    selection_evaluable = TRUE,
    true_positive = as.double(true_positive),
    false_positive = as.double(false_positive),
    false_negative = as.double(false_negative),
    true_negative = as.double(true_negative),
    true_positive_rate = recall,
    false_positive_rate = gwrs_selection_safe_ratio(
      false_positive, false_positive + true_negative
    ),
    specificity = gwrs_selection_safe_ratio(
      true_negative, true_negative + false_positive
    ),
    precision = precision,
    false_discovery_rate = gwrs_selection_safe_ratio(
      false_positive, true_positive + false_positive
    ),
    f1 = f1,
    mcc = mcc,
    support_iou = gwrs_selection_safe_ratio(
      true_positive, true_positive + false_positive + false_negative
    ),
    selection_accuracy = mean(selected == active_truth),
    mean_selected_predictors = mean(rowSums(selected))
  )
}
