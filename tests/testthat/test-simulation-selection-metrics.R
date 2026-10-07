selection_metric_helpers <- system.file(
  "simulations", "selection-metrics.R", package = "gwrs"
)
if (!nzchar(selection_metric_helpers)) {
  candidates <- c(
    file.path("inst", "simulations", "selection-metrics.R"),
    file.path("..", "..", "inst", "simulations", "selection-metrics.R")
  )
  selection_metric_helpers <- candidates[file.exists(candidates)][[1L]]
}
sys.source(selection_metric_helpers, envir = environment())

test_that("MCC covers perfect, inverse, and mixed support recovery", {
  truth <- matrix(c(1, 0, 1, 0), nrow = 2L, byrow = TRUE)

  perfect <- gwrs_selection_metrics(truth, truth)
  expect_equal(perfect$mcc, 1)
  expect_equal(perfect$f1, 1)
  expect_equal(perfect$support_iou, 1)
  expect_equal(perfect$specificity, 1)
  expect_equal(perfect$false_discovery_rate, 0)

  inverse <- gwrs_selection_metrics(1 - truth, truth)
  expect_equal(inverse$mcc, -1)

  mixed <- gwrs_selection_metrics(
    matrix(c(1, 1, 0, 0), nrow = 2L, byrow = TRUE), truth
  )
  expect_equal(mixed$true_positive, 1)
  expect_equal(mixed$false_positive, 1)
  expect_equal(mixed$false_negative, 1)
  expect_equal(mixed$true_negative, 1)
  expect_equal(mixed$mcc, 0)
})

test_that("degenerate MCC remains undefined while FPR is auditable", {
  truth <- matrix(0, nrow = 3L, ncol = 2L)
  estimated <- truth
  estimated[1L, 1L] <- 1
  metrics <- gwrs_selection_metrics(estimated, truth)

  expect_true(is.na(metrics$mcc))
  expect_equal(metrics$false_positive_rate, 1 / 6)
  expect_equal(metrics$specificity, 5 / 6)
  expect_equal(metrics$false_discovery_rate, 1)
  expect_equal(metrics$true_positive, 0)
  expect_equal(metrics$false_positive, 1)
})

test_that("constant selection on a two-class truth has MCC zero", {
  truth <- matrix(c(1, 0, 0, 0), nrow = 2L, byrow = TRUE)
  missed <- gwrs_selection_metrics(matrix(0, 2L, 2L), truth)
  all_selected <- gwrs_selection_metrics(matrix(1, 2L, 2L), truth)

  expect_equal(missed$mcc, 0)
  expect_equal(all_selected$mcc, 0)
})

test_that("F1 uses confusion counts at support boundary cases", {
  truth <- matrix(c(1, 0, 0, 0), nrow = 2L, byrow = TRUE)
  missed <- gwrs_selection_metrics(matrix(0, 2L, 2L), truth)
  empty <- gwrs_selection_metrics(
    matrix(0, 2L, 2L), matrix(0, 2L, 2L)
  )

  expect_equal(missed$f1, 0)
  expect_true(is.na(empty$f1))
})

test_that("selection metrics are NA for methods without exact selection", {
  metrics <- gwrs_selection_metrics(
    matrix(rnorm(12), nrow = 4L),
    matrix(0, nrow = 4L, ncol = 3L),
    evaluate = FALSE
  )

  expect_false(metrics$selection_evaluable)
  expect_true(all(is.na(unlist(metrics[-1L]))))
})

test_that("selection metric inputs are validated", {
  expect_error(
    gwrs_selection_metrics(matrix(0, 2L, 2L), matrix(0, 2L, 3L)),
    "equal size"
  )
  expect_error(
    gwrs_selection_metrics(matrix(0, 2L, 2L), matrix(0, 2L, 2L),
                           selection_tolerance = -1),
    "nonnegative"
  )
  expect_error(
    gwrs_selection_metrics(
      matrix(0, 2L, 2L), matrix(0, 2L, 2L),
      truth_support = matrix(TRUE, 2L, 3L)
    ),
    "truth_support"
  )
})

test_that("explicit support is not inferred from tapered coefficients", {
  truth <- matrix(c(1, 1e-16, 0, 0), nrow = 2L, byrow = TRUE)
  support <- matrix(c(TRUE, TRUE, FALSE, FALSE), nrow = 2L, byrow = TRUE)
  estimated <- matrix(c(1, 0, 0, 0), nrow = 2L, byrow = TRUE)

  inferred <- gwrs_selection_metrics(estimated, truth)
  explicit <- gwrs_selection_metrics(
    estimated, truth, truth_support = support
  )

  expect_equal(inferred$false_negative, 0)
  expect_equal(explicit$false_negative, 1)
})
