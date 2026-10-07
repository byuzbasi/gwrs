test_that("dashboard composes existing results with shared visual settings", {
  skip_if_no_heatmap_packages()
  skip_if_not_installed("patchwork")
  fit <- make_heatmap_fit()$fit
  before <- serialize(fit, NULL)
  result <- plot_gwrs_dashboard(fit, "income", resolution = 25, return_data = TRUE)
  expect_s3_class(result$plot, "patchwork")
  expect_named(result$panels, c("observed", "bivariate", "beta", "contour", "local_r2", "residuals"))
  expect_identical(serialize(fit, NULL), before)
  expect_equal(result$panels$bivariate$theme$legend.position, "bottom")
  expect_equal(result$settings$observed_rows, nrow(coef(fit)))
  expect_equal(result$surfaces$beta$points$value_1, coef(fit)[, "income"])
  for (surface in result$surfaces) {
    expect_equal(surface$bandwidth, result$settings$bandwidth)
    expect_true(surface$settings$visual_interpolation_only)
  }
  expect_equal(result$panels$observed$scales$get_scales("colour")$limits,
               result$panels$beta$scales$get_scales("fill")$limits)
  direct <- plot_heatmap(fit, "income", resolution = 25, return_data = TRUE)
  expect_equal(result$surfaces$beta$grid, direct$grid)
  expect_equal(result$surfaces$beta$points, direct$points)
  expect_equal(result$surfaces$local_r2$points$value_1, fit$local_r2)
  expect_equal(result$surfaces$residuals$points$value_1, fit$residuals)
  expect_silent(patchwork::patchworkGrob(result$plot))
})

test_that("atlas supports per-coefficient or common scales and explicit labels", {
  skip_if_no_heatmap_packages()
  skip_if_not_installed("patchwork")
  fit <- make_heatmap_fit()$fit
  free <- plot_gwrs_atlas(fit, c("income", "education"), resolution = 25,
    labels = c(income = "Income"), return_data = TRUE)
  fixed <- plot_gwrs_atlas(fit, c("income", "education"), resolution = 25,
    scales = "fixed", return_data = TRUE)
  expect_equal(free$panels$income$labels$title, "Income")
  expect_equal(free$panels$education$labels$title, "education")
  expect_equal(fixed$panels$income$scales$get_scales("fill")$limits,
               fixed$panels$education$scales$get_scales("fill")$limits)
  expect_equal(free$surfaces$income$points$value_1, coef(fit)[, "income"])
  expect_equal(free$surfaces$education$bandwidth, free$surfaces$income$bandwidth)
  expect_silent(patchwork::patchworkGrob(free$plot))
  expect_error(plot_gwrs_atlas(fit, "missing"), "coefficient names")
  expect_error(plot_gwrs_atlas(fit, c("income", "income")), "distinct")
  expect_error(plot_gwrs_atlas(fit, labels = "unnamed"), "named character")
})

test_that("unselected support, optional inference and constant maps render honestly", {
  skip_if_no_heatmap_packages()
  skip_if_not_installed("patchwork")
  fit <- make_heatmap_fit(diagnostics = "full", lambda = 0.35)$fit
  result <- plot_gwrs_dashboard(fit, "income", resolution = 25,
    significance = TRUE, return_data = TRUE)
  expect_true(any(!result$surfaces$beta$grid$selected))
  expect_true(result$surfaces$beta$settings$significance)
  expect_false(result$surfaces$local_r2$settings$significance)
  expect_silent(patchwork::patchworkGrob(result$plot))
  zero <- make_heatmap_fit(lambda = 100, alpha = 1)$fit
  zero_result <- plot_gwrs_dashboard(zero, "income", resolution = 25,
    return_data = TRUE)
  expect_true(all(coef(zero)[, "income"] == 0))
  expect_true(all(!zero_result$surfaces$beta$grid$selected))
  expect_silent(patchwork::patchworkGrob(zero_result$plot))
})

test_that("overlay CRS and invalid options fail clearly", {
  skip_if_no_heatmap_packages()
  skip_if_not_installed("patchwork")
  fit <- make_heatmap_fit()$fit
  boundary <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
  ))), crs = 3857)
  overlay <- sf::st_transform(boundary, 4326)
  result <- plot_gwrs_dashboard(fit, "income", boundary = boundary,
    overlays = overlay, crs = 3857, resolution = 25, return_data = TRUE)
  expect_silent(patchwork::patchworkGrob(result$plot))
  sf_layers <- Filter(function(layer) inherits(layer$geom, "GeomSf"), result$panels$beta$layers)
  expect_equal(sf::st_crs(sf_layers[[length(sf_layers)]]$data)$epsg, 3857)
  expect_error(plot_gwrs_dashboard(fit, "residuals"), "coefficient")
  expect_error(plot_gwrs_dashboard(fit, "income", resoultion = 25), "named plot_heatmap")
  expect_error(plot_gwrs_dashboard(fit, "income", ncol = 0), "positive integer")
  expect_error(plot_gwrs_dashboard(fit, "income", return_data = NA), "TRUE or FALSE")
  expect_error(plot_gwrs_dashboard(fit, "income", coefficient_limits = c(2, 1),
                                  resolution = 25), "must increase")
  expect_error(plot_gwrs_dashboard(fit, "income", overlays = sf::st_set_crs(boundary, NA),
                                  resolution = 25), "known.*CRS")
})
test_that("documentation additions preserve existing S3 registrations", {
  namespace <- asNamespace("gwrs")
  for (method in c("print.gwrs_fit", "print.gwrs_path", "print.gwrs_cv",
                   "summary.gwrs_fit", "print.summary.gwrs_fit")) {
    generic <- sub("\\..*$", "", method)
    class <- sub("^[^.]+\\.", "", method)
    expect_identical(getS3method(generic, class), get(method, namespace))
  }
})
