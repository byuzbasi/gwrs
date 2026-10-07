test_that("heat maps return documented data and ggplot objects", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  result <- plot_heatmap(
    fit, "income", resolution = 35L, show_points = TRUE,
    return_data = TRUE
  )

  expect_named(
    result,
    c(
      "plot", "grid", "points", "boundary", "bandwidth", "settings",
      "interpolation_points"
    )
  )
  expect_s3_class(result$plot, "ggplot")
  expect_true(all(c(
    "x", "y", "value", "selected", "interpolation_valid"
  ) %in% names(result$grid)))
  expect_true(inherits(result$boundary, "sfc"))
  expect_true(isTRUE(result$settings$visual_interpolation_only))
  expect_silent(ggplot2::ggplot_build(result$plot))

  alias <- plot_heatmap(
    fit, "beta_1", resolution = 30L, return_data = TRUE
  )
  expect_equal(alias$points$value, result$points$value)
  expect_error(
    plot_heatmap(fit, "not_a_variable", resolution = 30L),
    "Valid choices"
  )
})

test_that("value type selects diverging or sequential scales", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  residual <- plot_heatmap(fit, "residuals", resolution = 30L)
  local_r2 <- plot_heatmap(fit, "local_r2", resolution = 30L)

  residual_scale <- residual$scales$get_scales("fill")
  r2_scale <- local_r2$scales$get_scales("fill")
  expect_equal(r2_scale$limits, c(0, 1))
  residual_colors <- residual_scale$palette(c(0, 0.5, 1))
  expect_length(residual_colors, 3L)
  expect_equal(toupper(residual_colors[2L]), "#F7F7F7")
})

test_that("contour and combined map types build", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  contour <- plot_heatmap(
    fit, "fitted", type = "contour", resolution = 35L,
    contour_bins = 7L
  )
  combined <- plot_heatmap(
    fit, "response", type = "heatmap_contour", resolution = 35L
  )
  expect_true(any(vapply(
    contour$layers,
    function(layer) inherits(layer$geom, "GeomContour"),
    logical(1L)
  )))
  expect_silent(ggplot2::ggplot_build(contour))
  expect_silent(ggplot2::ggplot_build(combined))
})

test_that("sparse selection and significance are separate overlays", {
  skip_if_no_heatmap_packages()
  sparse <- make_heatmap_fit(lambda = 0.35)$fit
  selected <- plot_heatmap(
    sparse, "income", resolution = 40L, return_data = TRUE
  )
  expect_true(any(
    !selected$grid$selected & selected$grid$interpolation_valid,
    na.rm = TRUE
  ))
  expect_gte(length(selected$plot$layers), 3L)

  no_inference <- make_heatmap_fit(diagnostics = "standard")$fit
  expect_error(
    plot_heatmap(
      no_inference, "income", significance = TRUE, resolution = 30L
    ),
    "diagnostics = \"full\""
  )

  full <- make_heatmap_fit(diagnostics = "full")$fit
  significant <- plot_heatmap(
    full, "income", significance = TRUE, resolution = 35L,
    return_data = TRUE
  )
  expect_true(all(c(
    "significant", "significance_support"
  ) %in% names(significant$grid)))
  expect_true(isTRUE(significant$settings$significance))
  expect_silent(ggplot2::ggplot_build(significant$plot))
})

test_that("bivariate maps create a common quantile surface", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit(diagnostics = "full")$fit
  result <- plot_bivariate(
    fit, "income", "local_r2", bins = 3L, resolution = 35L,
    return_data = TRUE
  )
  significance <- plot_bivariate(
    fit, "income", "significance", bins = 3L, resolution = 35L,
    return_data = TRUE
  )

  expect_s3_class(result$plot, "ggplot")
  expect_true(all(c(
    "value_1", "value_2", "bin_1", "bin_2", "bivariate_class"
  ) %in% names(result$grid)))
  expect_lte(length(unique(stats::na.omit(result$grid$bivariate_class))), 9L)
  fill_scale <- result$plot$scales$get_scales("fill")
  expect_true(is.na(fill_scale$na.value))
  expect_false(fill_scale$na.translate)
  expect_true(result$plot$layers[[1L]]$show.legend)
  expect_true(significance$settings$significance_second)
  expect_silent(ggplot2::ggplot_build(result$plot))
  expect_silent(ggplot2::ggplot_build(significance$plot))
})

test_that("autoplot routes heat maps without altering base plot", {
  skip_if_no_heatmap_packages()
  fit <- make_heatmap_fit()$fit
  automatic <- ggplot2::autoplot(
    fit, type = "heatmap", variable = "income", resolution = 30L
  )
  bivariate <- ggplot2::autoplot(
    fit, type = "bivariate", variable = "income",
    variable2 = "local_r2", resolution = 30L
  )
  expect_s3_class(automatic, "ggplot")
  expect_s3_class(bivariate, "ggplot")
  expect_error(
    ggplot2::autoplot(fit, type = "choropleth", variable = "income"),
    "stores point coordinates"
  )
  grDevices::pdf(tempfile(fileext = ".pdf"))
  expect_invisible(plot(fit, type = "coefficient", coefficient = "income"))
  grDevices::dev.off()
})

test_that("KTDD process fields use a zero-centred diverging heat map", {
  skip_if_no_heatmap_packages()
  example <- make_heatmap_example(n = 36L, seed = 913L)
  neighbors <- gwr_neighbors(example$coords, k = 20L)
  laplacian <- gwr_graph_laplacian(neighbors, kernel = "gaussian")
  process_truth <- sin(2 * pi * example$coords[, 1L]) *
    cos(pi * example$coords[, 2L])
  process_truth <- process_truth - mean(process_truth)

  fit <- gwr_ktdd_fit(
    example$x,
    example$y,
    neighbors = neighbors,
    kernel = "gaussian",
    diffusion_operator = laplacian,
    source = as.double(laplacian %*% process_truth),
    coefficient_laplacian = laplacian,
    lambda_diffusion = 1,
    gamma = 0.05,
    standardize = FALSE,
    linear_tolerance = 1e-8,
    linear_max_iterations = 1000L,
    control = serial_control(tolerance = 1e-7, max_iterations = 1000L)
  )
  result <- plot_heatmap(
    fit, "process", resolution = 30L, return_data = TRUE
  )

  expect_equal(result$points$value, fit$process)
  expect_true(all(result$points$selected == 1))
  expect_identical(result$settings$value_type, "process")
  expect_s3_class(result$plot, "ggplot")
  process_scale <- result$plot$scales$get_scales("fill")
  expect_equal(
    toupper(process_scale$palette(c(0, 0.5, 1))[2L]),
    "#F7F7F7"
  )
  expect_silent(ggplot2::ggplot_build(result$plot))
})
