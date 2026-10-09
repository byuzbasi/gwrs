test_that("the county data retain geographic keys, units, and MOE flags", {
  env <- new.env()
  data("acs2020_2024_counties", package = "gwrs", envir = env)
  x <- env$acs2020_2024_counties
  expect_identical(dim(x), c(3107L, 91L))
  expect_false(anyNA(x))
  expect_false(anyDuplicated(x$fips) > 0L)
  expect_true(all(nchar(x$fips) == 5L))
  expect_false(any(c("35011", "48301") %in% x$fips))
  expect_equal(x$log_median_household_income_2024,
               log(x$median_household_income_2024), tolerance = 1e-12)
  expect_setequal(x$spatial_fold, 1:5)
  dictionary <- read.csv(system.file("extdata", "acs2020_2024_counties_dictionary.csv", package = "gwrs"))
  expect_identical(dictionary$column, names(x))
  flags <- names(x)[grepl("__moe90_controlled$", names(x))]
  expect_length(flags, 27L)
  for (flag in flags) {
    expect_type(x[[flag]], "logical")
    moe <- sub("_controlled$", "", flag)
    expect_true(all(x[[moe]][x[[flag]]] == 0))
  }
})

test_that("period-named data preserve the legacy table exactly", {
  env <- new.env()
  data("acs2020_2024_counties", package = "gwrs", envir = env)
  data("acs2024_counties", package = "gwrs", envir = env)
  expect_identical(env$acs2020_2024_counties, env$acs2024_counties)
})
