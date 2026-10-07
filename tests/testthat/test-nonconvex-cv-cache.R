cache_example <- function() {
  e <- make_spatial_example(n = 42, p = 4, seed = 403)
  e$fold <- rep(1:3, each = 14)
  e
}

test_that("shared preparation retains separate fold scaling and neighbors", {
  e <- cache_example()
  cache <- gwr_nonconvex_cv_cache(e$x,e$coords,27,e$fold)
  expect_true(environmentIsLocked(cache))
  expect_error(cache$inputs$x[1,1] <- 100, "locked")
  for (i in 1:3) {
    train <- e$fold != i
    expected <- gwrs:::standardize_design(e$x[train,,drop=FALSE],TRUE)
    expect_identical(cache$folds[[i]]$prepared,expected)
    expect_identical(cache$folds[[i]]$neighbors,
      gwr_neighbors(e$coords[train,,drop=FALSE],e$coords[!train,,drop=FALSE],k=27,include_self=FALSE))
  }
})

test_that("cached CV retains scores, selection, refits and method-specific anchors exactly", {
  e <- cache_example()
  ctl <- gwrs_control(nonconvex_solver="guarded_block",n_threads=1,
    max_iterations=4000,diagnostics="none")
  for (scaling in c(TRUE,FALSE)) {
    cache <- gwr_nonconvex_cv_cache(e$x,e$coords,27,e$fold,standardize=scaling)
    for (method in c("scad","mcp")) for (kernel in c("bisquare","gaussian")) {
      a <- list(x=e$x,y=e$y,coords=e$coords,k=27,fold=e$fold,
        penalty=method,kernel=kernel,standardize=scaling,n_lambda=4,
        lambda_min_ratio=.1,control=ctl)
      plain <- do.call(cv_gwr_nonconvex,a)
      shared <- do.call(cv_gwr_nonconvex,c(a,list(cv_cache=cache)))
      plain$call <- shared$call <- NULL
      expect_identical(plain,shared)
    }
  }
})

test_that("stale preparation is rejected and responses are never reused", {
  e <- cache_example()
  cache <- gwr_nonconvex_cv_cache(e$x,e$coords,27,e$fold)
  a <- list(x=e$x,y=e$y,coords=e$coords,k=27,fold=e$fold,
    penalty="mcp",lambda=c(100,50),cv_cache=cache,refit=FALSE)
  for (field in c("x","coords","k","fold","standardize")) {
    b <- a
    if (field == "x") b$x[1,1] <- b$x[1,1]+.01
    if (field == "coords") b$coords[1,1] <- b$coords[1,1]+.01
    if (field == "k") b$k <- 26
    if (field == "fold") b$fold <- rev(e$fold)
    if (field == "standardize") b$standardize <- FALSE
    expect_error(do.call(cv_gwr_nonconvex,b),"does not match")
  }
  a$y <- e$y+2
  reused <- do.call(cv_gwr_nonconvex,a)
  a$cv_cache <- NULL
  fresh <- do.call(cv_gwr_nonconvex,a)
  reused$call <- fresh$call <- NULL
  expect_identical(reused,fresh)
  expect_error(gwr_nonconvex_cv_cache(e$x,e$coords,29,e$fold),"training size")
  expect_error(gwr_nonconvex_cv_cache(e$x,e$coords,27,rep(1,42)),"two spatial folds")
})

test_that("guarded failure validity and observation counts are unchanged by caching", {
  e <- cache_example()
  cache <- gwr_nonconvex_cv_cache(e$x,e$coords,27,e$fold)
  a <- list(x=e$x,y=e$y,coords=e$coords,k=27,fold=e$fold,
    penalty="scad",lambda=c(100,.001,.0001),refit=FALSE,
    control=gwrs_control(nonconvex_solver="guarded_block",n_threads=1,
      max_iterations=1,diagnostics="none"))
  plain <- do.call(cv_gwr_nonconvex,a)
  shared <- do.call(cv_gwr_nonconvex,c(a,list(cv_cache=cache)))
  plain$call <- shared$call <- NULL
  expect_identical(plain,shared)
  expect_true(any(!shared$cv$valid))
  expect_true(all(shared$cv$observation_count == nrow(e$x)))
})
