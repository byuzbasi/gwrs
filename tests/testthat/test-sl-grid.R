sl_grid_fixture <- function(rho=0, p=6L) {
  e <- make_spatial_example(n=65,p=p,seed=20261004)
  e$x[,2] <- rho*e$x[,1]+sqrt(1-rho^2)*e$x[,2]
  # Include an exactly constant column to exercise zero local curvature.
  e$x[,p] <- 2
  e$y <- 1+e$x[,1]-e$x[,2]+.3*e$x[,3]+sin(e$coords[,1]*6)
  e$nn <- gwr_neighbors(e$coords,query_coords=e$coords[1:5,,drop=FALSE],
    k=35,include_self=FALSE)
  e
}

sl_grid_call <- function(e,alpha,d,lambda,screen=FALSE,limit=10000L) {
  gwrs:::run_cross_gwr_sl_grid(e$x,e$y,e$x[1:5,,drop=FALSE],e$nn,
    gwrs:::kernel_code("bisquare"),gwrs:::resolve_bandwidth(NULL),
    lambda,alpha,d,screen,TRUE,serial_control(tolerance=1e-7,max_iterations=limit))
}

test_that("Gram grid matches independent residual paths including heavy collinearity", {
  alpha <- c(.1,.9,.1,.9,1,1); d <- c(0,0,1,1,0,1)
  lambda <- c(1,.1,.001,0)
  for(rho in c(0,.998,1)) {
    e <- sl_grid_fixture(rho)
    grid <- sl_grid_call(e,alpha,d,lambda)
    for(b in seq_along(alpha)) {
      old <- gwrs:::run_cross_gwr_sl_path(e$x,e$y,e$x[1:5,,drop=FALSE],e$nn,
        gwrs:::kernel_code("bisquare"),gwrs:::resolve_bandwidth(NULL),lambda,
        alpha[b],d[b],FALSE,TRUE,serial_control(tolerance=1e-7,max_iterations=10000L))
      cols <- (b-1L)*length(lambda)+seq_along(lambda)
      expect_equal(grid$predictions[,cols],old$predictions,tolerance=2e-7)
      expect_equal(grid$coefficients[,,cols],old$coefficients,tolerance=2e-6)
      expect_identical(grid$converged[,cols],old$converged)
      expect_equal(grid$kkt[,cols],old$kkt,tolerance=2e-7)
    }
    expect_true(all(is.finite(grid$predictions)))
    expect_true(all(grid$path_seconds>=0))
    expect_lt(max(grid$kkt),2e-6)
  }
})

test_that("screening, block order and independent KKT checks preserve solutions", {
  e <- sl_grid_fixture(.998)
  lambda <- c(2,.8,.3,.1)
  alpha <- c(.1,.9,.5); d <- c(1,0,.5)
  full <- sl_grid_call(e,alpha,d,lambda)
  screened <- sl_grid_call(e,alpha,d,lambda,TRUE)
  reverse <- sl_grid_call(e,rev(alpha),rev(d),lambda)
  expect_equal(screened$predictions,full$predictions,tolerance=2e-6)
  for(b in seq_along(alpha)) {
    cols <- (b-1L)*length(lambda)+seq_along(lambda)
    reverse_cols <- (length(alpha)-b)*length(lambda)+seq_along(lambda)
    expect_identical(full$predictions[,cols],reverse$predictions[,reverse_cols])
    for(t in 1:5) {
      ids <- e$nn$index[,t]; dist <- e$nn$distance[,t]
      w <- pmax(0,1-(dist/(max(dist)*(1+1e-10)))^2)^2; w <- w/sum(w)
      xm <- colSums(e$x[ids,,drop=FALSE]*w); ym <- sum(w*e$y[ids])
      xc <- sweep(e$x[ids,,drop=FALSE],2,xm); yc <- e$y[ids]-ym
      g <- crossprod(xc*sqrt(w)); s <- crossprod(xc,w*yc)
      eig <- eigen(g,symmetric=TRUE)
      cutoff <- max(1,norm(g,"I"))*1e-10
      pilot <- drop(eig$vectors %*% (ifelse(eig$values>cutoff,1/eig$values,0)*
        drop(crossprod(eig$vectors,s))))
      for(l in seq_along(lambda)) {
        beta <- full$coefficients[-1,t,cols[l]]
        score <- drop(crossprod(xc,w*(yc-drop(xc%*%beta))))+
          lambda[l]*(1-alpha[b])*(d[b]*pilot-beta)
        violation <- ifelse(abs(beta)<=1e-12,pmax(0,abs(score)-lambda[l]*alpha[b]),
          abs(score-lambda[l]*alpha[b]*sign(beta)))
        expect_lt(abs(full$kkt[t,cols[l]]-max(violation)),1e-8)
      }
    }
  }
})

test_that("grid rejects invalid inputs and reports exhausted iteration budgets", {
  e <- sl_grid_fixture(.998)
  expect_error(sl_grid_call(e,c(.1,.9),0,c(1,.1)),"Invalid SL grid")
  expect_error(sl_grid_call(e,0,0,c(1,.1)),"Invalid SL grid")
  expect_error(sl_grid_call(e,.5,0,c(.1,1)),"decreasing")
  low <- sl_grid_call(e,.9,1,c(.1,.001),limit=1L)
  expect_true(any(low$converged==0))
  bad <- e; bad$nn$index[1,1] <- 0L
  invalid <- sl_grid_call(bad,.5,0,c(1,.1))
  expect_equal(invalid$status[1],4L)
  expect_true(all(is.na(invalid$predictions[1,])))
})

test_that("thirty-predictor Gram paths preserve the legacy objective and predictions", {
  e <- sl_grid_fixture(.998,p=30L)
  lambda <- c(1,.1,.001)
  grid <- sl_grid_call(e,c(.1,.9),c(0,1),lambda)
  for(b in 1:2) {
    old <- gwrs:::run_cross_gwr_sl_path(e$x,e$y,e$x[1:5,,drop=FALSE],e$nn,
      gwrs:::kernel_code("bisquare"),gwrs:::resolve_bandwidth(NULL),lambda,
      c(.1,.9)[b],c(0,1)[b],FALSE,TRUE,
      serial_control(tolerance=1e-7,max_iterations=10000L))
    cols <- (b-1L)*length(lambda)+seq_along(lambda)
    expect_lt(max(abs(grid$predictions[,cols]-old$predictions)),2e-6)
    expect_lt(max(abs(grid$coefficients[,,cols]-old$coefficients)),2e-5)
    expect_identical(grid$converged[,cols],old$converged)
    # Compare the frozen centered quadratic + l1 + shifted l2 objective,
    # using a shared center obtained from the unpenalized reference path.
    pilot <- gwrs:::run_cross_gwr_sl_path(e$x,e$y,e$x[1:5,,drop=FALSE],e$nn,
      gwrs:::kernel_code("bisquare"),gwrs:::resolve_bandwidth(NULL),0,
      .5,0,FALSE,TRUE,serial_control())$coefficients
    for(t in 1:5) {
      ids <- e$nn$index[,t]; dist <- e$nn$distance[,t]
      w <- pmax(0,1-(dist/(max(dist)*(1+1e-10)))^2)^2; w <- w/sum(w)
      objective <- function(coef,l) {
        res <- e$y[ids]-drop(cbind(1,e$x[ids,,drop=FALSE])%*%coef)
        .5*sum(w*res^2)+lambda[l]*c(.1,.9)[b]*sum(abs(coef[-1]))+
          .5*lambda[l]*(1-c(.1,.9)[b])*sum((coef[-1]-c(0,1)[b]*pilot[-1,t,1])^2)
      }
      for(l in seq_along(lambda)) expect_lt(abs(
        objective(grid$coefficients[,t,cols[l]],l)-
        objective(old$coefficients[,t,l],l)),1e-8)
    }
  }
})
