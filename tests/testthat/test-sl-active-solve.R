active_fixture <- function(p=10L,rho=.998) {
  e<-make_spatial_example(n=70,p=p,seed=20261005)
  for(j in 2:4) e$x[,j]<-sqrt(rho)*e$x[,1]+sqrt(1-rho)*e$x[,j]
  e$x[,p]<-1
  e$y<-1+e$x[,1]-.8*e$x[,2]+.4*e$x[,3]+.1*sin(6*e$coords[,1])
  e$nn<-gwr_neighbors(e$coords,query_coords=e$coords[1:4,,drop=FALSE],k=45,include_self=FALSE)
  e
}

active_grid <- function(e,enable=TRUE,alpha=.9,d=1,limit=10000L,tol=1e-7) {
  gwrs:::run_cross_gwr_sl_grid(e$x,e$y,e$x[1:4,,drop=FALSE],e$nn,
    gwrs:::kernel_code('bisquare'),gwrs:::resolve_bandwidth(NULL),
    c(1,.1,.001,.0001,0),alpha,d,FALSE,TRUE,
    serial_control(tolerance=tol,max_iterations=limit),active_solve=enable)
}

active_objective <- function(e,fit,center,alpha,d) {
  lambda<-c(1,.1,.001,.0001,0)
  value<-matrix(NA_real_,4,length(lambda))
  for(t in 1:4) {
    ids<-e$nn$index[,t];w<-adaptive_weights(e$nn,t)
    coef<-fit$coefficients[,t,]
    res<-e$y[ids]-cbind(1,e$x[ids,,drop=FALSE])%*%coef
    shifted<-coef[-1,,drop=FALSE]-d*center[-1,t]
    value[t,]<-.5*colSums(w*res^2)+lambda*alpha*colSums(abs(coef[-1,,drop=FALSE]))+
      .5*lambda*(1-alpha)*colSums(shifted^2)
  }
  value
}

test_that('accepted active solves preserve the objective and improve precision', {
  for(p in c(10L,30L)) {
    e<-active_fixture(p)
    for(d in c(0,1)) {
      old<-active_grid(e,FALSE,d=d);fast<-active_grid(e,TRUE,d=d)
      ref<-active_grid(e,FALSE,d=d,limit=100000L,tol=1e-11)
      # Some rank-poor p=30 fixtures exhaust the unchanged budget in both
      # engines. Acceleration must not introduce additional nonconvergence.
      expect_true(all(fast$converged>=old$converged))
      expect_true(all(ref$converged==1))
      expect_gt(sum(fast$active_solve_accepted),0)
      expect_lt(sum(fast$iterations),sum(old$iterations))
      pilot<-old$coefficients[,,5]
      expect_lt(max(active_objective(e,fast,pilot,.9,d)-active_objective(e,old,pilot,.9,d)),1e-9)
      expect_lte(max(abs(fast$predictions-ref$predictions)),max(abs(old$predictions-ref$predictions))+1e-8)
      expect_lte(max(abs(fast$coefficients-ref$coefficients)),max(abs(old$coefficients-ref$coefficients))+1e-8)
      accepted<-fast$active_solve_accepted>0
      limits<-matrix(rep(1e-7*(1+.9*c(1,.1,.001,.0001,0)),each=4),4)
      expect_true(all(fast$kkt[accepted]<=limits[accepted]+1e-12))
    }
  }
})

test_that('rank deficiency, pure lasso, lambda zero and exhausted budgets retain fallbacks', {
  e<-active_fixture(rho=1)
  fast<-active_grid(e,TRUE);old<-active_grid(e,FALSE)
  expect_true(all(is.finite(fast$predictions)))
  expect_equal(fast$predictions[,5],old$predictions[,5],tolerance=1e-12)
  expect_true(all(fast$active_solve_attempts[,5]==0))
  for(alpha in c(.9,1)) {
    for(budget in c(1L,16L)) {
      a<-active_grid(e,TRUE,alpha=alpha,limit=budget)
      b<-active_grid(e,FALSE,alpha=alpha,limit=budget)
      expect_identical(a$predictions,b$predictions)
      expect_identical(a$converged,b$converged)
      expect_true(all(a$active_solve_attempts==0))
    }
  }
  a<-active_grid(e,TRUE,alpha=1);b<-active_grid(e,FALSE,alpha=1)
  expect_identical(a$predictions,b$predictions)
  expect_true(all(a$active_solve_attempts==0))
})

test_that('rejected active faces fall back without claiming an accepted solution', {
  e<-active_fixture(30L)
  a<-active_grid(e,TRUE,alpha=.5,d=.5)
  expect_gte(sum(a$active_solve_attempts),sum(a$active_solve_accepted))
  expect_true(all(a$active_solve_accepted<=a$active_solve_attempts))
  expect_true(all(is.finite(a$predictions)))
})
