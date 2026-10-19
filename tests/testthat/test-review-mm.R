review_turn <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)

test_that("cell-weighted SO(2) fit agrees with a scalar global reference", {
  set.seed(19)
  A <- matrix(rnorm(16), 8, 2)
  B <- A %*% review_turn(.6) + matrix(rnorm(16, sd = .15), 8, 2)
  WA <- matrix(seq(.7, 2, length.out = 16), 8, 2)
  WB <- matrix(seq(4, .4, length.out = 16), 8, 2)
  # Eliminate the free consensus and translation analytically, then search the circle.
  W <- WA * WB / (WA + WB)
  reference <- function(theta) {
    E <- B %*% review_turn(theta) - A
    shift <- colSums(W * E) / colSums(W)
    sum(W * sweep(E, 2, shift, "-")^2)
  }
  grid <- seq(-pi, pi, length.out = 101)
  optimum <- min(vapply(seq_len(100), function(j) {
    optimize(reference, grid[c(j, j + 1)], tol = 1e-12)$objective
  }, numeric(1)))
  f <- gpa(list(A = A, B = B), transform = "SO", anchor = "A",
           metric = proc_metric(cell = list(A = WA, B = WB)),
           control = gpa_control(tolerance = 1e-10, max_iterations = 500))
  expect_equal(f$objective, optimum, tolerance = 1e-7)
  Y <- fitted(f)
  expected_M <- (WA * Y$A + WB * Y$B) / (WA + WB)
  expect_equal(consensus(f), expected_M, tolerance = 2e-6)
  expect_equal(f$optimality_status, "first_order_stationary")
  expect_lte(f$rank$stationarity, 1e-5)
})

test_that("Huber cell fitting solves the vector loss, rather than the L2 mean", {
  # SO(1) is {1}; translation and scale are prohibited. Only the consensus is free.
  v <- list(A = matrix(c(0, 1, -2), 3, 1), B = matrix(c(0, 1, -2), 3, 1),
            C = matrix(c(10, 11, 8), 3, 1))
  cells <- lapply(v, function(Y) matrix(1, 3, 1))
  f <- gpa(v, transform = proc_similarity("SO", translation = FALSE, scaling = "none"),
           metric = proc_metric(cell = cells), loss = proc_huber(k = 1),
           control = gpa_control(tolerance = 1e-12, max_iterations = 500))
  expect_equal(as.numeric(consensus(f)), c(.5, 1.5, -1.5), tolerance = 1e-6)
  expect_equal(f$objective, 27.75, tolerance = 1e-8)
  expect_equal(f$optimality_status, "first_order_stationary")
  expect_lte(f$rank$stationarity, 1e-6)
})

test_that("MM never certifies stagnation above the true first-order threshold", {
  set.seed(9)
  X <- matrix(rnorm(24), 12, 2)
  v <- list(A = X, B = X %*% review_turn(.5) + matrix(rnorm(24, sd=.2),12,2),
            C = X %*% review_turn(-.7) + matrix(rnorm(24, sd=.2),12,2))
  W <- lapply(v, function(Y) matrix(1, 12, 2)); W$C[1,1] <- 100
  tol <- 1e-8
  f <- gpa(v, transform = "SO", metric = proc_metric(cell = W),
           control = gpa_control(tolerance = tol, max_iterations = 500))
  Y <- fitted(f); M <- consensus(f)
  manual <- function(Y, M) sum(vapply(names(Y), function(nm) sum(W[[nm]]*(Y[[nm]]-M)^2), numeric(1)))
  expect_equal(f$objective, manual(Y,M), tolerance = 1e-10)
  expect_true(all(diff(f$history$objective) <= 1e-10))
  if (identical(f$optimality_status, "first_order_stationary")) {
    expect_lte(f$rank$stationarity, sqrt(tol))
    h <- 1e-6
    for (nm in names(Y)) {
      objective <- function(a) { Z<-Y; Z[[nm]]<-Y[[nm]]%*%review_turn(a); manual(Z,M) }
      expect_lt(abs((objective(h)-objective(-h))/(2*h)), 3e-4)
    }
    for (i in seq_along(M)) {
      Mp<-Mm<-M;Mp[i]<-Mp[i]+h;Mm[i]<-Mm[i]-h
      expect_lt(abs((manual(Y,Mp)-manual(Y,Mm))/(2*h)), 3e-4)
    }
  }
  short <- gpa(v, transform = "SO", metric = proc_metric(cell = W),
                control = gpa_control(tolerance = 1e-16, max_iterations = 1))
  expect_equal(short$optimality_status, "not_converged")
})

test_that("observed-loss chart gradients agree with independent finite differences", {
  set.seed(52)
  n <- 5; d <- 3
  v <- list(A = matrix(rnorm(n*d),n,d), B = matrix(rnorm(n*d),n,d), C = matrix(rnorm(n*d),n,d))
  dat <- proc_data(v)
  Q <- crossprod(matrix(rnorm(d*d),d,d)) + diag(d)
  W <- lapply(v, function(X) matrix(runif(n*d, .3, 2),n,d)); W$B[2,1]<-0
  alpha <- c(A=.5,B=2,C=1)
  for (mode in c("none", "gower", "fixed_consensus", "preshape")) {
    for (family in c("squared_l2", "huber", "tukey")) {
      loss <- switch(family, squared_l2=proc_squared_l2(), huber=proc_huber(k=2), tukey=proc_tukey(c=5))
      scale_on <- mode %in% c("gower", "fixed_consensus")
      spec <- proc_similarity("SO", translation=TRUE, scaling=if(scale_on) "isotropic" else "none")
      model <- list(data=dat,Q=Q,weights=W,alpha=alpha,loss=loss,spec=spec,
        mode=mode,anchor="A",translate=mode!="preshape",scale=scale_on,
        landmark=list(A=seq(.7,1.3,length.out=n),B=rep(1,n),C=rep(2,n)))
      model$D <- alpha * vapply(v, function(X) sum(scale(X,scale=FALSE)^2), numeric(1))
      model$energy <- sum(model$D)
      state <- list(transforms=lapply(v,function(X)proc_identity(spec,d)), M=matrix(rnorm(n*d),n,d))
      chart <- gprocrustes:::.mm_chart(state,model)
      par <- chart$par + rnorm(length(chart$par),sd=.1)
      dec <- gprocrustes:::.mm_decode(par,chart,model)
      val <- gprocrustes:::.mm_evaluate(dec$state,model)
      analytic <- gprocrustes:::.mm_chart_gradient(val,dec,chart,model)
      # Literal row-vector definition; no production loss, score or moments helpers.
      objective <- function(p) {
        st <- gprocrustes:::.mm_decode(p,chart,model)$state
        total<-0
        for(nm in names(v)) {
          tr<-st$transforms[[nm]]
          E<-sweep(tr$s*v[[nm]]%*%tr$R,2,tr$t,"+")-st$M
          for(j in seq_len(n)) {
            e<-sqrt(W[[nm]][j,])*E[j,]
            r<-sqrt(drop(e%*%Q%*%e))
            rho<-switch(family,squared_l2=r*r,huber=if(r<=2) .5*r*r else 2*r-2,
                        tukey=if(r<5)25/6*(1-(1-(r/5)^2)^3) else 25/6)
            total<-total+alpha[[nm]]*model$landmark[[nm]][j]*rho
          }
        }
        total
      }
      numeric_grad<-vapply(seq_along(par),function(j) {
        h<-1e-6*max(1,abs(par[j]));plus<-minus<-par
        plus[j]<-plus[j]+h;minus[j]<-minus[j]-h
        (objective(plus)-objective(minus))/(2*h)
      },numeric(1))
      expect_equal(val$value,objective(par),tolerance=1e-12,info=paste(mode,family))
      expect_lt(max(abs(numeric_grad-analytic)),1e-5,label=paste(mode,family,"gradient error"))
      for(tr in dec$state$transforms) {
        expect_equal(crossprod(tr$R),diag(d),tolerance=1e-12)
        expect_gt(det(tr$R),0)
      }
      if(mode=="gower") expect_equal(sum(model$D*vapply(dec$state$transforms,function(tr)tr$s^2,numeric(1))),model$energy,tolerance=1e-12)
      if(mode=="fixed_consensus") {
        expect_equal(sum(dec$state$M^2),1,tolerance=1e-12)
        expect_equal(colSums(dec$state$M),rep(0,d),tolerance=1e-12)
      }
    }
  }
})

test_that("public observed-data fits satisfy their scale constraints and true tangent derivatives", {
  set.seed(83)
  n <- 6
  X <- scale(matrix(rnorm(12),n,2),scale=FALSE)
  X <- X / sqrt(sum(X^2))
  v <- list(A=X, B=X%*%review_turn(.4)+matrix(rnorm(12,sd=.03),n,2),
            C=X%*%review_turn(-.6)+matrix(rnorm(12,sd=.03),n,2))
  v$C[1,]<-v$C[1,]+c(.4,-.2)
  W<-lapply(v,function(X)matrix(1,n,2));W$B[2,1]<-0
  alpha<-c(A=1,B=.7,C=1.3)
  landmark<-list(A=rep(1,n),B=seq(.5,1.2,length.out=n),C=seq(.7,1.3,length.out=n))
  Q<-matrix(c(1,.3,.3,2),2)
  metric<-proc_metric(configuration=alpha,landmark=landmark,cell=W,
                      precision=proc_covariance(coordinate=Q,as_precision=TRUE))
  energy<-alpha*vapply(v,function(X)sum(scale(X,scale=FALSE)^2),numeric(1))
  for(mode in c("gower","fixed_consensus","preshape")) {
    f<-gpa(v,transform=proc_similarity("O"),gauge=proc_gauge(scale=mode),metric=metric,
           loss=proc_huber(k=.1),anchor="A",control=gpa_control(tolerance=1e-10,max_iterations=500))
    Y<-fitted(f);M<-consensus(f);trs<-transformations(f)
    objective<-function(Y,M) {
      total<-0
      for(nm in names(Y)) for(j in seq_len(n)) {
        e<-sqrt(W[[nm]][j,])*(Y[[nm]][j,]-M[j,])
        r<-sqrt(drop(e%*%Q%*%e))
        total<-total+alpha[[nm]]*landmark[[nm]][j]*(if(r<=.1).5*r*r else .1*r-.005)
      }
      total
    }
    expect_equal(f$objective,objective(Y,M),tolerance=1e-12)
    expect_equal(f$optimality_status,"first_order_stationary",info=mode)
    expect_lte(f$rank$stationarity,1e-5)
    h<-1e-5
    derivative<-function(fun)(fun(h)-fun(-h))/(2*h)
    for(nm in c("B","C")) {
      expect_lt(abs(derivative(function(a){Z<-Y;Z[[nm]]<-Y[[nm]]%*%review_turn(a);objective(Z,M)})),1e-4)
      if(mode!="preshape") {
        for(k in 1:2) expect_lt(abs(derivative(function(a){Z<-Y;Z[[nm]][,k]<-Z[[nm]][,k]+a;objective(Z,M)})),1e-4)
      }
      if(mode=="fixed_consensus") {
        expect_lt(abs(derivative(function(a){Z<-Y;Z[[nm]]<-sweep(exp(a)*sweep(Y[[nm]],2,trs[[nm]]$t,"-"),2,trs[[nm]]$t,"+");objective(Z,M)})),1e-4)
      }
    }
    if(mode=="gower") {
      s<-vapply(trs,`[[`,numeric(1),"s")
      expect_equal(sum(energy*s^2),sum(energy),tolerance=1e-10)
      q<-sqrt(energy[c("B","C")])*s[c("B","C")]
      expect_lt(abs(derivative(function(a){
        qnew<-as.numeric(q%*%review_turn(a));Z<-Y
        for(i in 1:2){nm<-c("B","C")[[i]];Z[[nm]]<-sweep(sweep(Y[[nm]],2,trs[[nm]]$t,"-")*(qnew[i]/q[i]),2,trs[[nm]]$t,"+")}
        objective(Z,M)
      })),1e-4)
    }
    if(mode=="fixed_consensus") {
      expect_equal(sum(M^2),1,tolerance=1e-10)
      expect_equal(colSums(M),c(0,0),tolerance=1e-10)
    }
    if(mode=="preshape") expect_equal(vapply(Y,function(X)sum(X^2),numeric(1)),c(A=1,B=1,C=1),tolerance=1e-10)
    for(j in seq_along(M)) {
      expect_lt(abs(derivative(function(a){
        N<-M;N[j]<-N[j]+a
        if(mode=="fixed_consensus"){N<-scale(N,scale=FALSE);N<-N/sqrt(sum(N^2))}
        objective(Y,N)
      })),1e-4)
    }
  }
})

test_that("observed fitting retains O(d) reflection components and handles no free scale", {
  set.seed(4);X<-matrix(rnorm(16),8,2)
  W<-list(A=matrix(1,8,2),B=matrix(1,8,2))
  f<-gpa(list(A=X,B=X%*%diag(c(1,-1))),transform="O",anchor="A",metric=proc_metric(cell=W))
  expect_lt(det(transformations(f)$B$R),0)
  expect_lt(f$objective,1e-20)
  expect_no_warning(z<-gpa(list(A=X,B=X),transform=proc_similarity("SO"),anchor="A",
                          metric=proc_metric(configuration=c(A=1,B=0),cell=W)))
  expect_lt(z$objective,1e-20)
  expect_equal(z$optimality_status,"first_order_stationary")
})
