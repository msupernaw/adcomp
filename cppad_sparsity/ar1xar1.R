require(TMB)
compile("ar1xar1.cpp")
set.seed(123)
n <- 200 ## Size of problem = n*n

## ======================= Simulate separable 2D GMRF 
## - With exponential correlation in both directions
## - phi1 = 1-lag correlation in 1st direction
## - phi2 = 1-lag correlation in 2nd direction
ar1corr <- function(n,phi){
  phi^abs(outer(1:n,1:n,"-"))
}
simgmrf <- function(n1,n2,phi1,phi2){
  u <- matrix(rnorm(n1*n2),n1,n2)
  L1 <- t(chol(ar1corr(n1,phi1)))
  L2 <- t(chol(ar1corr(n2,phi2)))
  x <- L1%*%u         ## phi1 in 1st direction (fastest)
  x <- t(L2%*%t(x))   ## phi2 in 2nd direction
  x
}

## ======================= Simulate data
phi1=exp(-1/(.1*n)) ## Correlation range=10% of grid size first dimension
phi2=exp(-1/(.2*n)) ## Correlation range=20% of grid size second dimension
eta <- simgmrf(n,n,phi1,phi2)
N <- rpois(length(eta),exp(eta))
d <- expand.grid(x=factor(1:n),y=factor(1:n))
d$N <- N

## ======================= Parameterization of phi
f <- function(x) 2/(1 + exp(-2 * x)) - 1
invf <- function(y) -0.5 * log(2/(y + 1) - 1)

## ======================= Fit model
dyn.load(dynlib("ar1xar1"))
obj <- MakeADFun(data=list(N=N),
                 parameters=list(
                   eta=matrix(0,n,n),
                   transf_phi1=invf(0.5),
                   transf_phi2=invf(0.5)),
                 map=list(transf_phi1=factor(NA), transf_phi2=factor(NA)),
                 random=c("eta"))

## Time the current TMB pattern detection (and some additional stuff)
tim1 <- system.time(obj$env$spHess)

## Time the pattern detection using CppAD
DLL <- obj$env$DLL
ADFun <- obj$env$ADFun
theta <- obj$env$par
control <- list()
tim2 <- system.time(res <- .Call("EvalADFunTest", ADFun$ptr, theta, control, PACKAGE = DLL))

## We can check that the pattern is correct with this:
list2spmat <- function(res){
    n <- length(res)
    j <- rep(1:n,sapply(res,length))
    i <- unlist(res) + 1
    spMatrix(n,n,i,j,0*i+1)
}

c(tmb = tim1["elapsed"], cppad = tim2["elapsed"])
#}

## n <- c(20,50,80,90,100,110,120,130)
## y <- sapply(n,f)
## matplot(n*n,t(y),type="b",ylab="Time (sec)", xlab="Dimension")
## lm(I(log(y[1,]))~I(log(n*n)))
## lm(I(log(y[2,]))~I(log(n*n)))

