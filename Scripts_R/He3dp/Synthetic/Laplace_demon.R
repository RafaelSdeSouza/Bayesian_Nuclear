
require(LaplacesDemon)
N <- 10000
J <- 5
X <- matrix(1,N,J)
for (j in 2:J) {X[,j] <- rnorm(N,runif(1,-3,3),runif(1,0.1,1))}
beta <- runif(J,-3,3)
e <- rnorm(N,0,0.1)
y <- tcrossprod(X, t(beta)) + e
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
  beta <- rnorm(Data$J)
  sigma <- runif(1)
  return(c(beta, sigma))
}
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
               parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)


Model <- function(parm, Data)
{
  ### Parameters
  beta <- parm[Data$pos.beta]
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  ### Log-Prior
  beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
  sigma.prior <- dgamma(sigma, 25, log=TRUE)
  ### Log-Likelihood
  mu <- tcrossprod(Data$X, t(beta))
  LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
                   yhat=rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)
}


Initial.Values <- c(rep(0,J), 1)



Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
                     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))





