library(LaplacesDemon)
require(gsl)
source("Sfactor3.R")
N <- 150

#obsx1 <- runif(N,0,0.7)
obsx1 <- exp(runif(N,log(1e-3),log(1)))
sd <- 1
# Artificial values, here we will just fit 3 parameters Er, gi, gf
y <- rnorm(N, Sfactor3(obsx1,0.0912,0.0912,2.93,0.0794,6,5,0),sd = sd)

mon.names <- c("LP","Er","gi","gf")
parm.names <- c("Er","gi","gf","sigma")
pos.Er <- 1
pos.gi <- 2
pos.gf <- 3
pos.sigma <- 4

PGF <- function(Data) {
  Er <- runif(1)
  gi <- runif(1)
  gi <- runif(1)
  sigma <- runif(1)
  return(c(Er,gi,gi, sigma))
}
MyData <- list(J=4, PGF=PGF, X=obsx1, mon.names=mon.names,
               parm.names=parm.names, pos.Er=pos.Er,pos.gi = pos.gi,pos.gf=pos.gf, pos.sigma=pos.sigma, y=y)

Model <- function(parm, Data)
{
  ### Parameters
  Er <- parm[Data$pos.Er]
  gi <- parm[Data$pos.gi]
  gf <- parm[Data$pos.gf]
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma 
  ### Log(Prior Densities)
  Er.prior <- dunif(Er,0, 10, log=TRUE)
  gi.prior <- dunif(gi,0, 10, log=TRUE)
  gf.prior <- dunif(gf,0, 5, log=TRUE)
  sigma.prior <- dhalfcauchy(sigma, 5, log=TRUE)
  ### Log-Likelihood
  mu <- Sfactor3(obsx1,Er,0.0912,gi,gf,6,5,0)
  LL <- sum(dnorm(y, mu, sigma, log=TRUE))
  LP <- LL + Er.prior + gi.prior + gf.prior + sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor = c(LP, Er,gi,gf),
                   yhat=0,#rnorm(length(mu), mu, sigma),
                   parm=parm)
  return(Modelout)
}


Initial.Values <- rep(0.5,4)
########################  Laplace Approximation  ##########################

Fit <- LaplaceApproximation(Model, Initial.Values, Data=MyData,
                            Iterations=10000, Method="NM", CPUs=1)

print(Fit)


Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                     Iterations=5000, Status=100, Thinning=1,
                     Algorithm="HMC", Specs=list(epsilon=0.001, L=2, m=NULL))

Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=3000, Status=100, Thinning=1,
     Algorithm="RAM", Specs=list(alpha.star=0.234, B=NULL, Dist="N",
     gamma=0.66, n=0))
