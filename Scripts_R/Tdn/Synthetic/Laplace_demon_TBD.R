library(LaplacesDemon)
require(nuclear)
require(gsl)
N <- 200

#obsx1 <- runif(N,0,0.7)
obsx1 <- exp(runif(N,log(1e-3),log(1)))
sd <- 1
# Artificial values, here we will just fit 3 parameters Er, gi, gf
y <- rnorm(N, sfactorTdn_5p(obsx1,0.0912,2.93,0.0794,6,5),sd = sd)

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
  Er.prior <- dhalfnorm(Er, scale=1, log=TRUE)
  gi.prior <- dhalfnorm(gi, scale=1, log=TRUE)
  gf.prior <- dhalfnorm(gf, scale=1, log=TRUE)
  sigma.prior <- dhalfcauchy(sigma, 5, log=TRUE)
  ### Log-Likelihood
  mu <- sfactorTdn_5p(obsx1,Er,gi,gf,6,5)
  LL <- sum(dnorm(y, mu, sigma, log=TRUE))
  LP <- LL + Er.prior + gi.prior + gf.prior + sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor = c(LP, Er,gi,gf),
                   yhat=0,#rnorm(length(mu), mu, sigma),
                   parm=parm)
  return(Modelout)
}


Initial.Values <- rep(runif(4,0.001,1))
########################  Laplace Approximation  ##########################


FitAFSS <- LaplacesDemon(Model, Data=MyData, Initial.Values,Chains=2, 
                    Iterations=50000, Status=100, Thinning=10,
                     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))

plot(FitAFSS, BurnIn=25000, MyData, PDF=FALSE, Parms=c("Er","gi","gf","sigma"))


FitMWG <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                          Covar=NULL, Iterations=20000, Status=100, Thinning=1,
                          Algorithm="MWG", Specs=list(B=NULL))
                     
plot(FitMWG, BurnIn=1000, MyData, PDF=FALSE, Parms=c("Er","gi","gf","sigma"))

