
require(LaplacesDemon)
require(nuclear)

## ARTIFICIAL DATA GENERATION

N <- 100

#obsx1 <- runif(N,0,0.7)
tau <- 1
obsx1 <- exp(runif(N,log(1e-3),log(1)))
errobsy1 <- runif(N,0.1,0.2)
y1 <- rnorm(N, sfactor3Hedp(obsx1 ,0.35779,1.0085,0.025425),errobsy1)
obsy1 <- rnorm(N,y1,tau)

M <- 150
xx <- seq(min(obsx1),max(obsx1),length.out = M)
model.data <- list(obsy = obsy1,    # Response variable
                   obsx =  obsx1,   # Predictors
                   erry = errobsy1,
                   N = N, # Sample size
                   M = M,
                   xx = xx
)



y <- obsy1
X <- obsx1
mon.names <- "LP"
parm.names <- as.parm.names(list(Er=0, gin=0, gout=0, sigma=0))
pos.Er <- grep("Er", parm.names)
pos.gin <- grep("gin", parm.names)
pos.gout <- grep("gout", parm.names)
pos.sigma <- grep("sigma", parm.names)


MyData <- list(X=X, mon.names=mon.names,
               parm.names=parm.names, pos.Er=pos.Er, pos.gin = pos.gin, pos.gout = pos.gout, pos.sigma=pos.sigma, y=y)





Model <- function(parm, Data)
{
  ### Parameters
  Er <- parm[Data$pos.Er]
  gin <- parm[Data$pos.gin]
  gout <- parm[Data$pos.gout]
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  ### Log-Prior
  Er.prior <- dhalfcauchy(25)
  gin.prior <- dhalfcauchy(25)
  gout.prior <- dhalfcauchy(1)
  sigma.prior <- dgamma(sigma, 25, log=TRUE)
  ### Log-Likelihood
  mu <- sfactor3Hedp(obsx1 ,Er,gin,gout)
  LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
  ### Log-Posterior
  LP <- LL + Er.prior + gin.prior + gout.prior + sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
                   yhat=rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)
}


Initial.Values <- c(rep(1,3), 1)



Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                     Covar=NULL, Iterations=10000, Status=100, Thinning=1,
                     Algorithm="RJ")





