library(LaplacesDemon)
require(nuclear)
require(gsl)

## Read DATA 
ensamble <- read.csv("ensamble_Tdn.csv",header = T)

# data
N <- nrow(ensamble)
obsy <- ensamble$S    # Response variable
obsx <- ensamble$E   # Predictors
erry <- ensamble$Stat
set <- ensamble$dat
lab <- ensamble$invK
syst = c(unique(ensamble$Syst))


mon.names <- c("LP","Er","gi","gf","ri","rf")
parm.names <- c("Er","gi","gf","ri","rf","sigma")
pos.Er <- 1
pos.gi <- 2
pos.gf <- 3
pos.ri <- 4
pos.rf <- 5
pos.sigma <- 6

PGF <- function(Data) {
  Er <- runif(1,0.01,1)
  gi <- runif(1,0.01,4)
  gf <- runif(1,0.01,4)
  ri <- runif(1,3,6)
  rf <- runif(1,3,6)
  sigma <- runif(1)
  return(c(Er,gi,gf,ri,rf, sigma))
}
MyData <- list(J=6, PGF=PGF, X = obsx, mon.names=mon.names,
               parm.names=parm.names,
               pos.Er=pos.Er,
               pos.gi = pos.gi,
               pos.gf=pos.gf,
               pos.ri=pos.ri,
               pos.rf=pos.rf,
               pos.sigma=pos.sigma,
               y=obsy)

Model <- function(parm, Data)
{
  ### Parameters
  Er <- interval(parm[Data$pos.Er], 1e-100, Inf)
  gi <- parm[Data$pos.gi]
  gf <- parm[Data$pos.gf]
  ri <- parm[Data$pos.ri]
  rf <- parm[Data$pos.rf]
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma 
  ### Log(Prior Densities)
  Er.prior <- dhalfcauchy(sigma, 2.5, log=TRUE)
  gi.prior <- dhalfcauchy(sigma, 2.5, log=TRUE)
  gf.prior <- dhalfcauchy(sigma, 2.5, log=TRUE)
  ri.prior <- dhalfcauchy(sigma, 2.5, log=TRUE)
  rf.prior <- dhalfcauchy(sigma, 2.5, log=TRUE)
  sigma.prior <- dhalfcauchy(sigma, 2.5, log = TRUE)
  ### Log-Likelihood
  mu <- sfactorTdn_5p(obsx,Er,gi,gf,ri,rf)
  LL <- sum(dnorm(obsy, mu, sigma, log=TRUE))
  LP <- LL + Er.prior + gi.prior + gf.prior + sigma.prior + ri.prior + rf.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor = c(LP, Er,gi,gf,ri,rf),
                   yhat=0,#rnorm(length(mu), mu, sigma),
                   parm=parm)
  return(Modelout)
}

library(compiler)
Model <- cmpfun(Model) 

Initial.Values <- rep(runif(6,1,5))
########################  Laplace Approximation  ##########################


FitLA <- LaplaceApproximation(Model, Data=MyData, Initial.Values, 
                         Iterations=2000)

plot(FitLA, MyData, PDF=FALSE)
caterpillar.plot(FitLA, Parms=c("Er","gi","gf"))

FitAFSS <- LaplacesDemon(Model, Data=MyData, Initial.Values, 
                    Iterations=10000, Status=100, Thinning=10,
                     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))



FitAFSS <- LaplacesDemon(Model, Data = MyData, Initial.Values,
                         Covar = FitAFSS$Covar, Iterations = 130000, Status = 1449,
                         Thinning = 260, Algorithm = "AFSS", 
                         Specs = list(A = 5000, B = NULL, m = FitAFSS$Specs$m,
                         n = 5000, w = FitAFSS$CovarDHis[nrow(FitAFSS$CovarDHis),]))                                                                           +                          n=5000, w=FitAFSS$CovarDHis[nrow(FitAFSS$CovarDHis),]))

caterpillar.plot(FitAFSS, Parms=c("Er","gi","gf"))
LaplacesDemon(FitAFSS)

Consort(FitAFSS)
plot(FitAFSS, BurnIn=10000, MyData,  Parms=c("Er","gi","gf","sigma"))


FitMWG <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                          Covar=NULL, Iterations=20000, Status=100, Thinning=1,
                          Algorithm="MWG", Specs=list(B=NULL))
                     
plot(FitMWG, BurnIn=1000, MyData, PDF=FALSE, Parms=c("Er","gi","gf","sigma"))


FitHARM <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                        Covar=NULL, Iterations=20000, Status=100, Thinning=1,
                        Algorithm="HARM")

plot(FitHARM, BurnIn=1000, MyData, PDF=FALSE, Parms=c("Er","gi","gf","sigma"))


