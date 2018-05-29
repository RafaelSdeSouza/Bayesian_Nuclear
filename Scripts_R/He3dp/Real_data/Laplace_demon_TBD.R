library(LaplacesDemon)
require(nuclear)
require(gsl)
require(nuclear);library(magrittr)
require(dplyr)
source("..//..//auxiliar_functions/jagsresults.R")
source("..//..//auxiliar_functions/theme_rafa.R")
source("..//..//auxiliar_functions/pair_wise_plot.R")
source("..//..//auxiliar_functions/Gamma3Hedp.R")
source("..//..//auxiliar_functions/table_reaction.R")
## for block updating [we do not need to center predictor variables]



######################################################################
## Read DATA 
ensamble <- read.csv("ensamble.csv",header = T) %>%
  mutate(Syst=replace(Syst,Syst==0.06,0.078))  %>% filter(E <= 0.5)  


re <- as.numeric(ensamble$dat)
Nre <- length(unique(ensamble$dat))
ik <- as.numeric(ensamble$invK)
Nik <- length(unique(ensamble$invK))
# Radius
# r_i = 6
# r_f = 5

# Literature
#  0.35779   # resonance energy
#  1.0085    # reduced width incoming
#  0.025425   # reduced width outgoing


N <- nrow(ensamble)
obsy <- ensamble$S    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat
set <- ensamble$dat
lab <- ensamble$invK
syst = c(0.03,unique(ensamble$Syst))
#syst <- syst[-3]

m <- as.numeric(lab)

mon.names <- c("LP","Er","gi","gf","ri","rf","ue[1]","ue[2]")
parm.names <- as.parm.names(list(Er="Er",gi="gi",gf="gf",ri="ri",rf="rf",ue = c(1,2),sigma="sigma"))
pos.Er <- 1
pos.gi <- 2
pos.gf <- 3
pos.ri <- 4
pos.rf <- 5
pos.ue <- grep("ue", parm.names)
pos.sigma <- 8

PGF <- function(Data) {
  Er <- runif(1,0.01,1)
  gi <- runif(1,0.01,4)
  gf <- runif(1,0.01,4)
  ri <- runif(1,4,6)
  rf <- runif(1,4,6)
  ue <- as.vector(runif(2,0,0.001))
  sigma <- runif(1)
  return(c(Er,gi,gf,ri,rf,ue, sigma))
}
MyData <- list(J=6, PGF=PGF, X = obsx, mon.names=mon.names,
               parm.names=parm.names,
               pos.Er=pos.Er,
               pos.gi = pos.gi,
               pos.gf=pos.gf,
               pos.ri=pos.ri,
               pos.rf=pos.rf,
               pos.sigma=pos.sigma,
               y = obsy,
               m = m,
               erry = erry)

Model <- function(parm, Data)
{
  ### Parameters
  Er <- interval(parm[Data$pos.Er], 1e-10, Inf)
  gi <- interval(parm[Data$pos.gi], 1e-10, Inf)
  gf <- interval(parm[Data$pos.gf], 1e-10, Inf)
  ri <- interval(parm[Data$pos.ri], 1e-10, Inf)
  rf <- interval(parm[Data$pos.rf], 1e-10, Inf)
  ue <- interval(parm[grep("ue", Data$parm.names)], 1e-10, Inf)
  sigma <- interval(parm[Data$pos.sigma], 1e-10, Inf)
  parm[Data$pos.Er] <- Er 
  parm[Data$pos.gi] <- gi
  parm[Data$pos.gf] <- gf 
  parm[Data$pos.ri] <- ri 
  parm[Data$pos.rf] <- rf 
  parm[Data$pos.ue] <- ue
  parm[Data$pos.sigma] <- sigma 
  ### Log(Prior Densities)
  Er.prior <- dhalfnorm(Er, 0.5, log=TRUE)
  gi.prior <- dhalfnorm(gi, 0.5, log=TRUE)
  gf.prior <- dhalfnorm(gf, 0.5, log=TRUE)
#  ri.prior <- dhalfcauchy(ri, 1, log=TRUE)
#  rf.prior <- dhalfcauchy(rf, 1, log=TRUE)
  ri.prior <- dnorm(ri, 5, 0.01,log=TRUE)
  rf.prior <- dnorm(rf,  5, 0.01, log=TRUE)
  ue.prior <- sum(dhalfnorm(ue, 0.01, log=TRUE))
  sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)

  
  ### Log-Likelihood
  mu <- sfactor3Hedp_5p(obsx,Er,gi,gf,ri,rf,ue[Data$m])
 # epsilon <- dnorm(obsy, 0, erry, log=TRUE)
  LL <- sum(dnorm(obsy, mu, sigma, log=TRUE)) 

  
  LP <- LL + Er.prior + gi.prior + gf.prior + sigma.prior + ri.prior + rf.prior +
    ue.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor = c(LP, Er,gi,gf,ri,rf,ue),
                   yhat=0,#rnorm(length(mu), mu, sigma),
                   parm=parm)
  return(Modelout)
}

library(compiler)
Model <- cmpfun(Model) 

Initial.Values <- rep(runif(8,0.01,1))
########################  Laplace Approximation  ##########################


FitAFSS <- LaplacesDemon(Model, Data=MyData, Initial.Values, 
                   Covar=NULL, Iterations=5000, Status=100, Thinning=10,
                   Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))



FitAFSS <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                         Covar=FitAFSS$Covar, Iterations=25000, Status=1262, Thinning=50,
                         Algorithm="AFSS", Specs=list(A=5000, B=NULL, m=FitAFSS$Specs$m,
                                                      n=5000, w=FitAFSS$CovarDHis[nrow(FitAFSS$CovarDHis),]))

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


