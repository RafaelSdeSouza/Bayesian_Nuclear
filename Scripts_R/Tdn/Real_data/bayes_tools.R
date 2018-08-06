# preparation: remove all variables from the work space
rm(list=ls())
set.seed(123)
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb;

######################################################################
# import packages
require(RcppGSL);require(ggplot2);require(ggthemes)
require(nuclear);library(magrittr);
<<<<<<< HEAD
library(dplyr);library(BayesianTools)
## for block updating [we do not need to center predictor variables]

=======
library(dplyr);require(lessR);library(BayesianTools)
>>>>>>> c6bf149046631a4132f8a3b51e0cc79e250d3e9a


######################################################################
## Read DATA 
ensamble <- read.csv("ensamble_Tdn_extra.csv",header = T) %>%  filter(E <= 0.5)
#filter(dat!= "Mag75")
#%>% filter(E <= 0.5) %>%   filter(dat!= "Arn53") %>%
# droplevels(ensamble$dat)


re <- as.numeric(ensamble$dat) 
Nre <- length(unique(ensamble$dat))
Nik <- length(unique(ensamble$invK))
<<<<<<< HEAD
# Radius
# r_i = 6
# r_f = 5



=======
>>>>>>> c6bf149046631a4132f8a3b51e0cc79e250d3e9a

N <- nrow(ensamble)
obsy <- ensamble$S    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat
set <- ensamble$dat
lab <- ensamble$invK
syst = c(unique(ensamble$Syst))
<<<<<<< HEAD
i = seq(1:N)
=======
>>>>>>> c6bf149046631a4132f8a3b51e0cc79e250d3e9a


likelihood <- function(par){
  e1 = par[1]
  gin = par[2]
  gout = par[3]
  sigmax = par[4]
  scale = par[5:9]
<<<<<<< HEAD
#  y = par[10:133]
  
  llscale = sum(dlnorm(scale,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  
  lly <- sum(dnorm(obsy,scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5), sd = sigmax,  log = T))
    return(llscale   + lly)
=======
  y = par[10:(N + 9)]
  
  llRandom = sum(dlnorm(scale,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  lly <- sum(dnorm(y,mean = scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5), sd = sigmax,  log = T))
  llobs = sum(dnorm(obsy,mean = y,sd = erry,log = T))
  return(llRandom + llobs + lly)
>>>>>>> c6bf149046631a4132f8a3b51e0cc79e250d3e9a
  
}



setup <- createBayesianSetup(likelihood = likelihood,
<<<<<<< HEAD
lower = c(0.001,0.001,0.001,0.001,rep(0.5,5)),
upper = c(1,2,2,5,rep(1.5,5)))
settings <- list(iterations = 100000,adaptation=5000,
                 burnin = 20000, message=T)
=======
lower = c(0.001,0.001,0.001,0.001,rep(0.5,5),obsy - 5*abs(erry)),
upper = c(1,2,2,5,rep(1.5,5),obsy + 5*abs(erry)),
names = c("e0","gd2","gn2","sigma",to("scale", 5),to("y", N)))
settings <- list(iterations = 200000,
                 burnin = 50000, message = T)
>>>>>>> c6bf149046631a4132f8a3b51e0cc79e250d3e9a

res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")
summary(res)
<<<<<<< HEAD
tracePlot(sampler = res, thin = 20, start = 10000, whichParameters = c(1,2,3))

correlationPlot(res )


correlationPlot(out)
=======
tracePlot(sampler = res, thin = 10, start = 5000, whichParameters = c(1,2,3,4,5,6,7,8,9))



>>>>>>> c6bf149046631a4132f8a3b51e0cc79e250d3e9a




<<<<<<< HEAD



density = function(par){
  d1 = dnorm(par[1], 0.001,1, log =TRUE)
  d2 = dunif(par[2], 0.001,3, log =TRUE)
  d3 = dunif(par[3], 0.001,3, log =TRUE)
  return(d1 + d2 + d3)
}

sampler = function(n=1){
  d1 = rnorm(n, 0.001,1)
  d2 = rnorm(n, 0.001,1)
  d3 = rnorm(n, 0.001,1)
  return(cbind(d1,d2,d3))
}
prior <- createPrior(density = density, sampler = sampler, 
                     lower = c(0.001,0.001,0.001), upper = c(10,10,10), best = NULL)

=======
codaObject = getSample(res, start = 500, coda = TRUE)

as.mcmc(codaObject)
getmcmc_var <- function(outjags=outjags,vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
getmcmc_var(codaObject,vars = c("par 1","par 2","par 3","par 4"))
>>>>>>> c6bf149046631a4132f8a3b51e0cc79e250d3e9a
