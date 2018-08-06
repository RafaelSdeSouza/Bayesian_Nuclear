# 3Hedp analysis
#
# purpose: Real  DATA
#
# - 5 parameters are assumed: Er, gamma_d^2, gamma_n^2 [e1, gin, gout]
#
# - uses the function sfactorHe3dp(obsx1[i], e1, gin, gout), which
#   is a C++ version of a Fortran code that includes Coulomb wave
#   function calculations; JAGS has been recompiled with this C++ function
#
######################################################################
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
library(dplyr);library(BayesianTools)
## for block updating [we do not need to center predictor variables]



######################################################################
## Read DATA 
ensamble <- read.csv("ensamble_Tdn_extra.csv",header = T) %>%  filter(E <= 0.3)
#filter(dat!= "Mag75")
#%>% filter(E <= 0.5) %>%   filter(dat!= "Arn53") %>%
# droplevels(ensamble$dat)


re <- as.numeric(ensamble$dat) 
Nre <- length(unique(ensamble$dat))
Nik <- length(unique(ensamble$invK))
# Radius
# r_i = 6
# r_f = 5




N <- nrow(ensamble)
obsy <- ensamble$S    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat
set <- ensamble$dat
lab <- ensamble$invK
syst = c(unique(ensamble$Syst))
i = seq(1:N)


likelihood <- function(par){
  e1 = par[1]
  gin = par[2]
  gout = par[3]
  sigmax = par[4]
  scale = par[5:9]
#  y = par[10:133]
  
  llscale = sum(dlnorm(scale,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  
  lly <- sum(dnorm(obsy,scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5), sd = sigmax,  log = T))
    return(llscale   + lly)
  
}



setup <- createBayesianSetup(likelihood = likelihood,
lower = c(0.001,0.001,0.001,0.001,rep(0.5,5)),
upper = c(1,2,2,5,rep(1.5,5)))
settings <- list(iterations = 100000,adaptation=5000,
                 burnin = 20000, message=T)

res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")
summary(res)
tracePlot(sampler = res, thin = 20, start = 10000, whichParameters = c(1,2,3))

correlationPlot(res )


correlationPlot(out)







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

