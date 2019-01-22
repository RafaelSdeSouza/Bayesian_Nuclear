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
library(rjags);library(R2jags);library(mcmcplots)
require(RcppGSL);require(ggplot2);require(ggthemes)
require(nuclear);library(magrittr);library(wesanderson)
library(dplyr);require(ggsci);require(ggmcmc);require(plyr);library(latex2exp)
source("..//..//auxiliar_functions/jagsresults.R")
source("..//..//auxiliar_functions/theme_rafa.R")
source("..//..//auxiliar_functions/pair_wise_plot.R")
source("..//..//auxiliar_functions/Gamma3Hedp.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")
load.module("nuclear")


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
syst = c(unique(ensamble$Syst))
#syst <- syst[-3]
library(BayesianTools)

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



likelihood <- function(par){
  e1 = par[1]
  gin = par[2]
  gout = par[3]
  sigmax = par[4]
  scale = par[5:9]
  y = par[10:]
  
  llRandom = sum(dlnorm(scale,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  lly <- sum(dnorm(y,mean = scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5), sd = sigmax,  log = T))
  llobs = sum(dnorm(obsy,mean = y,sd = erry,log = T))
  
  #llobs = sum(dnorm(obsy,scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5),sd = sigmax,log = T)) 
  
#  obsy[i] ~ dnorm(y[i], pow(erry[i], -2))
#  y[i] ~ dnorm(scale[re[i]]*sfactorTdn(obsx[i], e1, ex,gin, gout,ri,rf,0),pow(tau, -2))
  
#  llytrue = sum(dnorm(scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5),
#                      sd = sigmax, log = T))
#  llobs = sum(dnorm(obsy,mean = llytrue,sd = erry,log = T))

  
  #llobs = sum(dnorm(scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5) - obsy,sd = sigmax,log = T))  
  
  return(llRandom + llobs + lles)
  
}



setup <- createBayesianSetup(likelihood = likelihood,
lower = c(0.001,0.001,0.001,0.001,rep(0.5,5),0.001),
upper = c(1,2,2,5,rep(1.5,5),1))
settings <- list(iterations = 100000,
                 burnin = 15000, message=T)

res <- runMCMC(bayesianSetup = setup, settings = settings)
summary(res)
tracePlot(sampler = res, thin = 10, start = 20000, whichParameters = c(1,2,3))
