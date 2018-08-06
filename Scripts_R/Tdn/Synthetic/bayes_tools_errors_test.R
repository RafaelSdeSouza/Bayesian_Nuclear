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
library(BayesianTools)
# Data
set.seed(1056)                   # set seed to replicate example
nobs = 100                   # number of obs in model

sdobsx <- 1.25
truex <- rnorm(nobs,0,2.5)       # normal variable
#errx <- rnorm(nobs, 0, sdobsx)
#obsx <- truex + errx 
obsx <- truex

beta1 <- -3
beta2 <- 8
sdy <- 1.5
sdobsy <- 0.5 # Variance to simulate reported errors

erry <- rnorm(nobs, 0, sdobsy) # reported errors
truey <- rnorm(nobs,beta1 + beta2*truex,sdy)
obsy <- truey + erry

K <- 2
i <- seq(4,nobs+3)
# Level of  mistake in the reported errors
Lambda <- 0.5

# Create a general prior distribution by specifying an arbitrary density function and a
# corresponding sampling function
density = function(par){
  d1 = dnorm(par[1], mean= 0, sd = 10, log =TRUE)
  d2 = dnorm(par[2], mean= 0, sd = 10, log =TRUE)
  d3 = dunif(par[3], 0,100, log =TRUE)
  d4 = dunif(par[4:(nobs + 3)], obsy - 3*abs(erry), obsy + 3*abs(erry), log =TRUE)
  return(d1 + d2 + d3 + d4)
}

# The sampling is optional but recommended because the MCMCs can generate automatic starting
# conditions if this is provided
sampler = function(n=1){
  d1 = rnorm(n, mean= 0, sd = 10)
  d2 = rnorm(n, mean= 0, sd = 10)
  d3 = runif(n, 0,100)
  d4 = runif(nobs,obsy - 3*abs(erry), obsy + 3*abs(erry))
  return(cbind(d1,d2,d3,d4))
}

prior <- createPrior(density = density, sampler = sampler, lower = c(-100,-100,0,obsy - 3*abs(erry)),
                     upper = c(100,100,100,obsy + 3*abs(erry)), best = NULL)


likelihood <- function(par){
#  
  a = par[1]
  b = par[2]
  sy <- par[3]
  y = par[4:(nobs + 3)]
  
  lly <- sum(dnorm(y,(a + b*obsx), sd = sy,  log = T))
  llobs = sum(dnorm(obsy,y, sd = abs(erry),log = T))
  
  #llobs = sum(dnorm(obsy,scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5),sd = sigmax,log = T)) 
  
#  obsy[i] ~ dnorm(y[i], pow(erry[i], -2))
#  y[i] ~ dnorm(scale[re[i]]*sfactorTdn(obsx[i], e1, ex,gin, gout,ri,rf,0),pow(tau, -2))
  
#  llytrue = sum(dnorm(scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5),
#                      sd = sigmax, log = T))
#  llobs = sum(dnorm(obsy,mean = llytrue,sd = erry,log = T))

  
  #llobs = sum(dnorm(scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5) - obsy,sd = sigmax,log = T))  
  
  return( lly + llobs)
  
}



setup <- createBayesianSetup(likelihood = likelihood,lower = c(-100,-100,0,obsy - 3*abs(erry)),
                             upper = c(100,100,100,obsy + 3*abs(erry)))

settings <- list(iterations = 100000,
                 burnin = 50000, message = T)

out <- runMCMC(bayesianSetup = setup, settings = settings,sample="DREAMzs")
summary(out)
tracePlot(sampler = out, thin = 10, start = 1000, whichParameters = c(1,2,3))
marginalPlot(out,whichParameters = c(1,2,3))