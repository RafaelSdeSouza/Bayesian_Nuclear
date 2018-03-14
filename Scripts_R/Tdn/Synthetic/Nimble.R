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
#rm(list=ls())
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
library(dplyr);require(ggsci);require(ggmcmc);require(gsl);library(latex2exp)
require(runjags);require(gsl)
source("..//..//auxiliar_functions/jagsresults.R")
source("..//..//auxiliar_functions/theme_rafa.R")
source("..//..//auxiliar_functions/pair_wise_plot.R")
source("..//..//auxiliar_functions/Gamma3Hedp.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")
load.module("nuclear")


######################################################################
## ARTIFICIAL DATA GENERATION


N <- 150

#obsx1 <- runif(N,0,0.7)
obsx1 <- exp(runif(N,log(1e-3),log(1)))
sd <- 0.1
obsy1 <- rnorm(N, Sfactor3(obsx1,0.0912,0.0912,2.93,0.0794,6,5,0),sd=sd)

M <- 150
xx <- seq(min(obsx1),max(obsx1),length.out = M)
model.data <- list(obsy = obsy1,    # Response variable
                   obsx =  obsx1,   # Predictors
                   #                   erry = errobsy1,
                   N = N, # Sample size
                   M = M,
                   xx = xx
)
#

sfactorTdn <- nimbleRcall(function(obsx1 = double(1),
                              e1 = double(1),ex = double(1),gin = double(1),
                              gout = double(1),ri = double(1),rf = double(1),ue = double(1)){}, Rfun = 'Sfactor3',
                     returnType = double(1), envir = .GlobalEnv)


model <- nimbleCode({
   for (i in 1:150) {
    obsy[i] ~ dnorm(sfactorTdn(obsx[i], e1,0.0912, gin, gout,6,5,0),sd)
  }

  sd ~  dunif(0.01,1)
  e1 ~  dunif(0.01,1)
  gin ~  dunif(0.01,1)
  gout ~ dunif(1,4)  
})
inits <- list(e1 = runif(1,0.01,1),gout=runif(1,0.01,1),gin=runif(1,0.01,1)) 

Rmodel <- nimbleModel(model,data = model.data, inits = inits)
mcmcConf <- configureMCMC(Rmodel, monitors = c("e1", "gin", "gout","sd"))
mcmc_CL <- buildMCMC(mcmcConf)
samplesList <- runMCMC(mcmc_CL, niter = 1000, nchains = 3, inits = inits)


Rmcmc <- buildMCMC(model)
Rmodel$setData(model.data)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)



mcmc.output <- nimbleMCMC(model, data = model.data, inits = inits,
                          monitors = c("e1", "gin", "gout","sd"), thin = 10,
                          niter = 20000, nburnin = 1000, nchains = 3,
                          summary = TRUE, WAIC = TRUE)



