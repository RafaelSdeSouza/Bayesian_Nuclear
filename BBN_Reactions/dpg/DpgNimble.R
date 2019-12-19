##########################################
# Dpg Posterior Simulation
# Attempt using Nimble
# Done by HK
# Adapted from different parts of code from
# Christian and Rafael
##########################################
# preparation: remove all variables from the work space
rm(list=ls())
library(nimble)
library(magicaxis)
library(MASS)
library(ggplot2)
library(magrittr)
## For block updating [we do not need to center predictor variables]
load.module("glm")

#######################################
# Importing Data (Dpgdat.csv)
#######################################

#setwd("C:/Users/Hong Kiat/Desktop/UROPS Materials/UROPS 1 - Thermonuclear Reaction Rates/Working Folder/Nimble/Dpg")
dat <- read.csv("Dpgdat.csv")
re <- as.numeric(dat$lab) # Change the label to a numeric vector
# Note that the numbers are assigned by "sorting" the labels in alphabetical order
Nre <- length(unique(dat$lab))
# Unique removes duplicated vector, we want to know how many groups of 
# data are there
N <- nrow(dat) # Total No of data sets
obsy <- dat$S    # Response variable in MeV
obsx <-  dat$E   # Predictors
erry <- dat$Stat # Error in MeV
set <- dat$lab # Get the labels as a vector

# Inputting the quoted systematic uncertainty:
syst = c(log(1.08),log(1.045),log(1.09),log(1.08))
# In accordance to 
# bys08a = 1, <8% ~ 8% = log(1.08)
# cas02a = 2, 4.5% = log(1.045)
# ma97a = 3, 9% = log(1.09)
# sch97a = 4, 8% = log(1.08)

################
# Nimble Function, myapprox
###############
# E is in MeV, S is in eV b; convert to MeV and MeV b, respectively
theory <- read.table("Marcucci2005.dat", header=FALSE)

interp.x <- theory[,1] 
interp.y <- theory[,2]*1e-6


myapproxR <- nimbleFunction(
  run = function(x = double(1,100000), y = double(1,100000), xout = double(0)){
    ind <- max(2,min(which(x>xout)))
    ans <- (y[ind]-y[ind-1])/(x[ind]-x[ind-1])*(xout-x[ind-1]) + y[ind-1]
    returnType(double(0))
    return(ans)
  }
)

samplerCode <-  nimbleCode({
  for (i in 1:N) {
    obsy[i] ~ dnorm(ya[i], sd = erry[i])
    # Propagating the errors in observations of y
    ya[i] ~ dnorm(y.norm[re[i]]*mut[i],sd = y.scat[re[i]])
    # Temp mean (mut) multiplied by normalising constant, 
    # re[i] refers to the index that we would use to refer to the category
    mut[i] <- a.scale*(myapproxR(interp.x[1:100000],interp.y[1:100000],obsx[i]))
    # Interpolate and scale
  }
  
  ## PRIORS
  # Scaling factor of theory 
  a.scale ~ T(dnorm(0, sd = 5),0,Inf)  
  
  for (k in 1:Nre){
    # Systematic Uncertainty as a highly informative prior
    y.norm[k] ~ dlnorm(log(1.0),sd=syst[k])
    y.scat[k] ~  T(dnorm(mt, sd=5),0,Inf)
  }
  mt ~  T(dnorm(0, sd=5),0,Inf) # Hyperprior
})

samplerData <- list(obsy = obsy    # Response variable
)

samplerConst <- list(N = N, # Sample size
                     Nre = Nre, 
                     re = re, # This is used to "iterate"
                     erry = erry,
                     syst = syst,
                     interp.x = interp.x,
                     interp.y = interp.y,
                     obsx = obsx
)

samplerInits <- list(a.scale = 1,
                     y.norm = rep(1,Nre),
                     y.scat = rep(0.01,Nre),
                     mt = 0,
                     ya = approx(interp.x,interp.y,obsx)$y
)
# Here, we are using the R's version of approx <-> interp.lin in JAGS

# All the sampled parameters have to the FULLY SPECIFIED

n.chains = 1
n.iter = 5000
n.burnin = 1000


###############################################################
# (Alternative) If invoking Nimble MCMC stepwise (but more customisable)
##############################################################
ourmodel <- nimbleModel(code = samplerCode, constants = samplerConst,
                        data = samplerData, inits = samplerInits, check = FALSE
)
compileNimble(ourmodel)
# Always compile the model after you are done setting up with it

ourmodel$calculate()
# Calculate the log likelihood (logProb). If the value is not NA,
# we have successfully initialised the model (FULLY)
# One iteration: simulate -> new values -> calculate

conf <- configureMCMC(ourmodel,print=TRUE)
# print = TRUE tells you what kind of samplers are being used for each stochastic node

conf$addMonitors(c('a.scale',
                   'y.norm','y.scat'))
# Add the parameters to monitor throughout the MCMC process

samplerMCMC <- buildMCMC(conf)
# Note that if no configuration is required, we can directly jump towards buildMCMC.
# But you know, then theres no point in using this entire chunk of code
# where you can just run nimbleMCMC

compiledMCMC <- compileNimble(samplerMCMC,project = ourmodel)
# Complie the configured model again once we are done with configuring it;
# Just before we perform the MCMC runs
# Can include project = ourmodel; to include all the nimbleFunctions, nimbleModels that
# you would want to inherit the functions and models from 
# resetFunctions = TRUE; if you would want to reset all the previously created functions
# in order to addd the new MCMC

system.time(
  mcmcChain <- runMCMC(compiledMCMC,niter = n.iter, nchain = n.chains, nburnin = n.burnin,
                       setSeed=15,samplesAsCodaMCMC = TRUE)
)