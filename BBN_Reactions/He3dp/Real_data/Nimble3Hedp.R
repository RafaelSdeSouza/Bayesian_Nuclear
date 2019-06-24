##########################################
# Dpg Posterior Simulation
# Attempt using Nimble
# Done by HK
# Adapted from different parts of code from
# Christian and Rafael
##########################################
# preparation: remove all variables from the work space
rm(list=ls())
library(nuclear)
library(nimble)
library(magicaxis)
library(MASS)
library(ggplot2)
library(magrittr)

#######################################
# Importing Data (Dpgdat.csv)
#######################################

setwd("C:/Users/Hong Kiat/Desktop/UROPS Materials/UROPS 1 - Thermonuclear Reaction Rates/Working Folder/Nimble R Matrix")
dat <- read.csv("3Hedp.csv")
# Trying to reconcile naming conflicts; Order of data = Number
re <- as.numeric(dat$no) # Change the label to a numeric vector
# Note that the numbers are assigned by "sorting" the labels in alphabetical order
Nre <- length(unique(dat$no))
# Unique removes duplicated vector, we want to know how many groups of 
# data are there
N <- nrow(dat) # Total No of data sets
obsy <- dat$S    # Response variable in MeV
obsx <-  dat$E   # Predictors
erry <- dat$Stat # Error in MeV
set <- dat$lab # Get the labels as a vector

sepab <- c(1,1,1,1,2,1,2)
# 1 - a, 2 - b

# Inputting the quoted systematic uncertainty:
syst = c(log(1.043),log(1.078),log(1.039),log(1.034),log(1.055),log(1.030),log(1.030))
# In accordance to 
# Gei99 - 4.3%
# Kra87 - 7.8%
# Moe80 - 3.9%
# Zhi77 - 3.4%
# Cos00 - 5.5%
# Ali01_a - 3.0%
# Ali01_b - 3.0%

################
# Nimble Function, sfact (using nuclear package)
###############

sfact <- nimbleRcall(function(ecm = double(1,212),e0 = double(0),
                                  gi = double(0), gf = double(0),
                                  ri = double(0), rf = double(0),
                                  ue = double(1,212)){}, 
                         Rfun = 'sfactorHe3dp',
                         returnType = double(1,212))

samplerCode <-  nimbleCode({
  for (i in 1:N) {
    obsy[i] ~ dnorm(re[i]*mut[i], sd = (y.scat[re[i]]^2 + erry[i]^2)^(0.5) 
    # Propagating the errors in observations of y
    # ya[i] ~ dnorm(y.norm[re[i]]*mut[i],sd = y.scat[re[i]])
    ## Skipping ya to reduce stochastic nodes
    # Temp mean (mut) multiplied by normalising constant, 
    # re[i] refers to the index that we would use to refer to the category
    # Interpolate and scale
  }
  
  mut <- sfact(obsx,e0,ga,gb,ra,rb,ue[sepab[re]])
  
  ## PRIORS - Get it to work first, understand later
  # e0:       eigenenergy
  # ga,gb:    initial reduced width, final reduced width [ga,gb = gamma^2]
  
  # energy eigenvalue
  e0 ~ T(dnorm(0, sd=1),0,Inf)       # positive since we see sigma peak 
  
  # reduced widths
  ga ~ T(dnorm(0, sd = 1),0,Inf)      
  gb ~ T(dnorm(0, sd = 2),0,Inf)     
  ##  ga ~ dnorm(0.0, pow(wl_d, -2))T(0,)      
  ##  gb ~ dnorm(0.0, pow(wl_p, -2))T(0,)      
  
  # channel radii 
  ##  ra ~ dunif(2.5, 8.0)
  ##  rb ~ dunif(2.5, 8.0)
  ###  ra ~ dnorm(3.25, pow(1.0, -2))T(0,)
  ###  rb ~ dnorm(5.77, pow(2.0, -2))T(0,)
  ra ~ T(dnorm(6, sd=0.00001),0,Inf)
  rb ~ T(dnorm(5, sd=0.00001),0,Inf)
  
  # screening potential:
  ue[1] ~ T(dnorm(0.0,sd=0.001),0,Inf)     # certainly less than 10 keV
  ue[2] ~ T(dnorm(0.0,sd=0.001),0,Inf)     # certainly less than 10 keV
  ###  uea ~ dnorm(200e-6, pow(80e-6, -2))T(0,)     # certainly less than 10 keV
  ####  uea ~ dnorm(180e-6, pow(80e-10, -2))T(0,)     # certainly less than 10 keV
  ####  ueb ~ dnorm(121e-6, pow(50e-6, -2))T(0,)     # certainly less than 10 keV
  
  for (k in 1:Nre){
    # Systematic Uncertainty as a highly informative prior
    y.norm[k] ~ dlnorm(log(1.0),sd=syst[k])
    y.scat[k] ~  T(dnorm(mt, sd=5),0,Inf)
  }
  mt ~ T(dnorm(0, sd=5),0,Inf) # Hyperprior
})

#  Wigner limit: wl = hbar^2/(m_red a_c^2) = 41.80159/(M_red a_c^2)
#  Not used in MCMC chain, brought outside of nimble code
#  deuteron channel: wl_d = 41.80159/(1.207266 a_c^2) = 34.6250/a_c^2
#  proton channel:   wl_p = 41.80159/(0.804710 a_c^2) = 51.9461/a_c^2
#  wl_d <- 34.6250*ra^(-2)
#  wl_p <- 51.9461*rb^(-2)

samplerData <- list(obsy = obsy    # Response variable
)

samplerConst <- list(N = N, # Sample size
                     Nre = Nre, 
                     re = re, # This is used to "iterate"
                     erry = erry,
                     syst = syst,
                     obsx = obsx,
                     sepab = sepab
)

samplerInits <- list(y.norm = rep(1,Nre),
                     y.scat = rep(0.01,Nre),
                     # ya = sfact(obsx,1,1,2,6,5,c(0.1,0.1)),
                     mt = 0, e0 = 0.5, ga = 4, gb = 0.1, ra = 1, rb = 1,
                     ue = c(0.001,0.001)
)
# Here, we are using the R's version of approx <-> interp.lin in JAGS

# All the sampled parameters have to the FULLY SPECIFIED

n.chains = 3
n.iter = 20000
n.burnin = 2000
# Note to self: Dropping stochastic nodes do speed things up!


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

conf$addMonitors(c('e0','ga','gb','ra','rb','ue',
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


pdf("try.pdf")
plot(mcmcChain)
dev.off()