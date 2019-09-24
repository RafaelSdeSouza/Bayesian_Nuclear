# preparation: remove all variables from the work space
rm(list=ls())
# import jags package
library(nimble)
require(gsl)
require(RcppGSL)
require(nuclear)
require(forcats)
require(utils)
require(stats)
require(methods)

path <- getwd()
setwd(path)
sigma7Benp7mod <- function(ecm,
                           e0_1, ga_1, gb_1, ra_1, rb_1,
                           e0_2, ga_2, gb_2, ra_2, rb_2,
                           e0_3, ga_3, gb_3, ra_3, rb_3,
                           e0_4, ga_4, gb_4, ra_4, rb_4,
                           e0_5, ga_5, gb_5, ra_5, rb_5,
                           e0_6, ga_6, gb_6, ra_6, rb_6,
                           e0_7, ga_7, gb_7, ra_7, rb_7){

  SF1 <-  sigma7Benp(ecm, e0_1, ga_1, gb_1, ra_1, rb_1, jr = 2, la = 0, lb = 0)
  SF2 <-  sigma7Benp(ecm, e0_2, ga_2, gb_2, ra_2, rb_2, jr = 3, la = 1, lb = 1)
  SF3 <-  sigma7Benp(ecm, e0_3, ga_3, gb_3, ra_3, rb_3, jr = 3, la = 1, lb = 1)
  SF4 <-  sigma7Benp(ecm, e0_4, ga_4, gb_4, ra_4, rb_4, jr = 1, la = 0, lb = 0)
  SF5 <-  sigma7Benp(ecm, e0_5, ga_5, gb_5, ra_5, rb_5, jr = 4, la = 3, lb = 3)
  SF6 <-  sigma7Benp(ecm, e0_6, ga_6, gb_6, ra_6, rb_6, jr = 2, la = 1, lb = 1)
  SF7 <-   sigma7Benp(ecm, e0_7, ga_7, gb_7, ra_7, rb_7, jr = 0, la = 1, lb = 1)
  SF <- SF1 + SF2 + SF3 + SF4 + SF5 + SF6 + SF7
  return(SF = SF)
}

# Vectorised sigma function
sigmaBe7mod <- nimbleRcall(function(ecm = double(1),
                                    e0_1 = double(0), ga_1 = double(0), gb_1 = double(0),
                                    ra_1 = double(0), rb_1 = double(0),
                                    e0_2 = double(0), ga_2 = double(0), gb_2 = double(0),
                                    ra_2 = double(0), rb_2 = double(0),
                                    e0_3 = double(0), ga_3 = double(0), gb_3 = double(0),
                                    ra_3 = double(0), rb_3 = double(0),
                                    e0_4 = double(0), ga_4 = double(0), gb_4 = double(0),
                                    ra_4 = double(0), rb_4 = double(0),
                                    e0_5 = double(0), ga_5 = double(0), gb_5 = double(0),
                                    ra_5 = double(0), rb_5 = double(0),
                                    e0_6 = double(0), ga_6 = double(0), gb_6 = double(0),
                                    ra_6 = double(0), rb_6 = double(0),
                                    e0_7 = double(0), ga_7 = double(0), gb_7 = double(0),
                                    ra_7 = double(0), rb_7 = double(0)
){},
Rfun = "sigma7Benp7mod", returnType = double(1))

######################################################################
## DATA SETS
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the
# latter is the individual statistical error of each datum [i];
# energy is in MeV, sqrt(Ecm)*sigma is in sqrt(MeV)b


Be7np <- read.csv("Be7np.csv")


Be7np$dat <- as.factor(Be7np$dat)
Be7np$dat <- fct_relevel(Be7np$dat, "Dam18","Gib59","Mar19","Koe88","Koe88b","Dam18b","Gib59b",
                                    "Her19","Cer89","Tom19")  

Be7np$type <- as.factor(Be7np$type)

re <- as.numeric(Be7np$dat) 
Nre <- length(unique(Be7np$dat))
# Unique removes duplicated vector, we want to know how many groups of
# data are there
N <- nrow(Be7np) # Total No of data sets
obsy <- Be7np$S    # Response variable in MeV
obsx <-  Be7np$E   # Predictors
erry <- Be7np$Stat # Error in MeV
set <- Be7np$dat # Get the labels as a vector
fu <- log(c(1.020,1.10,1.050,1.051,1.085,1.032))

samplerCode <- nimbleCode({

  ###################
  # LIKELIHOODS
  ###################
  # - careful: dnorm is differently defined in R and JAGS!
  # - precision=sigma^(-2)
  # - in a for loop, make sure **all** variables on the LEFT of an
  #   expression has the index [i]
  # - systematic error as normalization factor y.norm...

  ## Calling sigmaBe7mod once!
  sigmaBe7modT[1:N] <- sigmaBe7mod(obsx[1:N],
                                 e0_1, ga_1, gb_1, ra, rb,
                                 e0_2, ga_2, gb_2, ra, rb,
                                 e0_3, ga_3, gb_3, ra, rb,
                                 e0_4, ga_4, gb_4, ra, rb,
                                 e0_5, ga_5, gb_5, ra, rb,
                                 e0_6, ga_6, gb_6, ra, rb,
                                 e0_7, ga_7, gb_7, ra, rb)

  for (i in 1:N){
       obsy[i] ~ dnorm(yt[i], pow(sqrt(erry[i]^2+y.scat[re[i]]^2),-2))
        yt[i] <- y.norm[re[i]]*(sqrt(obsx[i])*sigmaBe7modT[i] + hbg )
       }
  
 # for (i in 1:N){
#    obsy[i] ~ dnorm(yt[i], sd = sqrt(erry[i]^2 + y.scat[re[i]]^2))
#     obsy[i] ~ dnorm(yt[i], sd = sqrt(erry[i]^2 + y.scat^2))
#     yt[i] <- y.norm[re[i]]*sqrt(obsx[i])*sigmaBe7modT[i]
#  }
  


  ###################
  # PRIORS
  ###################
  # parameters: ecm, e0, ga, gb, ra, rb, xj, xla, xlb

  # channel radii
  ra ~ T(dnorm(4.0, pow(0.5, -2)),0,Inf)
  rb ~ T(dnorm(4.0, pow(0.5, -2)),0,Inf)

  # Wigner limit: wl = hbar^2/(m_red a_c^2) = 41.80159/(M_red a_c^2)
  #
  # neutron channel: wl_n = 41.80159/(0.88186 a_c^2) = 47.40160/a_c^2
  # proton channel:  wl_p = 41.80159/(0.88079 a_c^2) = 47.45920/a_c^2
  #
  wl_n <- 47.40160*pow(ra, -2)
  wl_p <- 47.45920*pow(rb, -2)


  ###################################################################
  # RESONANCE 1: e0=0 MeV
  # resonance spin, orbital angular momenta   # 2-
  # energy eigenvalue
  e0_1 ~ T(dnorm(0.0, pow(0.1, -2)),0,Inf)

  # reduced widths
  ga_1 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_1 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 2: e0=0.150 MeV
  # resonance spin, orbital angular momenta # 3+


  # energy eigenvalue
  e0_2 ~ T(dnorm(0.15, pow(0.025, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_2 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_2 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)


  ##############
  # RESONANCE 3: e0=0.336 MeV
  # resonance spin, orbital angular momenta  # 3+


  # energy eigenvalue
  e0_3 ~ T(dnorm(0.336, pow(0.010, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_3 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_3 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 4: e0=0.510 MeV
  # resonance spin, orbital angular momenta # 1-

  # energy eigenvalue
  e0_4 ~ T(dnorm(0.51, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_4 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_4 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 5: e0=0.960 MeV
  # resonance spin, orbital angular momenta # 4+

  # energy eigenvalue
  e0_5 ~ T(dnorm(0.96, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_5 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_5 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 6: e0=1.23 MeV
  # resonance spin, orbital angular momenta # 2+


  # energy eigenvalue
  e0_6 ~ T(dnorm(1.23, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_6 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_6 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 7: e0 = 1.32 MeV
  # resonance spin, orbital angular momenta # 0+


  # energy eigenvalue
  e0_7 ~ T(dnorm(1.32, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak
  # reduced widths
  ga_7 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_7 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)
  ##################################################################


 
  hbg ~ dgamma(0.1,0.1)

  for (k in 1:4) {
    y.scat[k]  <- 0
  }
  for (k in 5:10) {
   y.scat[k] ~ dgamma(1,20)
  }

  
  
  for (k in 1:4) {
    # Systematic Uncertainty as a weakly informative prior
    nf[k] ~ dunif(-1,1)
    y.norm[k] <- 10^nf[k]
  }
  
  for (k in 5:10) {
    # Systematic Uncertainty as a highly informative prior
    y.norm[k] ~ dlnorm(log(1.0), pow(fu[k-4], -2))
  }
  


#  k ~ T(dnorm(0,sd=0.15),0,Inf)
  r_1 <- ga_1/gb_1
  r_4 <- ga_4/gb_4
})

samplerData <- list(obsy = obsy,
                    re = re, # Sample size
                    erry = erry,
                    obsx = obsx,
                    fu = fu)


samplerConst <- list(N = N,
                     Nre = Nre)


samplerInits <- list(y.norm = rep(1,Nre),
                     y.scat = rep(0.01,Nre),
                     nf =rep(0,Nre), 
                     y.norm = rep(1,Nre),
                     e0_1 = 0.05, gb_1 = 0.1, ga_1 = 0.1,
                     e0_2 = 0.15, ga_2 = 0.1, gb_2 = 0.1,
                     e0_3 = 0.336, ga_3 = 0.1, gb_3 = 0.1,
                     e0_4 = 0.51, ga_4 = 0.1, gb_4 = 0.1,
                     e0_5 = 0.96, ga_5 = 0.1,  gb_5 = 0.6,
                     e0_6 = 1.23, ga_6 = 0.1, gb_6 = 0.1,
                     e0_7 = 1.32, ga_7 = 0.1, gb_7 = 0.1,
                     ra = 4, rb = 4,
                     yt = rep(1,N),
                     hbg = 0.01
                     )



#nimbleOptions(oldConjugacyChecking = FALSE) 
#nimbleOptions(useNewConfigureMCMC = TRUE)

###############################################################
# (Alternative) If invoking Nimble MCMC stepwise (but more customisable)
##############################################################
ourmodel <- nimbleModel(code = samplerCode, constants = samplerConst,
                        data = samplerData, inits = samplerInits, check = FALSE
)

nimbleOptions(oldConjugacyChecking = FALSE)
nimbleOptions(useNewConfigureMCMC = TRUE)

compileNimble(ourmodel)
# Always compile the model after you are done setting up with it

ourmodel$calculate()
# Calculate the log likelihood (logProb). If the value is not NA,
# we have successfully initialised the model (FULLY)
# One iteration: simulate -> new values -> calculate

conf <- configureMCMC(ourmodel,print=TRUE,thin = 5)
# print = TRUE tells you what kind of samplers are being used for each stochastic node

conf$addMonitors(c('e0_1','ga_1','gb_1','e0_2','ga_2','gb_2',
                   'e0_3', 'ga_3', 'gb_3', 'e0_4', 'ga_4', 'gb_4',
                   'e0_5', 'ga_5', 'gb_5', 'e0_6', 'ga_6', 'gb_6',
                   'e0_7', 'ga_7', 'gb_7',
                   'r_1', 'r_4', 'ra', 'rb',
                   'y.norm',
                   'y.scat','hbg'))


conf$removeSampler(c('ga_1','gb_1','e0_4', 'ga_4', 'gb_4',
                     'e0_5', 'e0_6', 
                     'e0_7', 'ga_7', 
                    'ra', 'rb'))

conf$addSampler(target = c('ga_1','gb_1','e0_4', 'ga_4', 'gb_4',
                           'e0_5', 'e0_6', 
                           'e0_7', 'ga_7', 
                           'ra', 'rb'),
               type = "AF_slice")


# Add the parameters to monitor throughout the MCMC process

samplerMCMC <- buildMCMC(conf)
# Note that if no configuration is required, we can directly jump towards buildMCMC.
# But you know, then theres no point in using this entire chunk of code
# where you can just run nimbleMCMC

compiledMCMC <- compileNimble(samplerMCMC,project = ourmodel,showCompilerOutput = TRUE)
# Complie the configured model again once we are done with configuring it;
# Just before we perform the MCMC runs
# Can include project = ourmodel; to include all the nimbleFunctions, nimbleModels that
# you would want to inherit the functions and models from
# resetFunctions = TRUE; if you would want to reset all the previously created functions
# in order to addd the new MCMC

n.chains = 3
n.iter =   150000
n.burnin = 100000

system.time(
  mcmcChain <- runMCMC(compiledMCMC,niter = n.iter, nchains = n.chains, nburnin = n.burnin,
                       samplesAsCodaMCMC = TRUE)
)

samp <- as.matrix(mcmcChain)


write.csv(samp,"Be7MCMC_case2.csv")



