######################################################################
# test.R
#
# purpose: 
# - regression model for x and y data
#
# conditions:
# - x and y data have no error
# - linear relationship between x and y [could be changed]
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package
library(rjags)
require(MASS)
## for block updating [we do not need to center predictor variables]
load.module("glm") 
#set.seed(531)

######################################################################
## ARTIFICIAL DATA GENERATION
######################################################################
N <- 300

# generate uniformly random x values [without noise]; compute y at all
# x, then add Gaussian noise on y variable, with mean=0 and sd=5

obsx <- runif(N, 20, 100)
epsilon <- rnorm(N, 0, 15)      # last value is the standard deviation
obsy <- 1 * obsx + 0 + epsilon

######################################################################
## Output Observed data (for later ploting)
######################################################################
tempdata <- data.frame(obsx,obsy)
write.matrix(tempdata,"data")



######################################################################
## FUNCTIONS
######################################################################

## no functions


######################################################################
######################################################################
## BAYESIAN MODEL
######################################################################
cat('model {

###################
# LIKELIHOODS
###################
# - careful: dnorm is differently defined in R and JAGS! 
# - precision=sigma^(-2)
# - in a for loop, make sure **all** variables on the LEFT of an 
#   expression has the index [i]

for (i in 1:length(obsx)){
  # subject to statistical uncertainties
  obsy[i] ~ dnorm(y[i], pow(sigmay, -2))
  # true relationship
  y[i] <- alpha + beta * obsx[i]
} 

###################
# PRIORS
###################
alpha ~ dnorm(0.0,pow(100, -2))
beta ~ dnorm(0.0,pow(100, -2))
sigmay ~ dunif(0,200)

}', file={f <- tempfile()})

######################################################################
# n.adapt:  number of iterations in the chain for adaptation  
#           [JAGS will use to choose the sampler and to assure optimum 
#           mixing of the MCMC chain; will be discarded] 
# update(): performs the burn-in on each chain by running the MCMC for 
#           n.burn iterations without saving any of the posterior samples
# coda.samples(): runs each MCMC chain for the number of iterations 
#           specified by n.iter, but it does not save every iteration; 
#           instead, it saves only ever nth iteration, where n is given 
#           by thin
# n.chains: number of mcmc chains

n.chains = 3
n.adapt = 1000 
n.burn = 10000     
n.iter = 20000  
thin = 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model
ourmodel <- jags.model(f,
                data = list('obsx' = obsx, 'obsy' = obsy),
##              inits = list(alpha = 0, beta = 1, sigmax = 10, sigmay = 10),
                n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn)
    
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
               variable.names=c( 
                   'alpha', 
                   'beta', 
                   'sigmay'
                               ), 
               n.iter=n.iter, n.thin=n.thin)
               
# <---- rjags
######################################################################
######################################################################
# OUTPUT RESULTS TO SCREEN
######################################################################
cat("", "\n")    # output empty line

# sample size adjusted for autocorrelation
effectiveChainLength = effectiveSize(mcmcChain) 
show(effectiveChainLength)

cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n") 

######################################################################
# OUTPUT RESULTS TO FILES
######################################################################
# !!! make sure to check the order of the parameters in the MCMC output !!!
samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)

# output all results
write.matrix(samplesmat,"test")













