#tdn analysis
#
# obsy1[i] ~ dlnorm(yl1[i],pow(corrl.err1[i], -2))                
# yl1[i] <- log(y1[i])-0.5*log(1+pow(corr.err1[i],2)/pow(y1[i],2)) ## lognormal mu
# corrl.err1[i] <- sqrt(log(1+pow(corr.err1[i],2)/pow(y1[i],2)))   ## lognormal sigma
#
#
# conditions:
# - relationship between x and y given by theory [table] 
# - theory relation from calculations with sfactorTdn function and tdn.f
#
# purpose: 
# - fit a theoretical model to data
#
# details:
# - x data have no error; only y data have errors
# - different data sets are statistically independent
# - includes systematic errors in all sets [n.norm]
# - includes robust regression [outliers]
# - true relationship between variables is given by nuclear theory,
#   assuming only a single [scaling] parameter [a.scale]
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
#set.seed(123)
library(emdbook)
library(pracma)
library(R2jags)
library(coda)
library(rjags)
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i]

##We will work with energy is in MeV and S-factor in MeVb 
##(need to convert keVb to MeVb of data from the tables)

########################
########################

## DATA SET 1: AR52
tab1 <- read.table("tdn_data_AR52.dat", header=FALSE)

obsx1 <- c(tab1[,1])

obsy1 <- c(tab1[,2])

errobsy1 <- c(tab1[,3])

######################################################################
# convert latter to MeV b ???


######################################################################
### Model
# import jags package 
library('rjags')
## for block updating [we do not need to center predictor variables]
load.module("glm")  
##
cat('model {
  # 1. LIKELIHOOD
    for (i in 1:length(obsx1)) {
      obsy1[i] ~ dlnorm(yl1[i],pow(corrl.err1[i], -2))                 
      yl1[i] <- log(y1[i])-0.5*log(1+(pow(corr.err1[i],2)/pow(y1[i],2))) ## lognormal mu
      corrl.err1[i] <- sqrt(log(1+(pow(corr.err1[i],2)/pow(y1[i],2))))   ## lognormal sigma
      #p.alt1[i] ~ dcat(p1[])
      corr.err1[i] <- errobsy1[i] #*phi[p.alt1[i]]
      y1[i] <- n.norm1 * z1[i]    # systematic error as normalization factor
    }
  # 2. REALTIONSHIP BETWEEN TRUE VARIABLES
      #a 1000x2 table of Elab(not used), Ecm, cross-section(not used) and S(E)
      table1 = sfactorTdn(e1, gi, gf)
      for (j in 1:length(obsx1)) {
        z1[j] = interp.lin(obsx1[j], table1[1:1000,1], table1[1:1000,2])
      }
    
  # 3. PRIORS
  #e1, gi, gf are defined as in tdn.f (by Coc): resonance energy, initial width and final width
  #We assume normal distributions centered at Coc´s values, truncated at 0
      e1 ~ dnorm(0.0912, pow(0.5, -2))T(0,)
      gi ~ dnorm(2.93, pow(1, -2))T(0,)
      gf ~ dnorm(0.0794, pow(1, -2))T(0,)
      n.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
          logmu1 <- log(1.0)      # median of factor uncertainty is 1.0
          logsigma1 <- log(1.02)  # factor uncertainty is 1.02, i.e. 2% ****HAS TO BE CHANGED****
}', file={f <- tempfile()})





######################################################################
# n.adapt:  number of iterations in the chain for adaptation (n.adapt) 
#           [JAGS will use to choose the sampler and to assure optimum 
#           mixing of the MCMC chain; will be discarded] 
# n.udpate: number of iterations for burnin; these will be discarded to 
#           allow the chain to converge before iterations are stored
# n.iter:   number of iterations to store in the final chain as samples 
#           from the posterior distribution 
# n.chains: number of mcmc chains
# n.thin:   store every n.thin element [=1 keeps all samples]

n.adapt  <- 500   
n.update <- 2000  
n.iter   <- 7500  
n.chains <- 3
n.thin   <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model;
#
ourmodel <- jags.model(f,
                       data = list('obsx1' = obsx1, ## jags wants all data in a list
                                   'obsy1' = obsy1,
                                   'errobsy1' = errobsy1),
                       n.chains = n.chains,
                       n.adapt = n.adapt)
# burnin 
update(ourmodel, n.update, progress.bar="none") 

# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                          variable.names=c('e1','gi','gf',
                                           'n.norm1'), 
                          n.iter=n.iter, n.thin=n.thin)
######################################################################
# output results on screen
cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n") 





