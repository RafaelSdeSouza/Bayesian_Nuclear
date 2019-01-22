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
library(dplyr);require(lessR);library(BayesianTools)


######################################################################
## Read DATA
ensamble <- read.csv("ensamble.csv",header = T) %>%
  mutate(Syst=replace(Syst,Syst==0.06,0.078))  %>% filter(E <= 0.5)


re <- as.numeric(ensamble$dat)
Nre <- length(unique(ensamble$dat))
ik <- as.numeric(ensamble$invK)
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
syst = c(0.03,unique(ensamble$Syst))

likelihood <- function(par){
  e0 = par[1]
  gd2 = par[2]
  gp2 = par[3]
  sigmax = par[4]
  scale = par[5:11]
  ue = par[12:13]
  y = par[14:(N + 13)]
  
  llRandom = sum(dlnorm(scale,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  lly <- sum(dnorm(y,mean = scale[re]*sfactor3Hedp_5p(obsx, e0,e0,gd2, gp2,6,5,ue = ue[ik]), sd = sigmax,  log = T))
  llobs = sum(dnorm(obsy,mean = y,sd = erry,log = T))
  return(llRandom + llobs + lly)
  
}


low <- c(1e-3,1e-5,1e-5,1e-2,rep(0.8,7),rep(0,2),obsy - 2*abs(erry))
up <- c(1,3,1,5,rep(1.2,7),rep(300,2),obsy + 2*abs(erry))

density = function(par){
  d1 = dnorm(par[1], mean = 0, sd = 1, log = TRUE)
  d2 = sum(dnorm(par[2:3], mean = 0, sd = 3, log = TRUE))
  d3 = dnorm(par[4], mean = 0, sd = 1, log = TRUE)
  d4 = sum(dlnorm(par[5:11],1,log(1 + syst^2),log = TRUE))
  d5 = sum(dnorm(par[12:13], mean = 0, sd = 100, log = TRUE))
  d6 = sum(dnorm(par[14:(N + 13)],mean=obsy,sd=erry))
  return(d1 + d2 + d3 + d4 + d5 + d6)
}


prior <- createPrior(density = density, 
                     lower = low, upper = up, best = NULL)


#prior <- createUniformPrior(lower = low,
#                            upper = up)

#setup <- createBayesianSetup(likelihood = likelihood,prior = prior,
#        names = c("e0","gd2","gp2","sigma",to("scale", 7),to("ue", 2),to("y", N)))


#createBayesianSetup(likelihood = likelihood,prior = prior,
#                          names = c("e0","gd2","gp2","sigma",to("scale", 7),to("ue", 2),to("y", N)))
                    
                    
setup <- createBayesianSetup(likelihood = likelihood,lower = low,upper = up)

settings <- list(iterations = 500000,
                 burnin = 50000, message = T)


system.time(
res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")
)
tracePlot(sampler = res, thin = 10, start = 25000, whichParameters = c(1,2,3,12,13))





summary(res)




codaObject = getSample(res, start = 500, coda = TRUE)

as.mcmc(codaObject)
getmcmc_var <- function(outjags=outjags,vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
getmcmc_var(codaObject,vars = c("par 1","par 2","par 3","par 4"))
