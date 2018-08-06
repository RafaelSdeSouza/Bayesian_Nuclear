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
  er = par[2]
  gd2 = par[3]
  gp2 = par[4]
  ad   = par[5]
  ap =  par[6]
  sigmax = par[7]
  scale = par[8:14]
  ue = par[15:16]
  y = par[17:(N + 16)]
  
  llRandom = sum(dlnorm(scale,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  lly <- sum(dnorm(y,mean = scale[re]*sfactor3Hedp_5p(obsx, e0,er,gd2, gp2,ad,ap,ue = ue[ik]), sd = sigmax,  log = T))
  llobs = sum(dnorm(obsy,mean = y,sd = erry,log = T))
  return(llRandom + llobs + lly)
  
}


low <- c(1e-3,1e-3,1e-5,1e-5,1,1,1e-2,rep(0.8,7),rep(0,2),obsy - 2*abs(erry))
up <- c(1,2,10,10,15,15,5,rep(1.2,7),rep(300,2),obsy + 2*abs(erry))

density = function(par){
  d1 = dnorm(par[1], mean = 0, sd = 1, log = TRUE)
  d2 = dnorm(par[2], mean = 0, sd = 1, log = TRUE)
  d3 = sum(dnorm(par[3:4], mean = 0, sd = 3, log = TRUE))
  d4 = sum(dunif(par[5:6], 1, 15, log = TRUE))
  d5 = dnorm(par[7], mean = 0, sd = 1, log = TRUE)
  d6 = sum(dlnorm(par[5:11],1,log(1 + syst^2),log = TRUE))
  d7 = sum(dnorm(par[12:13], mean = 0, sd = 100, log = TRUE))
  d8 = sum(dnorm(par[14:(N + 13)],mean=obsy,sd=erry))
  return(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d7)
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
                 burnin = 100000, message = T)


system.time(
res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")
)
tracePlot(sampler = res, thin = 10, start = 50000, whichParameters = c(1,2,3,4,5,6,15,16))





summary(res)




codaObject = getSample(res, start = 500, coda = TRUE)

as.mcmc(codaObject)
getmcmc_var <- function(outjags=outjags,vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
getmcmc_var(codaObject,vars = c("par 1","par 2","par 3","par 4"))
