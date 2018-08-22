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
require(msm)

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


low <- c(0,0,rep(1e-4,2), 1,1,1e-4,rep(0.5,7),rep(0,2),obsy - 2*erry)
up <- c(1,1,rep(10,2),10,10,5,rep(1.5,7),rep(300,2),obsy + 2*erry)


createHedPrior <- function(lower, upper, best = NULL){
density = function(par){
  d1 = dtnorm(par[2], mean = 0, sd = 1,log = TRUE)
  d2 = dtnorm(par[2], mean = 0, sd = 1,log = TRUE)
  d3 = sum(dnorm(par[3:4], mean = 0, sd = 3, log = TRUE))
  d4 = sum(dunif(par[5:6], 2, 10, log = TRUE))
  d5 = dtnorm(par[7], mean = 0, sd = 5, log = TRUE)
  d6 = sum(dlnorm(par[8:14],log(1),log(1 + syst^2),log = TRUE))
  d7 = sum(dtnorm(par[15:16], mean = 0, sd = 100, log = TRUE))
  d8 = sum(dunif(par[17:(N + 16)],obsy - 2*erry,obsy + 2*erry,log = TRUE))
  return(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8)
}

sampler = function(){
   c(runif(1, 0, 1),
   runif(1, 0, 1),
    exp(runif(2, log(1e-3), log(10))),
    runif(2, 2, 10),
    runif(1, 0, 5),
    rlnorm(7, log(1), log(1 + syst^2)), #ynorm
    runif(2, 0, 300),
    runif(N, obsy - 2*erry,obsy + 2*erry)
)
}

out <- createPrior(density = density, sampler = sampler, lower = lower, upper = upper, best = best)
return(out)
}

prior <- createHedPrior(lower = low, upper = up)


#prior <- createUniformPrior(lower = low,
#                            upper = up)



setup <- createBayesianSetup(likelihood = likelihood,prior = prior,
names = c("e0","er","gd2","gp2","ad","ap","sigma",to("scale", 7),to("ue", 2),to("y", N)))

#setup <- createBayesianSetup(likelihood = likelihood,lower = low,upper = up,
#names = c("e0","er","gd2","gp2","ad","ap","sigma",to("scale", 7),to("ue", 2),to("y", N)))




settings <- list(iterations = 1e6,adaptation = 0.5,
                 burnin = 5e5, message = T,nrChains = 2)


system.time(
res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")
)


tracePlot(sampler = res, thin = 1, start = 5E4, whichParameters = c(1,2,3,4,5,6,15,16))

summary(res)

codaObject = getSample(res, start = 500, coda = TRUE)

as.mcmc(codaObject)
getmcmc_var <- function(outjags=outjags, vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
getmcmc_var(codaObject,vars = c("par 1","par 2","par 3","par 4"))
