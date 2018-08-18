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
require(nuclear);library(magrittr);require(gsl)
library(dplyr);require(lessR);library(BayesianTools);
require(truncnorm);require(msm)
source("SfacTdn.R")


######################################################################
## Read DATA
ensamble <- read.csv("ensamble_Tdn_extra.csv",header = T)
#%>%  filter(E <= 0.5)
#filter(dat!= "Mag75")
#%>% filter(E <= 0.5) %>%   filter(dat!= "Arn53") %>%
# droplevels(ensamble$dat)


re <- as.numeric(ensamble$dat)
Nre <- length(unique(ensamble$dat))
Nik <- length(unique(ensamble$invK))

# Radius
# r_i = 6
# r_f = 5


N <- nrow(ensamble)
obsy <- ensamble$S    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat
errx <- ensamble$E_stat
set <- ensamble$dat
lab <- ensamble$invK
syst = c(unique(ensamble$Syst))
systx <- c(0.000075,0.000009,0.0002,0.000006,0.0032)


likelihood <- function(par){
  e0 = par[1]
  er = par[2]
  gin = par[3]
  gout = par[4]
  ad   = par[5]
  ap =  par[6]
  yscat = par[7:11]
  ynorm = par[12:16]
  ue = par[17]
  y = par[18:(N + 17)]

  llynorm = sum(dlnorm(ynorm,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  lly <- sum(dnorm(y,mean = ynorm[re]*SfacTdn(obsx, e0 ,er,gin, gout,ad,ap,ue), sd = yscat[re],  log = T))
  llobs = sum(dnorm(obsy,mean = y,sd = erry,log = T))
  return(llynorm + llobs + lly)
  
}





low <- c(rep(1e-3,2),1e-4,1e-4, 2,2,rep(1e-3,5), rep(0.5,5), 1e-3,obsy - 2*abs(erry))
up <- c(1,0.3, rep(10,2), 10,10, rep(1,5),  rep(1.5,5), 100,obsy + 2*abs(erry))



createTdnPrior <- function(lower, upper, best = NULL){
density = function(par){
  d1 = dnorm(par[1], mean = 0, sd = 1, log = TRUE)
  d2 = sum(dtnorm(par[2], mean = 0, sd = 3,log = TRUE))
  d3 = sum(dnorm(par[3:4], mean = 0, sd = 3, log = TRUE))
  d4 = sum(dunif(par[5:6], 2, 10, log = TRUE))
  d5 = sum(dunif(par[7:11], 0, 1,log = TRUE))
  d6 = sum(dlnorm(par[12:16],log(1),log(1 + syst^2),log = TRUE))
  d7 = dgamma(par[17], 1, 0.05,log = TRUE)
  d8 = sum(dunif(par[18:(N + 17)],obsy - 2*erry,obsy + 2*erry,log = TRUE))
  return(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8)
}

sampler = function(){
  c(runif(1, 0,  1),
    exp(runif(1, log(0.001), log(0.3))),
    exp(runif(2, log(1e-3), log(10))),
    runif(2, 2, 10),
    exp(runif(5,log(1e-4), log(1))),
    rlnorm(5, log(1), log(1 + syst^2)),
    runif(5, 1e-5, 0.01),
    rgamma(1, 1, 0.05),
    runif(N, obsy - erry,obsy + erry)
  )
}

out <- createPrior(density = density, sampler = sampler, lower = lower, upper = upper, best = best)
return(out)
}


#prior <- createPrior(density = density,
#                    lower = low, upper = up, best = NULL)


setup <- createBayesianSetup(likelihood = likelihood, lower = low, upper = up,
names = c("e0","er","gd2","gn2","ad","an",to("yscat", 5),to("ynorm", 5),"ue", to("y", N)))

  settings <- list(iterations = 500000,
                   burnin = 200000, message = T,nrChains = 1,adaptation = 0.3)


  res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")

summary(res)
tracePlot(sampler = res,  start = 30000, whichParameters = c(1,2,3,4,5,6))

correlationPlot(res )


tracePlot(sampler = res, thin = 10, start = 50000, whichParameters = c(1,2,3,4,5,6,7,8,9))





density = function(par){
  d1 = dnorm(par[1], 0.001,1, log =TRUE)
  d2 = dunif(par[2], 0.001,3, log =TRUE)
  d3 = dunif(par[3], 0.001,3, log =TRUE)
  return(d1 + d2 + d3)
}

sampler = function(n=1){
  d1 = rnorm(n, 0.001,1)
  d2 = rnorm(n, 0.001,1)
  d3 = rnorm(n, 0.001,1)
  return(cbind(d1,d2,d3))
}
prior <- createPrior(density = density, sampler = sampler,
                     lower = c(0.001,0.001,0.001), upper = c(10,10,10), best = NULL)


codaObject = getSample(res, start = 500, coda = TRUE)

as.mcmc(codaObject)
getmcmc_var <- function(outjags=outjags,vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
getmcmc_var(codaObject,vars = c("par 1","par 2","par 3","par 4"))

