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
library(dplyr);require(lessR);library(BayesianTools)
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
set <- ensamble$dat
lab <- ensamble$invK
syst = c(unique(ensamble$Syst))




likelihood <- function(par){
  e0 = par[1]
  er = par[2]
  gin = par[3]
  gout = par[4]
  ad   = par[5]
  ap =  par[6]
  sigmax = par[7]
  scale = par[8:12]
  y = par[13:(N + 12)]
#  y = par[10:133]


  llRandom = sum(dlnorm(scale,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  lly <- sum(dnorm(y,mean = scale[re]*SfacTdn(obsx, e0 ,er,gin, gout,ad,ap,ue=0), sd = sigmax,  log = T))
  llobs = sum(dnorm(obsy,mean = y,sd = erry,log = T))
  return(llRandom + llobs + lly)


}



low <- c(1e-4,1e-3,1e-4,1e-4,2.5,2.5,1e-2,rep(0.8,5),obsy - 2*abs(erry))
up <- c(1,0.3,50,50,10,10,5,rep(1.2,5),obsy + 2*abs(erry))


density = function(par){
  d1 = dnorm(par[1], mean = 0, sd = 1, log = TRUE)
  d2 = dunif(par[2], 0.001,0.3, log = TRUE)
  d3 = sum(dnorm(par[3:4], mean = 0, sd = 3, log = TRUE))
  d4 = sum(dunif(par[5:6], 2.5, 10, log = TRUE))
  d5 = dnorm(par[7], mean = 0, sd = 1, log = TRUE)
  d6 = sum(dlnorm(par[8:12],1,log(1 + syst^2),log = TRUE))
  d7 = sum(dnorm(par[13:(N + 12)],mean=obsy,sd=erry))
  return(d1 + d2 + d3 + d4 + d5 + d6 + d7)
}


sampler = function(){
  d1 = runif(1, 0,  1)
  d2 = runif(1, 0.001, 0.3)
  d3 = runif(2, 0, 10)
  d4 = runif(2, 2.5, 10)
  d5 = runif(1, 0, 1)
  d6 = runif(5, 0.8, 1.2)
  d7 = rnorm(N,mean=obsy,sd=erry)
  return(cbind(d1,d2,d3,d4,d5,d6,d7))
}


#prior <- createPrior(density = density,
#                    lower = low, upper = up, best = NULL)


setup <- createBayesianSetup(likelihood = likelihood, lower = low, upper = up,
names = c("e0","er","gd2","gn2","ad","an","sigma",to("scale", 5),to("y", N)))

  settings <- list(iterations = 800000,
                   burnin = 75000, message = T,nrChains = 1,adaptation = 0.5)


  res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")

summary(res)
tracePlot(sampler = res,  start = 1000, whichParameters = c(1,2,3,4,5,6))

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

