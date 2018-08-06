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
ensamble <- read.csv("ensamble_Tdn_extra.csv",header = T) %>%  filter(E <= 0.5)
#filter(dat!= "Mag75")
#%>% filter(E <= 0.5) %>%   filter(dat!= "Arn53") %>%
# droplevels(ensamble$dat)


re <- as.numeric(ensamble$dat) 
Nre <- length(unique(ensamble$dat))
Nik <- length(unique(ensamble$invK))

N <- nrow(ensamble)
obsy <- ensamble$S    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat
set <- ensamble$dat
lab <- ensamble$invK
syst = c(unique(ensamble$Syst))


likelihood <- function(par){
  e1 = par[1]
  gin = par[2]
  gout = par[3]
  sigmax = par[4]
  scale = par[5:9]
  y = par[10:(N + 9)]
  
  llRandom = sum(dlnorm(scale,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  lly <- sum(dnorm(y,mean = scale[re]*sfactorTdn_5p(obsx, e1,gin, gout,6,5), sd = sigmax,  log = T))
  llobs = sum(dnorm(obsy,mean = y,sd = erry,log = T))
  return(llRandom + llobs + lly)
  
}



setup <- createBayesianSetup(likelihood = likelihood,
lower = c(0.001,0.001,0.001,0.001,rep(0.5,5),obsy - 5*abs(erry)),
upper = c(1,2,2,5,rep(1.5,5),obsy + 5*abs(erry)),
names = c("e0","gd2","gn2","sigma",to("scale", 5),to("y", N)))
settings <- list(iterations = 200000,
                 burnin = 50000, message = T)

res <- runMCMC(bayesianSetup = setup, settings = settings)
summary(res)
tracePlot(sampler = res, thin = 10, start = 5000, whichParameters = c(1,2,3,4,5,6,7,8,9))







codaObject = getSample(res, start = 500, coda = TRUE)

as.mcmc(codaObject)
getmcmc_var <- function(outjags=outjags,vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
getmcmc_var(codaObject,vars = c("par 1","par 2","par 3","par 4"))
