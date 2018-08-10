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
require(truncnorm)
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
  gin = par[2]
  gout = par[3]
  yscat = par[4:8]
  ynorm = par[9:13]
  xscat = par[14:18]
  xnorm = par[19:23]
  ue = par[24]
  y = par[25:(N + 24)]
  xx = par[(N + 25):(2*N + 24)]
  xtrue = par[(2*N + 25):(3*N + 24)]

  llxnorm = sum(dnorm(xnorm,mean=0,sd=systx,log=T))
  llynorm = sum(dlnorm(ynorm,meanlog = log(1), sdlog = log(1 + syst^2), log = T))
  
  llxt <- sum(dnorm(xx,mean = xtrue,sd=xscat[re],log=T))
  llxx <- sum(dnorm(obsx,mean = xx + xnorm[re],sd=errx,log=T))
  lly <- sum(dnorm(y,mean = ynorm[re]*SfacTdn(xtrue, e0,e0,gin, gout,6,5,ue), sd = yscat[re],  log = T))
  llobs = sum(dnorm(obsy,mean = y,sd = erry,log = T))


  return(llynorm + llobs + lly + llxx + llxt + llxnorm)
}

low <- c(1e-3,1e-3,1e-3,rep(1e-4,5),rep(0.8,5),rep(1e-5,5), rep(-0.8,5),0,obsy - 2*abs(erry),obsx - abs(errx),obsx - abs(errx))
up <- c(1,20,20,rep(1,5),rep(1.2,5),rep(0.03,5),rep(0.2,5),100,obsy + 2*abs(erry),obsx + abs(errx),obsx + abs(errx))

#best <- c(0.1,1,1,rep(1,5),rep(1,5),rep(1,5),rep(0.1,5),20,(obsy + 2*abs(erry))/2,(obsx + abs(errx))/2,(obsx + abs(errx))/2)


density = function(par){
  d1 = dnorm(par[1], mean = 0, sd = 1, log = TRUE)
  d2 = sum(dtnorm(par[2:3], mean = 0, sd = 2,log = TRUE))
  d3 = sum(dtnorm(par[4:8], 0, 1,log = TRUE))
  d4 = sum(dlnorm(par[9:13],1,log(1 + syst^2),log = TRUE))
  d5 = sum(dtnorm(par[14:18],0, 0.01,log = TRUE))
  d6 = sum(dnorm(par[19:23],0,sd=systx,log = TRUE))
  d7 = dtnorm(par[24], 0, 1000,log = TRUE)
  d8 = sum(dnorm(par[25:(N + 24)],mean=obsy,sd=erry,log = TRUE))
  d9 = sum(dunif(par[(N + 25):(2*N + 24)],0.001,0.3,log = TRUE))
  d10 = sum(dunif(par[(2*N + 25):(3*N + 24)],0.001,0.3,log = TRUE))
  return(d1 + d2 + d3 + d4 + d5 + d6 + d7 +d8 + d9 + d10)
}
sampler = function(){
  d <- list()
  d[1] = runif(1, 0,  1)
  d[2:3] = rtnorm(2, 0,  3,low=0)
  d[4:8] = rtnorm(5, 0, 1,low=0)
  d[9:13] = rtnorm(5, 1, log(1 + syst^2))
  d[14:18] = rtnorm(5, 0, 0.01,low=0)
  d[19:23] = rnorm(5, 0, systx)
  d[24] = rtnorm(1, 0, 100,low=0)
  d[25:(N + 24)] = rnorm(N, obsy, erry)
  d[(N + 25):(2*N + 24)] = runif(N,0.001,0.3)
  d[(2*N + 25):(3*N + 24)] = runif(N,0.001,0.3)
  return(as.numeric(d))
}


prior <- createPrior(density = density, sampler = sampler,
                    lower = low, upper = up, best = best)



setup <- createBayesianSetup(likelihood = likelihood,  priorSampler = sampler,lower = low, upper = up, best = best,
names = c("e0","gd2","gn2",to("yscat", 5),to("ynorm", 5),to("xscat", 5),
          to("xnorm", 5),"ue", to("y", N),to("x", N),to("xt", N)))

settings <- list(iterations = 200000,adaptation = 0.5,
                 burnin = 60000, message = T,nrChains = 1)

res <- runMCMC(bayesianSetup = setup, settings = settings,sampler = "DREAMzs")

tracePlot(sampler = res, start = 20000, whichParameters = c(1,2,3,24))


mcmc_out <- function(out = out, vars = vars){
  as.data.frame(do.call(rbind, out$chain[,vars]))
}


mo <- mcmc_out(res,vars=c("e0","gd2","gn2",to("yscat", 5),to("ynorm", 5),to("xscat", 5),
                          to("xnorm", 5),"ue"))[100000:900003,]

index <- sample(1:nrow(mo),50000,replace=F)
short_mo <- mo[index ,]

write.matrix(short_mo,"BT_Tdn_caseI.dat")









bayesianSetup <- createBayesianSetup(likelihood = ll, prior = newPrior)

settings = list(iterations = 1000,  message = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)


out2 <- DREAMzs(res,settings = settings)






plot(mo$gd2,mo$gn2,type="l")


correlationPlot(res )


tracePlot(sampler = res, thin = 10, start = 5000, whichParameters = c(1,2,3,4,5,6,7,8,9))





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

