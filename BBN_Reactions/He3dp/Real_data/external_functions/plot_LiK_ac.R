# Plot likelihood chanel radii

# preparation: remove all variables from the work space
rm(list=ls())
set.seed(27)
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
require(msm);require(LaplacesDemon);require(mcmcplots);require(ggmcmc);
require(ggridges)
source("..//..//auxiliar_functions/pair_wise_plot.R")
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
syst = 1 + c(0.03,unique(ensamble$Syst))

obsx1 <- exp(runif(N,log(1e-3),log(1)))
errobsy1 <- runif(N,0.1,0.5)
obsy1 <- rnorm(N, sfactor3Hedp_5p(obsx1 ,0.35779,1.0085,0.025425,5,6),errobsy1)

likelihood <- function(par){
  e0 = par[1]
  er = par[2]
  gd2 = par[3]
  gp2 = par[4]
  ad   = par[5]
  ap =  par[6]
  lly <- sum(dnorm(y,mean = sfactor3Hedp_5p(obsx, e0,er,gd2, gp2,ad,ap,ue = ue[ik]),sd = 0.4, log = T))
  return(lly)
  
}


best <- c(0.1,1e-3,rep(1e-4,2), 1,1,1e-4,rep(0.5,7),rep(0,2),obsy - 2*erry,1)



up <- c(0.4,1,rep(3,2),10,10,5,rep(1.5,7),rep(350,2),obsy + 2*erry,100)
