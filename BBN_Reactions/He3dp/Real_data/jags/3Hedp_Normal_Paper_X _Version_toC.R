# 3Hedp analysis
#
# purpose: Real  DATA
#
# - 5 parameters are assumed: Er, gamma_d^2, gamma_n^2 [e1, gin, gout]
#
# - uses the function sfactorHe3dp(obsx1[i], e1, gin, gout), which
#   is a C++ version of a Fortran code that includes Coulomb wave
#   function calculations; JAGS has been recompiled with this C++ function
#
######################################################################
# preparation: remove all variables from the work space
#rm(list=ls())
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb;

######################################################################
# import packages
library(rjags);library(R2jags);library(mcmcplots)
require(RcppGSL);require(ggplot2);require(ggthemes)
require(nuclear);library(magrittr);library(wesanderson)
library(dplyr);require(ggsci);require(ggmcmc);require(plyr);library(latex2exp)
require(MCMCvis);require(ggridges);require(cowplot)
source("jagsresults.R")
source("theme_rafa.R")

source("plot_Sfactor.R")
source("plot_normfactors.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")
load.module("nuclear")


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
#syst <- syst[-3]

M <- 500
xx <- seq(min(obsx),max(obsx),length.out = M)

model.data <- list(obsy = obsy,    # Response variable
                   obsy2 = obsy,    # Response variable
                   obsx =  obsx,   # Predictors
                   erry = erry,
                   N = nrow(ensamble), # Sample size
                   syst = syst,
                   Nre = Nre,
                   re = re,
                   Nik = Nik,
                   ik  = ik,
                   M = M,
                   xx = xx
 #                  ap  = 5,
 #                  ad = 6

)


# Conservative case
######################################################################
Model <- "model{
# LIKELIHOOD informative
for (i in 1:N) {
obsy[i] ~ dnorm(y[i], pow(erry[i], -2))
y[i] ~ dt(scale[re[i]]*sfactor3Hedpx(obsx[i], E0,  gd2, gp2, ad, ap, ue[ik[i]]), pow(tau, -2),nu)

}


# Residuals
for (i in 1:N) {
yy[i] ~ dnorm(scale[re[i]]*sfactor3Hedpx(obsx[i], E0,  gd2, gp2, ad, ap, ue[ik[i]]), pow(tau, -2))
res[i] <- obsy[i]-yy[i] 
#res[i] <- obsy[i]-sfactor3Hedpx(obsx[i], E0,  gd2, gp2, ad, ap,0)
nres[i] <- res[i]/obsy[i]

}


# LIKELIHOOD broad
for (i in 1:N) {
obsy2[i] ~ dnorm(y_2[i], pow(erry[i], -2))
y_2[i] ~ dnorm(scale[re[i]]*sfactor3Hedpx(obsx[i],  E0_b,  gd2_b, gp2_b, ad_b, ap_b, ue[ik[i]]),pow(tau_2, -2))
}
RSS <- sum(res^2)

# Predicted values
for (j in 1:M){

# Bare...

mux0[j] <- sfactor3Hedpx(xx[j], E0,  gd2, gp2, ad, ap,0)


mux0_2[j] <- sfactor3Hedpx(xx[j], E0_b,  gd2_b, gp2_b, ad_b, ap_b,0)

DeltaM[j] <- mux0[j]/mux0_2[j]



# No inverse Kinematics

mux1[j] <- sfactor3Hedpx(xx[j], E0,  gd2, gp2, ad, ap,ue[1])
yx1[j] ~ dt(mux1[j],pow(tau,-2),nu)

# With inverse Kinematics
mux2[j] <- sfactor3Hedpx(xx[j], E0,  gd2, gp2, ad, ap,ue[2])
yx2[j] ~ dt(mux1[j],pow(tau,-2),nu)

}

for (k in 1:Nre){
scale[k] ~ dlnorm(log(1.0),pow(log(syst[k]),-2))
}

for (z in 1:Nik){
ue[z] ~ dnorm(0,pow(0.1,-2))T(0,)
}


# PRIORS

# Case I
tau ~  dnorm(0, pow(1,-2))T(0,)
E0  ~  dnorm(0, pow(1,-2))T(0,)
gd2 ~  dnorm(0, pow(1,-2))T(0,)
gp2 ~  dnorm(0, pow(2,-2))T(0,)
nu <- nuMinusOne + 1
nuMinusOne ~ dexp( 1/29 )
ad  <- 6
ap  <- 5



# Case II
tau_2  ~  dnorm(0, pow(1,-2))T(0,)
E0_b  ~   dnorm(0.2,pow(0.02,-2))T(0,)

gd2_b  ~ dnorm(0.3, pow(0.05,-2))T(0,)
gp2_b ~  dnorm(0.04, pow(0.01,-2))T(0,)


ad_b  ~  dnorm(3.5,pow(0.5,-2))T(0,)
ap_b  ~  dnorm(5.5,pow(1,-2))T(0,)


ue_ev[1] <-1e6*ue[1]
ue_ev[2] <-1e6*ue[2]

S_0   <- sfactor3Hedpx(1e-4, E0,  gd2, gp2, ad, ap,0)
S_0b  <- sfactor3Hedpx(1e-4, E0_b,  gd2_b, gp2_b, ad_b, ap_b,0)


}"

inits <- function(){ list(E0 = runif(1,0.3,0.35),E0_b = 0.2,gd2 = 1,
                        gp2 = runif(1,0.01,0.06),gd2_b = 1) }


set.seed(24)
# JAGS model with R2Jags;
Normfit <- jags(data = model.data,
                inits = inits,
                parameters.to.save  = c("E0","gd2", "gp2","ue_ev","tau", "ad","ap",
                                        "RSS","mux0","mux1","mux2","scale","DeltaM","S_0",
                                        "S_0b","E0_b","gd2_b",
                                        "gp2_b","tau_2","ad_b","ap_b","res","nres","nu"),
                model.file  = textConnection(Model),
                n.thin = 5,
                n.chains = 3,
                n.burnin = 1000,
                n.iter = 3500)
jagsresults(x = Normfit, params = c("E0","gd2", "gp2","ue","tau", "ad","ap","ue_ev","S_0"),probs = c(0.16, 0.5, 0.84))



temp <- Normfit
temp <- update(temp, n.thin = 5, n.iter = 1000)


plot_Sfactor(Normfit)




