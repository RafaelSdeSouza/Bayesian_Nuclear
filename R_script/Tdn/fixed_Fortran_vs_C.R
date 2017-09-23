
# He3dp analysis
#
# purpose: ARTIFICIAL DATA
#
# - 3 parameters are assumed: Er, gamma_d^2, gamma_n^2 [e1, gin, gout]
#
# - uses the function sfactorTdn_fast(obsx1[i], e1, gin, gout), which
#   is a C++ version of a Fortran code that includes Coulomb wave
#   function calculations; JAGS has been recompiled with this C++ function 
#
# Fortran code tdn_plot.f is needed in this R script for plotting the
# S-factor only
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
#set.seed(123)
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical 
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb; 
######################################################################
## ARTIFICIAL DATA GENERATION 

N <- 5000

obsx1 <- exp(runif(N,log(0.001),log(1)))

res <- vector()
obsy1 <- vector()
errobsy1 <- vector()

# Barker values:
# Er  = 0.0912 MeV
# g^2_in = 2.93 MeV         ! reduced width of deuteron
# g^2_out = 0.0794 MeV      ! reduced width of neutron

#res[2] <-  0.35779  # resonance energy
res[2] <- 0.0912   # resonance energy
res[3] <-  2.93    # reduced width incoming
res[4] <- 0.0794   # reduced width outgoing

for (i in 1:length(obsx1)){
  res[1] <- obsx1[i]
  write.table(res, file="tdn_AD.in", quote=TRUE, 
              row.names=FALSE, col.names=FALSE)
  
  # Load the fortran code needed to calculate S-factor curve
  if(!is.loaded("tdn_AD_Sub"))
    dyn.load("tdn_AD.so") 
  .Fortran("tdn_AD_Sub")
  
  tab1 <- read.table("tdn_AD.out", header=FALSE)
  
  errobsy1[i] <- 0.001
  obsy1[i] <- rnorm( 1, tab1[1,2], errobsy1[i] )
}

######################################################################                 
# import jags package 
library(rjags)
library(runjags)
library(R2jags)
library(mcmcplots)
## for block updating [we do not need to center predictor variables]
load.module("glm")  
require(RcppGSL)
load.module("nuclear")  

# 
######################################################################                 
cat('model {
    
    # LIKELIHOOD
    for (i in 1:length(obsx1)) {
    obsy1[i] ~ dnorm(obsx1[i],tau)

    yt[i] <- sfactorBW(obsx1[i], 0.0912, 2.93, 0.0794)
    yk[i] <- sfactorTdn(obsx1[i], 0.0912 , 2.93, 0.0794)
    }    
   tau ~ dgamma(0.1,0.1)
   e1 ~ dgamma(0.1,0.1)
    }', file={f <- tempfile()})
######################################################################
inits <- function () { list(e1 = runif(1,0,0.3),gin=runif(1,2,10),gout=runif(1,0.01,1)) }
# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 



# JAGS model with R2Jags;
out <- jags(data = list('obsx1' = obsx1,'obsy1' = obsy1),
            parameters = c("e1","yt","yk"),
            model.file = f,
            n.thin = 10,
            n.chains = 3,
            n.burnin = 50,
            n.iter = 100,jags.module = c("nuclear"))

plot(obsx1,obsy1,col="red",xlim=c(0.001,0.3),ylim=c(0,30),cex=0.4,log="x")
par(new=TRUE)
plot(obsx1,out$BUGSoutput$mean$yk,col="green",xlim=c(0.001,0.3),ylim=c(0,30),cex=0.4,log="x")

hist((out$BUGSoutput$mean$yk-obsy1)/obsy1)

par(new=TRUE)
plot(obsx1,out$BUGSoutput$mean$yt,col="blue",xlim=c(0.001,0.3),ylim=c(0,25),cex=1.3,log="x")




