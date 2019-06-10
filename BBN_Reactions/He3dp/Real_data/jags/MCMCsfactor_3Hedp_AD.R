######################################################################
# Author: Christian Iliadis (09/26/2017)
#
# MCMCsfactor3Hedp_AD.R
#
# purpose: testing with ARTIFICIAL DATA
#
# five parameters are assumed: E0, gamma_i^2, gamma_f^2, R_i, R_f, U_e 
# [e0, ga, gb, ra, rout, ue]
#
# uses the C++ function sfactor3Hedpx(ecm, E0, ga, gb, ra, rb, ue), 
# which calls gsl libraries to compute Coulomb wave functions;
#
# the function "Sfactor" in this script is used only for plotting
# believable S-factors and for generating the artificial data.
# 
# MAKE SURE YOU ARE USING THE SAME VALUES FOR MASSES, ENERGIES, ETC.,
# IN BOTH THIS SCRIPT (FUNCTION Sfactor) AND THE C++ SFACTOR FUNCTION
# IN JAGS (sfactor3Hedp.cc)
#
# uses:
# - normal likelihoods
# - systematic uncertainties
# - extrinsic scatter
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package 
library("rjags")
## for block updating [we do not need to center predictor variables]
load.module("glm") 

require(gsl) 
require(RcppGSL)

## load external R-matrix functions
load.module("nuclear")   
# random number seed
# set.seed(3000)

######################################################################
## FUNCTIONS
######################################################################
# error bars
# w is the width of perpendicular end bars on errors; set to zero
add.error.bars <- function(X,Y,dX,dY,w,col){
  X0 = X 
  Y0 = (Y-dY) 
  X1 =X 
  Y1 = (Y+dY)
  arrows(X0, Y0, X1, Y1, lwd=2, code=3, angle=90, length=w, col=col)
  Y0 = Y 
  X0 = (X-dX) 
  Y1 =Y 
  X1 = (X+dX)
  arrows(X0, Y0, X1, Y1, lwd=2, code=3, angle=90, length=w, col=col)
}
######################################################################
# compute S-factor using gsl Coulomb wave functions;
# for plotting and artificial data generation purpose only
#
# ECM   :   center-of-mass energy
# E0    :   eigenenergy
# ga,gb :   reduced widths for incoming/outgoing channels
# ra,rb :   channel radii (fm) for incoming/outgoing channels
# ue    :   laboratory electron screening potential

Sfactor <- function(ecm, e0, ga, gb, ra, rb, ue){
  # input masses, charges, angular momenta
  m1_i = 3.014932162 
  m2_i = 2.013553198    # masses (amu) of 3He and d
  m1_f = 4.001506094 
  m2_f = 1.007276452    # masses (amu) of p and 4He
  z1_i = 2 
  z2_i = 1              # charges of t and d
  z1_f = 2 
  z2_f = 1          # charges of n and 4He
  jt = 0.5              # spins of target, projectile, resonance
  jp = 1.0 
  jr = 1.5               
  Q = 18.353053      # reaction Q-value (MeV)
  la = 0 
  lb = 2            # orbital angular momenta of d and p

  # reduced masses
  mue_i <- (m1_i*m2_i)/(m1_i+m2_i)
  mue_f <- (m1_f*m2_f)/(m1_f+m2_f)
  
  # constants
  pek <- 6.56618216e-1/mue_i
  omega <- (2*jr+1)/((2*jt+1)*(2*jp+1))

  ### CALCULATE S-FACTOR
  ## incoming channel 
  etpe_i=exp(0.98951013*z1_i*z2_i*sqrt(mue_i/ecm))
  eta_a=0.1574854*z2_i*z1_i*sqrt(mue_i)
  rho_a=0.218735*ra*sqrt(mue_i)    
  eta_i=eta_a/(sqrt(ecm))
  rho_i=rho_a*(sqrt(ecm))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k=0)
  # penetration and shift factor
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  s_i <- rho_i*(P3$val_F*P3$val_Fp + P3$val_G*P3$val_Gp)/(P3$val_F^2 + P3$val_G^2)
  # shift factor at E0 [eigenvalue]
  xeta_i=eta_a/(sqrt(e0))
  xrho_i=rho_a*(sqrt(e0))
  PX1 <- coulomb_wave_FG(xeta_i, xrho_i, la, k=0)
  b_i <- xrho_i*(PX1$val_F*PX1$val_Fp + PX1$val_G*PX1$val_Gp)/(PX1$val_F^2 + PX1$val_G^2)
  # partial width
  Ga <- 2*ga*p_i

  ## outgoing channel
  eta_b=0.1574854*z2_f*z1_f*sqrt(mue_f)
  rho_b=0.218735*rb*sqrt(mue_f)      
  eta_f=eta_b/(sqrt(ecm+Q))
  rho_f=rho_b*(sqrt(ecm+Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k=0)
  # penetration and shift factor
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  s_f <- rho_f*(P4$val_F*P4$val_Fp + P4$val_G*P4$val_Gp)/(P4$val_F^2 + P4$val_G^2)
  # shift factor at energy E0+Q
  xeta_f=eta_b/(sqrt(e0+Q))
  xrho_f=rho_b*(sqrt(e0+Q)) 
  PX2 <- coulomb_wave_FG(xeta_f, xrho_f, lb, k=0)
  b_f <- xrho_f*(PX2$val_F*PX2$val_Fp + PX2$val_G*PX2$val_Gp)/(PX2$val_F^2 + PX2$val_G^2)
  # partial width
  Gb <- 2*gb*p_f

  tapp <- (s_i-b_i)*ga+(s_f-b_f)*gb

  s1=pek*etpe_i*omega*Ga*Gb
  s2=((e0-ecm-tapp)^2)+0.25*((Ga+Gb)^2)
  SF <- exp( 0.5*0.98951013e0*z1_i*z2_i*sqrt(mue_i)*ue*(ecm^(-1.5)) )*s1/s2  

  return(SF = SF)
}

######################################################################
## ARTIFICIAL DATA GENERATION 
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical 
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb; 
N <-200

obsx1 <- runif(N, 0.001, 1.0)
errobsy1 <- runif(N, 0.1, 0.5)

bar <- vector()
obsy1 <- vector()

# Coc values:
# Er  = 0.35779 MeV
# g^2_in = 1.0085 MeV         ! reduced width of deuteron
# g^2_out = 0.025425 MeV      ! reduced width of neutron
# ra = 6.0 fm
# rb = 5.0 fm

bar[1] <- 0.35779   # resonance energy
bar[2] <- 1.0085    # reduced width incoming
bar[3] <- 0.025425   # reduced width outgoing
bar[4] <- 6.0
bar[5] <- 5.0
bar[6] <- 2e-4

for (i in 1:length(obsx1)){
  obsy1[i] <- rnorm( 1, Sfactor(obsx1[i], bar[1], bar[2], bar[3], bar[4], 
     bar[5], bar[6]), errobsy1[i] ) 
}

######################################################################                 
######################################################################                 
# rjags
cat('model {

###################
# LIKELIHOODS
###################
# - careful: dnorm is differently defined in R and JAGS! 
# - precision=sigma^(-2)
# - in a for loop, make sure **all** variables on the LEFT of an 
#   expression has the index [i]
# - systematic error as normalization factor y.norm...
# - ue is the electron screening potential for normal kinematics

for (i in 1:length(obsx1)) {
#
# S-FACTOR
  # measured number of decays including exp stat uncertainties
  obsy1[i] ~ dnorm(ya1[i], pow(errobsy1[i], -2))
  # actual number of decays subject to stat fluctuations
  ya1[i] ~ dnorm(yt1[i], pow(yscat, -2))    
  # true number of decays 
  yt1[i] <- sfactor3Hedpx(obsx1[i], e0, ga, gb, ra, rb, ue)    
}    

###################
# PRIORS
###################

# e0, ga, gb are:
# resonance energy, initial reduced width, final reduced width

# energy eigenvalue
#  e0 ~ dunif(0.0, 10.0)
  e0 ~ dnorm(0.0, pow(1, -2))T(0,)        # positive since we see sigma peak 

# reduced widths
  ga ~ dnorm(0.0, pow(10, -2))T(0,)
  gb ~ dnorm(0.0, pow(10, -2))T(0,)

# channel radii 
##  ra ~ dnorm(5.5, pow(1.0, -2))T(0,)
##  rb ~ dnorm(5.5, pow(1.0, -2))T(0,)
  ra ~ dnorm(6, pow(0.0001,-2))T(0,)
  rb ~ dnorm(5, pow(0.0001,-2))T(0,)

# screening potential:
  ue ~ dnorm(0.0, pow(0.1, -2))T(0,)

# extrinsic scatter
  yscat ~ dnorm(0.0, pow(1, -2))T(0,)

}', file={f <- tempfile()})

######################################################################
# n.adapt:  number of iterations in the chain for adaptation  
#           [JAGS will use to choose the sampler and to assure optimum 
#           mixing of the MCMC chain; will be discarded] 
# update(): performs the burn-in on each chain by running the MCMC for 
#           n.burn iterations without saving any of the posterior samples
# coda.samples(): runs each MCMC chain for the number of iterations 
#           specified by n.iter, but it does not save every iteration; 
#           instead, it saves only ever nth iteration, where n is given 
#           by thin
# n.chains: number of mcmc chains

n.chains <- 3
n.adapt  <- 5000   
n.burn   <- 5000  
n.iter   <- 5000 
thin     <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model;
#
ourmodel <- jags.model(f,                     ## jags wants all data in a list
              data = list('obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1
                         ),
#                    inits = list(e0 = 0.091, ga = 2.9, gb = 0.079),
                    n.chains = n.chains, n.adapt = n.adapt)
# burnin 
update(ourmodel, n.burn) 

# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                    variable.names=c('e0', 'ga', 'gb', 
                                     'ra', 'rb', 
                                     'ue', 'yscat' 
                                    ), 
                    n.iter=n.iter, thin=thin)

# <---- rjags
######################################################################
######################################################################
# OUTPUT RESULTS TO SCREEN
######################################################################
cat("", "\n")    # output empty line

# sample size adjusted for autocorrelation
effectiveChainLength = effectiveSize(mcmcChain) 
show(effectiveChainLength)

cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n") 

######################################################################
# PLOTTING
######################################################################
## set margins in order south, west, north, east
## oma is "outer margin" of entire figure
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
## mfcol=c(nrows, ncols) fills in the matrix by columns 
## tck sets tick mark lengths; negative value makes them point
##    outward
## las=1 shows tick mark labels in horizontal orientation
## mgp sets axis label locations relative to edge of inner plot window;
##   first value represents label locations [xlab, ylab], the second
##   the tick mark labels, the third the tick marks; default is c(3,1,0)
## cex controls symbol size
## cex.yy controls label sizes
## pch is the symbol id
## xlim describes plot limits
## main adds a title
library(sfsmisc) 
library(plotrix)
library(emdbook)
library(magicaxis)

######################################################################
# TRACES AND DENSITIES
######################################################################
pdf("MCMC_3Hedp_a_AD.pdf")
plot(mcmcChain)
dev.off()

######################################################################
# S-FACTOR FIT + DATA + true curve
######################################################################
pdf("MCMC_3Hedp_b_AD.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(2e-3,1)
yLim = c(0,25)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
       axes=FALSE, main="", xlab = "", ylab = "",
       cex=1.5 , cex.lab=1.3, cex.axis=1.0,
       cex.main=1.0, log="x" )
# control distance between axis and label [line=...]
title(xlab="Energy (MeV)", line=2, cex.lab=1.3)
title(ylab="S-Factor (MeV b)", line=2, cex.lab=1.3)
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.2,0))
box()

# plot legend
##legend(0.12, 3e-7, legend=c("MA97", "SCH97", "CAS02", "BYS08"), pch=c(1, 0, 6, 5))

text(0.4, 22, labels=expression(paste(NULL^"3","He(d,p)",NULL^"4","He")), cex=1.3)

# PLOT BELIEVABLE S-FACTORS

# !!! make sure to check the order of the parameters in the MCMC output !!!
samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)

# calculate for these energies and the set of Bayesian samples for
# Er, g^2_i, g^2_f the S_factor curve using Fortran code

# define grid of energy values in MeV on a log scale for plotting of credible 
# solutions; lseq is appropriate for log scale
x1 = lseq(0.001, 1, length=300)

for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1,Sfactor(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,3],
    samplesmat[i,4],samplesmat[i,5], 0),
    col=adjustcolor("red", alpha=0.02))
}

for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1,Sfactor(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,3],
    samplesmat[i,4],samplesmat[i,5], samplesmat[i,6]),
    col=adjustcolor("blue", alpha=0.02))
}

# add data - circles
points( obsx1, obsy1, col="black", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="black" )  

# PLOT INPUT S-FACTOR LINE 

lines(x1,Sfactor(x1,bar[1],bar[2],bar[3], bar[4], bar[5], bar[6]),
        col=adjustcolor("green"))

dev.off()

######################################################################
# output samples to file for rate calculation
#write.table(samplesmat, file="MCMCsamples3Hedp", quote=TRUE, 
#                        row.names=FALSE, col.names=FALSE, sep="   ")




