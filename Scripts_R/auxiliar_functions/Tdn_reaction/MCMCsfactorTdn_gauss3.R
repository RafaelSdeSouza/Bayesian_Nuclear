######################################################################
# Author: Christian Iliadis (02/12/2018)
#
# MCMCsfactorTdn_gauss3.R
#
# PURPOSE: analyzing real Tdn data
# - seven parameters are assumed: E0, Er, gamma_i^2, gamma_f^2, R_i, R_f, U_e 
#   [e0, er, gi, gf, ri, rf, ue]
# - uses the C++ function sfactorTdn3(E, Er, Er, gi, gf, ri, rf, ue), 
#   which calls gsl libraries to compute Coulomb wave functions
# - uses: normal likelihoods
#         systematic uncertainties
#         intrinsic scatter
#
#   E0: energy eigenvalue
#   Er: energy for Bc=Sc(Er), i.e., the energy at which we would like to set
#       the level shift equal to zero
#
# FUNCTIONS:
# - SfacTdn:  only used for plotting believable S-factors 
# - GammaTdn: only used for calculating partial widths 
# 
#   MAKE SURE YOU ARE USING THE SAME VALUES FOR MASSES, ENERGIES, ETC.,
#   IN BOTH THIS SCRIPT (FUNCTION Sfactor) AND THE C++ SFACTOR FUNCTION
#   IN JAGS (sfactorTdn3.cc)
#
#
# OUTPUT:
# - MCMCsfactor_a.pdf: chains and marginalized posteriors 
# - MCMCsfactor_b.pdf: S-factor vs energy
# - MCMCsfactor_c.pdf: posteriors for energy and reduced widths
# - MCMCsfactor_d.pdf: posteriors for energy and partial widths
# - MCMCsfactor_e.pdf: posterior for S-factor at given energy
# - MCMCsfactor_f.pdf: posteriors for normalization factors 
# - MCMCsfactor_g.pdf: posterior for electron screening potential
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package 
library(rjags)
## for block updating [we do not need to center predictor variables]
require(gsl) 
require(RcppGSL)
load.module("nuclear")   
load.module("glm") 
# random number seed
# set.seed(3000)

######################################################################
######################################################################
## FUNCTIONS
######################################################################
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
# compute S-factor using gsl Coulomb wave functions for plotting 
# purposes only
#
# ECM   :   center-of-mass energy
# Er    :   energy for Bc=Sc(Er), i.e., the energy at which we would 
#           like to set the level shift equal to zero
# E0    :   eigenenergy
# gi,gf :   reduced widths for incoming/outgoing channels
# ri,rf :   channel radii (fm) for incoming/outgoing channels
# ue    :   laboratory electron screening potential

SfacTdn <- function(ecm, e0, er, gi, gf, ri, rf, ue){
  # input masses, charges, angular momenta
  m1_i = 3.01550 
  m2_i = 2.01355        # masses (amu) of t and d
  m1_f = 4.00151 
  m2_f = 1.008664       # masses (amu) of n and 4He
  z1_i = 1 
  z2_i = 1              # charges of t and d
  z1_f = 2 
  z2_f = 0				# charges of n and 4He
  jt=0.5                # spins of target, projectile, resonance
  jp=1.0 
  jr=1.5                
  Q = 17.589293			# reaction Q-value (MeV)
  la = 0 
  lb = 2				# orbital angular momenta of d and n	

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
  rho_a=0.218735*ri*sqrt(mue_i)      
  eta_i=eta_a/(sqrt(ecm))
  rho_i=rho_a*(sqrt(ecm))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k=0)
  # penetration and shift factor 
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  s_i <- rho_i*(P3$val_F*P3$val_Fp + P3$val_G*P3$val_Gp)/(P3$val_F^2 + P3$val_G^2)
  # shift factor at energy Er
  xeta_i=eta_a/(sqrt(er))
  xrho_i=rho_a*(sqrt(er))
  PX1 <- coulomb_wave_FG(xeta_i, xrho_i, la, k=0)
  b_i <- xrho_i*(PX1$val_F*PX1$val_Fp + PX1$val_G*PX1$val_Gp)/(PX1$val_F^2 + PX1$val_G^2)
  # partial width
  Ga <- 2*gi*p_i

  ## outgoing channel
  eta_b=0.1574854*z2_f*z1_f*sqrt(mue_f)
  rho_b=0.218735*rf*sqrt(mue_f)      
  eta_f=eta_b/(sqrt(ecm+Q))
  rho_f=rho_b*(sqrt(ecm+Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k=0)
  # penetration and shift factor
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  s_f <- rho_f*(P4$val_F*P4$val_Fp + P4$val_G*P4$val_Gp)/(P4$val_F^2 + P4$val_G^2)
  # shift factor at energy Er+Q
  xeta_f=eta_b/(sqrt(er+Q))
  xrho_f=rho_b*(sqrt(er+Q)) 
  PX2 <- coulomb_wave_FG(xeta_f, xrho_f, lb, k=0)
  b_f <- xrho_f*(PX2$val_F*PX2$val_Fp + PX2$val_G*PX2$val_Gp)/(PX2$val_F^2 + PX2$val_G^2)
  # partial width
  Gb <- 2*gf*p_f
  
  tapp <- (s_i-b_i)*gi+(s_f-b_f)*gf
  
  s1=pek*etpe_i*omega*Ga*Gb
  s2=((e0-ecm-tapp)^2)+0.25*((Ga+Gb)^2)
  SF <- exp( 0.5*0.98951013e0*z1_i*z2_i*sqrt(mue_i)*ue*ecm^(-1.5) )*s1/s2

  return(SF = SF)
}
######################################################################
GammaTdn <- function(e0,gi,gf,ri,rf){
  # input masses, charges, angular momenta
  m1_i = 3.01550 
  m2_i = 2.01355        # masses (amu) of t and d
  m1_f = 4.00151 
  m2_f = 1.008664       # masses (amu) of n and 4He
  z1_i = 1 
  z2_i = 1              # charges of t and d
  z1_f = 2 
  z2_f = 0				# charges of n and 4He
  jt=0.5                # spins of target, projectile, resonance
  jp=1.0 
  jr=1.5                
  Q = 17.589293			# reaction Q-value (MeV)
  la = 0 
  lb = 2				# orbital angular momenta of d and n	

  # reduced masses
  mue_i <- (m1_i*m2_i)/(m1_i+m2_i)
  mue_f <- (m1_f*m2_f)/(m1_f+m2_f)

  # constants
  pek <- 6.56618216e-1/mue_i
  omega <- (2*jr+1)/((2*jt+1)*(2*jp+1))

  ## incoming channel   
  etpe_i=exp(0.98951013*z1_i*z2_i*sqrt(mue_i/e0))
  eta_a=0.1574854*z2_i*z1_i*sqrt(mue_i)
  rho_a=0.218735*ri*sqrt(mue_i)      
  eta_i=eta_a/(sqrt(e0))
  rho_i=rho_a*(sqrt(e0))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k=0)
  # penetration and shift factor 
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  # partial width
  Ga <- 2*gi*p_i

  ## outgoing channel
  eta_b=0.1574854*z2_f*z1_f*sqrt(mue_f)
  rho_b=0.218735*rf*sqrt(mue_f)      
  eta_f=eta_b/(sqrt(e0+Q))
  rho_f=rho_b*(sqrt(e0+Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k=0)
  # penetration and shift factor
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  # partial width
  Gb <- 2*gf*p_f

  return(list(Ga = Ga, Gb = Gb))
}
######################################################################
######################################################################
## DATA INPUT 
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical 
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb; 
## DATA SETS
tab1 <- read.table("data_jarmie.dat", header=FALSE)
obsx1 <- c(tab1[,1])
obsy1 <- c(tab1[,2])
errobsy1 <- c(tab1[,3])

tab2 <- read.table("data_brown.dat", header=FALSE)
obsx2 <- c(tab2[,1])
obsy2 <- c(tab2[,2])
errobsy2 <- c(tab2[,3])

tab3 <- read.table("data_kobzev.dat", header=FALSE)
obsx3 <- c(tab3[,1])
obsy3 <- c(tab3[,2])
errobsy3 <- c(tab3[,3])

tab4 <- read.table("data_arnold.dat", header=FALSE)
obsx4 <- c(tab4[,1])
obsy4 <- c(tab4[,2])
errobsy4 <- c(tab4[,3])

######################################################################                 
# rjags ----->
######################################################################                 
cat('model {

# LIKELIHOOD
# systematic error as normalization factor n.norm...
# ue is the electron screening potential for normal kinematics

for (i in 1:length(obsx1)) {
  # measured number of decays including exp stat uncertainties
  obsy1[i] ~ dnorm(ym1[i], pow(errobsy1[i], -2))
  # measured number of decays including exp syst uncertainties
  ym1[i] <- n.norm1 * ya1[i]
  # actual number of decays subject to stat fluctuations
  ya1[i] ~ dnorm(yt1[i], pow(intrscat, -2))    
  # true number of decays 
  yt1[i] <- sfactorTdn(obsx1[i], e0, er, gi, gf, ri, rf, ue)    
}    

for (i in 1:length(obsx2)) {
  # measured number of decays including exp stat uncertainties
  obsy2[i] ~ dnorm(ym2[i], pow(errobsy2[i], -2))
  # measured number of decays including exp syst uncertainties
  ym2[i] <- n.norm2 * ya2[i]
  # actual number of decays subject to stat fluctuations
  ya2[i] ~ dnorm(yt2[i], pow(intrscat, -2))    
  # true number of decays 
  yt2[i] <- sfactorTdn(obsx2[i], e0, er, gi, gf, ri, rf, ue)    
}    

for (i in 1:length(obsx3)) {
  # measured number of decays including exp stat uncertainties
  obsy3[i] ~ dnorm(ym3[i], pow(errobsy3[i], -2))
  # measured number of decays including exp syst uncertainties
  ym3[i] <- n.norm3 * ya3[i]
  # actual number of decays subject to stat fluctuations
  ya3[i] ~ dnorm(yt3[i], pow(intrscat, -2))    
  # true number of decays 
  yt3[i] <- sfactorTdn(obsx3[i], e0, er, gi, gf, ri, rf, ue)    
}    

for (i in 1:length(obsx4)) {
  # measured number of decays including exp stat uncertainties
  obsy4[i] ~ dnorm(ym4[i], pow(errobsy4[i], -2))
  # measured number of decays including exp syst uncertainties
  ym4[i] <- n.norm4 * ya4[i]
  # actual number of decays subject to stat fluctuations
  ya4[i] ~ dnorm(yt4[i], pow(intrscat, -2))    
  # true number of decays 
  yt4[i] <- sfactorTdn(obsx4[i], e0, er, gi, gf, ri, rf, ue)    
}    

# PRIORS
#
# Barker values:
# E0  = 0.0912 MeV
# g^2_in = 2.93 MeV         ! reduced width of deuteron
# g^2_out = 0.0794 MeV      ! reduced width of neutron
#
# Er = 0.0912 MeV
# ri = 6.0 fm
# rf = 5.0 fm
#
# e0:       eigenenergy
# er:       level shift is set zero at Er
# gi,gf:    initial reduced width, final reduced width [gi,gf = gamma^2]
#
# Wigner limit: wl = hbar^2/(m_red a_c^2) = 41.80159/(M_red a_c^2)
#
# deuteron channel: wl_d = 41.80159/(1.207357 a_c^2) = 34.6224/a_c^2
# neutron channel:  wl_n = 41.80159/(0.805597 a_c^2) = 51.8889/a_c^2
  wl_d <- 34.6224*pow(ri, -2)
  wl_n <- 51.8889*pow(rf, -2)
# gives for Barkers choices of radii:
# wl_d=0.962 MeV
# wl_n=2.076 MeV

###################
## first set of priors for paper [fixed a_c and B_c]
## Barker values for radii and boundary condition parameters
##  ri <- 6.0
##  rf <- 5.0
##  er <- 0.0912                  # level shift zero at er
##  ue <- 0.0                     # zero electron screening potential
#  ue ~ dunif(0, 0.01)           # certainly less than 10 keV
##  intrscat ~ dunif(0, 5)        # certainly less than 5 MeVb
##  e0 ~ dunif(0.0, 1.0)          # postive since we see sigma peak 
##  gi ~ dunif(0.0, 10*wl_d)       # x times Wigner limit
##  gf ~ dunif(0.0, 10*wl_n)       # x times Wigner limit
###################
# second set of priors by including a_c and B_c in sampling
#  ri ~ dnorm(5.5, pow(1.0, -2))T(0,)
#  rf ~ dnorm(5.5, pow(1.0, -2))T(0,)
  ri ~ dunif(2.5, 10.0)
  rf ~ dunif(2.5, 10.0)
  er ~ dunif(0.0, 1.0)          # level shift zero at er
  ue <- 0.0                     # zero electron screening potential
  intrscat ~ dunif(0, 5)        # certainly less than 5 MeVb
  e0 ~ dunif(0.02, 0.08)          # postive since we see sigma peak 
  gi ~ dunif(0.0, 10*wl_d)       # x times Wigner limit
  gf ~ dunif(0.0, 10*wl_n)       # x times Wigner limit

################
# systematic errors as normalization factors

### lognormal density:
# log(): natural logarithm
n.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
logmu1 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma1 <- log(1.0126) # factor uncertainty is 1.0126, i.e. 1.26% for Jar84

n.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
logmu2 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma2 <- log(10.0)   # factor uncertainty very large since Bro87 used
                         # the normalization of Jar84
n.norm3 ~ dlnorm(logmu3, pow(logsigma3, -2))
logmu3 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma3 <- log(1.025)  # factor uncertainty is 1.025, i.e., 2.5% for Kob66

n.norm4 ~ dlnorm(logmu4, pow(logsigma4, -2))
logmu4 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma4 <- log(1.020)  # factor uncertainty is 1.020, i.e., 2.0% for Arn53 

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
n.chains = 1
n.adapt = 200 
n.burn = 200     
n.iter = 1000   
thin = 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model
ourmodel <- jags.model(f, data = list(
                            'obsx1' = obsx1, ## jags wants all data in a list
                            'obsy1' = obsy1,
                            'errobsy1' = errobsy1,
                            'obsx2' = obsx2, 
                            'obsy2' = obsy2,
                            'errobsy2' = errobsy2,
                            'obsx3' = obsx3, 
                            'obsy3' = obsy3,
                            'errobsy3' = errobsy3,
                            'obsx4' = obsx4, 
                            'obsy4' = obsy4,
                            'errobsy4' = errobsy4
                            ),
#    inits = list(e0 = 0.090, gi = 4.0, gf = 0.1),
                    n.chains = n.chains,
                    n.adapt = n.adapt
                    )
update(ourmodel, n.burn)
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                    variable.names=c('e0', 'er', 'gi', 'gf', 'intrscat', 
                    'n.norm1', 'n.norm2', 'n.norm3', 'n.norm4',
                    'ri', 'rf', 'ue'), 
                    n.iter=n.iter, thin=thin)


source("SfacTdn.R")
source("nuclear_rate_Tdn.R")
source("reaction_table_Tdn.R")
Tgrid <- c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,
           0.013,0.014,0.015,0.016,0.018,0.020,0.025,0.030,0.040,0.050,0.060,0.070,
           0.080,0.090,0.100,0.110,0.120,0.130,0.140,0.150,0.160,0.180,0.200,0.250,0.300,
           0.350,0.400,0.450,0.500,0.600,0.700,0.800,0.900,1.000,1.250,1.500,1.750,2.000,2.500,3.000,3.500,4.000,5.000,
           6.000,7.000,8.000,9.000,10.000)

reaction_table_Tdn(mcmcChain,vars = c('e0', 'er', 'gi', 'gf','ri', 'rf','ue'),N=100,T9=Tgrid[1:4])




######################################################################
# <---- rjags
######################################################################
# sample size adjusted for autocorrelation
effectiveChainLength = effectiveSize(mcmcChain) 
show(effectiveChainLength)

# output results on screen
cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n") 
# PLOT TRACES AND DENSITIES
pdf("MCMCsfactor_a.pdf")
plot(mcmcChain)
dev.off()

######################################################################
## PLOTTING
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
require("logspline")

# !!! make sure to check the order of the parameters in the MCMC output !!!
samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)

######################################################################
# S-FACTOR
######################################################################
# first determine plot ranges
pdf("MCMCsfactor_b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(4e-3,0.3)
yLim = c(0,30)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
       axes=FALSE, main="", xlab = "", ylab = "", log="x", yaxs='i', xaxs='i' )
# control distance between axis and label [line=...]
title(xlab="Energy (MeV)", line=2.5, cex.lab=2.0)
title(ylab="S-Factor (MeV b)", line=2.5, cex.lab=2.0)

# control distance tick mark labels and axis; don't touch first number
# in mgp; second number controls distance tick mark labels and axis
# don't touch third number
# cex.axis controls size of tick mark labels
magaxis(mgp=c(0,0.4,0), cex.axis=1.3)
box()

# plot legend
legend(0.15, 27.0, legend=c("Jar84", "Bro87", "Kob66", "Arn53"), 
        pch=c(1, 5, 0, 6), cex=1.2,
        col=c("gray40", "gray40", "gray40", "gray40"))

# or plot in color:
#legend(0.15, 27.0, legend=c("Jar84", "Bro87", "Kob66", "Arn53"), 
#        pch=c(1, 18, 0, 6), 
#        col=c("red", "black", "green4", "blue"))

# plot reaction label
text(0.008, 27, labels=expression(paste(NULL^"3","H(d,n)",NULL^"4","He")), cex=2.0)

#################
# plot credible S-factors

# calculate for these energies and the set of Bayesian samples for
# Er, g^2_i, g^2_f the S_factor curve using R code above

# define grid of energy values in MeV on a log scale for plotting of credible 
# solutions; lseq is appropriate for log scale
x1 = lseq(0.001, 0.5, length=300)

## plot all credible lines between fixed limits
#for ( i in round(seq(from=1,to=nsamp,length=5000)) ) {
## randomly pick a subset of credible lines
# ECM, E0, Er, gi, gf, ri, rf, ue

# for inclusion of electron screening:
#for ( i in round(runif(500, min=1, max=nsamp)) ) {
#   lines(x1,SfacTdn(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,4],
#    samplesmat[i,3],samplesmat[i,11],samplesmat[i,10], samplesmat[i,12]),
#    col=adjustcolor("black", alpha=0.02))
#}

# 500 samples are randomly selected for plotting of S-factor
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1,SfacTdn(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,4],
    samplesmat[i,3],samplesmat[i,11],samplesmat[i,10], 0.0),
    col=adjustcolor("red", alpha=0.02), lw=0.1)
}
#################
# add data - circles
points( obsx1, obsy1, col="gray40", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="gray40" )  

points( obsx2, obsy2, col="gray40", pch=5, cex=1.0 )
add.error.bars(obsx2, obsy2, 0.0, errobsy2, 0.0, col="gray40" )  

points( obsx3, obsy3, col="gray40", pch=0, cex=1.0 )
add.error.bars(obsx3, obsy3, 0.0, errobsy3, 0.0, col="gray40" )  

points( obsx4, obsy4, col="gray40", pch=6, cex=1.0 )
add.error.bars(obsx4, obsy4, 0.0, errobsy4, 0.0, col="gray40" )  

# for plotting in color:
#points( obsx1, obsy1, col="red", pch=1, cex=1.2 )
#add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="red" )  

#points( obsx2, obsy2, col="black", pch=18, cex=1.8 )
#add.error.bars(obsx2, obsy2, 0.0, errobsy2, 0.0, col="black" )  

#points( obsx3, obsy3, col="green4", pch=0, cex=1.0 )
#add.error.bars(obsx3, obsy3, 0.0, errobsy3, 0.0, col="green4" )  

#points( obsx4, obsy4, col="blue", pch=6, cex=1.0 )
#add.error.bars(obsx4, obsy4, 0.0, errobsy4, 0.0, col="blue" )  

dev.off()
######################################################################
# POSTERIORS OF RESONANCE ENERGY AND REDUCED WIDTHS
######################################################################
# first determine plot ranges
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
pdf("MCMCsfactor_c.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,3), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)
   
# plot eigenenergy in first panel           
plot(density(samplesmat[,1]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(ylab="Probability (arb. units)", line=4.5, cex.lab=2.3)
title(xlab=expression(paste(E [0], " (MeV)")), line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,1]), col=adjustcolor("blue", alpha=0.5))

# E0,gi,gf,ri,rf
 
# plot deuteron partial width                   
plot(density(samplesmat[,4]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(gamma [d]^2, " (MeV)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')

polygon(density(samplesmat[,4]), col=adjustcolor("blue", alpha=0.5))

# plot neutron partial width
plot(density(samplesmat[,3]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(gamma [n]^2, " (MeV)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
    
polygon(density(samplesmat[,3]), col=adjustcolor("blue", alpha=0.5))
    
dev.off()

######################################################################
# POSTERIORS OF RESONANCE ENERGY AND PARTIAL WIDTHS
######################################################################
# first determine plot ranges
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
pdf("MCMCsfactor_d.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,3), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)
   
# plot eigenenergy in first panel           
plot(density(samplesmat[,1]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(ylab="Probability (arb. units)", line=4.5, cex.lab=2.3)
title(xlab=expression(paste(E [0], " (MeV)")), line=4.0, cex.lab=2.3)

# E0,gi,gf,ri,rf
 
# plot deuteron partial width                   
plot(density(
GammaTdn(samplesmat[,1], samplesmat[,4], samplesmat[,3], samplesmat[,11], 
   samplesmat[,10])$Ga), main="", xlab="", ylab="",
   cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5)
title(xlab=expression(paste(Gamma [d], " (MeV)")), line=4.0, cex.lab=2.3, 
   yaxs='i', xaxs='i')

# plot neutron partial width
plot(density(
GammaTdn(samplesmat[,1], samplesmat[,4], samplesmat[,3], samplesmat[,11], samplesmat[,10])$Gb
), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5)
title(xlab=expression(paste(Gamma [n], " (MeV)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
    
dev.off()
######################################################################
# POSTERIOR OF S-FACTOR AT GIVEN ENERGY
######################################################################
# prediction: plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

# define x value for which we would like to predict y
xchoice <- 0.04

# declare vector with y values
fitvec <- vector(mode = "numeric", length = nsamp)   

# calculate y values from all samples at given x value

for(i in 1:nsamp) fitvec[i] <- SfacTdn(xchoice,samplesmat[i,1],samplesmat[i,2],samplesmat[i,4],
    samplesmat[i,3],samplesmat[i,11],samplesmat[i,10], 0.0)

# define quantiles of y at xchoice
a16 <- quantile(fitvec, prob = 0.16)
a50 <- quantile(fitvec, prob = 0.50)
a84 <- quantile(fitvec, prob = 0.84)

# output
cat("", "\n") 
cat("PREDICTION FOR x.choice=", xchoice, "\n") 
cat("  16%        50%        84%", "\n")
cat(quantile(fitvec, probs = 0.16), quantile(fitvec, probs = 0.50),
             quantile(fitvec, probs = 0.84), "\n")

# plot density at xchoice
pdf("MCMCsfactor_e.pdf")
par(mfcol=c(1,1), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)

plot(density(fitvec), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, yaxs='i', xaxs='i')
###
polygon(density(fitvec), col=adjustcolor("blue", alpha=0.5))
###
title(xlab=expression(paste(S-factor, "  (MeV b)")), line=4.0, cex.lab=2.0)
title(ylab="Probability (arb. units)", line=4.5, cex.lab=2.0)
dev.off()
######################################################################
# DENSITIES OF NORMALIZATION FACTORS
######################################################################
pdf("MCMCsfactor_f.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(0.95, 1.04)
#yLim = c(0, 120)

# plot density in first panel           
plot(density(samplesmat[,9]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, 
     xlim=xLim, yaxs='i', xaxs='i')
title(ylab="Probability (arb. units)", line=4.5, cex.lab=2.3)
title(xlab="Normalization factor", line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,9]),
     col=adjustcolor("purple", alpha=0.5))

polygon(density(samplesmat[,7]), 
     col=adjustcolor("red", alpha=0.5))

polygon(density(samplesmat[,8]),  
     col=adjustcolor("green", alpha=0.5))

polygon(density(samplesmat[,6]),  
     col=adjustcolor("blue", alpha=0.5))

legend("topleft", inset=.02, 
   c("Jar84", "Bro87", "Kob66", "Arn53"), 
   fill=adjustcolor(c("blue", "red", "green", "purple"), alpha=0.5), 
   horiz=FALSE, cex=1.5, box.lty=0)

dev.off()

######################################################################
# POSTERIOR OF ELECTRON SCREENING POTENTIAL
######################################################################
pdf("MCMCsfactor_g.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.5,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(0, 30)
#yLim = c(0, 500)

dens <- density(1e6*samplesmat[,12])
a98 <- quantile(1e6*samplesmat[,12], prob = 0.975)
#m <- logspline(1e6*samplesmat[,12], lbound = 0)

# plot density in first panel  
# cex.axis controls tick mark label size
# xaxs and yaxs = "i" plots with exact limits         
plot(dens, main="", xlab="", ylab="",
     cex.axis=1.5, yaxs='i', xaxs='i', xlim=xLim)

title(ylab="Probability (arb. units)", line=4.5, cex.lab=2.3)
title(xlab="Electron screening potential (eV)", line=4.0, cex.lab=2.3)

polygon(dens, col=adjustcolor("blue", alpha=0.5))

dev.off()
######################################################################
# OUTPUT TABLES
######################################################################
# output samples to file for rate calculation; 
# chose only 5000 samples at random from full results
samplesmatx = samplesmat[sample(nrow(samplesmat), size=5000, replace=FALSE),]
write.table(samplesmatx, file="MCMCsamplesTdn", quote=TRUE, 
                        row.names=FALSE, col.names=FALSE, sep="   ")

# output S-factor density at given energy 
## write.table(fitvec, file="MCMCsdensityTdn", quote=TRUE, 
##                        row.names=FALSE, col.names=FALSE, sep="   ")



