######################################################################
# Author: Christian Iliadis (05/13/2019)
######################################################################
#
# MCMCsfactor_Tdn.R
#
# PURPOSE: analyzing real Tdn data
#
# - uses the C++ function sfactorTdn(E, E0, Er, gi, gf, ri, rf, ue),
#   which calls gsl libraries to compute Coulomb wave functions
#
#   E0: energy eigenvalue
#   Er: energy for Bc=Sc(Er), i.e., the energy at which we would like to
#       set the level shift equal to zero
#
# FEATURES: 
#   -- statistical uncertainties on S-factor
#   -- systematic uncertainties on S-factor
#   -- extrinsic scatter on S-factor
#
# FUNCTIONS:
# - SfacTdn:  only used for plotting believable S-factors 
# - GammaTdn: only used for calculating and plotting partial widths 
# 
#   MAKE SURE YOU ARE USING THE SAME VALUES FOR MASSES, ENERGIES, ETC.,
#   IN BOTH THIS SCRIPT (FUNCTION Sfactor) AND THE C++ SFACTOR FUNCTION
#   IN JAGS (sfactorTdn.cc)
#
# OUTPUT:
# - MCMCsamplesTdn: 5,000 parameter set samples, chosen randomly from
#                   all samples
#   [can be used to calculate rates, or re-plot figures] 
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package
library("rjags")
library(magicaxis)
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
# for plotting purposes only
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
  z2_f = 0				          # charges of n and 4He
  jt = 0.5              # spins of target, projectile, resonance
  jp = 1.0 
  jr = 1.5                
  Q = 17.589293			      # reaction Q-value (MeV)
  la = 0 
  lb = 2				            # orbital angular momenta of d and n	

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
# calculate partial widths;
# for plotting purposes only

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
## DATA INPUT 
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, S-factor in MeVb

## DATA SET 1: jarmie
obsx1    <- c(4.6809E-02, 4.6009E-02, 4.4007E-02, 4.2005E-02, 4.0004E-02,
              3.6001E-02, 3.1998E-02, 2.7996E-02, 2.3994E-02, 1.9992E-02,
              1.5990E-02, 1.1989E-02,  9.989E-03,  7.990E-03,  6.990E-03,
              5.990E-03,  4.992E-03)
obsy1    <- c(2.6640E+01, 2.6740E+01, 2.6300E+01, 2.6000E+01, 2.5280E+01,
              2.4020E+01, 2.2190E+01, 2.0700E+01, 1.8870E+01, 1.7350E+01,
              1.5810E+01, 1.4320E+01, 1.3920E+01, 1.3430E+01, 1.2830E+01,
              1.3480E+01, 1.2630E+01)
errobsy1 <- c(1.3853E-01, 1.2835E-01, 1.4202E-01, 1.2220E-01, 1.4157E-01,
              1.1289E-01, 1.1095E-01, 9.3150E-02, 7.9254E-02, 8.8485E-02,
              1.2964E-01, 9.7376E-02, 1.4059E-01, 2.6994E-01, 4.0414E-01,
              3.8553E-01, 5.7972E-01)

## DATA SET 2: brown
obsx2    <- c(4.7948E-02, 5.0947E-02, 5.3942E-02, 5.6942E-02, 5.9941E-02,
              6.2941E-02, 6.5941E-02, 6.9541E-02)
obsy2    <- c(2.6480E+01, 2.6840E+01, 2.5890E+01, 2.5500E+01, 2.4330E+01,
              2.3440E+01, 2.2020E+01, 2.0340E+01)
errobsy2 <- c(2.1184E-01, 2.1472E-01, 2.0712E-01, 2.0400E-01, 1.9464E-01,
              1.8752E-01, 1.7616E-01, 1.6272E-01)

## DATA SET 3: kobzev
obsx3    <- c( 46.0E-03,  48.0E-03,  52.0E-03,  56.0E-03,  60.0E-03,  64.0E-03,
               66.0E-03,  68.0E-03,  72.0E-03,  76.0E-03,  80.0E-03,  84.0E-03,
               88.0E-03,  92.0E-03,  96.0E-03, 100.0E-03, 104.0E-03, 108.0E-03,
              112.0E-03, 116.0E-03, 120.0E-03, 124.0E-03, 128.0E-03, 132.0E-03,
              136.0E-03, 140.0E-03, 144.0E-03, 148.0E-03, 152.0E-03, 156.0E-03,
              160.0E-03, 164.0E-03, 168.0E-03, 176.0E-03, 184.0E-03, 192.0E-03,
              200.0E-03, 208.0E-03, 216.0E-03, 224.0E-03, 232.0E-03, 240.0E-03,
              248.2E-03, 256.2E-03, 264.3E-03)
obsy3    <- c(    25.93,     25.96,     25.76,     25.28,     24.77,     23.66, 
                  22.85,     21.89,     19.98,     18.14,     16.53,     15.01, 
                  13.65,     12.50,     11.41,     10.45,      9.59,      8.76, 
                   7.98,      7.28,      6.65,      6.08,      5.61,      5.23, 
                   4.89,      4.60,      4.32,      4.11,      3.88,      3.69, 
                   3.50,      3.32,      3.15,      2.84,      2.62,      2.42, 
                   2.26,      2.13,      2.00,      1.89,      1.79,      1.69, 
                   1.60,      1.51,      1.44)
errobsy3 <- c(     0.52,      0.52,      0.52,      0.51,      0.50,      0.47,
                   0.46,      0.44,      0.40,      0.36,      0.33,      0.30,
                   0.27,      0.25,      0.23,      0.21,      0.19,      0.18,
                   0.16,      0.15,      0.13,      0.12,      0.11,      0.10,
                   0.10,      0.09,      0.09,      0.08,      0.08,      0.07,
                   0.07,      0.08,      0.08,      0.07,      0.07,      0.06,
                   0.06,      0.05,      0.05,      0.05,      0.04,      0.04,
                   0.04,      0.04,      0.04)

## DATA SET 4: arnold
obsx4    <- c(0.029951,  0.026383,  0.026095,  0.043292, 0.0093180,  0.020273,  
              0.031522,  0.014679,  0.037164,  0.025724,  0.025664,  0.037002, 
              0.014481,  0.031162,  0.019925, 0.0089763,  0.042495,  0.025322, 
             0.0095159, 0.0094739,  0.014894,  0.012826,  0.012808,  0.047250, 
              0.047220,  0.018348,  0.018330,  0.041254,  0.041230,  0.023973, 
              0.023955,  0.035383,  0.035359,  0.052826,  0.052796,  0.058678, 
              0.058660,  0.064543,  0.064513,  0.070443,  0.070395,  0.067391, 
              0.067373,  0.061431,  0.061395,  0.046650,  0.046638,  0.046614,    
              0.012034,  0.011992,  0.011950,  0.025262,  0.025166)
obsy4    <- c(  21.766,    20.525,    20.277,    27.067,    13.703,    17.721,
                22.695,    15.753,    25.184,    20.596,    19.920,    24.967,
                14.939,    22.749,    17.249,    13.340,    26.847,    19.969,
                13.600,    13.508,    15.448,    14.957,    14.302,    27.542,
                27.505,    16.989,    16.921,    26.514,    26.600,    18.366,
                18.969,    24.589,    24.314,    27.085,    26.975,    25.669,
                25.621,    23.157,    23.071,    20.227,    20.445,    21.951,
                22.002,    24.492,    24.593,    27.489,    27.365,    27.466,
                13.680,    13.849,    14.068,    20.755,    20.718)
errobsy4 <- c(0.021766,  0.020525,  0.020277,  0.027067,  0.027400,  0.035400,
              0.022695,  0.031400,  0.025184,  0.020596,  0.019920,  0.024967,
              0.030000,  0.022749,  0.034400,  0.026700,  0.026847,  0.019969,
              0.027200,  0.027000,  0.030000,  0.030000,  0.029000,  0.027542,
              0.027505,  0.032000,  0.034000,  0.026514,  0.026600,  0.036000,
              0.038000,  0.024589,  0.024314,  0.027085,  0.026975,  0.025669,
              0.025621,  0.023157,  0.023071,  0.020227,  0.020445,  0.021951,
              0.022002,  0.024492,  0.024593,  0.027489,  0.027365,  0.027466,
              0.027000,  0.028000,  0.028000,  0.020755,  0.020718)

## DATA SET 5: conner
obsx5    <- c(1.2420E-02, 1.5480E-02, 1.8600E-02, 2.0700E-02, 2.1780E-02, 2.4900E-02,
              2.8020E-02, 2.9100E-02, 3.1200E-02, 3.3240E-02, 3.4260E-02, 3.7380E-02,
              4.0500E-02, 4.1580E-02, 4.3680E-02, 4.5720E-02, 4.6800E-02, 4.9980E-02,
              5.4180E-02, 5.6220E-02, 5.8260E-02, 6.2400E-02, 6.5400E-02, 6.6600E-02,
              6.9000E-02, 7.5000E-02, 8.0400E-02, 8.1600E-02, 8.5800E-02, 8.7600E-02,
              9.1800E-02, 9.3600E-02, 9.7200E-02, 1.0020E-01, 1.0380E-01, 1.0980E-01,
              1.2300E-01, 1.3620E-01, 1.5060E-01, 1.6560E-01, 1.8120E-01, 1.9740E-01,
              2.1420E-01)
obsy5    <- c(1.3227E+01, 1.5174E+01, 1.5790E+01, 1.7328E+01, 1.7378E+01, 1.8234E+01,
              1.9699E+01, 2.0129E+01, 2.1804E+01, 2.2911E+01, 2.1594E+01, 2.3803E+01,
              2.5309E+01, 2.5725E+01, 2.5930E+01, 2.5899E+01, 2.5439E+01, 2.6833E+01,
              2.5526E+01, 2.6602E+01, 2.5886E+01, 2.4608E+01, 2.3429E+01, 2.2901E+01,
              2.1818E+01, 1.9229E+01, 1.6970E+01, 1.6604E+01, 1.4961E+01, 1.4266E+01,
              1.2897E+01, 1.2332E+01, 1.1024E+01, 1.0626E+01, 9.9085E+00, 8.9950E+00,
              6.7947E+00, 5.4398E+00, 4.4269E+00, 3.5523E+00, 2.8963E+00, 2.5082E+00,
              2.1603E+00)
errobsy5 <- c(      0.13,       0.15,       0.16,       0.17,       0.17,       0.18,
                    0.20,       0.20,       0.22,       0.23,       0.21,       0.24,
                    0.25,       0.26,       0.26,       0.26,       0.25,       0.27,
                    0.26,       0.27,       0.26,       0.25,       0.23,       0.23,
                    0.22,       0.20,       0.17,       0.17,       0.15,       0.14,
                    0.13,       0.12,       0.11,       0.11,       0.10,       0.09,
                    0.07,       0.05,       0.04,       0.04,       0.03,       0.03,
                    0.02)

######################################################################                 
######################################################################                 
# rjags ----->
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
  # ...subject to extrinsic scatter, if any:
  obsy1[i] ~ dnorm(ya1[i], pow(yscat1, -2))    
  # ...subject to stat uncertainties:
  ya1[i] ~ dnorm(ym1[i], pow(errobsy1[i], -2))
  # ...subject to syst uncertainties: 
  ym1[i] <- y.norm1 * yt1[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt1[i] <- sfactorTdn(obsx1[i], e0, er, gi, gf, ri, rf, ue)
}    

for (i in 1:length(obsx2)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter, if any:
  obsy2[i] ~ dnorm(ya2[i], pow(yscat2, -2))    
  # ...subject to stat uncertainties:
  ya2[i] ~ dnorm(ym2[i], pow(errobsy2[i], -2))
  # ...subject to syst uncertainties: 
  ym2[i] <- y.norm2 * yt2[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt2[i] <- sfactorTdn(obsx2[i], e0, er, gi, gf, ri, rf, ue)
}    

for (i in 1:length(obsx3)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter, if any:
  obsy3[i] ~ dnorm(ya3[i], pow(yscat3, -2))    
  # ...subject to stat uncertainties:
  ya3[i] ~ dnorm(ym3[i], pow(errobsy3[i], -2))
  # ...subject to syst uncertainties: 
  ym3[i] <- y.norm3 * yt3[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt3[i] <- sfactorTdn(obsx3[i], e0, er, gi, gf, ri, rf, ue)
}    

for (i in 1:length(obsx4)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter, if any:
  obsy4[i] ~ dnorm(ya4[i], pow(yscat4, -2))    
  # ...subject to stat uncertainties:
  ya4[i] ~ dnorm(ym4[i], pow(errobsy4[i], -2))
  # ...subject to syst uncertainties: 
  ym4[i] <- y.norm4 * yt4[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt4[i] <- sfactorTdn(obsx4[i], e0, er, gi, gf, ri, rf, ue)
}    

for (i in 1:length(obsx5)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter, if any:
  obsy5[i] ~ dnorm(ya5[i], pow(yscat5, -2))    
  # ...subject to stat uncertainties:
  ya5[i] ~ dnorm(ym5[i], pow(errobsy5[i], -2))
  # ...subject to syst uncertainties: 
  ym5[i] <- y.norm5 * yt5[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt5[i] <- sfactorTdn(obsx5[i], e0, er, gi, gf, ri, rf, ue)
}    

###################
# PRIORS
###################
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
#
  wl_d <- 34.6224*pow(ri, -2)
  wl_n <- 51.8889*pow(rf, -2)
#
# gives for Barkers choices of radii:
# wl_d=0.962 MeV
# wl_n=2.076 MeV

# R-matrix reduced widths:
  gi ~ dnorm(0.0, pow(wl_d, -2))T(0,)      # 1 MeV is WL for Barker choice of radius
  gf ~ dnorm(0.0, pow(wl_n, -2))T(0,)      # 2 MeV is WL for Barker choice of radius

# channel radius included in sampling; 
# B_c is not sampled; we chose Er = E0
  e0 ~ dnorm(0.0, pow(1, -2))T(0,)         # positive since we see sigma peak 
  er <- e0                                 # set Er = E0
  ri ~ dunif(2.5, 8.0)
  rf ~ dunif(2.5, 8.0)

# screening potential:
  ue ~ dnorm(0.0, pow(0.001, -2))T(0,)     # certainly less than 1000 eV
#  ue <- 0.0                               # zero electron screening potential

# extrinsic scatter
  yscat1 ~ dnorm(0.0, pow(5, -2))T(0,)     # sigma is 5 MeVb
  yscat2 ~ dnorm(0.0, pow(5, -2))T(0,)
  yscat3 ~ dnorm(0.0, pow(5, -2))T(0,)
  yscat4 ~ dnorm(0.0, pow(5, -2))T(0,)
  yscat5 ~ dnorm(0.0, pow(5, -2))T(0,)

# systematic normalization factor for S-factor:
# log(): natural logarithm
y.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
logmu1 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma1 <- log(1.0126) # factor uncertainty is 1.0126, i.e. 1.26% for Jar84

y.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
logmu2 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma2 <- log(10.0)   # factor uncertainty very large since Bro87 used
                         # the normalization of Jar84
y.norm3 ~ dlnorm(logmu3, pow(logsigma3, -2))
logmu3 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma3 <- log(1.025)  # factor uncertainty is 1.025, i.e., 2.5% for Kob66

y.norm4 ~ dlnorm(logmu4, pow(logsigma4, -2))
logmu4 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma4 <- log(1.020)  # factor uncertainty is 1.020, i.e., 2.0% for Arn53 

y.norm5 ~ dlnorm(logmu5, pow(logsigma5, -2))
logmu5 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma5 <- log(1.018)  # factor uncertainty is 1.018, i.e., 1.8% for Con52 

}', file={f <- tempfile()})

# dunif(min=0, max=1): gives uniform density   
# pow(a,b) = a^b
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

n.chains = 3
n.adapt = 1000 
n.burn = 1000     
n.iter = 1000  
thin = 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model
ourmodel <- jags.model(f, data = list(          ## jags wants all data in a list
                'obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1,
                'obsx2' = obsx2, 'obsy2' = obsy2, 'errobsy2' = errobsy2,
                'obsx3' = obsx3, 'obsy3' = obsy3, 'errobsy3' = errobsy3,
                'obsx4' = obsx4, 'obsy4' = obsy4, 'errobsy4' = errobsy4,
                'obsx5' = obsx5, 'obsy5' = obsy5, 'errobsy5' = errobsy5
                                    ),
                inits = list(e0 = 0.08, gi = 3.2, gf = 0.13, ri = 5.5, rf = 3.6), 
                n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn)

# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                variable.names=c(
                'e0', 'er', 'gi', 'gf', 'ri', 'rf', 'ue',
                'y.norm1', 'y.norm2', 'y.norm3', 'y.norm4', 'y.norm5',
                'yscat1',  'yscat2',   'yscat3',  'yscat4', 'yscat5'                    
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
# OUTPUT RESULTS TO FILES
######################################################################
# !!! make sure to check the order of the parameters in the MCMC output !!!
samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)
# output all results
capture.output(print(samplesmat, print.gap=3), file="MCMCresultsTdn")

# select first six colums [e0, er, gf, gi, rf, ri] from matrix
samplesmat2 <- samplesmat[,1:6]
# output samples to file for rate calculation
capture.output( print(samplesmat2[sample(nrow(samplesmat2), size=1000, 
                  replace=FALSE),], 
            print.gap=3), file="MCMCsamplesTdn" )

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
pdf("MCMC_Tdn_a.pdf")
plot(mcmcChain)
dev.off()

######################################################################
# S-FACTOR FIT [red] + DATA [black & white]
######################################################################
pdf("MCMC_Tdn_b.pdf",width=10,height=5,onefile=F)
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
legend(0.15, 27.0, legend=c("Jar84", "Bro87", "Kob66", "Arn53", "Con52"), 
        pch=c(1, 5, 0, 6, 2),
        col=c("gray40", "gray40", "gray40", "gray40", "gray40"))

# plot reaction label
text(0.008, 27, labels=expression(paste(NULL^"3","H(d,n)",NULL^"4","He")), cex=2.0)

# plot credible S-factors:
# calculate for these energies and the set of Bayesian samples for
# Er, g^2_i, g^2_f the S_factor curve using R code above

# define grid of energy values in MeV on a log scale for plotting of credible 
# solutions; lseq is appropriate for log scale
x1 = lseq(0.001, 0.5, length=300)

# ECM, E0, Er, gi, gf, ri, rf, ue

# 500 samples are randomly selected for plotting of S-factor including
# electron screening
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1,SfacTdn(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,4],
    samplesmat[i,3],samplesmat[i,6],samplesmat[i,5], samplesmat[i,7]),
    col=adjustcolor("red", alpha=0.02), lw=0.1)
}

# add data 
points( obsx1, obsy1, col="gray40", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="gray40" )  

points( obsx2, obsy2, col="gray40", pch=5, cex=1.0 )
add.error.bars(obsx2, obsy2, 0.0, errobsy2, 0.0, col="gray40" )  

points( obsx3, obsy3, col="gray40", pch=0, cex=1.0 )
add.error.bars(obsx3, obsy3, 0.0, errobsy3, 0.0, col="gray40" )  

points( obsx4, obsy4, col="gray40", pch=6, cex=1.0 )
add.error.bars(obsx4, obsy4, 0.0, errobsy4, 0.0, col="gray40" )  

points( obsx5, obsy5, col="gray40", pch=2, cex=1.0 )
add.error.bars(obsx5, obsy5, 0.0, errobsy5, 0.0, col="gray40" )  

dev.off()

######################################################################
# S-FACTOR DATA [color] ONLY
######################################################################
pdf("MCMC_Tdn_c.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(4e-3,0.3)
yLim = c(0,30)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
       axes=FALSE, main="", xlab = "", ylab = "",
       cex=1.5 , cex.lab=1.3, cex.axis=1.0,
       cex.main=1.0, log="x" )
# control distance between axis and label [line=...]
title(xlab="Energy (MeV)", line=2, cex.lab=2.0)
title(ylab="S-Factor (MeV b)", line=2, cex.lab=2.0)
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

# shade temperature region [low left corner to high right corner]
#rect(0.004, -2, 0.12, 32, col="snow2", border = NA)

magaxis(mgp=c(0,0.2,0))
box()

# plot legend
legend(0.15, 27.0, legend=c("Jar84", "Bro87", "Kob66", "Arn53", "Con52"), 
        pch=c(1, 18, 0, 6, 2), 
        col=c("red", "black", "green4", "blue", "purple"))

# plot reaction label
text(0.007, 27, labels=expression(paste(NULL^"3","H(d,n)",NULL^"4","He")), cex=2.0)

# add data 
points( obsx1, obsy1, col="red", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="red" )  

points( obsx2, obsy2, col="black", pch=18, cex=1.8 )
add.error.bars(obsx2, obsy2, 0.0, errobsy2, 0.0, col="black" )  

points( obsx3, obsy3, col="green4", pch=0, cex=1.0 )
add.error.bars(obsx3, obsy3, 0.0, errobsy3, 0.0, col="green4" )  

points( obsx4, obsy4, col="blue", pch=6, cex=1.0 )
add.error.bars(obsx4, obsy4, 0.0, errobsy4, 0.0, col="blue" )  

points( obsx5, obsy5, col="purple", pch=2, cex=1.0 )
add.error.bars(obsx5, obsy5, 0.0, errobsy5, 0.0, col="purple" )  

dev.off()

######################################################################
# POSTERIORS OF RESONANCE ENERGY AND REDUCED WIDTHS
######################################################################
# first determine plot ranges
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
pdf("MCMC_Tdn_d.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,3), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)
   
# plot eigenenergy in first panel           
plot(density(samplesmat[,1]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(ylab="Probability density", line=4.5, cex.lab=2.3)
title(xlab=expression(paste(E [0], " (MeV)")), line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,1]), col=adjustcolor("blue", alpha=0.5))

# order: E0, [Er,] gf, gi ,rf, ri
 
# plot deuteron reduced width                   
plot(density(samplesmat[,4]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(gamma [d]^2, " (MeV)")), line=4.0, cex.lab=2.3, 
      yaxs='i', xaxs='i')

polygon(density(samplesmat[,4]), col=adjustcolor("blue", alpha=0.5))

# plot neutron reduced width
plot(density(samplesmat[,3]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(gamma [n]^2, " (MeV)")), line=4.0, cex.lab=2.3, 
      yaxs='i', xaxs='i')
    
polygon(density(samplesmat[,3]), col=adjustcolor("blue", alpha=0.5))
    
dev.off()

######################################################################
# POSTERIORS OF CHANNEL RADII
######################################################################
# first determine plot ranges
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
pdf("MCMC_Tdn_e.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,2), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)
   
# plot a_d in first panel           
plot(density(samplesmat[,6]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(ylab="Probability density", line=4.5, cex.lab=2.3)
title(xlab=expression(paste(a [d], " (fm)")), line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,6]), col=adjustcolor("blue", alpha=0.5))
 
# plot a_n in second panel                   
plot(density(samplesmat[,5]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(a [n], " (fm)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')

polygon(density(samplesmat[,5]), col=adjustcolor("blue", alpha=0.5))
    
dev.off()

######################################################################
# POSTERIOR OF ELECTRON SCREENING POTENTIAL
######################################################################
pdf("MCMC_Tdn_f.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.5,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(0, 100)

dens <- density(1e6*samplesmat[,7])
a98 <- quantile(1e6*samplesmat[,7], prob = 0.975)

# plot density in first panel  
# cex.axis controls tick mark label size
# xaxs and yaxs = "i" plots with exact limits         
plot(dens, main="", xlab="", ylab="",
     cex.axis=1.5, yaxs='i', xaxs='i', xlim=xLim)

title(ylab="Probability density", line=4.5, cex.lab=2.3)
title(xlab="Electron screening potential (eV)", line=4.0, cex.lab=2.3)

polygon(dens, col=adjustcolor("blue", alpha=0.5))

dev.off()

######################################################################
# POSTERIOR OF S-FACTOR AT GIVEN ENERGY
######################################################################
# prediction: plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

# plot density at xchoice
pdf("MCMC_Tdn_g.pdf")
par(mfcol=c(1,1), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)

# define x value for which we would like to predict y
xchoice <- 0.04

# declare vector with y values
fitvec <- vector(mode = "numeric", length = nsamp)   

# calculate y values from all samples at given x value; 
# bare S-factor

for(i in 1:nsamp) fitvec[i] <- SfacTdn(xchoice,samplesmat[i,1],samplesmat[i,2],
    samplesmat[i,4], samplesmat[i,3],samplesmat[i,6],samplesmat[i,5], 
    0.0)

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

plot(density(fitvec), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, yaxs='i', xaxs='i')

polygon(density(fitvec), col=adjustcolor("blue", alpha=0.5))

title(xlab=expression(paste(S-factor, "  (MeV b)")), line=4.0, cex.lab=2.0)
title(ylab="Probability density", line=4.5, cex.lab=2.0)
dev.off()

######################################################################
# DENSITIES OF S-FACTOR NORMALIZATION FACTORS
######################################################################
pdf("MCMC_Tdn_h.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(0.93, 1.05)

# plot density in first panel           
plot(density(samplesmat[,12]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, 
     xlim=xLim, yaxs='i', xaxs='i')
title(ylab="Probability density", line=3.5, cex.lab=2.3)
title(xlab="S-factor normalization", line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,8]),  
     col=adjustcolor("red", alpha=0.5))
polygon(density(samplesmat[,9]), 
     col=adjustcolor("gray", alpha=0.5))
polygon(density(samplesmat[,10]),  
     col=adjustcolor("green", alpha=0.5))
polygon(density(samplesmat[,11]),  
     col=adjustcolor("blue", alpha=0.5))
polygon(density(samplesmat[,12]),  
     col=adjustcolor("purple", alpha=0.5))

legend("topleft", inset=.01, 
   c("Jar84", "Bro87", "Kob66", "Arn53", "Con52"), 
   fill=adjustcolor(c("red", "gray", "green", "blue", "purple"), alpha=0.5), 
   horiz=FALSE, cex=1.5, box.lty=0)

dev.off()

######################################################################
# POSTERIOR EXTRINSIC S-FACTOR SCATTER
######################################################################
pdf("MCMC_Tdn_i.pdf",width=10, height=6, onefile=F)
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
par(mfcol=c(2,3), mar=c(5.5,7.0,1.0,2.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# plot #1          
plot(density(samplesmat[,13]), main="", xlab="", ylab="", xlim=c(0,1.0),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3)
title(ylab="density", line=3.0, cex.lab=2.3)
legend("topright", legend="Jar84", pch=NA, cex=1.5)
polygon(density(samplesmat[,13]), col=adjustcolor("blue", alpha=0.5))
 
 # plot #2                  
plot(density(samplesmat[,14]), main="", xlab="", ylab="", xlim=c(0,1.2),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
title(ylab="Probability", line=3.0, cex.lab=2.3)
legend("topright", legend="Bro87", pch=NA, cex=1.5)
polygon(density(samplesmat[,14]), col=adjustcolor("blue", alpha=0.5))

# plot #3                  
plot(density(samplesmat[,15]), main="", xlab="", ylab="", xlim=c(0,0.4),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
legend("topright", legend="Kob66", pch=NA, cex=1.5)
polygon(density(samplesmat[,15]), col=adjustcolor("blue", alpha=0.5))

# plot #4                  
plot(density(samplesmat[,16]), main="", xlab="", ylab="", xlim=c(0,1.2),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
legend("topright", legend="Arn53", pch=NA, cex=1.5)
polygon(density(samplesmat[,16]), col=adjustcolor("blue", alpha=0.5))

# plot #5                  
plot(density(samplesmat[,17]), main="", xlab="", ylab="", xlim=c(0,2.0),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
legend("topright", legend="Con52", pch=NA, cex=1.5)
polygon(density(samplesmat[,17]), col=adjustcolor("blue", alpha=0.5))

dev.off()

######################################################################
# CORRELATION PLOT
######################################################################
pdf("MCMC_Tdn_j.pdf", width=6, height=6, onefile=F)
pairs(~e0+gf+gi+rf+ri, col=adjustcolor("red", alpha=0.5),  
   data=samplesmat2[sample(nrow(samplesmat2), size=1000, replace=FALSE),], 
   main="Simple Scatterplot Matrix")

dev.off()


