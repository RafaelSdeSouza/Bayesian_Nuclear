######################################################################
# Author: Christian Iliadis (05/13/2019)
######################################################################
#
# MCMCsfactor_3Hedp.R
#
# PURPOSE: analyzing real 3Hedp data
#
# - uses the C++ function sfactor3Hedp(E, E0, ga, gb, ra, rb, ue),
#   which calls gsl libraries to compute Coulomb wave functions
#
#   E0: energy eigenvalue
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
#   IN BOTH THIS SCRIPT (FUNCTION SfacHedp) AND THE C++ SFACTOR FUNCTION
#   IN JAGS (sfactor3Hedpx.cc)
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
# for plotting purpose only
#
# ECM   :   center-of-mass energy
# E0    :   eigenenergy
# ga,gb :   reduced widths for incoming/outgoing channels
# ra,rb :   channel radii (fm) for incoming/outgoing channels
# ue    :   laboratory electron screening potential

SfacHedp <- function(ecm, e0, ga, gb, ra, rb, ue){
  # input masses, charges, angular momenta
  m1_i = 3.014932162 
  m2_i = 2.013553198    # masses (amu) of 3He and d
  m1_f = 4.001506094 
  m2_f = 1.007276452    # masses (amu) of p and 4He
  z1_i = 2 
  z2_i = 1              # charges of t and d
  z1_f = 2 
  z2_f = 1				          # charges of n and 4He
  jt = 0.5              # spins of target, projectile, resonance
  jp = 1.0 
  jr = 1.5               
  Q = 18.353053			      # reaction Q-value (MeV)
  la = 0 
  lb = 2				            # orbital angular momenta of d and p	

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
## DATA INPUT 
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, S-factor in MeVb

## DATA SET 1: geist
obsx1    <- c( 1.5270E-01, 1.6446E-01, 1.7622E-01, 1.8804E-01, 1.8990E-01,
               1.9986E-01, 2.0166E-01, 2.1174E-01, 2.2356E-01, 2.2548E-01,
               2.3736E-01, 2.4930E-01, 2.6124E-01, 2.7318E-01, 2.7462E-01,
               2.8512E-01, 2.8656E-01, 2.9706E-01, 2.9850E-01, 3.0894E-01,
               3.1050E-01, 3.2244E-01, 3.3438E-01, 3.4632E-01, 3.5826E-01,
               3.7020E-01, 3.8214E-01, 3.8802E-01, 1.4682E-01, 1.7016E-01,
               1.7796E-01, 2.0916E-01, 2.1696E-01, 2.4042E-01, 2.5218E-01,
               2.5422E-01, 2.7174E-01, 2.8554E-01, 3.0294E-01, 3.1686E-01,
               3.4818E-01, 3.7956E-01, 4.1082E-01 )
obsy1    <- c( 1.5887E+01, 1.6762E+01, 1.7067E+01, 1.7360E+01, 1.7418E+01,
               1.7395E+01, 1.7562E+01, 1.7379E+01, 1.6360E+01, 1.6635E+01,
               1.5864E+01, 1.5084E+01, 1.4119E+01, 1.3300E+01, 1.3399E+01,
               1.2181E+01, 1.2464E+01, 1.1520E+01, 1.1428E+01, 1.0474E+01,
               1.0598E+01, 9.6450E+00, 8.9647E+00, 8.2734E+00, 7.6396E+00,
               7.1339E+00, 6.5708E+00, 6.3086E+00, 1.4635E+01, 1.6944E+01,
               1.6635E+01, 1.7016E+01, 1.7019E+01, 1.5786E+01, 1.4903E+01,
               1.4747E+01, 1.3266E+01, 1.2676E+01, 1.0686E+01, 1.0173E+01,
               8.0123E+00, 6.3943E+00, 5.3507E+00 )
errobsy1 <- c( 1.1559E-01, 2.0681E-01, 1.3776E-01, 2.0387E-01, 6.6959E-02,
               1.5533E-01, 9.4571E-02, 2.3166E-01, 2.0661E-01, 7.0310E-02,
               1.2357E-01, 1.1842E-01, 1.0117E-01, 4.5523E-02, 5.2230E-02,
               7.6976E-02, 8.3235E-02, 7.8651E-02, 6.3897E-02, 8.9604E-02,
               5.5349E-02, 5.3435E-02, 3.3043E-02, 4.5994E-02, 4.0650E-02,
               4.4876E-02, 3.7349E-02, 3.5649E-02, 2.2239E-01, 2.2923E-01,
               2.3348E-01, 2.2260E-01, 2.2109E-01, 2.0954E-01, 8.1260E-02,
               5.9155E-02, 1.1105E-01, 1.1537E-01, 1.4289E-01, 1.0692E-01,
               8.2994E-02, 7.3778E-02, 6.0356E-02 )

## DATA SET 2: krauss
obsx2    <- c(0.04970,  0.05369,  0.0590,  0.05952,  0.05966,  0.0653,  
               0.0830,   0.1007,  0.1184,   0.1360,   0.1537,  0.1713 )
obsy2    <- c(   7.24,     7.67,     8.4,     7.96,     8.26,     8.6,   
                  9.1,      9.9,    11.5,     13.1,     16.0,    17.2 )
errobsy2 <- c(  0.072,    0.077,     0.4,    0.080,    0.083,     0.5,
                  0.4,      0.5,     0.6,      0.6,      0.8,     0.9 )

## DATA SET 3: moeller
obsx3    <- c( 8.5600E-02, 1.4600E-01, 1.7360E-01, 2.0560E-01, 2.1960E-01,
               2.3360E-01, 2.6680E-01, 2.9640E-01, 3.0880E-01, 3.2800E-01,
               3.8560E-01, 4.4600E-01 )
obsy3    <- c( 7.1730E+00, 1.4182E+01, 1.6676E+01, 1.8580E+01, 1.7239E+01,
               1.7646E+01, 1.4769E+01, 1.3223E+01, 1.1221E+01, 1.0157E+01,
               7.3423E+00, 5.0221E+00 )
errobsy3 <- c( 5.0211E-01, 4.2545E-01, 5.0027E-01, 5.5739E-01, 5.1717E-01,
               5.2937E-01, 4.4306E-01, 3.9670E-01, 3.3662E-01, 3.0472E-01,
               2.2027E-01, 1.5066E-01 )

## DATA SET 4: zhichang
obsx4    <- c( 1.3440E-01, 1.4220E-01, 1.5060E-01, 1.6500E-01, 1.6770E-01,
               1.6920E-01, 1.7460E-01, 1.8000E-01, 1.8540E-01, 1.9140E-01,
               1.9740E-01, 2.0400E-01, 2.1000E-01, 2.2470E-01, 2.4000E-01,
               2.5200E-01, 2.5590E-01, 2.6760E-01, 2.7360E-01, 2.7480E-01,
               2.9490E-01, 3.1860E-01, 3.4680E-01, 3.8220E-01, 4.1100E-01,
               4.3140E-01 )
obsy4    <- c( 1.0074E+01, 1.2080E+01, 1.2544E+01, 1.2375E+01, 1.4590E+01,
               1.5314E+01, 1.5760E+01, 1.5987E+01, 1.7415E+01, 1.8116E+01,
               1.7763E+01, 1.7728E+01, 1.7825E+01, 1.7502E+01, 1.6640E+01,
               1.5987E+01, 1.5950E+01, 1.4505E+01, 1.4194E+01, 1.4231E+01,
               1.2757E+01, 1.0505E+01, 9.1881E+00, 6.7735E+00, 5.8631E+00,
               5.0714E+00 )
errobsy4 <- c( 6.0443E-02, 7.2483E-02, 7.5261E-02, 7.4252E-02, 8.7540E-02,
               9.1881E-02, 9.4560E-02, 9.5924E-02, 1.0449E-01, 1.0870E-01,
               1.0658E-01, 1.0637E-01, 1.0695E-01, 1.0501E-01, 9.9841E-02,
               9.5923E-02, 9.5702E-02, 8.7030E-02, 8.5162E-02, 8.5388E-02,
               7.6542E-02, 6.3032E-02, 5.5129E-02, 4.0641E-02, 3.5178E-02,
               3.0429E-02 )

## DATA SET 5: costantini [inverse kinematics]
obsx5    <- c( 0.00422, 0.00459, 0.00471, 0.00500, 0.00509, 0.00537, 
               0.00544, 0.00551, 0.00554, 0.00577, 0.00583, 0.00590,
               0.00593, 0.00623, 0.00631, 0.00652, 0.00657, 0.00672,
               0.00700, 0.00707, 0.00734, 0.00740, 0.00746, 0.00752,
               0.00812, 0.00819, 0.00822, 0.00829, 0.00834, 0.00851,
               0.00860, 0.00871, 0.00898, 0.00908, 0.00914, 0.00929,
               0.00938, 0.00948, 0.00978, 0.00987, 0.00991, 0.01007,
               0.01018, 0.01029, 0.01087, 0.01096, 0.01105, 0.01165,
               0.01178, 0.01185, 0.01241, 0.01255, 0.01265, 0.01307,
               0.01383 )
obsy5    <- c(   9.700,   9.000,   8.090,   8.950,   9.820,   8.860,
                 8.010,   9.680,   8.870,   8.990,   8.930,   8.740,
                 8.310,   8.150,   8.100,   8.260,   8.030,   8.320,
                 7.900,   8.320,   7.700,   7.690,   7.900,   7.830,
                 7.310,   7.310,   7.430,   7.850,   7.510,   7.380,
                 7.650,   7.820,   7.360,   7.660,   7.600,   7.540,
                 7.470,   7.590,   7.320,   7.520,   7.240,   7.390,
                 7.350,   7.440,   7.380,   7.350,   7.410,   7.240,
                 7.210,   7.350,   7.580,   7.310,   7.310,   7.360,
                 7.230 )
errobsy5 <- c(   1.650,  0.6700,  0.6200,  0.7000,   1.010,  0.4800,
                0.3900,  0.7000,  0.5500,  0.3100,  0.3600,  0.2500,
                0.5500,  0.2900,  0.2300,  0.4800,  0.2600,  0.1600,
                0.2700,  0.2800,  0.1600,  0.1500,  0.1800,  0.2600,
                0.1100,  0.2500,  0.2400,  0.2400,  0.2200,  0.2200,
                0.2300,  0.2800,  0.1900,  0.1700,  0.2400,  0.2300,
                0.1400,  0.1300,  0.2200,  0.1900,  0.1800,  0.1900,
                0.1400,  0.1700,  0.1600,  0.1300,  0.1700,  0.1400,
                0.1200,  0.1100,  0.1400,  0.1300,  0.1400,  0.1300,
                0.1300 )

## DATA SET 6: aliotta_a
obsx6    <- c( 0.00501, 0.00550, 0.00601, 0.00602, 0.00645, 0.00690,
               0.00751, 0.00780, 0.00818, 0.00896, 0.00902, 0.00966,
               0.01072, 0.01144, 0.01199, 0.01318, 0.01439, 0.01499,
               0.01675, 0.01799, 0.01914, 0.02094, 0.02155, 0.02271,
               0.02393, 0.02508, 0.02662, 0.02873, 0.02991, 0.03110,
               0.03289, 0.03349, 0.03586, 0.03858, 0.04067, 0.04187, 
               0.04306, 0.04485, 0.04544, 0.04786, 0.05021, 0.05083,
               0.05265, 0.05382, 0.05504, 0.05681, 0.05741, 0.05980 )
obsy6    <- c(   10.76,   10.31,    9.41,   10.37,    9.36,    8.82,
                  9.12,    8.66,    8.31,    7.85,    7.88,    7.68,
                  7.48,    7.29,    7.35,    7.28,    7.04,    7.07,
                  7.16,    7.02,     7.2,    6.81,    7.09,    7.13,
                  6.91,    7.18,    6.93,    7.22,    6.87,    7.24,
                  7.04,    7.33,    7.22,    7.33,    7.44,    7.29,
                  7.48,     7.4,     7.7,    7.63,    7.66,     7.7,
                   7.7,    7.77,    7.77,    7.88,    7.94,    8.12 )
errobsy6 <- c(    1.47,    0.81,    0.48,    0.35,    0.35,    0.31,
                  0.29,    0.28,    0.24,    0.26,    0.22,    0.22,
                   0.2,     0.2,     0.2,     0.2,     0.2,     0.2,
                   0.2,     0.2,     0.2,    0.18,     0.2,     0.2,
                  0.18,     0.2,    0.18,     0.2,    0.18,     0.2,
                   0.2,     0.2,    0.18,    0.18,    0.22,     0.2,
                  0.22,     0.2,    0.22,     0.2,    0.22,     0.2,
                  0.22,    0.22,    0.22,     0.2,    0.22,    0.22 )

## DATA SET 7: aliotta_b [inverse kinematics]
obsx7    <- c( 0.01195, 0.01395, 0.01595, 0.01794, 0.01993, 0.02192,
               0.02392, 0.02592, 0.02788, 0.02988, 0.03190, 0.03389,
               0.03589, 0.03589, 0.03788, 0.03987 )
obsy7    <- c(    7.19,    6.95,    6.87,    6.91,    6.75,    6.96,
                  6.92,    7.17,    7.01,    7.19,    7.17,    7.43,
                  7.39,    7.58,    7.63,    7.65 )
errobsy7 <- c(    0.35,    0.23,    0.18,    0.19,    0.18,    0.19,
                  0.19,    0.19,    0.19,    0.19,    0.19,     0.2,
                  0.19,     0.2,     0.2,     0.2 )

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

# uea is the electron screening potential for normal kinematics
# [deuteron beam], ueb is for inverse kinematics

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
  yt1[i] <- sfactor3Hedpx(obsx1[i], e0, ga, gb, ra, rb, uea)
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
  yt2[i] <- sfactor3Hedpx(obsx2[i], e0, ga, gb, ra, rb, uea)
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
  yt3[i] <- sfactor3Hedpx(obsx3[i], e0, ga, gb, ra, rb, uea)
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
  yt4[i] <- sfactor3Hedpx(obsx4[i], e0, ga, gb, ra, rb, uea)
}    

# inverse kinematics
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
  yt5[i] <- sfactor3Hedpx(obsx5[i], e0, ga, gb, ra, rb, ueb)
}    

# inverse kinematics
for (i in 1:length(obsx6)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter, if any:
  obsy6[i] ~ dnorm(ya6[i], pow(yscat6, -2))    
  # ...subject to stat uncertainties:
  ya6[i] ~ dnorm(ym6[i], pow(errobsy6[i], -2))
  # ...subject to syst uncertainties: 
  ym6[i] <- y.norm6 * yt6[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt6[i] <- sfactor3Hedpx(obsx6[i], e0, ga, gb, ra, rb, uea)
}    

# inverse kinematics
for (i in 1:length(obsx7)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter, if any:
  obsy7[i] ~ dnorm(ya7[i], pow(yscat7, -2))    
  # ...subject to stat uncertainties:
  ya7[i] ~ dnorm(ym7[i], pow(errobsy7[i], -2))
  # ...subject to syst uncertainties: 
  ym7[i] <- y.norm7 * yt7[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt7[i] <- sfactor3Hedpx(obsx7[i], e0, ga, gb, ra, rb, ueb)
}    

###################
# PRIORS
###################
#
# e0:       eigenenergy
# ga,gb:    initial reduced width, final reduced width [ga,gb = gamma^2]
#
# Wigner limit: wl = hbar^2/(m_red a_c^2) = 41.80159/(M_red a_c^2)
#
# deuteron channel: wl_d = 41.80159/(1.207266 a_c^2) = 34.6250/a_c^2
# proton channel:   wl_p = 41.80159/(0.804710 a_c^2) = 51.9461/a_c^2
#
  wl_d <- 34.6250*pow(ra, -2)
  wl_p <- 51.9461*pow(rb, -2)

# energy eigenvalue
  e0 ~ dnorm(0.0, pow(1, -2))T(0,)        # positive since we see sigma peak 

# reduced widths
  ga ~ dnorm(0, pow(1,-2))T(0,)      
  gb ~ dnorm(0, pow(2,-2))T(0,)     
##  ga ~ dnorm(0.0, pow(wl_d, -2))T(0,)      
##  gb ~ dnorm(0.0, pow(wl_p, -2))T(0,)      

# channel radius included in sampling; 
##  ra ~ dunif(2.5, 8.0)
##  rb ~ dunif(2.5, 8.0)
###  ra ~ dnorm(3.25, pow(1.0, -2))T(0,)
###  rb ~ dnorm(5.77, pow(2.0, -2))T(0,)
  ra ~ dnorm(6, pow(0.0001,-2))T(0,)
  rb ~ dnorm(5, pow(0.0001,-2))T(0,)

# screening potential:
  uea ~ dnorm(0.0, pow(0.1, -2))T(0,)     # certainly less than 10 keV
  ueb ~ dnorm(0.0, pow(0.1, -2))T(0,)     # certainly less than 10 keV
###  uea ~ dnorm(200e-6, pow(80e-6, -2))T(0,)     # certainly less than 10 keV
####  uea ~ dnorm(180e-6, pow(80e-10, -2))T(0,)     # certainly less than 10 keV
####  ueb ~ dnorm(121e-6, pow(50e-6, -2))T(0,)     # certainly less than 10 keV

# extrinsic scatter
  yscat1 ~ dnorm(0.0, pow(1, -2))T(0,)     # sigma is 5 MeVb
  yscat2 ~ dnorm(0.0, pow(1, -2))T(0,)
  yscat3 ~ dnorm(0.0, pow(1, -2))T(0,)
  yscat4 ~ dnorm(0.0, pow(1, -2))T(0,)
  yscat5 ~ dnorm(0.0, pow(1, -2))T(0,)
  yscat6 ~ dnorm(0.0, pow(1, -2))T(0,)
  yscat7 ~ dnorm(0.0, pow(1, -2))T(0,)

# systematic normalization factor for S-factor:
# log(): natural logarithm
y.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
logmu1 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma1 <- log(1.043)  # factor uncertainty is 1.043, i.e. 4.3% for Gei99

y.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
logmu2 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma2 <- log(1.078)  # factor uncertainty is 1.078, i.e., 7.8% for Kra87

y.norm3 ~ dlnorm(logmu3, pow(logsigma3, -2))
logmu3 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma3 <- log(1.039)  # factor uncertainty is 1.039, i.e., 3.9% for Moe80

y.norm4 ~ dlnorm(logmu4, pow(logsigma4, -2))
logmu4 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma4 <- log(1.034)  # factor uncertainty is 1.034, i.e., 3.4% for Zhi77 

y.norm5 ~ dlnorm(logmu5, pow(logsigma5, -2))
logmu5 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma5 <- log(1.055)  # factor uncertainty is 1.055, i.e., 5.5% for Cos00 

y.norm6 ~ dlnorm(logmu6, pow(logsigma6, -2))
logmu6 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma6 <- log(1.030)  # factor uncertainty is 1.030, i.e., 3.0% for Ali01_a 

y.norm7 ~ dlnorm(logmu7, pow(logsigma7, -2))
logmu7 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma7 <- log(1.030)  # factor uncertainty is 1.030, i.e., 3.0% for Ali01_b 

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
n.adapt = 1500 
n.burn = 1500     
n.iter = 1500  
thin = 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model
ourmodel <- jags.model(f, data = list(          ## jags wants all data in a list
                'obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1,
                'obsx2' = obsx2, 'obsy2' = obsy2, 'errobsy2' = errobsy2,
                'obsx3' = obsx3, 'obsy3' = obsy3, 'errobsy3' = errobsy3,
                'obsx4' = obsx4, 'obsy4' = obsy4, 'errobsy4' = errobsy4,
                'obsx5' = obsx5, 'obsy5' = obsy5, 'errobsy5' = errobsy5,
                'obsx6' = obsx6, 'obsy6' = obsy6, 'errobsy6' = errobsy6,
                'obsx7' = obsx7, 'obsy7' = obsy7, 'errobsy7' = errobsy7
                                    ),
#                inits = list(e0 = 0.236, ga = 0.36, gb = 0.035, ra = 3.25, rb = 5.77, 
#                             uea=2e-4, ueb=1e-4), 
                n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn)

# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                variable.names=c(
                'e0', 'ga', 'gb', 
                'ra', 'rb', 
                'uea', 'ueb',
                'y.norm1', 'y.norm2', 'y.norm3', 'y.norm4', 'y.norm5', 
                'y.norm6', 'y.norm7',
                'yscat1',  'yscat2',   'yscat3',  'yscat4', 'yscat5',  
                'yscat6', 'yscat7'                    
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
capture.output(print(samplesmat, print.gap=3), file="MCMCresults3Hedp")

# select first six colums [e0, ga, gb, ra, rb] from matrix
samplesmat2 <- samplesmat[,1:5]
# output samples to file for rate calculation
capture.output( print(samplesmat2[sample(nrow(samplesmat2), size=1000, 
                  replace=FALSE),], 
            print.gap=3), file="MCMCsamples3Hedp" )

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
pdf("MCMC_3Hedp_a.pdf")
plot(mcmcChain)
dev.off()

######################################################################
# S-FACTOR FIT [red] + DATA [black & white]
######################################################################
pdf("MCMC_3Hedp_b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(3e-3,0.9)
yLim = c(0,20)

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
legend(0.4, 19.0, legend=c("Gei99", "Kra87", "Moe80", "Zhi77", "Cos00", 
        "Ali01a", "Ali01b"), pch=c(1, 5, 0, 6, 3, 2, 4),
        col=c("gray40", "gray40", "gray40", "gray40", "gray40", "gray40", 
        "gray40"))

# plot reaction label
text(0.015, 17, labels=expression(paste(NULL^"3","He(d,p)",NULL^"4","He")), cex=2.0)

text(0.0045, 5, labels=expression(paste("bare")), cex=1.3)

# plot credible S-factors:
# calculate for these energies and the set of Bayesian samples for
# Er, g^2_i, g^2_f the S_factor curve using R code above

# define grid of energy values in MeV on a log scale for plotting of credible 
# solutions; lseq is appropriate for log scale
x1 = lseq(0.001, 0.5, length=300)

# ECM, E0, ga, gb, ra, rb, uea, ueb

# 500 samples are randomly selected for plotting of S-factor including
# electron screening
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1,SfacHedp(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,3],
    samplesmat[i,4],samplesmat[i,5],samplesmat[i,6]),
    col=adjustcolor("royalblue1", alpha=0.02), lw=0.1)
}

for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1,SfacHedp(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,3],
    samplesmat[i,4],samplesmat[i,5], samplesmat[i,7]),
    col=adjustcolor("black", alpha=0.02), lw=0.1)
}

for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1,SfacHedp(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,3],
    samplesmat[i,4],samplesmat[i,5], 0.0),
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

points( obsx5, obsy5, col="gray40", pch=3, cex=1.0 )
add.error.bars(obsx5, obsy5, 0.0, errobsy5, 0.0, col="gray40" )  

points( obsx6, obsy6, col="gray40", pch=2, cex=1.0 )
add.error.bars(obsx6, obsy6, 0.0, errobsy6, 0.0, col="gray40" )  

points( obsx7, obsy7, col="gray40", pch=4, cex=1.0 )
add.error.bars(obsx7, obsy7, 0.0, errobsy7, 0.0, col="gray40" )  

dev.off()

######################################################################
# S-FACTOR DATA [color] ONLY
######################################################################
pdf("MCMC_3Hedp_c.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(3e-3,0.9)
yLim = c(0,20)

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
legend(0.4, 19.0, legend=c("Gei99", "Kra87", "Moe80", "Zhi77", "Cos00", 
        "Ali01a", "Ali01b"), pch=c(1, 5, 0, 6, 3, 2, 4),
        col=c("red", "black", "green4", "blue", "purple", "orange", 
        "gray40"))

# plot reaction label
text(0.015, 17, labels=expression(paste(NULL^"3","He(d,p)",NULL^"4","He")), cex=2.0)

# add data 
points( obsx1, obsy1, col="red", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="red" )  

points( obsx2, obsy2, col="black", pch=5, cex=1.0 )
add.error.bars(obsx2, obsy2, 0.0, errobsy2, 0.0, col="black" )  

points( obsx3, obsy3, col="green4", pch=0, cex=1.0 )
add.error.bars(obsx3, obsy3, 0.0, errobsy3, 0.0, col="green4" )  

points( obsx4, obsy4, col="blue", pch=6, cex=1.0 )
add.error.bars(obsx4, obsy4, 0.0, errobsy4, 0.0, col="blue" )  

points( obsx5, obsy5, col="purple", pch=3, cex=1.0 )
add.error.bars(obsx5, obsy5, 0.0, errobsy5, 0.0, col="purple" )  

points( obsx6, obsy6, col="orange", pch=2, cex=1.0 )
add.error.bars(obsx6, obsy6, 0.0, errobsy6, 0.0, col="orange" )  

points( obsx7, obsy7, col="gray40", pch=4, cex=1.0 )
add.error.bars(obsx7, obsy7, 0.0, errobsy7, 0.0, col="gray40" )  

dev.off()

######################################################################
# POSTERIORS OF RESONANCE ENERGY AND REDUCED WIDTHS
######################################################################
# first determine plot ranges
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
pdf("MCMC_3Hedp_d.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,3), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)
   
# plot eigenenergy in first panel           
plot(density(samplesmat[,1]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(ylab="Probability density", line=4.5, cex.lab=2.3)
title(xlab=expression(paste(E [0], " (MeV)")), line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,1]), col=adjustcolor("blue", alpha=0.5))

# order: E0, ga, gb ,ra, rb
 
# plot deuteron reduced width                   
plot(density(samplesmat[,2]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(gamma [d]^2, " (MeV)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')

polygon(density(samplesmat[,2]), col=adjustcolor("blue", alpha=0.5))

# plot proton reduced width
plot(density(samplesmat[,3]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(gamma [p]^2, " (MeV)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')
    
polygon(density(samplesmat[,3]), col=adjustcolor("blue", alpha=0.5))
    
dev.off()

######################################################################
# POSTERIORS OF CHANNEL RADII
######################################################################
# first determine plot ranges
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
pdf("MCMC_3Hedp_e.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,2), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)
   
# plot a_d in first panel           
plot(density(samplesmat[,4]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(ylab="Probability density", line=4.5, cex.lab=2.3)
title(xlab=expression(paste(a [d], " (fm)")), line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,4]), col=adjustcolor("blue", alpha=0.5))
 
# plot a_p in second panel                   
plot(density(samplesmat[,5]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(a [p], " (fm)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')

polygon(density(samplesmat[,5]), col=adjustcolor("blue", alpha=0.5))
    
dev.off()

######################################################################
# POSTERIOR OF ELECTRON SCREENING POTENTIAL
######################################################################
pdf("MCMC_3Hedp_f.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,2), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)

xLim = c(0, 2000)

dens_a <- density(1e6*samplesmat[,6])
a98_a <- quantile(1e6*samplesmat[,6], prob = 0.975)

dens_b <- density(1e6*samplesmat[,7])
a98_b <- quantile(1e6*samplesmat[,7], prob = 0.975)

# plot density in first panel  
# cex.axis controls tick mark label size
# xaxs and yaxs = "i" plots with exact limits         

## plot 1:
plot(dens_a, main="", xlab="", ylab="",
     cex.axis=1.5, yaxs='i', xaxs='i', xlim=xLim)

#text(400, 0.02, labels=expression(paste(NULL^"3","He(d,p)",NULL^"4","He")), 
#   cex=2.0)
legend("topleft", inset=.01, 
   c(expression(paste(NULL^"3","He(d,p)",NULL^"4","He"))), 
   horiz=FALSE, cex=1.5, box.lty=0)
title(ylab="Probability density", line=4.5, cex.lab=2.3)
title(xlab=expression(paste("U"[e], " (eV)")), line=4.0, cex.lab=2.3)

polygon(dens_a, col=adjustcolor("blue", alpha=0.5))

## plot 2:
plot(dens_b, main="", xlab="", ylab="",
     cex.axis=1.5, yaxs='i', xaxs='i', xlim=xLim)

legend("topleft", inset=.01, 
   c(expression(paste("d(",NULL^"3","He,p)",NULL^"4","He"))), 
   horiz=FALSE, cex=1.5, box.lty=0)
title(xlab=expression(paste("U"[e], " (eV)")), line=4.0, cex.lab=2.3)

polygon(dens_b, col=adjustcolor("blue", alpha=0.5))

dev.off()

######################################################################
# POSTERIOR OF S-FACTOR AT GIVEN ENERGY
######################################################################
# prediction: plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

# plot density at xchoice
pdf("MCMC_3Hedp_g.pdf")
par(mfcol=c(1,1), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)

# define x value for which we would like to predict y
xchoice <- 0.04

# declare vector with y values
fitvec <- vector(mode = "numeric", length = nsamp)   

# calculate y values from all samples at given x value; 
# bare S-factor

for(i in 1:nsamp) fitvec[i] <- SfacHedp(xchoice, samplesmat[i,1], samplesmat[i,2],
    samplesmat[i,3], samplesmat[i,4],samplesmat[i,5], 0.0)

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

title(xlab=expression(paste(S [0], "  (MeV b)")), line=4.0, cex.lab=2.0)
title(ylab="Probability density", line=4.5, cex.lab=2.0)
dev.off()

######################################################################
# DENSITIES OF S-FACTOR NORMALIZATION FACTORS
######################################################################
pdf("MCMC_3Hedp_h.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(0.90, 1.1)

# plot density in first panel           
plot(density(samplesmat[,10]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, 
     xlim=xLim, yaxs='i', xaxs='i')
title(ylab="Probability density", line=3.5, cex.lab=2.3)
title(xlab="S-factor normalization", line=4.0, cex.lab=2.3)

# order: "red", "black", "green4", "blue", "purple", "orange", "gray40"
polygon(density(samplesmat[,8]),  
     col=adjustcolor("red", alpha=0.5))
polygon(density(samplesmat[,9]), 
     col=adjustcolor("black", alpha=0.5))
polygon(density(samplesmat[,10]),  
     col=adjustcolor("green4", alpha=0.5))
polygon(density(samplesmat[,11]),  
     col=adjustcolor("blue", alpha=0.5))
polygon(density(samplesmat[,12]),  
     col=adjustcolor("purple", alpha=0.5))
polygon(density(samplesmat[,13]),  
     col=adjustcolor("orange", alpha=0.5))
polygon(density(samplesmat[,14]),  
     col=adjustcolor("yellow", alpha=0.5))

legend("topleft", inset=.01, 
   c("Gei99", "Kra87", "Moe80", "Zhi77", "Cos00", "Ali01a", "Ali01b"), 
   fill=adjustcolor(c("red", "black", "green4", "blue", "purple", "orange", 
       "yellow"), alpha=0.5), 
   horiz=FALSE, cex=1.5, box.lty=0)

dev.off()

######################################################################
# POSTERIOR EXTRINSIC S-FACTOR SCATTER
######################################################################
pdf("MCMC_3Hedp_i.pdf",width=10, height=6, onefile=F)
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
par(mfcol=c(2,4), mar=c(5.5,4.0,1.0,2.0), oma=c(0.5,7.0,0.5,1.0), tck=0.02, 
     las=1)

# plot #1          
plot(density(samplesmat[,15]), main="", xlab="", ylab="", xlim=c(0,1.0),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3)
title(ylab="density", line=3.0, cex.lab=2.3)
legend("topright", legend="Gei99", pch=NA, cex=1.5)
polygon(density(samplesmat[,15]), col=adjustcolor("blue", alpha=0.5))
 
 # plot #2                  
plot(density(samplesmat[,16]), main="", xlab="", ylab="", xlim=c(0,1.2),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
title(ylab="Probability", line=3.0, cex.lab=2.3)
legend("topright", legend="Kra87", pch=NA, cex=1.5)
polygon(density(samplesmat[,16]), col=adjustcolor("blue", alpha=0.5))

# plot #3                  
plot(density(samplesmat[,17]), main="", xlab="", ylab="", xlim=c(0,0.4),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
legend("topright", legend="Moe80", pch=NA, cex=1.5)
polygon(density(samplesmat[,17]), col=adjustcolor("blue", alpha=0.5))

# plot #4                  
plot(density(samplesmat[,18]), main="", xlab="", ylab="", xlim=c(0,1.2),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
legend("topright", legend="Zhi77", pch=NA, cex=1.5)
polygon(density(samplesmat[,18]), col=adjustcolor("blue", alpha=0.5))

# plot #5                  
plot(density(samplesmat[,19]), main="", xlab="", ylab="", xlim=c(0,2.0),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
legend("topright", legend="Cos00", pch=NA, cex=1.5)
polygon(density(samplesmat[,19]), col=adjustcolor("blue", alpha=0.5))

# plot #6                  
plot(density(samplesmat[,20]), main="", xlab="", ylab="", xlim=c(0,2.0),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
legend("topright", legend="Ali01a", pch=NA, cex=1.5)
polygon(density(samplesmat[,20]), col=adjustcolor("blue", alpha=0.5))

# plot #7                  
plot(density(samplesmat[,21]), main="", xlab="", ylab="", xlim=c(0,2.0),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, yaxs='i', xaxs='i')
legend("topright", legend="Ali01b", pch=NA, cex=1.5)
polygon(density(samplesmat[,21]), col=adjustcolor("blue", alpha=0.5))

dev.off()

######################################################################
# CORRELATION PLOT
######################################################################
pdf("MCMC_3Hedp_j.pdf", width=6, height=6, onefile=F)
pairs(~e0+ga+gb+ra+rb, col=adjustcolor("red", alpha=0.5),  
   data=samplesmat2[sample(nrow(samplesmat2), size=1000, replace=FALSE),], 
   main="Simple Scatterplot Matrix")

dev.off()

