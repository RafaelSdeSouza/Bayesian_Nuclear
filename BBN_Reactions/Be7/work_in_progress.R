######################################################################
# Author: Christian Iliadis (07/13/2019)
######################################################################
#
# this script splits each data sets into two parts: a "relative cross  
# section" part, and an "absolute cross section" part; the relative 
# part is modeled with a normalization factor [e.g., y.norm1] for which 
# the prior is sampled from a very broad density; the absolute part 
# consists of just one data point. It does have a statistical uncertainty. 
# It is modeled with anormalization factor [e.g., y.norm11] for which 
# the prior is sampled from a highly informative lognormal density, where 
# the lognormal parameters mu and sigma a chosen according to the 
# systematic uncertainty reported for a given experiment;
#
# In the plot of the reduced cross section, the red data points show
# the "absolute" cross section of each data set
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
# compute cross section using gsl Coulomb wave functions;
# for plotting purposes only
#
# ECM   :   center-of-mass energy
# E0    :   eigenenergy
# ga,gb :   reduced widths for incoming/outgoing channels
# ra,rb :   channel radii (fm) for incoming/outgoing channels

Sigma7Benp <- function(ecm, e0, ga, gb, ra, rb, jr, la, lb){
  # input masses, charges, angular momenta
  m1_i = 7.0147344 
  m2_i = 1.0086649158   # masses (amu) of 7Be and n
  m1_f = 7.014357697 
  m2_f = 1.007276452    # masses (amu) of 7Li and p
  z1_i = 4 
  z2_i = 0              # charges of 7Be and n
  z1_f = 3 
  z2_f = 1              # charges of 7Li and p
  jt = 1.5              # spins of target, projectile
  jp = 0.5 
  Q = 1.6442402         # reaction Q-value (MeV)

  # reduced masses
  mue_i <- (m1_i*m2_i)/(m1_i+m2_i)
  mue_f <- (m1_f*m2_f)/(m1_f+m2_f)

  # constants
  pek <- 6.56618216e-1/mue_i
  omega <- (2*jr+1)/((2*jt+1)*(2*jp+1))

  ### CALCULATE S-FACTOR
  ## incoming channel   
  eta_a=0.1574854*z2_i*z1_i*sqrt(mue_i)
  rho_a=0.218735*ra*sqrt(mue_i)      
  eta_i=eta_a/(sqrt(ecm))
  rho_i=rho_a*(sqrt(ecm))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k=0)
  # penetration and shift factor 
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  s_i <- rho_i*(P3$val_F*P3$val_Fp + P3$val_G*P3$val_Gp)/(P3$val_F^2 + P3$val_G^2)
  # shift factor at energy E0 [eigenvalue]
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
  
  s1=pek*omega*Ga*Gb
  s2=((e0-ecm-tapp)^2)+0.25*((Ga+Gb)^2)
  SF <- (s1/s2)*(1/ecm)

  return(SF = SF)
}

######################################################################
## DATA INPUT: ENERGY, sqrt(E) * SIGMA 
## [statistical uncertainties only!]
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, sqrt(Ecm)*sigma is in sqrt(MeV)b; 

### RELATIVE CROSS SECTIONS:

## DATA SET 1: damone 2018;
##             includes 7Li ground and first excited state
obsx1    <- c( 1.7909E-08, 2.8383E-08, 4.4985E-08, 7.1296E-08, 1.1300E-07,
               1.7909E-07, 2.8383E-07, 4.4985E-07, 7.1296E-07, 1.1300E-06,
               1.7909E-06, 2.8383E-06, 4.4985E-06, 7.1296E-06, 1.1300E-05,
               1.7909E-05, 2.8383E-05, 4.4985E-05, 7.1296E-05, 1.1300E-04,
               1.7909E-04, 2.8383E-04, 4.4985E-04, 7.1296E-04, 1.1300E-03,
               1.7909E-03, 2.8383E-03, 4.4985E-03, 7.1296E-03, 1.1300E-02,
               1.7909E-02, 2.8383E-02 )
obsy1    <- c( 7.8092E+00, 7.8195E+00, 7.7802E+00, 7.7661E+00, 7.6519E+00,
               7.5714E+00, 7.5733E+00, 7.6167E+00, 7.5710E+00, 7.4976E+00,
               7.5745E+00, 7.5757E+00, 7.6668E+00, 7.5095E+00, 7.7362E+00,
               7.4136E+00, 7.7286E+00, 7.6756E+00, 7.5196E+00, 7.8341E+00,
               7.5144E+00, 7.2445E+00, 6.9702E+00, 6.8301E+00, 7.0019E+00,
               5.6533E+00, 5.3962E+00, 4.9763E+00, 5.1530E+00, 3.5495E+00,
               3.4931E+00, 3.8298E+00 )
errobsy1 <- c( 2.7114E-02, 1.6001E-02, 1.5102E-02, 1.7252E-02, 2.4333E-02,
               3.6717E-02, 4.5931E-02, 5.2407E-02, 5.9494E-02, 6.6445E-02,
               7.4939E-02, 8.3775E-02, 9.3976E-02, 1.0263E-01, 1.1767E-01,
               1.2603E-01, 1.4380E-01, 1.5864E-01, 1.7379E-01, 1.9866E-01,
               2.1469E-01, 2.5093E-01, 2.5641E-01, 2.7691E-01, 3.1250E-01,
               2.9929E-01, 3.3142E-01, 3.3659E-01, 3.9149E-01, 3.2201E-01,
               3.2915E-01, 4.1614E-01 )

## DATA SET 2: gibbons & macklin 1959; data with Eplab>2371 keV excluded;
##             also excluded data below Encm=0.01 MeV; use relativistic 
##             results;
##             includes only 7Li ground state
obsx2    <- c( 1.331079E-02, 1.855621E-02, 2.467585E-02,
               3.516666E-02, 4.478322E-02, 5.614824E-02, 6.663900E-02,
               7.800398E-02, 8.762049E-02, 9.374008E-02, 1.051050E-01,
               1.173442E-01, 1.278348E-01, 1.365771E-01, 1.496904E-01,
               1.636779E-01, 1.750427E-01, 1.942755E-01, 2.152566E-01,
               2.327408E-01, 2.467282E-01, 2.598413E-01, 2.720802E-01,
               2.729544E-01, 2.825707E-01, 2.948095E-01, 3.079225E-01,
               3.201614E-01, 3.324001E-01, 3.463873E-01, 3.551293E-01,
               3.708648E-01, 3.831035E-01, 3.970906E-01, 4.206937E-01 )
obsy2    <- c( 3.403512E+00, 3.025938E+00, 2.718642E+00,
               2.372191E+00, 2.154046E+00, 1.929502E+00, 1.795300E+00,
               1.670390E+00, 1.579014E+00, 1.537671E+00, 1.445473E+00,
               1.393007E+00, 1.342572E+00, 1.348650E+00, 1.302343E+00,
               1.286430E+00, 1.282212E+00, 1.300835E+00, 1.370100E+00,
               1.470022E+00, 1.480242E+00, 1.616691E+00, 1.663510E+00,
               1.698281E+00, 1.815184E+00, 1.906191E+00, 1.965685E+00,
               2.009236E+00, 1.967049E+00, 1.927036E+00, 1.787446E+00,
               1.703345E+00, 1.581327E+00, 1.505659E+00, 1.320614E+00 )              
errobsy2 <- c( 3.403512E-02, 3.025938E-02, 2.718642E-02,
               2.372191E-02, 2.154046E-02, 1.929502E-02, 1.795300E-02,
               1.670390E-02, 1.579014E-02, 1.537671E-02, 1.445473E-02,
               1.393007E-02, 1.342572E-02, 1.348650E-02, 1.302343E-02,
               1.286430E-02, 1.282212E-02, 1.300835E-02, 1.370100E-02,
               1.470022E-02, 1.480242E-02, 1.616691E-02, 1.663510E-02,
               1.698281E-02, 1.815184E-02, 1.906191E-02, 1.965685E-02,
               2.009236E-02, 1.967049E-02, 1.927036E-02, 1.787446E-02,
               1.703345E-02, 1.581327E-02, 1.505659E-02, 1.320614E-02 )

## DATA SET 3: martin-hernandez 2018; absolute pn cross directly obtained 
##             from author; converted to np using relativistic kinematics; 
##             data below Encm=0.002 MeV are excluded because uncertainty 
##             in 7Be mass impacts results;
##             includes only 7Li ground state
obsx3    <- c( 2.183971E-03, 2.385885E-03, 2.592686E-03,
               2.804365E-03, 3.020310E-03, 3.240364E-03, 3.463872E-03,
               3.690902E-03, 3.921238E-03, 4.154371E-03, 4.390563E-03,
               4.628933E-03, 4.869689E-03, 5.112692E-03, 5.357338E-03,
               5.603375E-03, 5.851169E-03, 6.100352E-03, 6.350480E-03,
               6.601578E-03, 6.853464E-03, 7.105873E-03, 7.359096E-03,
               7.612747E-03, 7.866494E-03, 8.120285E-03, 8.374636E-03,
               8.629143E-03, 8.882742E-03, 9.389572E-03  )
obsy3    <- c( 5.595375E+00, 5.456809E+00, 5.305795E+00,
               5.207356E+00, 5.052560E+00, 5.175160E+00, 4.925496E+00,
               4.805613E+00, 4.939230E+00, 5.024024E+00, 4.971533E+00,
               5.005326E+00, 4.567042E+00, 4.611864E+00, 4.691854E+00,
               4.420279E+00, 4.222387E+00, 4.694943E+00, 4.354219E+00,
               4.726389E+00, 4.137051E+00, 4.456333E+00, 4.135836E+00,
               4.159703E+00, 3.942435E+00, 4.362248E+00, 4.143442E+00,
               4.424349E+00, 4.026532E+00, 3.879684E+00  )
errobsy3 <- c( 8.382935E-02, 8.384005E-02, 8.433243E-02,
               8.552457E-02, 8.472940E-02, 9.031378E-02, 9.017204E-02,
               9.020569E-02, 9.367794E-02, 9.723094E-02, 9.851541E-02,
               1.041874E-01, 9.891482E-02, 1.026130E-01, 1.093601E-01,
               1.108938E-01, 1.076231E-01, 1.255085E-01, 1.199439E-01,
               1.345929E-01, 1.259839E-01, 1.447373E-01, 1.365988E-01,
               1.435971E-01, 1.551009E-01, 1.652441E-01, 1.539148E-01,
               1.872132E-01, 1.951334E-01, 1.964327E-01  )

## DATA SET 4: koehler 1988; for energies below 8e-5 MeV [constant reduced
##             cross section] data have been rebinned by adding 5 data points
##             and calculating weighted mean;
##             includes 7Li ground and first excited state
obsx4    <- c( 2.608865E-08, 3.456921E-08, 4.784085E-08, 6.119992E-08,
               7.321259E-08, 8.924697E-08, 1.108593E-07, 1.416341E-07,
               1.872718E-07, 2.552911E-07, 3.168407E-07, 3.888818E-07,
               4.887251E-07, 6.331569E-07, 8.503292E-07, 1.206513E-06,
               1.550981E-06, 1.734581E-06, 1.949655E-06, 2.208443E-06,
               2.521437E-06, 2.907871E-06, 3.390476E-06, 4.002475E-06,
               4.794577E-06, 5.854210E-06, 7.303774E-06, 9.377577E-06,
               1.243233E-05, 1.533495E-05, 1.829003E-05, 2.224180E-05,
               2.757494E-05, 3.509378E-05, 4.619720E-05, 6.356049E-05,               
               9.279657E-05, 1.197770E-04, 1.320170E-04, 1.460055E-04,
               1.626169E-04, 1.827255E-04, 2.054569E-04, 2.343083E-04,
               2.684054E-04, 3.103710E-04, 3.637024E-04, 4.327709E-04,
               5.219479E-04, 6.425992E-04, 8.104618E-04, 1.057884E-03,
               1.425084E-03, 2.037083E-03, 3.138682E-03, 5.473022E-03,
               1.180284E-02 )
obsy4    <- c( 5.882701E+00, 5.757024E+00, 5.684058E+00, 5.719731E+00,
               5.695627E+00, 5.733079E+00, 5.664370E+00, 5.627097E+00,
               5.674665E+00, 5.712598E+00, 5.710666E+00, 5.800993E+00,
               5.654563E+00, 5.751555E+00, 5.819887E+00, 5.759546E+00,
               5.691495E+00, 5.839202E+00, 5.656230E+00, 5.707883E+00,
               5.675726E+00, 5.622623E+00, 5.676084E+00, 5.655390E+00,
               5.661571E+00, 5.659802E+00, 5.656776E+00, 5.761900E+00,
               5.723593E+00, 5.629112E+00, 5.681136E+00, 5.753711E+00,
               5.692271E+00, 5.612388E+00, 5.706071E+00, 5.671267E+00,               
               5.432789E+00, 5.395524E+00, 5.584074E+00, 5.509973E+00,
               5.355897E+00, 5.407040E+00, 5.461165E+00, 5.342189E+00,
               5.488332E+00, 5.267588E+00, 5.244520E+00, 5.200786E+00,
               5.300314E+00, 5.145957E+00, 4.925070E+00, 4.895033E+00,
               4.794286E+00, 4.508892E+00, 4.123362E+00, 3.713790E+00,
               3.139722E+00 )
errobsy4 <- c( 1.330698E-01, 9.837351E-02, 1.169144E-01, 5.394071E-02,               
               5.390683E-02, 5.417379E-02, 5.350609E-02, 5.307900E-02,
               5.327827E-02, 5.409253E-02, 5.380607E-02, 5.454583E-02,
               5.303551E-02, 5.467026E-02, 5.500018E-02, 5.470012E-02,
               5.386171E-02, 5.568360E-02, 5.435516E-02, 5.477597E-02,
               5.647405E-02, 7.140915E-02, 7.332296E-02, 7.231835E-02,
               7.209151E-02, 7.128335E-02, 7.257372E-02, 8.851526E-02,
               9.468445E-02, 1.589130E-01, 1.632504E-01, 1.682113E-01,
               1.643368E-01, 1.583425E-01, 1.784229E-01, 1.801364E-01,               
               1.570300E-01, 3.611608E-01, 3.102263E-01, 2.899986E-01,
               2.932991E-01, 3.109048E-01, 3.153429E-01, 3.061426E-01,
               3.604278E-01, 3.347297E-01, 3.432777E-01, 3.328503E-01,
               3.426927E-01, 3.041945E-01, 2.846861E-01, 2.732111E-01,
               2.491519E-01, 2.166435E-01, 1.960838E-01, 1.553577E-01,
               2.172818E-01 )

### ABSOLUTE NORMALIZATIONS:

## DATA SET 10: koehler 1988; thermal cross section
obsx10    <- c( 2.21e-8 )
obsy10    <- c( 5.708 )
errobsy10 <- c( 0.057 ) 

## DATA SET 11: damone 2018; thermal cross section
obsx11    <- c( 2.21e-8 )
obsy11    <- c( 7.775 )
errobsy11 <- c( 0.078 ) 

## DATA SET 12: gibbons 1959; lowest-energy data point
obsx12    <- c( 9.813847E-03 )
obsy12    <- c( 3.688387E+00 )
errobsy12 <- c( 3.688387E-02 ) 

## DATA SET 13: hernandez 2019; lowest-energy data point
obsx13    <- c( 1.987364E-03 )
obsy13    <- c( 5.710998E+00 )
errobsy13 <- c( 8.205219E-02 ) 

## DATA SET 14: cervena 1989; thermal cross section
obsx14    <- c( 2.21e-8 )
obsy14    <- c( 6.818 )
errobsy14 <- c( 0.068 ) 

## DATA SET 15: tomandl 2019; thermal cross section
obsx15    <- c( 2.21e-8 )
obsy15    <- c( 6.482 )
errobsy15 <- c( 0.091 ) 

######################################################################                  
# JAGS MODEL
######################################################################                 
# rjags ----->
cat('model {

######################################
#
# LIKELIHOODS
#
######################################
# - careful: dnorm is differently defined in R and JAGS! 
# - precision=sigma^(-2)
# - in a for loop, make sure **all** variables on the LEFT of an 
#   expression has the index [i]
# - systematic error as normalization factor y.norm...

## RELATIVE DATA

for (i in 1:length(obsx1)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy1[i] ~ dnorm(ya1[i], pow(errobsy1[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya1[i] ~ dnorm(ym1[i], pow(yscat1, -2))
  # ...subject to syst uncertainties: 
  ym1[i] <- yt1[i] * 10^y.norm1
  # true sigma [calculated from theory]: 
  yt1[i] <- sqrt(obsx1[i]) * 
            (
            sigma7Benpx(obsx1[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx1[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx1[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx1[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx1[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx1[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx1[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}    

for (i in 1:length(obsx2)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy2[i] ~ dnorm(ya2[i], pow(errobsy2[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya2[i] ~ dnorm(ym2[i], pow(yscat2, -2))
  # ...subject to syst uncertainties: 
  ym2[i] <- yt2[i] * 10^y.norm2
  # true sigma [calculated from theory]: 
  yt2[i] <- sqrt(obsx2[i]) * 
            (
            sigma7Benpx(obsx2[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx2[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx2[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx2[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx2[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx2[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx2[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}    

for (i in 1:length(obsx3)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy3[i] ~ dnorm(ya3[i], pow(errobsy3[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya3[i] ~ dnorm(ym3[i], pow(yscat3, -2))
  # ...subject to syst uncertainties: 
  ym3[i] <- yt3[i] * 10^y.norm3
  # true sigma [calculated from theory]: 
  yt3[i] <- sqrt(obsx3[i]) * 
            (
            sigma7Benpx(obsx3[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx3[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx3[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx3[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx3[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx3[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx3[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}    

for (i in 1:length(obsx4)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy4[i] ~ dnorm(ya4[i], pow(errobsy4[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya4[i] ~ dnorm(ym4[i], pow(yscat4, -2))
  # ...subject to syst uncertainties: 
  ym4[i] <- yt4[i] * 10^y.norm4
  # true sigma [calculated from theory]: 
  yt4[i] <- sqrt(obsx4[i]) * 
            (
            sigma7Benpx(obsx4[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx4[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx4[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx4[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx4[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx4[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx4[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}    

## ABSOLUTE DATA [no "extrinsic scatter" needed since there is only one data point]

for (i in 1:length(obsx10)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy10[i] ~ dnorm(ym10[i], pow(errobsy10[i], -2))    
  # ...subject to syst uncertainties: 
  ym10[i] <- yt10[i] * y.norm10
  # true sigma [calculated from theory]: 
  yt10[i] <- sqrt(obsx10[i]) * 
            (
            sigma7Benpx(obsx10[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx10[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx10[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx10[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx10[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx10[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx10[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx11)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy11[i] ~ dnorm(ym11[i], pow(errobsy11[i], -2))    
  # ...subject to syst uncertainties: 
  ym11[i] <- yt11[i] * y.norm11
  # true sigma [calculated from theory]: 
  yt11[i] <- sqrt(obsx11[i]) * 
            (
            sigma7Benpx(obsx11[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx11[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx11[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx11[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx11[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx11[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx11[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx12)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy12[i] ~ dnorm(ym12[i], pow(errobsy12[i], -2))    
  # ...subject to syst uncertainties: 
  ym12[i] <- yt12[i] * y.norm12
  # true sigma [calculated from theory]: 
  yt12[i] <- sqrt(obsx12[i]) * 
            (
            sigma7Benpx(obsx12[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx12[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx12[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx12[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx12[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx12[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx12[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx13)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy13[i] ~ dnorm(ym13[i], pow(errobsy13[i], -2))    
  # ...subject to syst uncertainties: 
  ym13[i] <- yt13[i] * y.norm13
  # true sigma [calculated from theory]: 
  yt13[i] <- sqrt(obsx13[i]) * 
            (
            sigma7Benpx(obsx13[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx13[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx13[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx13[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx13[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx13[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx13[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx14)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy14[i] ~ dnorm(ym14[i], pow(errobsy14[i], -2))    
  # ...subject to syst uncertainties: 
  ym14[i] <- yt14[i] * y.norm14
  # true sigma [calculated from theory]: 
  yt14[i] <- sqrt(obsx14[i]) * 
            (
            sigma7Benpx(obsx14[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx14[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx14[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx14[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx14[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx14[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx14[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

for (i in 1:length(obsx15)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy15[i] ~ dnorm(ym15[i], pow(errobsy15[i], -2))    
  # ...subject to syst uncertainties: 
  ym15[i] <- yt15[i] * y.norm15
  # true sigma [calculated from theory]: 
  yt15[i] <- sqrt(obsx15[i]) * 
            (
            sigma7Benpx(obsx15[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx15[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
          + sigma7Benpx(obsx15[i], e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3)
          + sigma7Benpx(obsx15[i], e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4)
          + sigma7Benpx(obsx15[i], e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5)
          + sigma7Benpx(obsx15[i], e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6)
          + sigma7Benpx(obsx15[i], e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
            )
}

######################################
#
# PRIORS
#
######################################
#
# NUCLEAR PRIORS
# parameters: ecm, e0, ga, gb, ra, rb, xj, xla, xlb

# channel radii 
##  ra ~ dnorm(5.0, pow(2.0, -2))T(0,)
##  rb ~ dnorm(5.0, pow(2.0, -2))T(0,)

  ra <- 5
  rb <- 5

# Wigner limit: wl = hbar^2/(m_red a_c^2) = 41.80159/(M_red a_c^2)
#
# for the reduced width priors, we will assume sigma = 0.5 * wl, so 
# that Wigner limit corresponds to 2 * sigma
#
# neutron channel: wl_n = 41.80159/(0.88186 a_c^2) = 47.40160/a_c^2
# proton channel:  wl_p = 41.80159/(0.88079 a_c^2) = 47.45920/a_c^2
#
  wl_n <- 47.40160*pow(ra, -2)
  wl_p <- 47.45920*pow(rb, -2)

###########
# RESONANCE 1: e0=0.01136 MeV                      
# resonance spin, orbital angular momenta
  xj_1  <- 2      # 2-
  xla_1 <- 0
  xlb_1 <- 0
# energy eigenvalue 
  e0_1 ~ dnorm(0.01136, pow(0.1, -2))T(0,)        
# reduced widths
  ga_1 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)      
  gb_1 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)  

###########
# RESONANCE 2: e0=0.170 MeV                      
# resonance spin, orbital angular momenta
  xj_2  <- 3      # 3+
  xla_2 <- 1
  xlb_2 <- 1

# energy eigenvalue 
  e0_2 ~ dnorm(0.170, pow(0.1, -2))T(0,)        
# reduced widths
  ga_2 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)      
  gb_2 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)  

###########
# RESONANCE 3: e0=0.336 MeV                      
# resonance spin, orbital angular momenta
  xj_3  <- 3      # 3+
  xla_3 <- 1
  xlb_3 <- 1
# energy eigenvalue 
  e0_3 ~ dnorm(0.336, pow(0.1, -2))T(0,)        
# reduced widths
  ga_3 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)      
  gb_3 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)  
  
###########
# RESONANCE 4: e0=0.501 MeV                      
# resonance spin, orbital angular momenta
  xj_4  <- 1      # 1-
  xla_4 <- 0
  xlb_4 <- 0
# energy eigenvalue 
  e0_4 ~ dnorm(0.501, pow(0.1, -2))T(0,)        
# reduced widths
  ga_4 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)      
  gb_4 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)  
  
###########
# RESONANCE 5: e0=0.961 MeV                      
# resonance spin, orbital angular momenta
  xj_5  <- 4      # 4+
  xla_5 <- 3
  xlb_5 <- 3
# energy eigenvalue 
  e0_5 ~ dnorm(0.961, pow(0.1, -2))T(0,)        
# reduced widths
  ga_5 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)      
  gb_5 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)  

###########
# RESONANCE 6: e0=1.201 MeV                      
# resonance spin, orbital angular momenta
  xj_6  <- 2      # 2+
  xla_6 <- 1
  xlb_6 <- 1
# energy eigenvalue 
  e0_6 ~ dnorm(1.201, pow(0.1, -2))T(0,)        
# reduced widths
  ga_6 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)      
  gb_6 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)  

###########
# RESONANCE 7: e0=1.301 MeV                      
# resonance spin, orbital angular momenta
  xj_7  <- 0      # 0+
  xla_7 <- 1
  xlb_7 <- 1
# energy eigenvalue 
  e0_7 ~ dnorm(1.301, pow(0.1, -2))T(0,)        
# reduced widths
  ga_7 ~ dnorm(0.0, pow(0.5 * wl_n, -2))T(0,)      
  gb_7 ~ dnorm(0.0, pow(0.5 * wl_p, -2))T(0,)        

######################################
#
# DATA SET PRIORS

###########
# RELATIVE DATA SETS
# extrinsic scatter
  yscat1 ~ dnorm(0.0, pow(2, -2))T(0,)     # 2 [sqrt(MeV)]b
  yscat2 ~ dnorm(0.0, pow(2, -2))T(0,)     
  yscat3 ~ dnorm(0.0, pow(2, -2))T(0,)     
  yscat4 ~ dnorm(0.0, pow(2, -2))T(0,)     
 
# systematic normalization factor for S-factor:
# log(): natural logarithm

# chose very broad [weakly informative] for the data sets...
# sampled values are logs of multiplication factors
  y.norm1 ~ dunif(-3, 3)
  y.norm2 ~ dunif(-3, 3)
  y.norm3 ~ dunif(-3, 3)
  y.norm4 ~ dunif(-3, 3)

###########
## ABSOLUTE DATA SETS
# ...and informative for the absolute cross section points
  y.norm10 ~ dlnorm(logmu10, pow(logsigma10, -2))
  logmu10 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma10 <- log(1.020)  # factor uncertainty is 1.02, i.e. 2.0% for Koe88

  y.norm11 ~ dlnorm(logmu11, pow(logsigma11, -2))
  logmu11 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma11 <- log(1.10)   # factor uncertainty is 1.10, i.e. 10% for Dam18

  y.norm12 ~ dlnorm(logmu12, pow(logsigma12, -2))
  logmu12 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma12 <- log(1.050)  # factor uncertainty is 1.050, i.e. 5.0% for Gib59

  y.norm13 ~ dlnorm(logmu13, pow(logsigma13, -2))
  logmu13 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma13 <- log(1.051)  # factor uncertainty is 1.051, i.e. 5.1% for Her19

  y.norm14 ~ dlnorm(logmu14, pow(logsigma14, -2))
  logmu14 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma14 <- log(1.085)  # factor uncertainty is 1.085, i.e. 8.5% for Cer89

  y.norm15 ~ dlnorm(logmu15, pow(logsigma15, -2))
  logmu15 <- log(1.0)       # median of factor uncertainty is 1.0
  logsigma15 <- log(1.032)  # factor uncertainty is 1.032, i.e. 3.2% for Tom19
  
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

n.chains <- 3
n.adapt  <- 5000   
n.burn   <- 5000 
n.iter   <- 5000  
thin     <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model;

ourmodel <- jags.model(f, data = list( 
                 'obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1
                ,'obsx2' = obsx2, 'obsy2' = obsy2, 'errobsy2' = errobsy2
                ,'obsx3' = obsx3, 'obsy3' = obsy3, 'errobsy3' = errobsy3
                ,'obsx4' = obsx4, 'obsy4' = obsy4, 'errobsy4' = errobsy4
                ,'obsx10' = obsx10, 'obsy10' = obsy10, 'errobsy10' = errobsy10
                ,'obsx11' = obsx11, 'obsy11' = obsy11, 'errobsy11' = errobsy11
                ,'obsx12' = obsx12, 'obsy12' = obsy12, 'errobsy12' = errobsy12
                ,'obsx13' = obsx13, 'obsy13' = obsy13, 'errobsy13' = errobsy13
                ,'obsx14' = obsx14, 'obsy14' = obsy14, 'errobsy14' = errobsy14
                ,'obsx15' = obsx15, 'obsy15' = obsy15, 'errobsy15' = errobsy15
                                     ),
                 n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn) 
    
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                 variable.names=c(
                  'e0_1', 'ga_1', 'gb_1'
                 ,'e0_2', 'ga_2', 'gb_2'
                 ,'e0_3', 'ga_3', 'gb_3'
                 ,'e0_4', 'ga_4', 'gb_4'
                 ,'e0_5', 'ga_5', 'gb_5'
                 ,'e0_6', 'ga_6', 'gb_6'
                 ,'e0_7', 'ga_7', 'gb_7'
                 ,'yscat1', 'y.norm1'
                 ,'yscat2', 'y.norm2'
                 ,'yscat3', 'y.norm3'
                 ,'yscat4', 'y.norm4'
                 ,'y.norm10', 'y.norm11', 'y.norm12' 
                 ,'y.norm13', 'y.norm14', 'y.norm15'
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
# matrix samplesmat contains the samples from all n.chains with rows:
# i, alpha, beta,...
# !!! make sure to check the order of the parameters in the MCMC output !!!
samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)
# output all results
capture.output(print(samplesmat, print.gap=3), file="MCMCresults7Benp")

# select first three colums [alpha, beta, gamma] from matrix
samplesmat2 <- samplesmat[,1:3]
capture.output( print(samplesmat2[sample(nrow(samplesmat2), size=1000, 
                  replace=FALSE),], 
            print.gap=3), file="MCMCsamples7Benp")

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
pdf("MCMC_7Benp_a.pdf")
plot(mcmcChain)
dev.off()

######################################################################
# FIT + DATA 
######################################################################
pdf("MCMC_7Benp_b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(4e-9, 10.0)
yLim = c(0,10)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
       axes=FALSE, main="", xlab = "", ylab = "", log="x", yaxs='i', xaxs='i' )
# control distance between axis and label [line=...]
title(xlab="Energy (MeV)", line=2, cex.lab=2.0)
title(ylab=expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)")), line=2, cex.lab=2.0)

# shade temperature region [low left corner to high right corner]
#rect(0.004, -2, 0.12, 32, col="snow2", border = NA)

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

magaxis(mgp=c(0,0.6,0), cex.axis=1.3)
box()

# plot reaction label
text(1e-7, 2, labels=expression(paste(NULL^"7","Be(n,p)",NULL^"7","Li")), cex=2.0)

# plot credible S-factors:
# calculate for these energies and the set of Bayesian samples for
# Er, g^2_i, g^2_f the S_factor curve using R code above

# define grid of energy values in MeV on a log scale for plotting of credible 
# solutions; lseq is appropriate for log scale
x1 = lseq(1e-9, 5.0, length=500)

# ecm, E0, ga, gb, ra, rb, jr, la, lb

# 500 samples are randomly selected for plotting... 

### ...individual resonances
##for ( i in round(runif(500, min=1, max=nsamp)) ) {
##   lines(x1, sqrt(x1)*(
##      Sigma7Benp(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,3], 5, 5, 2, 0, 0)
##                      ),
##    col=adjustcolor("blue", alpha=0.02), lw=0.1)
##}

### sum of all resonances
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15],  5, 5, 2, 0, 0)
    + Sigma7Benp(x1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16],  5, 5, 3, 1, 1)
    + Sigma7Benp(x1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17], 5, 5, 3, 1, 1)
    + Sigma7Benp(x1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18], 5, 5, 1, 0, 0)
    + Sigma7Benp(x1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19], 5, 5, 4, 3, 3)
    + Sigma7Benp(x1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20], 5, 5, 2, 1, 1)
    + Sigma7Benp(x1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21], 5, 5, 0, 1, 1)
                      ),
    col=adjustcolor("red", alpha=0.02), lw=0.1)
}

### plot legend
legend(0.5, 9.8, legend=c("Dam18", "Gib59", "Mar19", "Koe88", 
                          "Koe88", "Dam18", "Cer89", "Tom19"), 
        pch=c(1, 2, 0, 4, 15, 17, 2, 0, 16,  18), 
        col=c("gray40", "gray40", "gray40", "gray40", "gray40", 
              "gray40", "gray40", "gray40"), 
        pt.cex=c(1, 1, 1, 1, 1, 1, 1.1, 1.4))

### add data 

## DATA SET 1: damone 2018;
points( obsx1, obsy1, col="gray40", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="gray40" )  

## DATA SET 2: gibbons & macklin 1959;
points( obsx2, obsy2, col="gray40", pch=2, cex=1.2 )
add.error.bars(obsx2, obsy2, 0.0, errobsy2, 0.0, col="gray40" )  

## DATA SET 3: martin-hernandez 2018;
points( obsx3, obsy3, col="gray40", pch=0, cex=1.2 )
add.error.bars(obsx3, obsy3, 0.0, errobsy3, 0.0, col="gray40" )  

## DATA SET 4: koehler 1988;
points( obsx4, obsy4, col="gray40", pch=4, cex=1.2 )
add.error.bars(obsx4, obsy4, 0.0, errobsy4, 0.0, col="gray40" )  


## DATA SET 10: koehler 1988; 
points( obsx10, obsy10, col="gray40", pch=15, cex=1.2 )
add.error.bars(obsx10, obsy10, 0.0, errobsy10, 0.0, col="gray40" )  

## DATA SET 11: damone 2018;
points( obsx11, obsy11, col="gray40", pch=17, cex=1.2 )
add.error.bars(obsx11, obsy11, 0.0, errobsy11, 0.0, col="gray40" )  

## DATA SET 12: gibbons 1959;
points( obsx12, obsy12, col="gray40", pch=2, cex=1.2 )
add.error.bars(obsx12, obsy12, 0.0, errobsy12, 0.0, col="gray40" )  

## DATA SET 13: hernandez 2019;
points( obsx13, obsy13, col="gray40", pch=0, cex=1.2 )
add.error.bars(obsx13, obsy13, 0.0, errobsy13, 0.0, col="gray40" )  

## DATA SET 14: cervena 1989;
points( obsx14, obsy14, col="gray40", pch=16, cex=1.2 )
add.error.bars(obsx14, obsy14, 0.0, errobsy14, 0.0, col="gray40" )  

## DATA SET 15: tomandl 2019;
points( obsx15, obsy15, col="gray40", pch=18, cex=1.5 )
add.error.bars(obsx15, obsy15, 0.0, errobsy15, 0.0, col="gray40" )  


### plot Maxwell-Boltzmann "factor" at T=1 GK [arbitrary scale]
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, 3e2*(x1*(2.718^(-x1/(0.086173*1.0)))) ,
    col=adjustcolor("gray"), lw=1, lty=1)
}

dev.off()

