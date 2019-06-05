######################################################################
# Author: Christian Iliadis (05/13/2019)
######################################################################
#
# MCMCsfactor_7Benp.R
#
# PURPOSE: analyzing real 7Benp data
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
# E0    :   eigenenergy
# ga,gb :   reduced widths for incoming/outgoing channels
# ra,rb :   channel radii (fm) for incoming/outgoing channels
# ue    :   laboratory electron screening potential

Sigma7Benp <- function(ecm, e0, ga, gb, ra, rb, jr, la, lb){
  # input masses, charges, angular momenta
  m1_i = 7.0147344 
  m2_i = 1.0086649158   # masses (amu) of t and d
  m1_f = 7.014357697 
  m2_f = 1.007276452    # masses (amu) of n and 4He
  z1_i = 4 
  z2_i = 0              # charges of t and d
  z1_f = 3 
  z2_f = 1              # charges of n and 4He
  jt = 1.5              # spins of target, projectile, resonance
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
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, sqrt(Ecm)*sigma is in sqrt(MeV)b

## DATA SET 1: damone 2018
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

## DATA SET 2: koehler 1988
obsx2    <- c( 2.3518E-08, 2.4655E-08, 2.5966E-08, 2.7365E-08, 2.8939E-08,
               3.0600E-08, 3.2436E-08, 3.4359E-08, 3.6545E-08, 3.8906E-08,
               4.1441E-08, 4.4326E-08, 4.7474E-08, 5.1058E-08, 5.4905E-08,
               5.7091E-08, 5.9102E-08, 6.1112E-08, 6.3211E-08, 6.5484E-08,
               6.7844E-08, 7.0380E-08, 7.3090E-08, 7.5888E-08, 7.8860E-08,
               8.2008E-08, 8.5330E-08, 8.9177E-08, 9.2674E-08, 9.7046E-08,
               1.0142E-07, 1.0579E-07, 1.1016E-07, 1.1541E-07, 1.2153E-07,
               1.2765E-07, 1.3377E-07, 1.4076E-07, 1.4863E-07, 1.5737E-07,
               1.6611E-07, 1.7573E-07, 1.8622E-07, 1.9759E-07, 2.1070E-07,
               2.2469E-07, 2.4043E-07, 2.5704E-07, 2.7190E-07, 2.8239E-07,
               2.9289E-07, 3.0425E-07, 3.1649E-07, 3.2873E-07, 3.4185E-07,
               3.5671E-07, 3.7157E-07, 3.8818E-07, 4.0479E-07, 4.2315E-07,
               4.4326E-07, 4.6425E-07, 4.8698E-07, 5.1146E-07, 5.3769E-07,
               5.6654E-07, 5.9714E-07, 6.3036E-07, 6.6620E-07, 7.0555E-07,
               7.4839E-07, 7.9472E-07, 8.4631E-07, 9.0051E-07, 9.6171E-07,
               1.0317E-06, 1.1103E-06, 1.1978E-06, 1.2939E-06, 1.3989E-06,
               1.4863E-06, 1.5213E-06, 1.5475E-06, 1.5825E-06, 1.6174E-06,
               1.6524E-06, 1.6961E-06, 1.7311E-06, 1.7748E-06, 1.8185E-06,
               1.8535E-06, 1.9059E-06, 1.9497E-06, 1.9934E-06, 2.0458E-06,
               2.0983E-06, 2.1507E-06, 2.2032E-06, 2.2644E-06, 2.3256E-06,
               2.3868E-06, 2.4480E-06, 2.5179E-06, 2.5879E-06, 2.6666E-06,
               2.7365E-06, 2.8239E-06, 2.9026E-06, 2.9901E-06, 3.0862E-06,
               3.1824E-06, 3.2786E-06, 3.3835E-06, 3.4971E-06, 3.6108E-06,
               3.7332E-06, 3.8643E-06, 3.9955E-06, 4.1354E-06, 4.2840E-06,
               4.4414E-06, 4.6075E-06, 4.7823E-06, 4.9747E-06, 5.1670E-06,
               5.3768E-06, 5.6042E-06, 5.8402E-06, 6.0938E-06, 6.3560E-06,
               6.6446E-06, 6.9506E-06, 7.2828E-06, 7.6325E-06, 8.0084E-06,
               8.4194E-06, 8.8303E-06, 9.3548E-06, 9.8794E-06, 1.0404E-05,
               1.1016E-05, 1.1628E-05, 1.2415E-05, 1.3202E-05, 1.3901E-05,
               1.4338E-05, 1.4775E-05, 1.5300E-05, 1.5825E-05, 1.6437E-05,
               1.6961E-05, 1.7573E-05, 1.8273E-05, 1.8972E-05, 1.9671E-05,
               2.0458E-05, 2.1333E-05, 2.2207E-05, 2.3081E-05, 2.4130E-05,
               2.5179E-05, 2.6316E-05, 2.7453E-05, 2.8764E-05, 3.0163E-05,
               3.1649E-05, 3.3223E-05, 3.4971E-05, 3.6807E-05, 3.8818E-05,
               4.1004E-05, 4.3365E-05, 4.5987E-05, 4.8785E-05, 5.1845E-05,
               5.5255E-05, 5.9014E-05, 6.3123E-05, 6.7670E-05, 7.2740E-05,
               7.8423E-05, 8.4806E-05, 9.1800E-05, 9.9668E-05, 1.0929E-04,
               1.1978E-04, 1.3202E-04, 1.4601E-04, 1.6262E-04, 1.8273E-04,
               2.0546E-04, 2.3431E-04, 2.6841E-04, 3.1037E-04, 3.6370E-04,
               4.3277E-04, 5.2195E-04, 6.4260E-04, 8.1046E-04, 1.0579E-03,
               1.4251E-03, 2.0371E-03, 3.1387E-03, 5.4730E-03, 1.1803E-02 )
obsy2    <- c( 5.9042E+00, 5.7783E+00, 5.8494E+00, 5.9553E+00, 5.9200E+00,
               5.8426E+00, 5.7092E+00, 5.7833E+00, 5.7675E+00, 5.7004E+00,
               5.6959E+00, 5.6845E+00, 5.6432E+00, 5.6490E+00, 5.7174E+00,
               5.6843E+00, 5.7446E+00, 5.7526E+00, 5.7197E+00, 5.6989E+00,
               5.7069E+00, 5.7197E+00, 5.6639E+00, 5.7024E+00, 5.6866E+00,
               5.6759E+00, 5.7196E+00, 5.7635E+00, 5.7110E+00, 5.7974E+00,
               5.6877E+00, 5.6691E+00, 5.6623E+00, 5.6528E+00, 5.6509E+00,
               5.6592E+00, 5.6031E+00, 5.6352E+00, 5.6209E+00, 5.6173E+00,
               5.6123E+00, 5.6047E+00, 5.6574E+00, 5.7653E+00, 5.7378E+00,
               5.6787E+00, 5.7467E+00, 5.6428E+00, 5.7619E+00, 5.7339E+00,
               5.7582E+00, 5.6814E+00, 5.6089E+00, 5.7163E+00, 5.7941E+00,
               5.8351E+00, 5.8640E+00, 5.7507E+00, 5.7770E+00, 5.7830E+00,
               5.6924E+00, 5.6007E+00, 5.6176E+00, 5.7070E+00, 5.6608E+00,
               5.6451E+00, 5.7415E+00, 5.7323E+00, 5.7706E+00, 5.8798E+00,
               5.8740E+00, 5.7946E+00, 5.9337E+00, 5.8076E+00, 5.7075E+00,
               5.7692E+00, 6.0168E+00, 5.4831E+00, 5.8355E+00, 5.7244E+00,
               5.7787E+00, 5.8068E+00, 5.7111E+00, 5.5841E+00, 5.5920E+00,
               5.9439E+00, 5.8723E+00, 5.7641E+00, 5.7418E+00, 5.8809E+00,
               5.8255E+00, 6.0800E+00, 5.5349E+00, 5.5331E+00, 5.4138E+00,
               5.6797E+00, 6.0025E+00, 5.6582E+00, 5.7423E+00, 5.5052E+00,
               5.5988E+00, 5.6357E+00, 5.6696E+00, 5.8540E+00, 5.6125E+00,
               5.7237E+00, 5.5169E+00, 5.6699E+00, 5.5420E+00, 5.6744E+00,
               5.6747E+00, 5.6167E+00, 5.7648E+00, 5.6700E+00, 5.6626E+00,
               5.4274E+00, 5.6458E+00, 5.7328E+00, 5.6472E+00, 5.9051E+00,
               5.5405E+00, 5.8771E+00, 5.6137E+00, 5.7321E+00, 5.5918E+00,
               5.6556E+00, 5.6082E+00, 5.6961E+00, 5.7295E+00, 5.6019E+00,
               5.5833E+00, 5.8291E+00, 5.6132E+00, 5.6138E+00, 5.6627E+00,
               5.8206E+00, 5.8035E+00, 5.7471E+00, 5.6702E+00, 5.7285E+00,
               5.7984E+00, 5.6844E+00, 5.7433E+00, 5.7444E+00, 5.6262E+00,
               5.4489E+00, 5.9580E+00, 5.6169E+00, 5.8875E+00, 5.5056E+00,
               5.7081E+00, 5.6299E+00, 5.6297E+00, 5.8584E+00, 5.6372E+00,
               5.7624E+00, 5.8612E+00, 5.5889E+00, 5.7171E+00, 5.8701E+00,
               5.8509E+00, 5.8276E+00, 5.6639E+00, 5.6314E+00, 5.5305E+00,
               5.5020E+00, 5.5737E+00, 5.7895E+00, 5.5452E+00, 5.7008E+00,
               5.5710E+00, 5.6237E+00, 5.6760E+00, 6.0976E+00, 5.7243E+00,
               5.5750E+00, 5.6002E+00, 5.8316E+00, 5.7007E+00, 5.6972E+00,
               5.5702E+00, 5.5162E+00, 5.4613E+00, 5.3511E+00, 5.3002E+00,
               5.3955E+00, 5.5841E+00, 5.5100E+00, 5.3559E+00, 5.4070E+00,
               5.4612E+00, 5.3422E+00, 5.4883E+00, 5.2676E+00, 5.2445E+00,
               5.2008E+00, 5.3003E+00, 5.1460E+00, 4.9251E+00, 4.8950E+00,
               4.7943E+00, 4.5089E+00, 4.1234E+00, 3.7138E+00, 3.1397E+00 )
errobsy2 <- c( 2.7604E-01, 3.2974E-01, 2.9005E-01, 3.8048E-01, 2.5517E-01,
               2.6239E-01, 2.1612E-01, 2.4097E-01, 1.7587E-01, 2.3669E-01,
               1.9746E-01, 2.3159E-01, 3.9219E-01, 2.9375E-01, 3.0461E-01,
               1.1947E-01, 1.2155E-01, 1.2113E-01, 1.2068E-01, 1.2027E-01,
               1.1982E-01, 1.2203E-01, 1.1895E-01, 1.2121E-01, 1.2075E-01,
               1.2028E-01, 1.1977E-01, 1.2244E-01, 1.2177E-01, 1.2149E-01,
               1.2101E-01, 1.2034E-01, 1.1949E-01, 1.1890E-01, 1.1853E-01,
               1.1790E-01, 1.1704E-01, 1.2006E-01, 1.1951E-01, 1.1901E-01,
               1.1820E-01, 1.1738E-01, 1.2083E-01, 1.2002E-01, 1.1935E-01,
               1.1850E-01, 1.2258E-01, 1.2168E-01, 1.1993E-01, 1.2222E-01,
               1.1906E-01, 1.2135E-01, 1.1814E-01, 1.2040E-01, 1.2278E-01,
               1.2542E-01, 1.2191E-01, 1.1838E-01, 1.2088E-01, 1.2360E-01,
               1.1984E-01, 1.1583E-01, 1.1863E-01, 1.2158E-01, 1.1732E-01,
               1.2043E-01, 1.2364E-01, 1.1909E-01, 1.2243E-01, 1.2600E-01,
               1.2111E-01, 1.2481E-01, 1.2879E-01, 1.2336E-01, 1.1768E-01,
               1.2188E-01, 1.2645E-01, 1.2039E-01, 1.2513E-01, 1.1827E-01,
               1.2191E-01, 1.2334E-01, 1.2067E-01, 1.1825E-01, 1.1828E-01,
               1.2597E-01, 1.2372E-01, 1.2236E-01, 1.2390E-01, 1.2676E-01,
               1.2525E-01, 1.3115E-01, 1.1869E-01, 1.1860E-01, 1.1586E-01, 
               1.2168E-01, 1.2906E-01, 1.2171E-01, 1.2339E-01, 1.1742E-01,
               1.2050E-01, 1.2047E-01, 1.2218E-01, 1.2548E-01, 1.4860E-01,
               1.6046E-01, 1.5964E-01, 1.6185E-01, 1.5217E-01, 1.6514E-01,
               1.5699E-01, 1.5934E-01, 1.6923E-01, 1.7018E-01, 1.6532E-01,
               1.5071E-01, 1.6119E-01, 1.5991E-01, 1.6065E-01, 1.8007E-01,
               1.5174E-01, 1.7387E-01, 1.6839E-01, 1.6059E-01, 1.5457E-01,
               1.6000E-01, 1.5861E-01, 1.5467E-01, 1.5799E-01, 1.6639E-01,
               1.5982E-01, 1.6873E-01, 1.6192E-01, 1.5747E-01, 1.6414E-01,
               1.5959E-01, 2.1693E-01, 2.1410E-01, 2.0745E-01, 2.1288E-01,
               2.1906E-01, 2.1142E-01, 2.0084E-01, 1.9620E-01, 2.3862E-01,
               3.0293E-01, 4.2283E-01, 3.3639E-01, 4.3758E-01, 3.3244E-01,
               3.9537E-01, 3.4375E-01, 3.4625E-01, 4.0943E-01, 3.4595E-01,
               3.4375E-01, 4.0183E-01, 3.7228E-01, 3.6993E-01, 4.0280E-01,
               3.8136E-01, 3.7448E-01, 3.8772E-01, 3.5933E-01, 3.4051E-01,
               3.4317E-01, 3.4007E-01, 4.1396E-01, 3.5188E-01, 3.3644E-01,
               3.5219E-01, 3.6877E-01, 4.0010E-01, 4.6797E-01, 4.3922E-01,
               3.8653E-01, 3.8410E-01, 4.5287E-01, 4.1953E-01, 3.8380E-01,
               3.7194E-01, 3.5915E-01, 3.5451E-01, 3.2945E-01, 3.4498E-01,
               3.6116E-01, 3.1023E-01, 2.9000E-01, 2.9330E-01, 3.1090E-01,
               3.1534E-01, 3.0614E-01, 3.6043E-01, 3.3473E-01, 3.4328E-01,
               3.3285E-01, 3.4269E-01, 3.0419E-01, 2.8469E-01, 2.7321E-01,
               2.4915E-01, 2.1664E-01, 1.9608E-01, 1.5536E-01, 2.1728E-01 )

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

for (i in 1:length(obsx1)) {
  #
  # S-FACTOR
  # ...subject to stat uncertainties:
  obsy1[i] ~ dnorm(ya1[i], pow(errobsy1[i], -2))    
  # ...subject to extrinsic scatter, if any:
  ya1[i] ~ dnorm(ym1[i], pow(yscat1, -2))
  # ...subject to syst uncertainties: 
  ym1[i] <- y.norm1 * yt1[i]
  # true sigma [calculated from theory: 
  yt1[i] <- sqrt(obsx1[i]) * 
            (
            sigma7Benpx(obsx1[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx1[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
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
  ym2[i] <- y.norm2 * yt2[i]
  # true sigma [calculated from theory: 
  yt2[i] <- sqrt(obsx2[i]) * 
            (
            sigma7Benpx(obsx2[i], e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1)
          + sigma7Benpx(obsx2[i], e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2)
            )    
}
###################
# PRIORS
###################
# parameters: ecm, e0, ga, gb, ra, rb, xj, xla, xlb

# channel radii 
##  ra ~ dnorm(5.0, pow(2.0, -2))T(0,)
##  rb ~ dnorm(5.0, pow(2.0, -2))T(0,)

  ra <- 5
  rb <- 5

##############
# RESONANCE 1: e0=0.00267 MeV
# resonance spin, orbital angular momenta
  xj_1  <- 2      # 2-
  xla_1 <- 0
  xlb_1 <- 0

# energy eigenvalue 
  e0_1 ~ dnorm(0.0, pow(1, -2))T(0,)         # positive since we see sigma peak 

# reduced widths
  ga_1 ~ dnorm(0.0, pow(1, -2))T(0,)      
  gb_1 ~ dnorm(0.0, pow(1, -2))T(0,)      
##############
# RESONANCE 2: e0=0.330 MeV
# resonance spin, orbital angular momenta
  xj_2  <- 3      # 3+
  xla_2 <- 1
  xlb_2 <- 1

# energy eigenvalue 
  e0_2 ~ dnorm(0.0, pow(1, -2))T(0,)         # positive since we see sigma peak 

# reduced widths
  ga_2 ~ dnorm(0.0, pow(1, -2))T(0,)      
  gb_2 ~ dnorm(0.0, pow(1, -2))T(0,)      

##############
  
# extrinsic scatter
  yscat1 ~ dnorm(0.0, pow(2, -2))T(0,)     # 2 sqrt(MeV)b
  yscat2 ~ dnorm(0.0, pow(2, -2))T(0,)
  
# systematic normalization factor for S-factor:
# log(): natural logarithm
y.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
logmu1 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma1 <- log(1.15)   # factor uncertainty is 1.0xx, i.e. x.x% for xxx

y.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
logmu2 <- log(1.0)       # median of factor uncertainty is 1.0
logsigma2 <- log(10.0)   # factor uncertainty very large since koehler did
                         # relative measurement
  
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
n.adapt  <- 1000   
n.burn   <- 1000 
n.iter   <- 1000  
thin     <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model;

ourmodel <- jags.model(f, data = list( 
                 'obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1
                ,'obsx2' = obsx2, 'obsy2' = obsy2, 'errobsy2' = errobsy2 
                                     ),
                 n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn) 
    
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                 variable.names=c(
                 'e0_1', 'ga_1', 'gb_1', 'yscat1', 'y.norm1'
                ,'e0_2', 'ga_2', 'gb_2', 'yscat2', 'y.norm2'
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
# S-FACTOR DATA 
######################################################################
pdf("MCMC_7Benp_b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(4e-9, 1.0)
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

magaxis(mgp=c(0,0.4,0), cex.axis=1.3)
box()

# plot legend
legend(0.10, 10.0, legend=c("Dam18"), 
        pch=c(1), 
        col=c("black"))

# plot reaction label
text(1e-7, 2, labels=expression(paste(NULL^"7","Be(n,p)",NULL^"7","Li")), cex=2.0)

# plot credible S-factors:
# calculate for these energies and the set of Bayesian samples for
# Er, g^2_i, g^2_f the S_factor curve using R code above

# define grid of energy values in MeV on a log scale for plotting of credible 
# solutions; lseq is appropriate for log scale
x1 = lseq(1e-9, 1.0, length=300)

# ecm, E0, ga, gb, ra, rb, jr, la, lb

# 500 samples are randomly selected for plotting... 

# ...individual resonances
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,1],samplesmat[i,3],samplesmat[i,5], 5, 5, 2, 0, 0)
                      ),
    col=adjustcolor("blue", alpha=0.02), lw=0.1)
}
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, sqrt(x1)*(
     Sigma7Benp(x1,samplesmat[i,2],samplesmat[i,4],samplesmat[i,6], 5, 5, 3, 1, 1) 
                      ),
    col=adjustcolor("green4", alpha=0.02), lw=0.1)
}

# ...total S-factor 
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,1],samplesmat[i,3],samplesmat[i,5], 5, 5, 2, 0, 0)
    + Sigma7Benp(x1,samplesmat[i,2],samplesmat[i,4],samplesmat[i,6], 5, 5, 3, 1, 1) 
                      ),
    col=adjustcolor("red", alpha=0.02), lw=0.1)
}

# add data 
points( obsx1, obsy1, col="black", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="black" )  

dev.off()

