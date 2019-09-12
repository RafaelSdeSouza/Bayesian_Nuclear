######################################################################
# Author: Christian Iliadis (05/13/2019) and HK
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
library(nimble)
library(magicaxis)
require(gsl)
require(RcppGSL)
library(emdbook)
setwd("C:/Users/Hong Kiat/Desktop/UROPS Materials/UROPS 1 - Thermonuclear Reaction Rates/Working Folder/Nimble 7Benp/Resonances/New Set 2")
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

sigma7Benp <- function(ecm, e0, ga, gb, ra, rb, jr, la, lb){
  # input masses, charges, angular momenta
  m1_i = 7.01473482886 
  m2_i = 1.00866491582  # masses (amu) of 7Be and n
  m1_f = 7.01435791572
  m2_f = 1.00727646658  # masses (amu) of 7Li and p
  z1_i = 4 
  z2_i = 0              # charges of 7Be and n
  z1_f = 3 
  z2_f = 1              # charges of 7Li and p
  jt = 1.5              # spins of target, projectile
  jp = 0.5 
  Q = 1.644425          # reaction Q-value (MeV) [from nuclear masses]
  
  
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
## DATA SETS
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, sqrt(Ecm)*sigma is in sqrt(MeV)b

######################################################################
## ARTIFICIAL DATA GENERATION
## DATA INPUT: ENERGY, sqrt(E) * SIGMA 
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, sqrt(Ecm)*sigma is in sqrt(MeV)b
set.seed(10)
rn <- 4.7
rp <- 3.7
###### RESONANCE 1 PARAMETERS
p_1 <- vector()            # parameter vector
# ecm, e0, ga, gb, ra, rb, jr, la, lb

p_1[1] <- 0.013    # resonance energy
p_1[2] <- 0.7      # neutron reduced width
p_1[3] <- 0.2      # proton reduced widtrh
p_1[4] <- rn        # radius neutron channel
p_1[5] <- rp        # radius proton channel
p_1[6] <- 2        # resonance spin
p_1[7] <- 0        # orbital angular momentum neutron
p_1[8] <- 0        # orbital angular momentum proton

###### RESONANCE 2 PARAMETERS
p_2 <- vector()
# ecm, e0, ga, gb, ra, rb, jr, la, lb

p_2[1] <- 0.170    # resonance energy
p_2[2] <- 0.16      # neutron reduced width
p_2[3] <- 0.02    # proton reduced width
p_2[4] <- rn        # radius neutron channel
p_2[5] <- rp        # radius proton channel
p_2[6] <- 3        # resonance spin
p_2[7] <- 1        # orbital angular momentum neutron
p_2[8] <- 1        # orbital angular momentum proton


###### RESONANCE 3 PARAMETERS
p_3 <- vector()
# ecm, e0, ga, gb, ra, rb, jr, la, lb

p_3[1] <- 0.350    # resonance energy
p_3[2] <- 0.6     # neutron reduced width
p_3[3] <- 0.2      # proton reduced width
p_3[4] <- rn        # radius neutron channel
p_3[5] <- rp        # radius proton channel
p_3[6] <- 3        # resonance spin
p_3[7] <- 1        # orbital angular momentum neutron
p_3[8] <- 1        # orbital angular momentum proton

###### RESONANCE 4 PARAMETERS
p_4 <- vector()
# ecm, e0, ga, gb, ra, rb, jr, la, lb; 3+

p_4[1] <- 0.6      # resonance energy
p_4[2] <- 0.5     # neutron reduced width
p_4[3] <- 0.2     # proton reduced width
p_4[4] <- rn        # radius neutron channel
p_4[5] <- rp        # radius proton channel
p_4[6] <- 1        # resonance spin
p_4[7] <- 0        # orbital angular momentum neutron
p_4[8] <- 0        # orbital angular momentum proton

###### RESONANCE 5 PARAMETERS
p_5 <- vector()
# ecm, e0, ga, gb, ra, rb, jr, la, lb; 3+

p_5[1] <- 1      # resonance energy
p_5[2] <- 1.2     # neutron reduced width
p_5[3] <- 0.4     # proton reduced width
p_5[4] <- rn        # radius neutron channel
p_5[5] <- rp        # radius proton channel
p_5[6] <- 4        # resonance spin
p_5[7] <- 3        # orbital angular momentum neutron
p_5[8] <- 3        # orbital angular momentum proton

###### RESONANCE 6 PARAMETERS
p_6 <- vector()
# ecm, e0, ga, gb, ra, rb, jr, la, lb; 3+

p_6[1] <- 1.4      # resonance energy
p_6[2] <- 1     # neutron reduced width
p_6[3] <- 0.4     # proton reduced width
p_6[4] <- rn        # radius neutron channel
p_6[5] <- rp        # radius proton channel
p_6[6] <- 2        # resonance spin
p_6[7] <- 1        # orbital angular momentum neutron
p_6[8] <- 1        # orbital angular momentum proton

###### RESONANCE 7 PARAMETERS
p_7 <- vector()
# ecm, e0, ga, gb, ra, rb, jr, la, lb; 3+

p_7[1] <- 1.5      # resonance energy
p_7[2] <- 0.5     # neutron reduced width
p_7[3] <- 0.14     # proton reduced width
p_7[4] <- rn        # radius neutron channel
p_7[5] <- rp        # radius proton channel
p_7[6] <- 0        # resonance spin
p_7[7] <- 1        # orbital angular momentum neutron
p_7[8] <- 1        # orbital angular momentum proton

###### RESONANCE 8 PARAMETERS
p_8 <- vector()
# ecm, e0, ga, gb, ra, rb, jr, la, lb; 3+

p_8[1] <- 2.287      # resonance energy
p_8[2] <- 1     # neutron reduced width
p_8[3] <- 0.1     # proton reduced width
p_8[4] <- rn        # radius neutron channel
p_8[5] <- rp        # radius proton channel
p_8[6] <- 4        # resonance spin
p_8[7] <- 2        # orbital angular momentum neutron
p_8[8] <- 2        # orbital angular momentum proton

###### RESONANCE 9 PARAMETERS
p_9 <- vector()
# ecm, e0, ga, gb, ra, rb, jr, la, lb; 3+

p_9[1] <- 3.54      # resonance energy
p_9[2] <- 0.1       # neutron reduced width
p_9[3] <- 0.1       # proton reduced width
p_9[4] <- rn        # radius neutron channel
p_9[5] <- rp        # radius proton channel
p_9[6] <- 1        # resonance spin
p_9[7] <- 0        # orbital angular momentum neutron
p_9[8] <- 0        # orbital angular momentum proton

############ DATA SET 1 (Mimic Damone)
N1 <- 30

p1 <- vector()            # parameter vector
obsy1 <- vector()
errobsy1 <- vector()
y1 <- vector()            # "true" value
stat1 <- 5                # x% statistical uncertainty 

obsx1 <- lseq(10^(-7),10^(-3),N1)

# generate true value first
for (i in 1:length(obsx1)){
  y1[i] <- sqrt(obsx1[i]) * 
    (
      sigma7Benp(obsx1[i], p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
        sigma7Benp(obsx1[i], p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
        sigma7Benp(obsx1[i], p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
        sigma7Benp(obsx1[i], p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
        sigma7Benp(obsx1[i], p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
        sigma7Benp(obsx1[i], p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
        sigma7Benp(obsx1[i], p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
        sigma7Benp(obsx1[i], p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
        sigma7Benp(obsx1[i], p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
    )
  # statistical uncertainty given by x% of true value          
  errobsy1[i] <- abs(rnorm(1, 0, 0.01 * stat1 * y1[i]))          
}


# randomize to find observed value
for (i in 1:length(obsx1)){
  obsy1[i] <- 1.2 * rnorm( 1, y1[i], errobsy1[i] ) 
}

# reported statistical uncertainty should also be scaled
errobsy1 <- errobsy1 * 1.2

############ DATA SET 2 (Mimic Gibbons)
N2 <- 40

p2 <- vector()            # parameter vector
obsy2 <- vector()
errobsy2 <- vector()
y2 <- vector()            # "true" value
stat2 <- 5                # x% statistical uncertainty 

obsx2 <- lseq(10^(-7),10^(-3),N2)

# generate true value first
for (i in 1:length(obsx2)){
  y2[i] <- sqrt(obsx2[i]) * 
    (
      sigma7Benp(obsx2[i], p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
        sigma7Benp(obsx2[i], p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
        sigma7Benp(obsx2[i], p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
        sigma7Benp(obsx2[i], p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
        sigma7Benp(obsx2[i], p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
        sigma7Benp(obsx2[i], p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
        sigma7Benp(obsx2[i], p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
        sigma7Benp(obsx2[i], p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
        sigma7Benp(obsx2[i], p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
    )
  # statistical uncertainty given by x% of true value          
  errobsy2[i] <- abs(rnorm(1, 0, 0.01 * stat2 * y2[i]))          
}


# randomize to find observed value
for (i in 1:length(obsx2)){
  obsy2[i] <- 0.85 * rnorm( 1, y2[i], errobsy2[i] ) 
}

# reported statistical uncertainty should also be scaled
errobsy2 <- errobsy2 * 0.85

############ DATA SET 3 (Mimic Martin-Her)
N3 <- 20

p3 <- vector()            # parameter vector
obsy3 <- vector()
errobsy3 <- vector()
y3 <- vector()            # "true" value
stat3 <- 5                # x% statistical uncertainty 

obsx3 <- lseq(10^(-4),10^(-2),N3)

# generate true value first
for (i in 1:length(obsx3)){
  y3[i] <- sqrt(obsx3[i]) * 
    (
      sigma7Benp(obsx3[i], p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
        sigma7Benp(obsx3[i], p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
        sigma7Benp(obsx3[i], p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
        sigma7Benp(obsx3[i], p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
        sigma7Benp(obsx3[i], p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
        sigma7Benp(obsx3[i], p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
        sigma7Benp(obsx3[i], p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
        sigma7Benp(obsx3[i], p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
        sigma7Benp(obsx3[i], p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
    )
  # statistical uncertainty given by x% of true value          
  errobsy3[i] <- abs(rnorm(1, 0, 0.01 * stat3 * y3[i]))          
}


# randomize to find observed value
for (i in 1:length(obsx3)){
  obsy3[i] <- 1.08 * rnorm( 1, y3[i], errobsy3[i] ) 
}

# reported statistical uncertainty should also be scaled
errobsy3 <- errobsy3 * 1.08

############ DATA SET 4 (Mimic Koehler)
N4 <- 30

p4 <- vector()            # parameter vector
obsy4 <- vector()
errobsy4 <- vector()
y4 <- vector()            # "true" value
stat4 <- 5                # x% statistical uncertainty 


obsx4 <- c(lseq(10^(-2),10^(-1),N4/2),seq(10^(-1),10^(-0.4),length = N4/2))

# generate true value first
for (i in 1:length(obsx4)){
  y4[i] <- sqrt(obsx4[i]) * 
    (
      sigma7Benp(obsx4[i], p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
        sigma7Benp(obsx4[i], p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
        sigma7Benp(obsx4[i], p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
        sigma7Benp(obsx4[i], p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
        sigma7Benp(obsx4[i], p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
        sigma7Benp(obsx4[i], p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
        sigma7Benp(obsx4[i], p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
        sigma7Benp(obsx4[i], p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
        sigma7Benp(obsx4[i], p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
    )
  # statistical uncertainty given by x% of true value          
  errobsy4[i] <- abs(rnorm(1, 0, 0.01 * stat4 * y4[i]))          
}


# randomize to find observed value
for (i in 1:length(obsx4)){
  obsy4[i] <- 0.95 * rnorm( 1, y4[i], errobsy4[i] ) 
}

# reported statistical uncertainty should also be scaled
errobsy4 <- errobsy4 * 0.95



### the following "data sets" consist of a single data point each, which carries the
### absolute normalization, one for each data set;
### notation: "11" labels the absolute cross section for data set "1", etc.;
### we simply chose the first data point in each set for the absolute normalization;
### notice there is no statistical uncertainty on this single data point

## DATA SET 10: Damone; thermal cross section
obsx10    <- 2.21e-8 
y10 <- sqrt(obsx10) * 
  (
    sigma7Benp(obsx10, p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
      sigma7Benp(obsx10, p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
      sigma7Benp(obsx10, p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
      sigma7Benp(obsx10, p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
      sigma7Benp(obsx10, p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
      sigma7Benp(obsx10, p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
      sigma7Benp(obsx10, p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
      sigma7Benp(obsx10, p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
      sigma7Benp(obsx10, p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
  )
errobsy10 <- abs(rnorm(1, 0, 0.01 * 5 * y10)) # stat at 5%
obsy10 <- 1.2*rnorm(1,y10,errobsy10)
errobsy10 <- 1.2*errobsy10

## DATA SET 11: Koehler; thermal cross section
obsx11    <- 2.21e-8
y11 <- sqrt(obsx11) * 
  (
    sigma7Benp(obsx11, p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
      sigma7Benp(obsx11, p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
      sigma7Benp(obsx11, p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
      sigma7Benp(obsx11, p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
      sigma7Benp(obsx11, p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
      sigma7Benp(obsx11, p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
      sigma7Benp(obsx11, p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
      sigma7Benp(obsx11, p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
      sigma7Benp(obsx11, p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
  )
errobsy11 <- abs(rnorm(1, 0, 0.01 * 5 * y11)) # stat at 5%
obsy11 <- 0.92*rnorm(1,y11,errobsy11)
errobsy11 <- 0.92*errobsy11 

## DATA SET 12: Martin H; lowest-energy data point
obsx12    <- 1.987364E-03 
y12 <- sqrt(obsx12) * 
  (
    sigma7Benp(obsx12, p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
      sigma7Benp(obsx12, p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
      sigma7Benp(obsx12, p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
      sigma7Benp(obsx12, p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
      sigma7Benp(obsx12, p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
      sigma7Benp(obsx12, p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
      sigma7Benp(obsx12, p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
      sigma7Benp(obsx12, p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
      sigma7Benp(obsx12, p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
  )
errobsy12 <- abs(rnorm(1, 0, 0.01 * 5 * y12)) # stat at 5%
obsy12 <- 1.08*rnorm(1,y12,errobsy12)
errobsy12 <- 1.08*errobsy12

## DATA SET 13: Gibb; lowest-energy data point
obsx13    <- 9.813847E-03 
y13 <- sqrt(obsx13) * 
  (
    sigma7Benp(obsx13, p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
      sigma7Benp(obsx13, p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
      sigma7Benp(obsx13, p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
      sigma7Benp(obsx13, p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
      sigma7Benp(obsx13, p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
      sigma7Benp(obsx13, p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
      sigma7Benp(obsx13, p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
      sigma7Benp(obsx13, p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
      sigma7Benp(obsx13, p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
  )
errobsy13 <- abs(rnorm(1, 0, 0.01 * 5 * y13)) # stat at 5%
obsy13 <- 0.95*rnorm(1,y13,errobsy13)
errobsy13 <- 0.95*errobsy13

## DATA SET 14: cervena 1989; thermal cross section
obsx14    <-  2.21e-8 
y14 <- sqrt(obsx14) * 
  (
    sigma7Benp(obsx14, p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
      sigma7Benp(obsx14, p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
      sigma7Benp(obsx14, p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
      sigma7Benp(obsx14, p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
      sigma7Benp(obsx14, p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
      sigma7Benp(obsx14, p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
      sigma7Benp(obsx14, p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
      sigma7Benp(obsx14, p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
      sigma7Benp(obsx14, p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
  )
errobsy14 <- abs(rnorm(1, 0, 0.01 * 5 * y14)) # stat at 5%
obsy14 <- 1.11*rnorm(1,y14,errobsy14)
errobsy14 <- 1.11*errobsy14

## DATA SET 15: tomandl 2019; thermal cross section
obsx15    <- 2.21e-8
y15 <- sqrt(obsx15) * 
  (
    sigma7Benp(obsx15, p_1[1], p_1[2], p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8]) +
      sigma7Benp(obsx15, p_2[1], p_2[2], p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8]) +
      sigma7Benp(obsx15, p_3[1], p_3[2], p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8]) +
      sigma7Benp(obsx15, p_4[1], p_4[2], p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8]) +
      sigma7Benp(obsx15, p_5[1], p_5[2], p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8]) +
      sigma7Benp(obsx15, p_6[1], p_6[2], p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8]) +
      sigma7Benp(obsx15, p_7[1], p_7[2], p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8]) +
      sigma7Benp(obsx15, p_8[1], p_8[2], p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8]) +
      sigma7Benp(obsx15, p_9[1], p_9[2], p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])
  )
errobsy15 <- abs(rnorm(1, 0, 0.01 * 5 * y15)) # stat at 5%
obsy15 <- 1.03*rnorm(1,y15,errobsy15)
errobsy15 <- 1.03*errobsy15


sigma7Benp7 <- function(ecm, 
                        e0_1, ga_1, gb_1, ra_1, rb_1, jr_1, la_1, lb_1,
                        e0_2, ga_2, gb_2, ra_2, rb_2, jr_2, la_2, lb_2,
                        e0_3, ga_3, gb_3, ra_3, rb_3, jr_3, la_3, lb_3, 
                        e0_4, ga_4, gb_4, ra_4, rb_4, jr_4, la_4, lb_4, 
                        e0_5, ga_5, gb_5, ra_5, rb_5, jr_5, la_5, lb_5,
                        e0_6, ga_6, gb_6, ra_6, rb_6, jr_6, la_6, lb_6,
                        e0_7, ga_7, gb_7, ra_7, rb_7, jr_7, la_7, lb_7){
  
  SF1 <- sigma7Benp(ecm, e0_1, ga_1, gb_1, ra_1, rb_1, jr_1, la_1, lb_1)
  SF2 <- sigma7Benp(ecm, e0_2, ga_2, gb_2, ra_2, rb_2, jr_2, la_2, lb_2)
  SF3 <- sigma7Benp(ecm, e0_3, ga_3, gb_3, ra_3, rb_3, jr_3, la_3, lb_3)
  SF4 <- sigma7Benp(ecm, e0_4, ga_4, gb_4, ra_4, rb_4, jr_4, la_4, lb_4)
  SF5 <- sigma7Benp(ecm, e0_5, ga_5, gb_5, ra_5, rb_5, jr_5, la_5, lb_5)
  SF6 <- sigma7Benp(ecm, e0_6, ga_6, gb_6, ra_6, rb_6, jr_6, la_6, lb_6)
  SF7 <- sigma7Benp(ecm, e0_7, ga_7, gb_7, ra_7, rb_7, jr_7, la_7, lb_7)
  SF <- SF1 + SF2 + SF3 + SF4 + SF5 + SF6 + SF7
  return(SF = SF)
}

# Vectorised sigma function
sigmaBe7 <- nimbleRcall(function(ecm = double(1), 
                                 e0_1 = double(0), ga_1 = double(0), gb_1 = double(0), 
                                 ra_1 = double(0), rb_1 = double(0), jr_1 = double(0),
                                 la_1 = double(0), lb_1 = double(0),
                                 e0_2 = double(0), ga_2 = double(0), gb_2 = double(0), 
                                 ra_2 = double(0), rb_2 = double(0), jr_2 = double(0),
                                 la_2 = double(0), lb_2 = double(0),
                                 e0_3 = double(0), ga_3 = double(0), gb_3 = double(0), 
                                 ra_3 = double(0), rb_3 = double(0), jr_3 = double(0),
                                 la_3 = double(0), lb_3 = double(0),
                                 e0_4 = double(0), ga_4 = double(0), gb_4 = double(0), 
                                 ra_4 = double(0), rb_4 = double(0), jr_4 = double(0),
                                 la_4 = double(0), lb_4 = double(0),
                                 e0_5 = double(0), ga_5 = double(0), gb_5 = double(0), 
                                 ra_5 = double(0), rb_5 = double(0), jr_5 = double(0),
                                 la_5 = double(0), lb_5 = double(0),
                                 e0_6 = double(0), ga_6 = double(0), gb_6 = double(0), 
                                 ra_6 = double(0), rb_6 = double(0), jr_6 = double(0),
                                 la_6 = double(0), lb_6 = double(0),
                                 e0_7 = double(0), ga_7 = double(0), gb_7 = double(0), 
                                 ra_7 = double(0), rb_7 = double(0), jr_7 = double(0),
                                 la_7 = double(0), lb_7 = double(0)
){},
Rfun = "sigma7Benp7", returnType = double(1))



samplerCode <- nimbleCode({
  
  ###################
  # LIKELIHOODS
  ###################
  # - careful: dnorm is differently defined in R and JAGS! 
  # - precision=sigma^(-2)
  # - in a for loop, make sure **all** variables on the LEFT of an 
  #   expression has the index [i]
  # - systematic error as normalization factor y.norm...
  
  ## Calling sigmaBe7mod once!
  obsxall[1:tl] <- c(obsx1[1:l1],obsx2[1:l2],obsx3[1:l3],obsx4[1:l4],obsx10,obsx11,
                     obsx12,obsx13,obsx14,obsx15)
  sigmatemp[1:tl] <- sigmaBe7(obsxall[1:tl],
                              e0_1, ga_1, gb_1, ra, rb, xj_1, xla_1, xlb_1,
                              e0_2, ga_2, gb_2, ra, rb, xj_2, xla_2, xlb_2,
                              e0_3, ga_3, gb_3, ra, rb, xj_3, xla_3, xlb_3,
                              e0_4, ga_4, gb_4, ra, rb, xj_4, xla_4, xlb_4,
                              e0_5, ga_5, gb_5, ra, rb, xj_5, xla_5, xlb_5,
                              e0_6, ga_6, gb_6, ra, rb, xj_6, xla_6, xlb_6,
                              e0_7, ga_7, gb_7, ra, rb, xj_7, xla_7, xlb_7)
  
  # Damone19
  for (i in 1:l1) {
    #
    # S-FACTOR
    # ...subject to stat and syst uncertainties:
    obsy1[i] ~ dnorm(ym1[i], pow(sqrt(errobsy1[i]^2 + yscat1^2), -2)) 
  }    
  ym1[1:l1] <- y.norm1 * (sqrt(obsx1[1:l1]) * sigmatemp[1:l1] + k )
  
  # Gibbons
  for (i in 1:l2) {
    #
    # S-FACTOR
    # ...subject to stat and syst uncertainties:
    obsy2[i] ~ dnorm(ym2[i], pow(sqrt(errobsy2[i]^2+yscat2^2), -2))  
  }    
  ym2[1:l2] <- y.norm2 * (sqrt(obsx2[1:l2]) * sigmatemp[is2:ie2] + k)
  
  
  # Martin Hern
  for (i in 1:l3) {
    #
    # S-FACTOR
    # ...subject to stat and syst uncertainties:
    obsy3[i] ~ dnorm(ym3[i], pow(sqrt(errobsy3[i]^2+yscat3^2), -2))  
  }    
  ym3[1:l3] <- y.norm3 * (sqrt(obsx3[1:l3]) * sigmatemp[is3:ie3] + k)
  
  # Koehler
  
  for (i in 1:l4) {
    #
    # S-FACTOR
    # ...subject to stat and syst uncertainties:
    obsy4[i] ~ dnorm(ym4[i], pow(sqrt(errobsy4[i]^2+yscat4^2), -2)) 
  }    
  ym4[1:l4] <- y.norm4 * (sqrt(obsx4[1:l4]) * sigmatemp[is4:ie4] + k)
  
  # Single data point, no extrinsic scatter
  # Koehler 88
  obsy10 ~ dnorm(ym10, sd = errobsy10)
  ym10 <- y.norm10 * (sqrt(obsx10) * sigmatemp[is10] + k)
  
  # Damone 18
  obsy11 ~ dnorm(ym11, sd = errobsy11)
  ym11 <- y.norm11 * (sqrt(obsx11) * sigmatemp[is11] + k)
  
  # Gibbons 59
  obsy12 ~ dnorm(ym12, sd = errobsy12)
  ym12 <- y.norm12 * (sqrt(obsx12) * sigmatemp[is12] + k)
  
  # Hern 19
  obsy13 ~ dnorm(ym13, sd = errobsy13)
  ym13 <- y.norm13 * (sqrt(obsx13) * sigmatemp[is13] + k)
  
  # Cervena
  obsy14 ~ dnorm(ym14, sd = errobsy14)
  ym14 <- y.norm14 * (sqrt(obsx14) * sigmatemp[is14] + k)
  
  #Tomandl
  obsy15 ~ dnorm(ym15, sd = errobsy15)
  ym15 <- y.norm15 * (sqrt(obsx15) * sigmatemp[is15] + k)
  
  ###################
  # PRIORS
  ###################
  # parameters: ecm, e0, ga, gb, ra, rb, xj, xla, xlb
  
  # channel radii 
  ra ~ T(dnorm(4.0, pow(0.5, -2)),0,Inf)
  rb ~ T(dnorm(4.0, pow(0.5, -2)),0,Inf)
  
  # Wigner limit: wl = hbar^2/(m_red a_c^2) = 41.80159/(M_red a_c^2)
  #
  # neutron channel: wl_n = 41.80159/(0.88186 a_c^2) = 47.40160/a_c^2
  # proton channel:  wl_p = 41.80159/(0.88079 a_c^2) = 47.45920/a_c^2
  #
  wl_n <- 47.40160*pow(ra, -2)
  wl_p <- 47.45920*pow(rb, -2)
  
  
  ###################################################################
  # RESONANCE 1: e0=0 MeV                      
  # resonance spin, orbital angular momenta
  xj_1  <- 2      # 2-
  xla_1 <- 0
  xlb_1 <- 0
  
  # energy eigenvalue 
  e0_1 ~ T(dnorm(0.0, pow(0.1, -2)),0,Inf)        
  
  # reduced widths
  ga_1 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)      
  gb_1 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)
  
  ##############
  # RESONANCE 2: e0=0.150 MeV
  # resonance spin, orbital angular momenta
  xj_2  <- 3      # 3+
  xla_2 <- 1
  xlb_2 <- 1
  
  # energy eigenvalue 
  e0_2 ~ T(dnorm(0.15, pow(0.025, -2)),0,Inf)         # positive since we see sigma peak 
  
  # reduced widths
  ga_2 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)      
  gb_2 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)
  
  
  ##############
  # RESONANCE 3: e0=0.336 MeV 
  # resonance spin, orbital angular momenta
  xj_3  <- 3      # 3+
  xla_3 <- 1
  xlb_3 <- 1
  
  # energy eigenvalue 
  e0_3 ~ T(dnorm(0.336, pow(0.010, -2)),0,Inf)         # positive since we see sigma peak 
  
  # reduced widths    
  ga_3 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)      
  gb_3 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)
  
  ##############
  # RESONANCE 4: e0=0.510 MeV                        
  # resonance spin, orbital angular momenta
  xj_4  <- 1      # 1-
  xla_4 <- 0      
  xlb_4 <- 0      
  
  # energy eigenvalue 
  e0_4 ~ T(dnorm(0.51, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak 
  
  # reduced widths
  ga_4 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)      
  gb_4 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)
  
  ##############
  # RESONANCE 5: e0=0.960 MeV                        
  # resonance spin, orbital angular momenta
  xj_5  <- 4      # 4+
  xla_5 <- 3      
  xlb_5 <- 3      
  
  # energy eigenvalue 
  e0_5 ~ T(dnorm(0.96, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak 
  
  # reduced widths
  ga_5 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)      
  gb_5 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)
  
  ##############
  # RESONANCE 6: e0=1.23 MeV                        
  # resonance spin, orbital angular momenta
  xj_6  <- 2      # 2+
  xla_6 <- 1      
  xlb_6 <- 1      
  
  # energy eigenvalue 
  e0_6 ~ T(dnorm(1.23, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak 
  
  # reduced widths
  ga_6 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)      
  gb_6 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)
  
  ##############
  # RESONANCE 7: e0 = 1.32 MeV                        
  # resonance spin, orbital angular momenta
  xj_7  <- 0      # 0+
  xla_7 <- 1      
  xlb_7 <- 1      
  
  # energy eigenvalue 
  e0_7 ~ T(dnorm(1.32, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak 
  # reduced widths
  ga_7 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)      
  gb_7 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)
  ##################################################################
  
  # extrinsic scatter
  yscat1 ~ T(dnorm(0.0, pow(2, -2)),0,Inf)     # 2 [sqrt(MeV)]b
  yscat2 ~ T(dnorm(0.0, pow(2, -2)),0,Inf)
  yscat3 ~ T(dnorm(0.0, pow(2, -2)),0,Inf)  
  yscat4 ~ T(dnorm(0.0, pow(2, -2)),0,Inf)  
  
  # Relative data
  y.norm1f ~ dunif(-1,1)
  y.norm1 <- 10^y.norm1f
  y.norm2f ~ dunif(-1,1)
  y.norm2 <- 10^y.norm2f
  y.norm3f ~ dunif(-1,1)
  y.norm3 <- 10^y.norm3f
  y.norm4f ~ dunif(-1,1)
  y.norm4 <- 10^y.norm4f
  
  # Absolute Data
  # Damone
  y.norm10 ~ dlnorm(logmu10, pow(logsigma10, -2))   
  logmu10 <- log(1.0)       
  logsigma10 <- log(1.15)   
  
  #Gibbons
  y.norm11 ~ dlnorm(logmu11, pow(logsigma11, -2))   
  logmu11 <- log(1.0)       
  logsigma11 <- log(1.1)   
  
  #Martin
  y.norm12 ~ dlnorm(logmu12, pow(logsigma12, -2))   
  logmu12 <- log(1.0)       
  logsigma12 <- log(1.1)   
  
  #Koehler
  y.norm13 ~ dlnorm(logmu13, pow(logsigma13, -2))   
  logmu13 <- log(1.0)       
  logsigma13 <- log(1.05)   
  
  #Cervena
  y.norm14 ~ dlnorm(logmu14, pow(logsigma14, -2))   
  logmu14 <- log(1.0)       
  logsigma14 <- log(1.1)   
  
  #Tomandl
  y.norm15 ~ dlnorm(logmu15, pow(logsigma15, -2))   
  logmu15 <- log(1.0)       
  logsigma15 <- log(1.05)   
  
  
  ############################################
  # Other Priors
  ############################################
  k ~ T(dnorm(0,sd=0.15),0,Inf)
  r_1 <- ga_1/gb_1
  r_4 <- ga_4/gb_4
})

samplerData <- list(obsy1 = obsy1, obsy2 = obsy2, obsy3 = obsy3, obsy4 = obsy4,
                    obsy10 = obsy10, obsy11 = obsy11, obsy12 = obsy12, obsy13 = obsy13,
                    obsy14 = obsy14, obsy15 = obsy15,
                    obsx1 = obsx1, obsx2 = obsx2, obsx3 = obsx3, obsx4 = obsx4)

l1 <- length(obsx1)
l2 <- length(obsx2)
l3 <- length(obsx3)
l4 <- length(obsx4)

samplerConst <- list(obsx10 = obsx10, obsx11 = obsx11, obsx12 = obsx12, obsx13 = obsx13,
                     obsx14 = obsx14, obsx15 = obsx15, 
                     errobsy1 = errobsy1, errobsy2 = errobsy2, errobsy3 = errobsy3, errobsy4 = errobsy4,
                     errobsy10 = errobsy10, errobsy11 = errobsy11, errobsy12 = errobsy12, errobsy13 = errobsy13,
                     errobsy14 = errobsy14, errobsy15 = errobsy15,
                     l1 = l1, l2 = l2, l3 = l3, l4 = l4,
                     tl = l1 + l2 + l3 + l4 + 6,
                     is2 = l1+1, ie2 = l1+l2, is3 = l1+l2+1, ie3 = l1+l2+l3, is4 = l1+l2+l3+1,
                     ie4 = l1+l2+l3+l4, is10 = l1+l2+l3+l4+1, is11 = l1+l2+l3+l4+2, is12 = l1+l2+l3+l4+3,
                     is13 = l1+l2+l3+l4+4, is14 = l1+l2+l3+l4+5, is15 = l1+l2+l3+l4+6)

samplerInits <- list(y.norm1f = 0, y.norm2f = 0, y.norm3f = 0, y.norm4f = 0,
                     y.norm10 = 1, y.norm11 = 1, y.norm12 = 1, y.norm13 = 1,
                     y.norm14 = 1, y.norm15 = 1,
                     yscat1 = 0.1, yscat2 = 0.1, yscat3 = 0.1, yscat4 = 0.1,
                     e0_1 = 0.05, gb_1 = 0.6, ga_1 = 0.6,
                     e0_2 = 0.15, ga_2 = 0.6, gb_2 = 0.6,
                     e0_3 = 0.336, ga_3 = 0.6, gb_3 = 0.6,
                     e0_4 = 0.51, ga_4 = 0.6, gb_4 = 0.6,
                     e0_5 = 0.96, ga_5 = 0.6, gb_5 = 0.6,
                     e0_6 = 1.23, ga_6 = 0.6, gb_6 = 0.6,
                     e0_7 = 1.32, ga_7 = 0.6, gb_7 = 0.6,
                     ra = 4, rb = 4,
                     k = 0.1)


###############################################################
# (Alternative) If invoking Nimble MCMC stepwise (but more customisable)
##############################################################
ourmodel <- nimbleModel(code = samplerCode, constants = samplerConst,
                        data = samplerData, inits = samplerInits, check = FALSE
)
compileNimble(ourmodel)
# Always compile the model after you are done setting up with it

ourmodel$calculate()
# Calculate the log likelihood (logProb). If the value is not NA,
# we have successfully initialised the model (FULLY)
# One iteration: simulate -> new values -> calculate

conf <- configureMCMC(ourmodel,print=TRUE)
# print = TRUE tells you what kind of samplers are being used for each stochastic node

conf$addMonitors(c('e0_1','ga_1','gb_1','e0_2','ga_2','gb_2',
                   'e0_3', 'ga_3', 'gb_3', 'e0_4', 'ga_4', 'gb_4',
                   'e0_5', 'ga_5', 'gb_5', 'e0_6', 'ga_6', 'gb_6', 
                   'e0_7', 'ga_7', 'gb_7', 
                   'r_1', 'r_4', 'ra', 'rb', 'k',
                   'y.norm1','y.norm2','y.norm3', 'y.norm4',
                   'y.norm10', 'y.norm11', 'y.norm12', 'y.norm13', 'y.norm14',
                   'y.norm15',
                   'yscat1','yscat2','yscat3','yscat4'))
# Add the parameters to monitor throughout the MCMC process

conf$removeSampler(c("ga_1","gb_1","ga_4","gb_4","ra","rb","k"))

conf$addSampler(target = c("ga_1","gb_1","ga_4","gb_4","ra","rb","k"),
                type = "AF_slice")

samplerMCMC <- buildMCMC(conf)
# Note that if no configuration is required, we can directly jump towards buildMCMC.
# But you know, then theres no point in using this entire chunk of code
# where you can just run nimbleMCMC

compiledMCMC <- compileNimble(samplerMCMC,project = ourmodel)
# Complie the configured model again once we are done with configuring it;
# Just before we perform the MCMC runs
# Can include project = ourmodel; to include all the nimbleFunctions, nimbleModels that
# you would want to inherit the functions and models from 
# resetFunctions = TRUE; if you would want to reset all the previously created functions
# in order to addd the new MCMC

n.chains = 2
n.iter = 40000
n.burnin = 15000

system.time(
  mcmcChain <- runMCMC(compiledMCMC,niter = n.iter, nchain = n.chains, nburnin = n.burnin,
                       samplesAsCodaMCMC = TRUE)
)

save(mcmcChain, file = "Final32.RData")

load(file="Final32.RData")

pdf("Final32.pdf")
plot(mcmcChain)
dev.off()

samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)
samplesmat2 <- samplesmat[,1:3]
library(sfsmisc) 
library(plotrix)
library(emdbook)
library(magicaxis)

quan <- summary(mcmcChain, quantile=c(0.16,0.5,0.84))[[2]]
true <- c(p_1[1],p_2[1],p_3[1],p_4[1],p_5[1],p_6[1],p_7[1],
          p_1[2],p_2[2],p_3[2],p_4[2],p_5[2],p_6[2],p_7[2],
          p_1[3],p_2[3],p_3[3],p_4[3],p_5[3],p_6[3],p_7[3],
          0.1,p_1[2]/p_1[3],p_4[2]/p_4[3],rn,rp,1.2,1.2,0.92,1.08,
          0.95,1.11,1.03,"-",0.85,"-",1.08,"-",0.95,rep("-",5))
out <- cbind(quan,true)

pdf("Final32b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
    las=1)

xLim = c(4e-9, 10.0)
yLim = c(0,12)

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


# ecm, E0, ga, gb, ra, rb, jr, la, lb
# 500 samples are randomly selected for plotting... 
######### plot individual resonance curves
# resonance 1
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    sigma7Benp(x1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15], samplesmat[i,25], samplesmat[i,26], 2, 0, 0)
  ),
  col=adjustcolor("blue", alpha=0.02), lw=0.1)
}

# resonance 2
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    sigma7Benp(x1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16],  samplesmat[i,25], samplesmat[i,26], 3, 1, 1)
  ),
  col=adjustcolor("black", alpha=0.02), lw=0.1)
}

# resonance 3
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    sigma7Benp(x1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17],  samplesmat[i,25], samplesmat[i,26], 3, 1, 1)
  ),
  col=adjustcolor("green4", alpha=0.02), lw=0.1)
}

# resonance 4
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    sigma7Benp(x1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18],  samplesmat[i,25], samplesmat[i,26], 1, 0, 0)
  ),
  col=adjustcolor("chocolate1", alpha=0.02), lw=0.1)
}

# resonance 5
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    sigma7Benp(x1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19],  samplesmat[i,25], samplesmat[i,26], 4, 3, 3)
  ),
  col=adjustcolor("violet", alpha=0.02), lw=0.1)
}

# resonance 6
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    sigma7Benp(x1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20],  samplesmat[i,25], samplesmat[i,26], 2, 1, 1)
  ),
  col=adjustcolor("orange", alpha=0.02), lw=0.1)
}


# resonance 7
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    sigma7Benp(x1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21],  samplesmat[i,25], samplesmat[i,26], 0, 1, 1)
  ),
  col=adjustcolor("skyblue", alpha=0.02), lw=0.1)
}


### sum of all resonances
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    sigma7Benp(x1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15],  samplesmat[i,25], samplesmat[i,26], 2, 0, 0)
    + sigma7Benp(x1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16],  samplesmat[i,25], samplesmat[i,26], 3, 1, 1)
    + sigma7Benp(x1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17],  samplesmat[i,25], samplesmat[i,26], 3, 1, 1)
    + sigma7Benp(x1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18],  samplesmat[i,25], samplesmat[i,26], 1, 0, 0)
    + sigma7Benp(x1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19],  samplesmat[i,25], samplesmat[i,26], 4, 3, 3)
    + sigma7Benp(x1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20],  samplesmat[i,25], samplesmat[i,26], 2, 1, 1)
    + sigma7Benp(x1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21],  samplesmat[i,25], samplesmat[i,26], 0, 1, 1)
  ) + samplesmat[i,22],
  col=adjustcolor("red", alpha=0.02), lw=0.1)
}


### plot legend
legend(0.5, 9.8, legend=c("Dam18", "Gib59" , "Mar19", "Koe88", 
                          "Dam18", "Koe88", "Mar19", "Gib59", "Cer89", "Tom19"), 
       pch=c(1, 2, 0, 4, 18, 19, 15, 17, 21, 22),
       pt.cex=c(1, 1, 1, 1, 1, 1, 1.1, 1.4))



lines(x1, sqrt(x1)*(
  sigma7Benp(x1,0.6,0.5,0.2,5,5,1,0,0)),col="brown")
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


## DATA SET 10: Damone 2018; 
points( obsx10, obsy10, col="gray40", pch=18, cex=1.2 )
add.error.bars(obsx10, obsy10, 0.0, errobsy10, 0.0, col="gray40" )  

## DATA SET 11: Koehler;
points( obsx11, obsy11, col="gray40", pch=19, cex=1.2 )
add.error.bars(obsx11, obsy11, 0.0, errobsy11, 0.0, col="gray40" )  

## DATA SET 12: Hernandexs 2019;
points( obsx12, obsy12, col="gray40", pch=15, cex=1.2 )
add.error.bars(obsx12, obsy12, 0.0, errobsy12, 0.0, col="gray40" )  

## DATA SET 13: Gibbons 1959;
points( obsx13, obsy13, col="gray40", pch=17, cex=1.2 )
add.error.bars(obsx13, obsy13, 0.0, errobsy13, 0.0, col="gray40" )  

## DATA SET 14: cervena 1989;
points( obsx14, obsy14, col="gray40", pch=21, cex=1.2 )
add.error.bars(obsx14, obsy14, 0.0, errobsy14, 0.0, col="gray40" )  

## DATA SET 15: tomandl 2019;
points( obsx15, obsy15, col="gray40", pch=22, cex=1.5 )
add.error.bars(obsx15, obsy15, 0.0, errobsy15, 0.0, col="gray40" )  


### plot Maxwell-Boltzmann "factor" at T=1 GK [arbitrary scale]
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, 3e2*(x1*(2.718^(-x1/(0.086173*0.5))))) 
}

dev.off()

pdf("Final32c.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
    las=1)

xLim = c(4e-9, 10.0)
yLim = c(0,12)

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
    sigma7Benp(x1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15],  samplesmat[i,25], samplesmat[i,26], 2, 0, 0)
    + sigma7Benp(x1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16],  samplesmat[i,25], samplesmat[i,26], 3, 1, 1)
    + sigma7Benp(x1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17],  samplesmat[i,25], samplesmat[i,26], 3, 1, 1)
    + sigma7Benp(x1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18],  samplesmat[i,25], samplesmat[i,26], 1, 0, 0)
    + sigma7Benp(x1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19],  samplesmat[i,25], samplesmat[i,26], 4, 3, 3)
    + sigma7Benp(x1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20],  samplesmat[i,25], samplesmat[i,26], 2, 1, 1)
    + sigma7Benp(x1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21],  samplesmat[i,25], samplesmat[i,26], 0, 1, 1)
  ) + samplesmat[i,22],
  col=adjustcolor("red", alpha=0.02), lw=0.1)
}

### sum of all resonances
  lines(x1, sqrt(x1)*(
    sigma7Benp(x1,p_1[1],p_1[2],p_1[3], p_1[4], p_1[5], p_1[6], p_1[7], p_1[8])
    + sigma7Benp(x1,p_2[1],p_2[2],p_2[3], p_2[4], p_2[5], p_2[6], p_2[7], p_2[8])
    + sigma7Benp(x1,p_3[1],p_3[2],p_3[3], p_3[4], p_3[5], p_3[6], p_3[7], p_3[8])
    + sigma7Benp(x1,p_4[1],p_4[2],p_4[3], p_4[4], p_4[5], p_4[6], p_4[7], p_4[8])
    + sigma7Benp(x1,p_5[1],p_5[2],p_5[3], p_5[4], p_5[5], p_5[6], p_5[7], p_5[8])
    + sigma7Benp(x1,p_6[1],p_6[2],p_6[3], p_6[4], p_6[5], p_6[6], p_6[7], p_6[8])
    + sigma7Benp(x1,p_7[1],p_7[2],p_7[3], p_7[4], p_7[5], p_7[6], p_7[7], p_7[8])
    + sigma7Benp(x1,p_8[1],p_8[2],p_8[3], p_8[4], p_8[5], p_8[6], p_8[7], p_8[8])
    + sigma7Benp(x1,p_9[1],p_9[2],p_9[3], p_9[4], p_9[5], p_9[6], p_9[7], p_9[8])), col = "blue")



### plot legend
legend(0.5, 9.8, legend=c("Dam18", "Gib59" , "Mar19", "Koe88", 
                          "Dam18", "Koe88", "Mar19", "Gib59", "Cer89", "Tom19"), 
       pch=c(1, 2, 0, 4, 18, 19, 15, 17, 21, 22),
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


## DATA SET 10: Damone 2018; 
points( obsx10, obsy10, col="gray40", pch=18, cex=1.2 )
add.error.bars(obsx10, obsy10, 0.0, errobsy10, 0.0, col="gray40" )  

## DATA SET 11: Koehler;
points( obsx11, obsy11, col="gray40", pch=19, cex=1.2 )
add.error.bars(obsx11, obsy11, 0.0, errobsy11, 0.0, col="gray40" )  

## DATA SET 12: Hernandexs 2019;
points( obsx12, obsy12, col="gray40", pch=15, cex=1.2 )
add.error.bars(obsx12, obsy12, 0.0, errobsy12, 0.0, col="gray40" )  

## DATA SET 13: Gibbons 1959;
points( obsx13, obsy13, col="gray40", pch=17, cex=1.2 )
add.error.bars(obsx13, obsy13, 0.0, errobsy13, 0.0, col="gray40" )  

## DATA SET 14: cervena 1989;
points( obsx14, obsy14, col="gray40", pch=21, cex=1.2 )
add.error.bars(obsx14, obsy14, 0.0, errobsy14, 0.0, col="gray40" )  

## DATA SET 15: tomandl 2019;
points( obsx15, obsy15, col="gray40", pch=22, cex=1.5 )
add.error.bars(obsx15, obsy15, 0.0, errobsy15, 0.0, col="gray40" )  


### plot Maxwell-Boltzmann "factor" at T=1 GK [arbitrary scale]
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, 3e2*(x1*(2.718^(-x1/(0.086173*0.5))))) 
}

dev.off()