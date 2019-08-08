######################################################################
# Author: Christian Iliadis (07/13/2019)
######################################################################
# 
# ANALYSIS OF 7Be(n,p)7Li DATA
#
# REACTION MODEL: 
# - incoherent sum of single-level, two-channel R-matrix expressions 
#   plus a constant background term 
# - boundary condition parameter is not varied independently; we
#   chose E_lambda = E_r
#
# DATA TREATMENT:
# - we split each data sets into two parts: a "relative cross  
#   section" part, and an "absolute cross section" part; the relative 
#   part is modeled with a normalization factor [e.g., y.norm1] for  
#   which the prior is sampled from a very broad density; the absolute 
#   part consists of just one data point. It does have a statistical  
#   uncertainty. It is modeled with a normalization factor [e.g., 
#   y.norm11] for which the prior is sampled from a highly informative 
#   lognormal density, where the lognormal parameters mu and sigma 
#   a chosen according to the systematic uncertainty reported for a 
#   given experiment;
#
# OUTPUT:
# - in the plot of the reduced cross section, the red data points show
#   the "absolute" cross section of each data set
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())

require(gsl) 
require(RcppGSL)
library(sfsmisc) 
library(plotrix)
library(emdbook)
library(magicaxis)

######################################################################
## LOAD MCMC CHAINS
######################################################################
samplesmat <- read.table("7Benp_SAMP", header = T)
nsamp = nrow(samplesmat)

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
## DATA INPUT: ENERGY, sqrt(E) * SIGMA 
## [statistical uncertainties only!]
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, sqrt(Ecm)*sigma is in sqrt(MeV)b; 

### RELATIVE CROSS SECTIONS:

## DATA SET 1: damone 2018;
##             includes 7Li ground and first excited state
tab1     <- read.table("exp_damone_rel.dat", header=FALSE)
obsx1    <- c(tab1[,1])
obsy1    <- c(tab1[,2])
errobsy1 <- c(tab1[,3])

## DATA SET 2: gibbons & macklin 1959; data with Eplab>2371 keV excluded;
##             also excluded data below Encm=0.01 MeV; use relativistic 
##             results;
##             includes only 7Li ground state
tab2     <- read.table("exp_gibbons_rel.dat", header=FALSE)
obsx2    <- c(tab2[,1])
obsy2    <- c(tab2[,2])
errobsy2 <- c(tab2[,3])

## DATA SET 3: martin-hernandez 2018; absolute pn cross directly obtained 
##             from author; converted to np using relativistic kinematics; 
##             data below Encm=0.002 MeV are excluded because uncertainty 
##             in 7Be mass impacts results;
##             includes only 7Li ground state
tab3     <- read.table("exp_martin_rel.dat", header=FALSE)
obsx3    <- c(tab3[,1])
obsy3    <- c(tab3[,2])
errobsy3 <- c(tab3[,3])

## DATA SET 4: koehler 1988; for energies below 8e-5 MeV [constant reduced
##             cross section] data have been rebinned by adding 5 data points
##             and calculating weighted mean;
##             includes 7Li ground and first excited state
tab4     <- read.table("exp_koehler_rel.dat", header=FALSE)
obsx4    <- c(tab4[,1])
obsy4    <- c(tab4[,2])
errobsy4 <- c(tab4[,3])

### ABSOLUTE NORMALIZATIONS:

## DATA SET 10: koehler 1988; thermal cross section
tab10     <- read.table("exp_koehler_abs.dat", header=FALSE)
obsx10    <- c(tab10[,1])
obsy10    <- c(tab10[,2])
errobsy10 <- c(tab10[,3])

## DATA SET 11: damone 2018; thermal cross section
tab11     <- read.table("exp_damone_abs.dat", header=FALSE)
obsx11    <- c(tab11[,1])
obsy11    <- c(tab11[,2])
errobsy11 <- c(tab11[,3])

## DATA SET 12: gibbons 1959; lowest-energy data point
tab12     <- read.table("exp_gibbons_abs.dat", header=FALSE)
obsx12    <- c(tab12[,1])
obsy12    <- c(tab12[,2])
errobsy12 <- c(tab12[,3])

## DATA SET 13: hernandez 2019; lowest-energy data point
tab13     <- read.table("exp_martin_abs.dat", header=FALSE)
obsx13    <- c(tab13[,1])
obsy13    <- c(tab13[,2])
errobsy13 <- c(tab13[,3])

## DATA SET 14: cervena 1989; thermal cross section
tab14     <- read.table("exp_cervena_abs.dat", header=FALSE)
obsx14    <- c(tab14[,1])
obsy14    <- c(tab14[,2])
errobsy14 <- c(tab14[,3])

## DATA SET 15: tomandl 2019; thermal cross section
tab15     <- read.table("exp_tomandl_abs.dat", header=FALSE)
obsx15    <- c(tab15[,1])
obsy15    <- c(tab15[,2])
errobsy15 <- c(tab15[,3])

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

######################################################################
# TRACES AND DENSITIES
######################################################################
#pdf("MCMC_7Benp_a.pdf")
#traceplot(samplesmat)
#dev.off()

######################################################################
# DATA ONLY
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
title(xlab="Energy (MeV)", line=2.5, cex.lab=2.0)
title(ylab=expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)")), line=2.0, cex.lab=2.0)

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

### plot legend
legend(0.5, 9.8, legend=c("Dam18", "Gib59", "Mar19", "Koe88", "Koe88", "Dam18", 
                          "Cer89", "Tom19"), 
        pch=c(1, 2, 0, 4, 15, 16, 17, 18), 
        col=c("black", "blue", "purple", "darkgreen", "red", "red", "red", "red"), 
        pt.cex=c(1, 1, 1, 1, 1.0, 1.1, 1, 1.4))

## relative data

## DATA SET 1: damone 2018;
points( obsx1, obsy1, col="black", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="black" )  

## DATA SET 2: gibbons & macklin 1959;
points( obsx2, obsy2, col="blue", pch=2, cex=1.2 )
add.error.bars(obsx2, obsy2, 0.0, errobsy2, 0.0, col="blue" )  

## DATA SET 3: martin-hernandez 2018;
points( obsx3, obsy3, col="purple", pch=0, cex=1.2 )
add.error.bars(obsx3, obsy3, 0.0, errobsy3, 0.0, col="purple" )  

## DATA SET 4: koehler 1988;
points( obsx4, obsy4, col="darkgreen", pch=4, cex=1.2 )
add.error.bars(obsx4, obsy4, 0.0, errobsy4, 0.0, col="darkgreen" )  

## absolute data

## DATA SET 10: koehler 1988; 
points( 2.21e-8, 5.708, col="red", pch=15, cex=1.2 )
add.error.bars(2.21e-8, 5.708, 0.0, 0.127, 0.0, col="red" )  

## DATA SET 11: damone 2018;
points( 2.21e-8, 7.775, col="red", pch=16, cex=1.2 )
add.error.bars(2.21e-8, 7.775, 0.0, 0.78, 0.0, col="red" )  

## DATA SET 12: gibbons 1959;
points( obsx12, obsy12, col="gray40", pch=2, cex=1.2 )
add.error.bars(obsx12, obsy12, 0.0, errobsy12, 0.0, col="gray40" )  

## DATA SET 13: hernandez 2019;
points( obsx13, obsy13, col="gray40", pch=0, cex=1.2 )
add.error.bars(obsx13, obsy13, 0.0, errobsy13, 0.0, col="gray40" )  

## DATA SET 14: cervena 1989;
points( 1.9e-8, 6.818, col="red", pch=17, cex=1.2 )
add.error.bars(1.9e-8, 6.818, 0.0, 0.583, 0.0, col="red" )  

## DATA SET 15: tomandl 2019;
points( 2.4e-8, 6.482, col="red", pch=18, cex=1.5 )
add.error.bars(2.4e-8, 6.482, 0.0, 0.226, 0.0, col="red" )  

######### plot Maxwell-Boltzmann "factor" at T=1 GK [arbitrary scale]
x1 = lseq(1e-9, 5.0, length=500)
lines(x1, 6e2*(x1*(2.718^(-x1/(0.086173*0.5)))) ,
    col=adjustcolor("gray"), lw=1, lty=1)

dev.off()

######################################################################
# FIT + DATA [total fit only]
######################################################################
pdf("MCMC_7Benp_c.pdf",width=10,height=5,onefile=F)
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
# rect(0.004, -2, 0.12, 32, col="snow2", border = NA)

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

### sum of all resonances
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, samplesmat[i,22] + sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15],  samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
    + Sigma7Benp(x1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16],  samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(x1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17], samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(x1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18], samplesmat[i,23], samplesmat[i,24], 1, 0, 0)
    + Sigma7Benp(x1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19], samplesmat[i,23], samplesmat[i,24], 4, 3, 3)
    + Sigma7Benp(x1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20], samplesmat[i,23], samplesmat[i,24], 2, 1, 1)
    + Sigma7Benp(x1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21], samplesmat[i,23], samplesmat[i,24], 0, 1, 1)
                                         ),
    col=adjustcolor("red", alpha=0.02), lw=0.1)
}

## ADD DATA

## plot legend
legend(0.5, 9.8, legend=c("Dam18", "Gib59", "Mar19", "Koe88", "Koe88", "Dam18", 
                          "Cer89", "Tom19"), 
        pch=c(1, 2, 0, 4, 15, 16, 17, 18), 
        col=c("gray40", "gray40", "gray40", "gray40", "gray40", "gray40", 
              "gray40", "gray40"), 
        pt.cex=c(1, 1, 1, 1, 1.0, 1.1, 1, 1.4))

## relative data

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

## absolute data

## DATA SET 10: koehler 1988; 
points( 2.21e-8, 5.708, col="gray40", pch=15, cex=1.2 )
add.error.bars(2.21e-8, 5.708, 0.0, 0.127, 0.0, col="gray40" )  

## DATA SET 11: damone 2018;
points( 2.21e-8, 7.775, col="gray40", pch=16, cex=1.2 )
add.error.bars(2.21e-8, 7.775, 0.0, 0.78, 0.0, col="gray40" )  

## DATA SET 12: gibbons 1959;
points( obsx12, obsy12, col="gray40", pch=2, cex=1.2 )
add.error.bars(obsx12, obsy12, 0.0, errobsy12, 0.0, col="gray40" )  

## DATA SET 13: hernandez 2019;
points( obsx13, obsy13, col="gray40", pch=0, cex=1.2 )
add.error.bars(obsx13, obsy13, 0.0, errobsy13, 0.0, col="gray40" )  

## DATA SET 14: cervena 1989;
points( 1.9e-8, 6.818, col="gray40", pch=17, cex=1.2 )
add.error.bars(1.9e-8, 6.818, 0.0, 0.583, 0.0, col="gray40" )  

## DATA SET 15: tomandl 2019;
points( 2.4e-8, 6.482, col="gray40", pch=18, cex=1.5 )
add.error.bars(2.4e-8, 6.482, 0.0, 0.226, 0.0, col="gray40" )  


### plot Maxwell-Boltzmann "factor" at T=1 GK [arbitrary scale]
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, 6e2*(x1*(2.718^(-x1/(0.086173*0.5)))) ,
    col=adjustcolor("gray"), lw=1, lty=1)
}

dev.off()

######################################################################
# FIT + DATA [total fit and contributions]
######################################################################
pdf("MCMC_7Benp_d.pdf",width=10,height=5,onefile=F)
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
# rect(0.004, -2, 0.12, 32, col="snow2", border = NA)

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

magaxis(mgp=c(0,0.6,0), cex.axis=1.3)
box()

# plot reaction label
text(1e-7, 2, labels=expression(paste(NULL^"7","Be(n,p)",NULL^"7","Li")), cex=2.0)

## plot legend
legend(0.3, 9.8, legend=c("2- , 0.00", "3+ , 0.15", "3+ , 0.34", "1- , 0.51", 
                          "4+ , 0.96", "2+ , 1.23", "0+ , 1.32"), 
        pch=c("-", "-", "-", "-", "-", "-", "-"), 
        col=c("blue", "darkolivegreen4", "aquamarine3", 
              "cyan3", "darkgoldenrod2", "darkgray", 
              "blueviolet"), 
        pt.cex=c(1, 1, 1, 1, 1, 1, 1))

# plot credible S-factors:
# calculate for these energies and the set of Bayesian samples for
# Er, g^2_i, g^2_f the S_factor curve using R code above

# define grid of energy values in MeV on a log scale for plotting of credible 
# solutions; lseq is appropriate for log scale
x1 = lseq(1e-9, 5.0, length=500)

# ecm, E0, ga, gb, ra, rb, jr, la, lb

# 500 samples are randomly selected for plotting... 

# ...individual resonances
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, samplesmat[i,22] + sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15], samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
                      ),
      col=adjustcolor("blue", alpha=0.02), lw=0.1)
}
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, samplesmat[i,22] + sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16], samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
                      ),
      col=adjustcolor("darkolivegreen4", alpha=0.02), lw=0.1)
}
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, samplesmat[i,22] + sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17], samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
                      ),
      col=adjustcolor("aquamarine3", alpha=0.02), lw=0.1)
}
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, samplesmat[i,22] + sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18], samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
                      ),
      col=adjustcolor("cyan3", alpha=0.02), lw=0.1)
}
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, samplesmat[i,22] + sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19], samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
                      ),
      col=adjustcolor("darkgoldenrod2", alpha=0.02), lw=0.1)
}
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, samplesmat[i,22] + sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20], samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
                      ),
      col=adjustcolor("darkgray", alpha=0.02), lw=0.1)
}
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, samplesmat[i,22] + sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21], samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
                      ),
      col=adjustcolor("blueviolet", alpha=0.02), lw=0.1)
}

### sum of all resonances
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, samplesmat[i,22] + sqrt(x1)*(
      Sigma7Benp(x1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15],  samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
    + Sigma7Benp(x1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16],  samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(x1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17], samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(x1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18], samplesmat[i,23], samplesmat[i,24], 1, 0, 0)
    + Sigma7Benp(x1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19], samplesmat[i,23], samplesmat[i,24], 4, 3, 3)
    + Sigma7Benp(x1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20], samplesmat[i,23], samplesmat[i,24], 2, 1, 1)
    + Sigma7Benp(x1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21], samplesmat[i,23], samplesmat[i,24], 0, 1, 1)
                                         ),
    col=adjustcolor("red", alpha=0.02), lw=0.1)
}

## ADD DATA

## relative data

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

## absolute data

## DATA SET 10: koehler 1988; 
points( 2.21e-8, 5.708, col="gray40", pch=15, cex=1.2 )
add.error.bars(2.21e-8, 5.708, 0.0, 0.127, 0.0, col="gray40" )  

## DATA SET 11: damone 2018;
points( 2.21e-8, 7.775, col="gray40", pch=16, cex=1.2 )
add.error.bars(2.21e-8, 7.775, 0.0, 0.78, 0.0, col="gray40" )  

## DATA SET 12: gibbons 1959;
points( obsx12, obsy12, col="gray40", pch=2, cex=1.2 )
add.error.bars(obsx12, obsy12, 0.0, errobsy12, 0.0, col="gray40" )  

## DATA SET 13: hernandez 2019;
points( obsx13, obsy13, col="gray40", pch=0, cex=1.2 )
add.error.bars(obsx13, obsy13, 0.0, errobsy13, 0.0, col="gray40" )  

## DATA SET 14: cervena 1989;
points( 1.9e-8, 6.818, col="gray40", pch=17, cex=1.2 )
add.error.bars(1.9e-8, 6.818, 0.0, 0.583, 0.0, col="gray40" )  

## DATA SET 15: tomandl 2019;
points( 2.4e-8, 6.482, col="gray40", pch=18, cex=1.5 )
add.error.bars(2.4e-8, 6.482, 0.0, 0.226, 0.0, col="gray40" )  


### plot Maxwell-Boltzmann "factor" at T=1 GK [arbitrary scale]
for ( i in round(runif(500, min=1, max=nsamp)) ) {
   lines(x1, 6e2*(x1*(2.718^(-x1/(0.086173*0.5)))) ,
    col=adjustcolor("gray"), lw=1, lty=1)
}

dev.off()
######################################################################
# POSTERIOR OF REDUCED CROSS SECTION AT SELECTED ENERGIES
######################################################################
# prediction: plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

pdf("MCMC_7Benp_e.pdf", width=10, height=5, onefile=F)
par(mfcol=c(1,3), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)

# define x value for which we would like to predict y
xchoice1 <- 1e-6
xchoice2 <- 1e-4
xchoice3 <- 1e-2

# declare vector with y values
fitvec1 <- vector(mode = "numeric", length = nsamp)   
fitvec2 <- vector(mode = "numeric", length = nsamp)   
fitvec3 <- vector(mode = "numeric", length = nsamp)   

# calculate y values from all samples at given x value

for(i in 1:nsamp) fitvec1[i] <- samplesmat[i,22] + sqrt(xchoice1)*(
      Sigma7Benp(xchoice1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15],  samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
    + Sigma7Benp(xchoice1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16],  samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(xchoice1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17], samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(xchoice1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18], samplesmat[i,23], samplesmat[i,24], 1, 0, 0)
    + Sigma7Benp(xchoice1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19], samplesmat[i,23], samplesmat[i,24], 4, 3, 3)
    + Sigma7Benp(xchoice1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20], samplesmat[i,23], samplesmat[i,24], 2, 1, 1)
    + Sigma7Benp(xchoice1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21], samplesmat[i,23], samplesmat[i,24], 0, 1, 1)
                                                            )
for(i in 1:nsamp) fitvec2[i] <- samplesmat[i,22] + sqrt(xchoice2)*(
      Sigma7Benp(xchoice2,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15],  samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
    + Sigma7Benp(xchoice2,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16],  samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(xchoice2,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17], samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(xchoice2,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18], samplesmat[i,23], samplesmat[i,24], 1, 0, 0)
    + Sigma7Benp(xchoice2,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19], samplesmat[i,23], samplesmat[i,24], 4, 3, 3)
    + Sigma7Benp(xchoice2,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20], samplesmat[i,23], samplesmat[i,24], 2, 1, 1)
    + Sigma7Benp(xchoice2,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21], samplesmat[i,23], samplesmat[i,24], 0, 1, 1)
                                                            )
for(i in 1:nsamp) fitvec3[i] <- samplesmat[i,22] + sqrt(xchoice3)*(
      Sigma7Benp(xchoice3,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15],  samplesmat[i,23], samplesmat[i,24], 2, 0, 0)
    + Sigma7Benp(xchoice3,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16],  samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(xchoice3,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17], samplesmat[i,23], samplesmat[i,24], 3, 1, 1)
    + Sigma7Benp(xchoice3,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18], samplesmat[i,23], samplesmat[i,24], 1, 0, 0)
    + Sigma7Benp(xchoice3,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19], samplesmat[i,23], samplesmat[i,24], 4, 3, 3)
    + Sigma7Benp(xchoice3,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20], samplesmat[i,23], samplesmat[i,24], 2, 1, 1)
    + Sigma7Benp(xchoice3,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21], samplesmat[i,23], samplesmat[i,24], 0, 1, 1)
                                                            )

# define quantiles of y at xchoice
a16_1 <- quantile(fitvec1, prob = 0.16)
a50_1 <- quantile(fitvec1, prob = 0.50)
a84_1 <- quantile(fitvec1, prob = 0.84)

a16_2 <- quantile(fitvec2, prob = 0.16)
a50_2 <- quantile(fitvec2, prob = 0.50)
a84_2 <- quantile(fitvec2, prob = 0.84)

a16_3 <- quantile(fitvec3, prob = 0.16)
a50_3 <- quantile(fitvec3, prob = 0.50)
a84_3 <- quantile(fitvec3, prob = 0.84)

# output
cat("", "\n") 
cat("PREDICTION FOR x.choice=", xchoice1, "\n") 
cat("  16%        50%        84%", "\n")
cat(quantile(fitvec1, probs = 0.16), quantile(fitvec1, probs = 0.50),
             quantile(fitvec1, probs = 0.84), "\n")

cat("", "\n") 
cat("PREDICTION FOR x.choice=", xchoice2, "\n") 
cat("  16%        50%        84%", "\n")
cat(quantile(fitvec2, probs = 0.16), quantile(fitvec2, probs = 0.50),
             quantile(fitvec2, probs = 0.84), "\n")

cat("", "\n") 
cat("PREDICTION FOR x.choice=", xchoice3, "\n") 
cat("  16%        50%        84%", "\n")
cat(quantile(fitvec3, probs = 0.16), quantile(fitvec3, probs = 0.50),
             quantile(fitvec3, probs = 0.84), "\n")

# plot density at xchoice

plot(density(fitvec1), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.5, yaxs='i', xaxs='i')
polygon(density(fitvec1), col=adjustcolor("blue", alpha=0.5))
title(ylab="Probability density", line=4.5, cex.lab=2.5)
title(xlab=expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)")), 
      line=4.5, cex.lab=2.5)
title(main=expression(paste("E=10",NULL^"-6 ","MeV")), cex.main=2.0)

plot(density(fitvec2), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.5, yaxs='i', xaxs='i')
polygon(density(fitvec2), col=adjustcolor("blue", alpha=0.5))
title(xlab=expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)")), 
      line=4.5, cex.lab=2.5)
title(main=expression(paste("E=10",NULL^"-4 ","MeV")), cex.main=2.0)

plot(density(fitvec3), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.5, yaxs='i', xaxs='i')
polygon(density(fitvec3), col=adjustcolor("blue", alpha=0.5))
title(xlab=expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)")), 
      line=4.5, cex.lab=2.5)
title(main=expression(paste("E=10",NULL^"-2 ","MeV")), cex.main=2.0)

dev.off()

######################################################################
# DENSITIES OF SYSTEMATIC SHIFT FACTORS FOR REDUCED CROSS SECTION
######################################################################
pdf("MCMC_7Benp_f.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(0.9, 1.5)
        
# plot only results for "absolute data"   
plot(density(samplesmat$y.norm10), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, 
     xlim=xLim, yaxs='i', xaxs='i')
title(ylab="Probability density", line=3.5, cex.lab=2.3)
title(xlab="Normalization f", line=4.0, cex.lab=2.3)

polygon(density(samplesmat$y.norm10),  
     col=adjustcolor("orange", alpha=0.5))
polygon(density(samplesmat$y.norm11),  
     col=adjustcolor("red", alpha=0.5))
polygon(density(samplesmat$y.norm12), 
     col=adjustcolor("gray", alpha=0.5))
polygon(density(samplesmat$y.norm13),  
     col=adjustcolor("green", alpha=0.5))
polygon(density(samplesmat$y.norm14),  
     col=adjustcolor("blue", alpha=0.5))
polygon(density(samplesmat$y.norm15),  
     col=adjustcolor("purple", alpha=0.5))

legend("topleft", inset=.01, 
   c("Koe88", "Dam18", "Gib59", "Mar19", "Cer89", "Tom19"), 
   fill=adjustcolor(c("orange", "red", "gray", "green", "blue", "purple"), alpha=0.5), 
   horiz=FALSE, cex=1.5, box.lty=0)

dev.off()



