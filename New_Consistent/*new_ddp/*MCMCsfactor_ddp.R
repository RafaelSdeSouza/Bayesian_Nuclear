######################################################################
# Author: Christian Iliadis (05/13/2019)
######################################################################
#
# MCMCsfactor_ddp.R
#
# CONDITIONS:
# - relationship between x and y given by theory [table] 
# - theory relation from Arai et al., PRL 107, 132502 (2011)
#
# PURPOSE: 
# - fit a theoretical model to data
#
# DETAILS:
# - x data have no error; only y data have errors
# - no outlier identification
# - likelihoods have a Gaussian shape
# - true relationship between variables is given by nuclear theory,
#   assuming only a single [scaling] parameter [a.scale]
# - different data sets are statistically independent
# - includes systematic errors in all sets [y.norm...]
# - includes extrinsic scatter [unreported additional statistical scatter]
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package 
library("rjags")
library(magicaxis)
## for block updating [we do not need to center predictor variables]
load.module("glm")  
# random number seed
#set.seed(123)

######################################################################
## FUNCTIONS
######################################################################
# error bars
# w is the width of perpendicular end bars on errors; set to zero
add.error.bars <- function(X,Y,dX,dY,w,col){
   X0 = X; 
   Y0 = (Y-dY); 
   X1 =X; 
   Y1 = (Y+dY);
   arrows(X0, Y0, X1, Y1, lwd=2, code=3, angle=90, length=w, col=col);
   Y0 = Y; 
   X0 = (X-dX); 
   Y1 =Y; 
   X1 = (X+dX);
   arrows(X0, Y0, X1, Y1, lwd=2, code=3, angle=90, length=w, col=col);
}

######################################################################
# DATA INPUT
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, S-factor in keVb

## DATA SET 1: Leo06
obsx1    <- c(5.6100E-02, 8.6400E-02, 1.1640E-01, 1.5730E-01, 1.8960E-01,
              2.3510E-01, 2.7865E-01, 3.2305E-01)
obsy1    <- c(7.3361E+01, 7.4820E+01, 8.0457E+01, 8.8004E+01, 9.8428E+01,
              1.0653E+02, 1.2432E+02, 1.3815E+02) * 1e-3
errobsy1 <- c(3.1541E+00, 2.2281E+00, 1.7523E+00, 2.0381E+00, 1.7427E+00,
              1.4576E+00, 1.7546E+00, 1.9829E+00) * 1e-3

## DATA SET 2: Gre95
obsx2    <- c(2.4660E-02, 3.9620E-02, 4.1000E-02, 5.0000E-02, 5.2000E-02,
              7.6000E-02, 1.0200E-01, 1.2800E-01)
obsy2    <- c(6.4500E+01, 6.8500E+01, 7.0100E+01, 6.9000E+01, 6.7500E+01,
              7.4000E+01, 7.8500E+01, 8.6600E+01) * 1e-3
errobsy2 <- c(2.7000E+00, 2.9000E+00, 3.3000E+00, 2.5000E+00, 2.4000E+00,
              2.5000E+00, 2.7000E+00, 3.3000E+00) * 1e-3

## DATA SET 3: Bro90
obsx3    <- c(1.9964E-02, 2.4961E-02, 2.9959E-02, 3.4957E-02, 3.9956E-02,
              4.4956E-02, 4.9955E-02, 5.4955E-02, 5.8455E-02)
obsy3    <- c(5.9390E+01, 5.9960E+01, 6.1130E+01, 6.1980E+01, 6.2610E+01,
              6.4210E+01, 6.4990E+01, 6.5850E+01, 6.6770E+01) * 1e-3
errobsy3 <- c(5.2000E-01, 4.7000E-01, 4.5000E-01, 4.2000E-01, 4.0000E-01,
              2.7000E-01, 3.8000E-01, 3.7000E-01, 5.0000E-01) * 1e-3

## DATA SET 4: Kra87 (B)
obsx4    <- c(1.9880E-02, 2.2360E-02, 2.9820E-02, 3.4780E-02, 3.9790E-02,
              4.4720E-02, 4.9670E-02)
obsy4    <- c(5.2900E+01, 5.3700E+01, 5.5100E+01, 5.6200E+01, 5.6200E+01,
              5.8200E+01, 5.6400E+01) * 1e-3
errobsy4 <- c(3.4000E+00, 3.4000E+00, 3.5000E+00, 3.6000E+00, 3.6000E+00,
              3.7000E+00, 3.6000E+00) * 1e-3

## DATA SET 5: Kra87 (M)
obsx5    <- c(1.9600E-02, 2.4000E-02, 2.9000E-02, 3.4800E-02, 3.9700E-02,
              4.4600E-02, 4.9500E-02, 5.4300E-02, 5.9300E-02, 6.4300E-02,
              7.4000E-02, 8.3800E-02, 9.3600E-02, 1.0340E-01, 1.1320E-01,
              1.2300E-01, 1.3300E-01, 1.4250E-01, 1.5250E-01, 1.6250E-01)
obsy5    <- c(5.6800E+01, 5.5700E+01, 5.5300E+01, 5.4800E+01, 5.6200E+01,
              5.6800E+01, 5.6400E+01, 5.9200E+01, 6.3000E+01, 6.3100E+01,
              6.6800E+01, 7.0000E+01, 7.2600E+01, 7.6100E+01, 7.7800E+01,
              8.0800E+01, 8.4400E+01, 8.8300E+01, 8.5100E+01, 9.1400E+01) * 1e-3
errobsy5 <- c(8.1000E+00, 7.4000E+00, 5.8000E+00, 5.7000E+00, 4.7000E+00,
              4.6000E+00, 4.6000E+00, 4.8000E+00, 5.1000E+00, 5.1000E+00,
              5.4000E+00, 5.7000E+00, 5.9000E+00, 6.2000E+00, 6.3000E+00,
              6.6000E+00, 6.9000E+00, 7.2000E+00, 6.9000E+00, 7.4000E+00) * 1e-3

######################################################################
# INPUT OF THEORY MODEL [ARAI et al.]
######################################################################
# E is in MeV, S is in keVb; convert latter to MeV b
theory <- read.table("Arai_ddp_2011.dat", header=FALSE)

interp.x <- theory[,1]
interp.y <- theory[,2] * 1e-3

# we will use JAGS interp.lin function to use this theoretical S-factor:
#
# - the columns of this table define vectors x and y
# - a single point is given by x_i, y_i
# - interp.lin gives the y value for the x value provided as argument e,
#   interp.lin(e, x, y)

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
  # ...subject to extrinsic scatter, if any:
  obsy1[i] ~ dnorm(ya1[i], pow(yscat1, -2))    
  # ...subject to stat uncertainties:
  ya1[i] ~ dnorm(ym1[i], pow(errobsy1[i], -2))
  # ...subject to syst uncertainties: 
  ym1[i] <- y.norm1 * yt1[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt1[i] <- a.scale * interp.lin(obsx1[i], interp.x, interp.y)
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
  yt2[i] <- a.scale * interp.lin(obsx2[i], interp.x, interp.y)
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
  yt3[i] <- a.scale * interp.lin(obsx3[i], interp.x, interp.y)
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
  yt4[i] <- a.scale * interp.lin(obsx4[i], interp.x, interp.y)
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
  yt5[i] <- a.scale * interp.lin(obsx5[i], interp.x, interp.y)
}    

###################
# PRIORS
###################

### scaling factor of theory [a.scale] 
  a.scale ~ dnorm(0.0, pow(10, -2))T(0,)    #  "a.scale" cannot become negative!
  
# extrinsic scatter in S-factor
  yscat1 ~ dnorm(0.0, pow(0.2, -2))T(0,)       # sigma assumed to be 0.2 MeVb
  yscat2 ~ dnorm(0.0, pow(0.2, -2))T(0,)
  yscat3 ~ dnorm(0.0, pow(0.2, -2))T(0,)
  yscat4 ~ dnorm(0.0, pow(0.2, -2))T(0,)
  yscat5 ~ dnorm(0.0, pow(0.2, -2))T(0,)
  
# systematic normalization factor for S-factor:
# log(): natural logarithm
  y.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
  logmu1 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma1 <- log(1.02)  # factor uncertainty is 1.02, i.e. 2% for Leo06
  
  y.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
  logmu2 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma2 <- log(1.03)  # factor uncertainty is 1.03, i.e., 3% for Gre95

  y.norm3 ~ dlnorm(logmu3, pow(logsigma3, -2))
  logmu3 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma3 <- log(1.013)  # factor uncertainty is 1.013, i.e., 1.3% for Bro90
  
  y.norm4 ~ dlnorm(logmu4, pow(logsigma4, -2))
  logmu4 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma4 <- log(1.064) # factor uncertainty is 1.064, i.e., 6.4% for Kra87 (B) 

  y.norm5 ~ dlnorm(logmu5, pow(logsigma5, -2))
  logmu5 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma5 <- log(1.082) # factor uncertainty is 1.082, i.e., 8.2% for Kra87 (M) 
  
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
n.iter   <- 10000  
thin     <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model;
#
ourmodel <- jags.model(f, data = list(        ## jags wants all data in a list
                 'obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1,
                 'obsx2' = obsx2, 'obsy2' = obsy2, 'errobsy2' = errobsy2,
                 'obsx3' = obsx3, 'obsy3' = obsy3, 'errobsy3' = errobsy3,
                 'obsx4' = obsx4, 'obsy4' = obsy4, 'errobsy4' = errobsy4,
                 'obsx5' = obsx5, 'obsy5' = obsy5, 'errobsy5' = errobsy5,
                 'interp.x' = interp.x, 'interp.y' = interp.y
                                     ),
                 n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn) 
    
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                 variable.names=c(
                 'a.scale', 
                 'y.norm1', 'y.norm2', 'y.norm3', 'y.norm4', 'y.norm5',
                 'yscat1',  'yscat2',  'yscat3',  'yscat4',  'yscat5'
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
capture.output(print(samplesmat, print.gap=3), file="MCMCresultsDdp")

# select first colum [a.scale] from matrix
samplesmat2 <- samplesmat[,1]
# output samples to file for rate calculation
capture.output(print(samplesmat2, print.gap=3), file="MCMCsamplesDdp")

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
pdf("MCMC_Ddp_a.pdf")
plot(mcmcChain)
dev.off()

######################################################################
# S-FACTOR FIT + DATA 
######################################################################
pdf("MCMC_Ddp_b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# determine plot ranges
xLim = c(5e-3,0.5)
yLim = c(4e-2,2e-1)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
       axes=FALSE, main="", xlab = "", ylab = "", log="xy" )
# control distance between axis and label [line=...]
title(xlab="Energy (MeV)", line=2, cex.lab=2.0)
title(ylab="S-Factor (MeV b)", line=2.5, cex.lab=2.0)
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.2,0), cex.axis=1.3)
box()

# plot legend
legend(0.18, 0.08, legend=c("Leo06", "Gre95", "Bro90", "Kra87(B)", "Kra87(M)"), 
        pch=c(1, 5, 0, 6, 2))

text(0.01, 0.17, labels=expression(paste("D(d,p)",NULL^"3","H")), cex=2.0)

# PLOT BELIEVABLE S-FACTORS
# matrix samplesmat contains the samples from all n.chains with rows:
# i, a.scale,...

# define grid of 201 x values on a log scale for plotting of credible solutions;
# lseq is appropriate for log scale
xComb = lseq(xLim[1],xLim[2],length=201)

# find y values via interpolation; helpmat contains helpmat$x=xComb and 
# helpmat$y=theory S-factors [unscaled]
helpmat <- approx(interp.x, interp.y, xComb, method="linear")

# next, pick number of samples to take into account for calculation 
# of y at each x; here we take 1001 samples into account
for ( i in round(seq(from=1,to=nsamp,length=1001))) {
  lines( xComb, samplesmat[i, "a.scale"] * helpmat$y, 
         col=adjustcolor("black", alpha=0.01))
}

# PLOT CREDIBLE REGION [red solid lines]
# the credible region is obtained by computing credible lines using
# all sampled values of "a.scale";  
# from the sampled credible lines, we can compute y values at each point 
# of our x-grid; we can thus define suitable quantiles at each x value; 
# this defines the credible region for the regression;
#
# at any given x value, the corresponding true y value most likely is
# located between the two red solid lines

# calculate y values from all regression samples at each new xx value;
# rows denote x-values at which regression line values are being extracted;
# columns denote different mcmc samples
fitmat <- matrix(nrow=length(xComb), ncol=nsamp)

for(i in 1:nsamp) fitmat[,i] <- helpmat$y * samplesmat2[i]

# extract 16 and 84 percentiles for each x-value from fitmat and draw lines
# for the lower and upper limits of the credible interval;
#  
lines(xComb, apply(fitmat, 1, quantile, prob=0.16), lty=1, lwd=2, 
          col=adjustcolor("red", alpha=0.5))
lines(xComb, apply(fitmat, 1, quantile, prob=0.84), lty=1, lwd=2, 
          col=adjustcolor("red", alpha=0.5))
# interpretation: we are 68% sure that the true line is within the
# credible interval

# PLOT MEDIAN REGRESSION LINE [blue solid line]
lines(xComb, apply(fitmat, 1, quantile, prob=0.50), lty=1, lwd=2, 
          col=adjustcolor("blue", alpha=0.5))

# add data Leo06 - circles
points( obsx1, obsy1, col="gray40", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="gray40" )  

# add data Gre95 - diamonds
points( obsx2, obsy2, col="gray40", pch=5, cex=1.2 )
add.error.bars( obsx2, obsy2, 0.0, errobsy2, 0.0, col="gray40" )

# add data Bro90 - squares
points( obsx3, obsy3, col="gray40", pch=0, cex=1.2 )
add.error.bars( obsx3, obsy3, 0.0, errobsy3, 0.0, col="gray40" )

# add data Kra87(B) - triangles down
points( obsx4, obsy4, col="gray40", pch=6, cex=1.2 )
add.error.bars( obsx4, obsy4, 0.0, errobsy4, 0.0, col="gray40" )

# add data Kra87(M) - triangles up
points( obsx5, obsy5, col="gray40", pch=2, cex=1.2 )
add.error.bars( obsx5, obsy5, 0.0, errobsy5, 0.0, col="gray40" )

dev.off()

######################################################################
# S-FACTOR DATA 
######################################################################
pdf("MCMC_Ddp_c.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# determine plot ranges
xLim = c(5e-3,0.5)
yLim = c(4e-2,2e-1)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
       axes=FALSE, main="", xlab = "", ylab = "", log="xy" )
# control distance between axis and label [line=...]
title(xlab="Energy (MeV)", line=2, cex.lab=2.0)
title(ylab="S-Factor (MeV b)", line=2.5, cex.lab=2.0)
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.2,0), cex.axis=1.3)
box()

# plot legend
legend(0.18, 0.08, legend=c("Leo06", "Gre95", "Bro90", "Kra87(B)", "Kra87(M)"), 
        pch=c(1, 5, 0, 6, 2), col=c(68,91,258,"darkorchid","darkorange"))

text(0.01, 0.17, labels=expression(paste("D(d,p)",NULL^"3","H")), cex=2.0)

# add data Leo06 - circles
points( obsx1, obsy1, col=68, pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col=68 )  

# add data Gre95 - diamonds
points( obsx2, obsy2, col=91, pch=5, cex=1.2 )
add.error.bars( obsx2, obsy2, 0.0, errobsy2, 0.0, col=91 )

# add data Bro90 - squares
points( obsx3, obsy3, col=258, pch=0, cex=1.2 )
add.error.bars( obsx3, obsy3, 0.0, errobsy3, 0.0, col=258 )

# add data Kra87(B) - triangles down
points( obsx4, obsy4, col="darkorchid", pch=6, cex=1.2 )
add.error.bars( obsx4, obsy4, 0.0, errobsy4, 0.0, col="darkorchid" )

# add data Kra87(M) - triangles up
points( obsx5, obsy5, col="darkorange", pch=2, cex=1.2 )
add.error.bars( obsx5, obsy5, 0.0, errobsy5, 0.0, col="darkorange" )

dev.off()

######################################################################
# POSTERIOR OF S-FACTOR AT GIVEN ENERGY 
######################################################################
# plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

pdf("MCMC_Ddp_d.pdf")

# define x value for which we would like to predict y
xchoice <- 0.00

# declare vector with y values
fitvec <- vector(mode = "numeric", length = nsamp)   

# calculate y values from all regression samples at given x value
predmat <- approx(interp.x, interp.y, xchoice, method="linear")

for(i in 1:nsamp) fitvec[i] <- samplesmat2[i] * predmat$y

# define quantiles of y at xchoice
a16 <- quantile(fitvec, prob = 0.16)
a50 <- quantile(fitvec, prob = 0.50)
a84 <- quantile(fitvec, prob = 0.84)

# output
cat("", "\n") 
cat("PREDICTION FOR x.choice=",xchoice, "\n") 
cat("  16%        50%        84%", "\n")
cat(quantile(fitvec, probs = 0.16),quantile(fitvec, probs = 0.50),
             quantile(fitvec, probs = 0.84), "\n")

# plot density at xchoice
## plot(density(fitvec))
par(mfcol=c(1,1), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)

#plot(density(fitvec), main="", xlab="", ylab="",
#     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, yaxs='i', xaxs='i')
plot(1, type="n", lwd=2, col="black",  ylim=c(0, 8), xlim=c(5.0, 5.65),
     main="", xlab="", ylab="", axes=FALSE,
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, log="", yaxs='i', xaxs='i')

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.3,0), minorn=0, cex=2.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

polygon(density(fitvec*1e2), col=adjustcolor("blue", alpha=0.5))

title(xlab=expression(paste(S(0), "  (MeV b)")), line=4.0, cex.lab=2.0)
title(ylab="Probability density", line=4.5, cex.lab=2.0)

mtext(expression(paste("x10"^{-2})), side=1, line=2.0, adj=1.0, cex=1.5)

dev.off()

######################################################################
# THEORY SCALE FACTOR 
######################################################################
# plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

pdf("MCMC_Ddp_e.pdf")

# plot density at xchoice
## plot(density(fitvec))
par(mfcol=c(1,1), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)

plot(1, type="n", lwd=2, col="black",  
     main="", xlab="", ylab="", axes=FALSE, xlim=c(0.9, 1.02), ylim=c(1,50),
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, log="", yaxs='i', xaxs='i')

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.3,0), minorn=0, cex=2.0, cex.lab=1.3, cex.main=1.0, cex.axis=1.5)
box()

polygon(density(samplesmat[,1]), col=adjustcolor("blue", alpha=0.5))

title(xlab="Theory scale factor", line=3.0, cex.lab=2.0)
title(ylab="Probability density", line=3.0, cex.lab=2.0)

dev.off()

######################################################################
# DENSITIES OF S-FACTOR NORMALIZATION FACTORS
######################################################################
pdf("MCMC_Ddp_f.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(0.8, 1.15)

# plot density in first panel           
plot(density(samplesmat[,4]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, 
     xlim=xLim, yaxs='i', xaxs='i')
title(ylab="Probability density", line=3.5, cex.lab=2.3)
title(xlab="S-factor normalization", line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,2]),  
     col=adjustcolor("blue", alpha=0.5))
polygon(density(samplesmat[,3]), 
     col=adjustcolor("black", alpha=0.5))
polygon(density(samplesmat[,4]),  
     col=adjustcolor("red", alpha=0.5))
polygon(density(samplesmat[,5]),  
     col=adjustcolor("green", alpha=0.5))
polygon(density(samplesmat[,6]),  
     col=adjustcolor("violet", alpha=0.5))

legend("topleft", inset=.01, 
   c("Leo06", "Gre95", "Bro90", "Kra87(B)", "Kra87(M)"), 
   fill=adjustcolor(c("blue", "black", "red", "green", "violet"), alpha=0.5), 
   horiz=FALSE, cex=1.5, box.lty=0)

dev.off()

######################################################################
# POSTERIOR EXTRINSIC S-FACTOR SCATTER
######################################################################
pdf("MCMC_Ddp_g.pdf",width=10, height=6, onefile=F)
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
par(mfcol=c(2,3), mar=c(5.5,7.0,1.0,2.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# plot #1          
plot(density(samplesmat[,7]*1e2), main="", xlab="", ylab="", xlim=c(0, 2),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
#title(ylab="density", line=2.5, cex.lab=2.3)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3)
legend("topright", legend="Leo06", pch=NA, cex=1.5)
polygon(density(samplesmat[,7]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)
 
# plot #2                  
plot(density(samplesmat[,8]*1e2), main="", xlab="", ylab="", xlim=c(0, 2),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')
#title(ylab="Probability", line=2.5, cex.lab=2.3)
legend("topright", legend="Gre95", pch=NA, cex=1.5)
polygon(density(samplesmat[,8]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)

# plot #3                  
plot(density(samplesmat[,9]*1e2), main="", xlab="", ylab="", xlim=c(0, 0.4),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')
legend("topright", legend="Bro90", pch=NA, cex=1.5)
polygon(density(samplesmat[,9]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)

# plot #4                  
plot(density(samplesmat[,10]*1e2), main="", xlab="", ylab="", xlim=c(0, 1.5),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')
legend("topright", legend="Kra87(B)", pch=NA, cex=1.5)
polygon(density(samplesmat[,10]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)

# plot #5                  
plot(density(samplesmat[,11]*1e2), main="", xlab="", ylab="", xlim=c(0, 1.5),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')
legend("topright", legend="Kra87(M)", pch=NA, cex=1.5)
polygon(density(samplesmat[,11]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)

dev.off()

