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
## DATA INPUT 
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
#  obsy1[i] ~ dnorm(ya1[i], pow(yscat1, -2))    
  # ...subject to stat uncertainties:
  obsy1[i] ~ dnorm(yt1[i], pow(errobsy1[i], -2))
  # ...subject to syst uncertainties: 
#  ym1[i] <- y.norm1 * yt1[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt1[i] <- (alpha + beta * obsx1[i] + gamma * (obsx1[i]^2))
}    

###################
# PRIORS
###################

### polynomial parameters
  alpha ~ dnorm(0.0,pow(100, -2))
  beta ~ dnorm(0.0,pow(100, -2))
  gamma ~ dnorm(0.0,pow(100, -2))

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
n.iter   <- 10000  
thin     <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model;

ourmodel <- jags.model(f, data = list( 'obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1),
                 n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn) 
    
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                 variable.names=c('alpha', 'beta', 'gamma'
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
capture.output(print(samplesmat, print.gap=3), file="MCMCresults3He3He")

# select first three colums [alpha, beta, gamma] from matrix
samplesmat2 <- samplesmat[,1:3]
capture.output( print(samplesmat2[sample(nrow(samplesmat2), size=5000, 
                  replace=FALSE),], 
            print.gap=3), file="MCMCsamples3He3He")

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
pdf("MCMC_7Benp_c.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(4e-9, 1.0)
yLim = c(0,10)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
       axes=FALSE, main="", xlab = "", ylab = "",
       cex=1.5 , cex.lab=1.3, cex.axis=1.0,
       cex.main=1.0, log="x" )
# control distance between axis and label [line=...]
title(xlab="Energy (MeV)", line=2, cex.lab=2.0)
title(ylab=expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)")), line=2, cex.lab=2.0)
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

# shade temperature region [low left corner to high right corner]
#rect(0.004, -2, 0.12, 32, col="snow2", border = NA)

magaxis(mgp=c(0,0.2,0))
box()

# plot legend
legend(0.10, 10.0, legend=c("Dam18"), 
        pch=c(1), 
        col=c("red"))

# plot reaction label
text(1e-7, 2, labels=expression(paste(NULL^"7","Be(n,p)",NULL^"7","Li")), cex=2.0)

# PLOT BELIEVABLE REGRESSION LINES

# define grid of 201 x values on a log scale for plotting of credible solutions;
# lseq is appropriate for log scale
xComb = lseq(xLim[1],xLim[2],length=201)

# next, pick number of samples to take into account for calculation 
# of y at each x; here we take 1001 samples into account
for ( i in round(seq(from=1,to=nsamp,length=1001)) ) {
  lines( xComb, samplesmat[i,"alpha"] + samplesmat[i,"beta"] * xComb 
                                      + samplesmat[i,"gamma"] * (xComb^2), 
         col=adjustcolor("black", alpha=0.01))}

# PLOT CREDIBLE REGION [red solid lines]
# the credible region is obtained by computing regression lines using
# all sampled slopes and intercepts; since we use all samples, we are 
# effectively marginalizing over all other parameters; 
# from the sampled regression lines, we can compute y values at each point 
# of our x-grid; we can thus define suitable quantiles at each x value; 
# this defines the credible region for the regression;
#
# at any given x value, the corresponding true y value most likely is
# located between the two red solid lines

# define xx values for which we would like to have the fitted values; 
# save in new data frame newdat; again, we chose 201 xx values
newdat <- data.frame(xx=lseq(xLim[1],xLim[2],length=201))

# create new model matrix that contains new xx values; the rows of the
# matrix are: 1, xx, xx^2
newmodmat <- model.matrix(~ poly(xx, degree=2, raw=TRUE), data=newdat)

# calculate y values from all regression samples at each new xx value, 
# using matrix multiplication (%*%); save values in fitmat;
# rows denote x-values at which regression line values are being extracted;
# columns denote different mcmc samples
fitmat <- matrix(ncol=nsamp, nrow=nrow(newdat))
for(i in 1:nsamp) fitmat[,i] <- newmodmat %*% samplesmat2[i,]

# in the lines() command, "newdat$xx" is the vector with the x-values; the
# corresponding y-values are found using the apply command; the first argument
# in apply() is the matrix name, the "1" indicates to perform the action over
# rows, the next arguments specify to find quantiles etc.
#
# extract 16 and 84 percentiles for each x-value from fitmat and draw lines
# for the lower and upper limits of the credible interval;
#  
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.16), lty=1, lwd=2, 
          col=adjustcolor("red", alpha=0.5))
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.84), lty=1, lwd=2, 
          col=adjustcolor("red", alpha=0.5))
# interpretation: we are 68% sure that the true regression line is within the
# credible interval

# PLOT MEDIAN REGRESSION LINE [blue solid line]
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.50), lty=1, lwd=2, 
          col=adjustcolor("blue", alpha=0.5))

# add data 
points( obsx1, obsy1, col="red", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="red" )  

dev.off()

