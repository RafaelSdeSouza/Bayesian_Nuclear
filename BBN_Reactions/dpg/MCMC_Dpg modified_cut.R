######################################################################
# MCMC_Dpg.R
#






# conditions:
# - relationship between x and y given by theory [table] 
# - theory relation from Marcucci et al.
#
# purpose: 
# - fit a theoretical model to data
#
# details:
# - x data have no error; only y data have errors
# - different data sets are statistically independent
# - includes systematic errors in all sets [n.norm]
# - includes robust regression [outliers]
# - true relationship between variables is given by nuclear theory,
#   assuming only a single [scaling] parameter [a.scale]
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package 
library('rjags')
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
# energy is in MeV, S-factor in MeVb

## DATA SET 1: ma97a
obsx1 <-    c(7.490E-02, 1.070E-01, 1.330E-01, 1.730E-01)
obsy1    <- c(6.850E-07, 7.080E-07, 9.560E-07, 1.260E-06)
errobsy1 <- c(7.020E-08, 6.840E-08, 8.400E-08, 9.820E-08)

## DATA SET 2: bys08a
obsx2 <-    c(8.280E-03, 9.490E-03, 10.10E-03)
obsy2    <- c(2.370E-07, 2.770E-07, 2.980E-07)
errobsy2 <- c(7.100E-08, 6.400E-08, 6.500E-08)

## DATA SET 3: sch97a
#obsx3 <-    c(1.000E-02, 1.670E-02, 2.330E-02, 3.000E-02, 3.670E-02, 4.330E-02,
#              5.000E-02)
#obsy3    <- c(2.425E-07, 2.740E-07, 3.452E-07, 3.974E-07, 4.452E-07, 4.738E-07,
#              4.744E-07)
#errobsy3 <- c(1.250E-08, 7.500E-09, 6.500E-09, 6.100E-09, 5.700E-09, 7.200E-09,
#              6.400E-09)



## DATA SET 3: sch97a
obsx3 <-    c(1.000E-02, 1.670E-02, 2.330E-02)
obsy3    <- c(2.425E-07, 2.740E-07, 3.452E-07)
errobsy3 <- c(1.250E-08, 7.500E-09, 6.500E-09)


## DATA SET 4: cas02a
obsx4 <-    c(2.500E-03, 2.600E-03,  2.040E-02, 2.120E-02)
obsy4    <- c(2.300E-07, 1.590E-07,  3.090E-07, 3.280E-07)
errobsy4 <- c(1.000E-07, 4.900E-08,  1.700E-08, 1.200E-08)

######################################################################
# INPUT OF THEORY MODEL [MARCUCCI et al.]
######################################################################
# E is in MeV, S is in eV b; convert latter to MeV b
theory <- read.table("Marcucci2005.dat", header=FALSE)

interp.x <- theory[,1]
interp.y <- theory[,2] * 1e-6

# we will use JAGS interp.lin function to use this theoretical S-factor:
#
# - the columns of this table define vectors x and y
# - a single point is given by x_i, y_i
# - interp.lin gives the y value for the x value provided as argument e,
#   interp.lin(e, x, y)

######################################################################                 
# rjags ----->
######################################################################                 
cat('model {

###################
# LIKELIHOODS
###################
# systematic error as normalization factor y.norm...

for (i in 1:length(obsx1)) {
  #
  # S-FACTOR
  # ...subject to extrinsic scatter:
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
  # ...subject to extrinsic scatter:
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
  # ...subject to extrinsic scatter:
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
  # ...subject to extrinsic scatter:
  obsy4[i] ~ dnorm(ya4[i], pow(yscat4, -2))    
  # ...subject to stat uncertainties:
  ya4[i] ~ dnorm(ym4[i], pow(errobsy4[i], -2))
  # ...subject to syst uncertainties: 
  ym4[i] <- y.norm4 * yt4[i]
  # true S-factor [calculated from theory and then scaled]: 
  yt4[i] <- a.scale * interp.lin(obsx4[i], interp.x, interp.y)
}    


###################
# PRIORS
###################

### scaling factor of theory [a.scale] 
  a.scale ~ dnorm(0.0, pow(10, -2))T(0,)        #  "a.scale" cannot become negative!

# extrinsic scatter in S-factor

# what happens to likelihood if yscat = 4e-6 ??
  yscat1 ~ dnorm(0.0, pow(2e-6, -2))T(0,)       # sigma is 2e-6 MeVb
  yscat2 ~ dnorm(0.0, pow(2e-6, -2))T(0,)
  yscat3 ~ dnorm(0.0, pow(2e-6, -2))T(0,)
  yscat4 ~ dnorm(0.0, pow(2e-6, -2))T(0,)

# systematic normalization factor for S-factor:
# log(): natural logarithm
  y.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
  logmu1 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma1 <- log(1.09)  # factor uncertainty is 1.09, i.e. 9% for MA97

  y.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
  logmu2 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma2 <- log(1.08)  # factor uncertainty is <1.08, i.e., <8% for BYS08
  # we are using "8%", not "<8%"

  y.norm3 ~ dlnorm(logmu3, pow(logsigma3, -2))
  logmu3 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma3 <- log(1.09)  # factor uncertainty is 1.09, i.e., 9% for SCH97

  y.norm4 ~ dlnorm(logmu4, pow(logsigma4, -2))
  logmu4 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma4 <- log(1.045) # factor uncertainty is 1.045, i.e., 4.5% for CAS02 


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
n.adapt  <- 50000   
n.burn <- 50000 
n.iter   <- 100000  
thin   <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model;
#
ourmodel <- jags.model(f, data = list(        ## jags wants all data in a list
                 'obsx1' = obsx1, 'obsy1' = obsy1, 'errobsy1' = errobsy1,
                 'obsx2' = obsx2, 'obsy2' = obsy2, 'errobsy2' = errobsy2,
                 'obsx3' = obsx3, 'obsy3' = obsy3, 'errobsy3' = errobsy3,
                 'obsx4' = obsx4, 'obsy4' = obsy4, 'errobsy4' = errobsy4,
                 'interp.x' = interp.x, 'interp.y' = interp.y
                                     ),
                 n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn, progress.bar="none") 
    
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, variable.names=c(
                 'a.scale', 
                 'y.norm1', 'y.norm2', 'y.norm3', 'y.norm4',
                 'yscat1',  'yscat2',  'yscat3',  'yscat4'
                                                    ), 
                  n.iter=n.iter, thin=thin)

######################################################################
# <---- rjags
######################################################################
# output results on screen
cat("", "\n")    # output empty line

# sample size adjusted for autocorrelation
effectiveChainLength = effectiveSize(mcmcChain) 
show(effectiveChainLength)

cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n") 

# output
#cat("REJECTION RATES:", "\n")  
#cat("", "\n")
#show(rejectionRate(mcmcChain[,1:5]))
#cat("-------------------------------------------", "\n")
#cat("GELMAN DIAGNOSTICS:", "\n")
#cat("", "\n")
#show(gelman.diag(mcmcChain[,1:5]))
#cat("-------------------------------------------", "\n")
#cat("HEIDELBERGER DIAGNOSTIC:", "\n")
#cat("", "\n")
#show(heidel.diag(mcmcChain[,1:5]))
#cat("-------------------------------------------", "\n")

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

# !!! make sure to check the order of the parameters in the MCMC output !!!
samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)

######################################################################
# TRACES AND DENSITIES
######################################################################
pdf("MCMC_Dpg_a.pdf")
plot(mcmcChain)
dev.off()

#pdf("MCMCsfactorDpg_b.pdf")
## traceplot(mcmcChain)
## densplot(mcmcChain)
#gelman.plot(mcmcChain[,1:5])
#dev.off()

######################################################################
# S-FACTOR FIT + DATA 
######################################################################
pdf("MCMC_Dpg_b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# determine plot ranges
xLim = c(2e-3,0.3)
yLim = c(8e-8,3e-6)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
       axes=FALSE, main="", xlab = "", ylab = "",
       cex=1.5 , cex.lab=1.3, cex.axis=1.0,
       cex.main=1.0, log="xy" )
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
legend(0.12, 3e-7, legend=c("MA97", "SCH97", "CAS02", "BYS08"), pch=c(1, 0, 6, 5))

text(4e-3, 2.0e-6, labels=expression(paste("D(p,",gamma,")",NULL^"3","He")), cex=1.3)

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

# select first colum [a.scale] from matrix
samplesmat2 <- samplesmat[,1]
# output samples to file for rate calculation
capture.output( print(samplesmat2, print.gap=3), file="MCMCsamplesDpg")

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

# add data MA97 - circles
points( obsx1, obsy1, col="black", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="black" )  

# add data BYS08 - diamonds
points( obsx2, obsy2, col="black", pch=5, cex=1.2 )
add.error.bars( obsx2, obsy2, 0.0, errobsy2, 0.0, col="black" )

# add data SCH97 - squares
points( obsx3, obsy3, col="black", pch=0, cex=1.2 )
add.error.bars( obsx3, obsy3, 0.0, errobsy3, 0.0, col="black" )

# add data CAS02 - triangles
points( obsx4, obsy4, col="black", pch=6, cex=1.2 )
add.error.bars( obsx4, obsy4, 0.0, errobsy4, 0.0, col="black" )

dev.off()

######################################################################
# PREDICTION 
######################################################################
# plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

pdf("MCMC_Dpg_c.pdf")

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
plot(density(fitvec))
dev.off()
######################################################################
######################################################################


