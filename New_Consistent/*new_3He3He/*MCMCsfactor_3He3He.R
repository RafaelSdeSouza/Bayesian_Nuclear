######################################################################
# Author: Christian Iliadis (05/13/2019)
######################################################################
#
# MCMCsfactor_3He3He.R
#
# CONDITIONS:
# - relationship between x and y given by theory [table] 
# - theory relation from Arai et al., PRL 107, 132502 (2011)
# - quadratic relationship between x and y for the "bare" S-factor,
#   multiplied by an exponential: exp(pi * eta * U_e/E), where
#   U_e is the screening potential [see Adelberger et al. (2011)]
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
# energy is in MeV, S-factor in MeVb

## DATA SET 1: kudomi 2004
obsx1    <- c(0.0312, 0.0331, 0.0352, 0.0373,    
              0.0393, 0.0413, 0.0433, 0.0453 )
obsy1    <- c(6.40, 5.48, 5.62, 5.46,
              5.69, 5.51, 5.43, 5.39 )
errobsy1 <- c(0.39, 0.22, 0.21, 0.20,
              0.25, 0.18, 0.14, 0.09 )

## DATA SET 2: bonetti 1999
obsx2    <- c(0.01699, 0.01846, 0.01898, 0.01946, 0.01993, 0.02143,
              0.02337, 0.02436 )
obsy2    <- c(13.15, 7.86, 8.25, 7.67,  5.10, 4.72,
              7.31, 5.44 )
errobsy2 <- c(4.98, 2.97, 2.29, 2.22, 1.70, 0.65,
              0.63, 0.34 )

## DATA SET 3: junker 1998
obsx3    <- c(0.02076, 0.02123, 0.02175, 0.02228, 0.02233, 0.02278,    
              0.02282, 0.02315, 0.02321, 0.02370, 0.02425, 0.02430,    
              0.02452, 0.02470, 0.02480, 0.04582, 0.05064, 0.05594,    
              0.06106, 0.06606, 0.07122, 0.07629, 0.08150, 0.08651,    
              0.09170 )
obsy3    <- c(6.80, 7.15, 7.63, 5.85, 7.27, 5.97,
              7.21, 6.82, 7.50, 6.87, 6.66, 6.90,
              7.10, 6.23, 5.96, 6.14, 5.63, 5.50,
              5.41, 5.43, 5.43, 5.32, 5.33, 5.23,
              5.15 )
errobsy3 <- c(0.82, 1.06, 0.91, 0.89, 1.05, 0.64,
              0.84, 1.47, 1.02, 0.74, 0.74, 0.72,
              0.79, 0.37, 0.62, 0.23, 0.14, 0.16,
              0.14, 0.15, 0.14, 0.11, 0.12, 0.11,
              0.11 )

## DATA SET 4: krauss 1987
obsx4    <- c(0.02451, 0.02655, 0.02900, 0.03145, 0.03390, 0.03634, 0.03909,    
              0.04124, 0.04373, 0.04648, 0.04808, 0.04900, 0.04932, 0.05440,    
              0.05940, 0.06440, 0.06460, 0.06800, 0.06930, 0.07270, 0.07340,    
              0.07780, 0.07940, 0.08450, 0.08630, 0.08950, 0.09160, 0.09400,    
              0.09720, 0.10340, 0.10920, 0.11600, 0.12150, 0.13360, 0.14130,    
              0.14600, 0.15630, 0.15790, 0.16890, 0.17050, 0.19540, 0.21980,    
              0.24430, 0.26880, 0.29330, 0.31790, 0.34250 )
obsy4    <- c(5.07, 5.18, 5.23, 5.45, 5.26, 5.35, 5.77,
              5.03, 4.88, 4.98, 5.08, 5.06, 5.86, 5.71,
              5.10, 5.18, 5.56, 5.39, 5.93, 5.30, 5.55,
              5.27, 5.26, 5.12, 4.92, 5.31, 4.69, 4.86,
              4.97, 4.93, 4.77, 4.89, 4.67, 4.56, 4.62,
              4.97, 4.63, 4.56, 4.67, 4.73, 4.68, 4.35,
              4.57, 4.73, 5.09, 4.40, 4.41 )
errobsy4 <- c(1.34, 1.06, 0.58, 0.45, 0.52, 0.41, 0.35,
              0.43, 0.24, 0.26, 0.16, 0.19, 0.32, 0.32,
              0.36, 0.20, 0.23, 0.31, 0.17, 0.18, 0.25,
              0.20, 0.18, 0.18, 0.13, 0.30, 0.07, 0.08,
              0.08, 0.10, 0.16, 0.08, 0.08, 0.13, 0.09,
              0.10, 0.05, 0.08, 0.05, 0.05, 0.21, 0.22,
              0.19, 0.26, 0.28, 0.26, 0.24 )

## DATA SET 5: dawarakanath 1971
obsx5    <- c(0.088, 0.126,	0.155, 0.193, 0.234, 0.288, 0.338,	   
              0.379, 0.435,	0.488, 0.591, 0.691, 0.746, 0.792,	   
              0.895, 0.997, 1.081 )
obsy5    <- c(4.864, 4.882, 4.956, 4.513, 4.680, 4.449, 4.342,
              4.503, 4.214, 4.131, 3.901, 3.762, 3.699, 3.503,
              3.511, 3.661, 3.504 )
errobsy5 <- c(0.340, 0.340, 0.350, 0.316, 0.328, 0.311, 0.174,
              0.180, 0.169, 0.165, 0.156, 0.150, 0.148, 0.140,
              0.140, 0.146, 0.140 )

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
  yt1[i] <- (exp(2.429819509 * i.screen * (obsx1[i]^(-1.5)))) * 
            (alpha + beta * obsx1[i] + gamma * (obsx1[i]^2))
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
  yt2[i] <- (exp(2.429819509 * i.screen * (obsx2[i]^(-1.5)))) * 
            (alpha + beta * obsx2[i] + gamma * (obsx2[i]^2))
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
  yt3[i] <- (exp(2.429819509 * i.screen * (obsx3[i]^(-1.5)))) * 
            (alpha + beta * obsx3[i] + gamma * (obsx3[i]^2))
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
  yt4[i] <- (exp(2.429819509 * i.screen * (obsx4[i]^(-1.5)))) * 
            (alpha + beta * obsx4[i] + gamma * (obsx4[i]^2))
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
  yt5[i] <- (exp(2.429819509 * i.screen * (obsx5[i]^(-1.5)))) * 
          (alpha + beta * obsx5[i] + gamma * (obsx5[i]^2))
}    

###################
# PRIORS
###################

### polynomial parameters
  alpha ~ dnorm(0.0,pow(100, -2))T(0,)    #  alpha, i.e., S(E=0), cannot become negative!
  beta ~ dnorm(0.0,pow(100, -2))
  gamma ~ dnorm(0.0,pow(100, -2))
  i.screen ~ dnorm(0.0,pow(100, -2))T(0,)   # screening potential cannot become negative!

# extrinsic scatter in S-factor
  yscat1 ~ dnorm(0.0, pow(20, -2))T(0,)       # sigma assumed to be 20 MeVb
  yscat2 ~ dnorm(0.0, pow(20, -2))T(0,)
  yscat3 ~ dnorm(0.0, pow(20, -2))T(0,)
  yscat4 ~ dnorm(0.0, pow(20, -2))T(0,)
  yscat5 ~ dnorm(0.0, pow(20, -2))T(0,)

# systematic normalization factor for S-factor:
# log(): natural logarithm
y.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
logmu1 <- log(1.0)      # median of factor uncertainty is 1.0
logsigma1 <- log(1.038)  # factor uncertainty is 1.038, i.e. 3.8% for KUD04

y.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
logmu2 <- log(1.0)      # median of factor uncertainty is 1.0
logsigma2 <- log(1.057)   # factor uncertainty is <1.057, i.e., 5.7% for BON99
# we are using "8%", not "<8%"

y.norm3 ~ dlnorm(logmu3, pow(logsigma3, -2))
logmu3 <- log(1.0)      # median of factor uncertainty is 1.0
logsigma3 <- log(1.037)   # factor uncertainty is 1.037, i.e., 3.7% for JUN98

y.norm4 ~ dlnorm(logmu4, pow(logsigma4, -2))
logmu4 <- log(1.0)      # median of factor uncertainty is 1.0
logsigma4 <- log(1.045)   # factor uncertainty is 1.045, i.e., 4.5% for KRA87 

y.norm5 ~ dlnorm(logmu5, pow(logsigma5, -2))
logmu5 <- log(1.0)      # median of factor uncertainty is 1.0
logsigma5 <- log(1.082)   # factor uncertainty is 1.082, i.e., 8.2% for DAW71 

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
                 'obsx5' = obsx5, 'obsy5' = obsy5, 'errobsy5' = errobsy5),
                 inits = list(alpha = 1, beta = -1, gamma = 1, i.screen = 1e-6),
                 n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn) 
    
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                 variable.names=c('alpha', 'beta', 'gamma', 'i.screen',
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
pdf("MCMC_3He3He_a.pdf")
plot(mcmcChain)
dev.off()

######################################################################
# S-FACTOR FIT + DATA 
######################################################################
pdf("MCMC_3He3He_b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(8e-3,1.2)
yLim = c(2,30)

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
legend(0.5, 27, legend=c("Kud04", "Bon99", "Jun98", "Kra87", "Daw71"), 
         pch=c(1, 5, 0, 6, 2))
         
text(6e-2, 22, labels=expression(paste(NULL^"3","He(",NULL^"3","He,", "2p",")",NULL^"4","He")), 
                   cex=2.0)

text(1e-2, 20,  labels=expression(paste("total")), cex=1.4)
text(1e-2, 4.2, labels=expression(paste("bare")), cex=1.4)

# PLOT BELIEVABLE REGRESSION LINES

# define grid of 201 x values on a log scale for plotting of credible solutions;
# lseq is appropriate for log scale
xComb = lseq(xLim[1],xLim[2],length=201)

# next, pick number of samples to take into account for calculation 
# of y at each x; here we take 1001 samples into account
for ( i in round(seq(from=1,to=nsamp,length=1001)) ) {
  lines( xComb, samplesmat[i,"alpha"] + samplesmat[i,"beta"] * xComb 
                                      + samplesmat[i,"gamma"] * (xComb^2), 
         col=adjustcolor("black", alpha=0.01))
         
  lines( xComb, 
       exp(2.429819509 * samplesmat[i,"i.screen"] * (xComb^(-1.5))) * 
             (samplesmat[i,"alpha"] + samplesmat[i,"beta"] * xComb 
                                    + samplesmat[i,"gamma"] * (xComb^2)), 
         col=adjustcolor("black", alpha=0.01))
}

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

# add data KUD04 - circles
points( obsx1, obsy1, col="gray40", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="gray40" )  

# add data BON99 - diamonds
points( obsx2, obsy2, col="gray40", pch=5, cex=1.2 )
add.error.bars( obsx2, obsy2, 0.0, errobsy2, 0.0, col="gray40" )

# add data JUN98 - squares
points( obsx3, obsy3, col="gray40", pch=0, cex=1.2 )
add.error.bars( obsx3, obsy3, 0.0, errobsy3, 0.0, col="gray40" )

# add data KRA87 - triangles
points( obsx4, obsy4, col="gray40", pch=6, cex=1.2 )
add.error.bars( obsx4, obsy4, 0.0, errobsy4, 0.0, col="gray40" )

# add data DAW71 - triangles
points( obsx5, obsy5, col="gray40", pch=2, cex=1.2 )
add.error.bars( obsx5, obsy5, 0.0, errobsy5, 0.0, col="gray40" )

dev.off()

######################################################################
# S-FACTOR DATA 
######################################################################
pdf("MCMC_3He3He_c.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(8e-3,1.2)
yLim = c(2,30)

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
legend(0.5, 27, legend=c("Kud04", "Bon99", "Jun98", "Kra87", "Daw71"), 
      pch=c(1, 5, 0, 6, 2), col=c(68,91,258,"darkorchid","darkorange"))

text(6e-2, 20, labels=expression(paste(NULL^"3","He(",NULL^"3","He,", "2p",")",NULL^"4","He")), 
                   cex=2.0)

# add data KUD04 - circles
points( obsx1, obsy1, col=68, pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col=68 )  

# add data BON99 - diamonds
points( obsx2, obsy2, col=91, pch=5, cex=1.2 )
add.error.bars( obsx2, obsy2, 0.0, errobsy2, 0.0, col=91 )

# add data JUN98 - squares
points( obsx3, obsy3, col=258, pch=0, cex=1.2 )
add.error.bars( obsx3, obsy3, 0.0, errobsy3, 0.0, col=258 )

# add data KRA87 - triangles
points( obsx4, obsy4, col="darkorchid", pch=6, cex=1.2 )
add.error.bars( obsx4, obsy4, 0.0, errobsy4, 0.0, col="darkorchid" )

# add data DAW71 - triangles
points( obsx5, obsy5, col="darkorange", pch=2, cex=1.2 )
add.error.bars( obsx5, obsy5, 0.0, errobsy5, 0.0, col="darkorange" )

dev.off()

######################################################################
# POLYNOMIAL PARAMETERS 
######################################################################
# plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

pdf("MCMC_3He3He_d.pdf", width=10,height=5,onefile=F)
par(mfcol=c(1,3), mar=c(5.5,5.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
     las=1)

# plot S(0) -> alpha
plot(density(samplesmat[,1]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(S(0), "  (MeV b)")), line=3.0, cex.lab=2.0)
title(ylab="Probability density", line=3.5, cex.lab=2.0)
polygon(density(samplesmat[,1]), col=adjustcolor("blue", alpha=0.5))

# plot S'(0) -> beta
plot(density(samplesmat[,2]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(S,"'",(0), "  (b)")), line=3.0, cex.lab=2.0)
polygon(density(samplesmat[,2]), col=adjustcolor("blue", alpha=0.5))

# plot S''(0) -> 2 * gamma
plot(density(2*samplesmat[,3]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(S,"''",(0), "  (b/MeV)")), line=3.0, cex.lab=2.0)
polygon(density(2*samplesmat[,3]), col=adjustcolor("blue", alpha=0.5))

dev.off()

######################################################################
# POSTERIOR OF ELECTRON SCREENING POTENTIAL
######################################################################
pdf("MCMC_3He3He_e.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.5,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(100, 500)

dens <- density(1e6*samplesmat[,4])
a98 <- quantile(1e6*samplesmat[,4], prob = 0.975)

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
# DENSITIES OF S-FACTOR NORMALIZATION FACTORS
######################################################################
pdf("MCMC_3He3He_f.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(0.8, 1.15)

# plot density in first panel           
plot(density(samplesmat[,8]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, 
     xlim=xLim, yaxs='i', xaxs='i')
title(ylab="Probability density", line=3.5, cex.lab=2.3)
title(xlab="S-factor normalization", line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,5]),  
     col=adjustcolor("blue", alpha=0.5))
polygon(density(samplesmat[,6]), 
     col=adjustcolor("black", alpha=0.5))
polygon(density(samplesmat[,7]),  
     col=adjustcolor("red", alpha=0.5))
polygon(density(samplesmat[,8]),  
     col=adjustcolor("green", alpha=0.5))
polygon(density(samplesmat[,9]),  
     col=adjustcolor("violet", alpha=0.5))

legend("topleft", inset=.01, 
   c("Kud04", "Bon99", "Jun98", "Kra87", "Daw71"), 
   fill=adjustcolor(c("blue", "black", "red", "green", "violet"), alpha=0.5), 
   horiz=FALSE, cex=1.5, box.lty=0)

dev.off()

######################################################################
# POSTERIOR EXTRINSIC S-FACTOR SCATTER
######################################################################
pdf("MCMC_3He3He_g.pdf",width=10, height=6, onefile=F)
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
par(mfcol=c(2,3), mar=c(5.5,7.0,1.0,2.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# plot #1          
plot(density(samplesmat[,10]), main="", xlab="", ylab="", xlim=c(0, 1.0),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
#title(ylab="density", line=2.5, cex.lab=2.3)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3)
legend("topright", legend="Kud04", pch=NA, cex=1.5)
polygon(density(samplesmat[,10]), col=adjustcolor("blue", alpha=0.5))
 
# plot #2                  
plot(density(samplesmat[,11]), main="", xlab="", ylab="", xlim=c(0, 7),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')
#title(ylab="Probability", line=2.5, cex.lab=2.3)
legend("topright", legend="Bon99", pch=NA, cex=1.5)
polygon(density(samplesmat[,11]), col=adjustcolor("blue", alpha=0.5))

# plot #3                  
plot(density(samplesmat[,12]), main="", xlab="", ylab="", xlim=c(0, 0.5),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')
legend("topright", legend="Jun98", pch=NA, cex=1.5)
polygon(density(samplesmat[,12]), col=adjustcolor("blue", alpha=0.5))

# plot #4                  
plot(density(samplesmat[,13]), main="", xlab="", ylab="", xlim=c(0, 0.7),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')
legend("topright", legend="Kra87", pch=NA, cex=1.5)
polygon(density(samplesmat[,13]), col=adjustcolor("blue", alpha=0.5))

# plot #5                  
plot(density(samplesmat[,14]), main="", xlab="", ylab="", xlim=c(0, 0.5),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
     )
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
     yaxs='i', xaxs='i')
legend("topright", legend="Daw71", pch=NA, cex=1.5)
polygon(density(samplesmat[,14]), col=adjustcolor("blue", alpha=0.5))

dev.off()

######################################################################
# CORRELATION PLOT
######################################################################
pdf("MCMC_3He3He_h.pdf", width=6, height=6, onefile=F)
pairs(~alpha+beta+gamma, col=adjustcolor("red", alpha=0.5),  
   data=samplesmat2[sample(nrow(samplesmat2), size=1000, replace=FALSE),], 
   main="Simple Scatterplot Matrix")

dev.off()



