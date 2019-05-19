######################################################################
#
# test.R
#
# purpose: 
# - regression model for x and y data
#
# conditions:
# - x data have no error; only y data have errors
# - linear relationship between x and y [could be changed]
# - user simulated data
# - includes effect of systematic error
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package 
library("rjags")
## for block updating [we do not need to center predictor variables]
load.module("glm")  
# random number seed
#set.seed(123)

######################################################################
# GENERATE ARTIFICIAL DATA
######################################################################
# sample data input/generation; the input is of the form: obsx, obsy, 
# errobsy, where the latter is the individual error of each datum [i]

N <- 10
## DATA SET
# x values are sampled uniformly
obsx <- runif(N, 0, 10)

# y values:
# - mean=linear model of slope=0.5 and intercept=0
# - y errors are drawn from a uniform distribution with sd=individual 
#   errors
# - finally, option to shift all data by a normalization factor [for
#   offset, change '*' into '+'

errobsy <- runif(N, 0.2, 0.5)
obsy.hat <- rnorm(N, 0.5 * obsx + 0, errobsy)
# now shift everything by known amount
obsy <- obsy.hat * 1.1     # *true* systematic shift 

# suppose the *estimated* systematic error on each data point is 20%

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

for (i in 1:length(obsx)) {
 obsy[i] ~ dnorm(y[i], pow(errobsy[i], -2))
 # ...subject to syst uncertainties:
 y[i] <- y.norm * z[i]    
 # true value [only mother nature knows]
 z[i] <- alpha + beta * obsx[i]         
}

###################
# PRIORS
###################

alpha ~ dnorm(0.0, pow(100, -2))
beta ~ dnorm(0.0, pow(100, -2))

# systematic normalization factor for S-factor:
# log(): natural logarithm
y.norm ~ dlnorm(logmu, pow(logsigma, -2))
logmu <- log(1.0)      # median of factor uncertainty is 1.0
logsigma <- log(1.2)   # factor uncertainty is 1.2, i.e. 20%

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
n.iter   <- 25000   
thin   <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model;
#
ourmodel <- jags.model(f, data = list(         ## jags wants all data in a list
               'obsx' = obsx, 'obsy' = obsy, 'errobsy' = errobsy
                                     ),
               n.chains = n.chains, n.adapt = n.adapt)

update(ourmodel, n.burn) 
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                  variable.names=c('alpha', 'beta', 'y.norm'), 
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
pdf("test_a.pdf")
plot(mcmcChain)
dev.off()
######################################################################
# PLOT DATA
######################################################################
pdf("test_b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# determine plot ranges
xRang = max(obsx)-min(obsx)
yRang = max(obsy)-min(obsy)
limMult = 0.25
xLim= c( -5 , 15 )
yLim= c( -5 , 15 )

plot( obsx , obsy , cex=1.5 , lwd=2 , col="black" , xlim=xLim , 
      ylim=yLim, xlab="X" , ylab="Y" , cex.lab=1.5 ,
      main="Data with credible regression lines" , cex.main=1.33 )
      
# add error bars;
# w is the width of perpendicular end bars on errors; set to zero
add.error.bars <- function(X,Y,dX,dY,w,col=1){
X0 = X; 
Y0 = (Y-dY); 
X1 =X; 
Y1 = (Y+dY);
arrows(X0, Y0, X1, Y1, lwd=3, code=3, angle=90, length=w, col="black");
Y0 = Y; 
X0 = (X-dX); 
Y1 =Y; 
X1 = (X+dX);
arrows(X0, Y0, X1, Y1, lwd=3, code=3, angle=90, length=w, col="black");
}
add.error.bars(obsx, obsy, 0.0, errobsy, 0.0)   

# PLOT BELIEVABLE REGRESSION LINES
# first, set up vector with 201 x values; 
# seq(from, to, nlength) generates a sequence of "length" values between 
# "from" and "to"

# matrix samplesmat contains the samples from all n.chains with rows:
# i, alpha, beta,...
samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)

xComb = seq(xLim[1],xLim[2],length=201)
# next, determine number of samples to take into account [here 1000] for 
# calculation of y at each x; here we take 1000 samples into account
for ( i in round(seq(from=1,to=nsamp,length=1000)) ) {
  lines( xComb , 
# change_change_change_change_change_change_change_change_change_change
         samplesmat[i,"alpha"] + samplesmat[i,"beta"] * xComb, 
# change_change_change_change_change_change_change_change_change_change
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
 
# select first two colums [alpha and beta] from matrix
samplesmat2 <- samplesmat[,1:2]

# define xx values for which we would like to have the fitted values; 
# save in new data frame newdat; again, we chose 201 xx values
newdat <- data.frame(xx=seq(xLim[1],xLim[2],length=201))

# create new model matrix that contains new xx values
newmodmat <- model.matrix(~xx,data=newdat)

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
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.16), lty=1, lwd=1, col='red')
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.84), lty=1, lwd=1, col='red')
# interpretation: we are 68% sure that the true regression line is within the
# credible interval

# PLOT MEDIAN REGRESSION LINE [blue solid line]
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.50), lty=1, lwd=1, col='blue')
# PLOT MEAN INTRINSIC SCATTER LINES [dashed lines]
# if we aquire new data under the same conditions, they will most likely be 
# located between the two blue dashed lines

#lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.50) + mean(samplesmat[,3]), 
#       lty=2, lwd=2, col='black')
#lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.50) - mean(samplesmat[,3]), 
#       lty=2, lwd=2, col='black')

dev.off()
