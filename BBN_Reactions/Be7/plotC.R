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
# Read Dataset
######################################################################

temp <- read.table("data",header=T)
obsx <- temp[,1]
obsy <- temp[,2]  

######################################################################
# Read MCMC Chain
######################################################################

mcmcChain <- read.table("test",header = T)




######################################################################
# TRACES AND DENSITIES
######################################################################
#pdf("test_a.pdf")
#plot(mcmcChain)
#dev.off()

######################################################################
# PLOT DATA+ FIT
######################################################################
pdf("test_b.pdf")
par(mfcol=c(1,1), mar=c(7,7,7,7), oma=c(0,0,0,0), tck=0.02, 
    las=1)

# first determine plot ranges
xRang = max(obsx)-min(obsx)
yRang = max(obsy)-min(obsy)
limMult = 0.25
xLim= c( min(obsx)-limMult*xRang , max(obsx)+limMult*xRang )
yLim= c( min(obsy)-limMult*yRang , max(obsy)+limMult*yRang )

plot( obsx , obsy , cex=1.5 , lwd=2 , col="black" , xlim=xLim , 
      ylim=yLim, xlab="" , ylab="" , cex.lab=1.5 ,
      main= "" , cex.main=1.33, axes=FALSE )  
title(xlab="x", line=2, cex.lab=1.3)
title(ylab=expression("y'"), line=2, cex.lab=1.3)
magaxis(mgp=c(0,0.2,0))
box()

# PLOT BELIEVABLE REGRESSION LINES
# first, set up vector with 201 x values; 
# seq(from, to, nlength) generates a sequence of "length" values between 
# "from" and "to"

# matrix samplesmat contains all accepted samples from all n.chains with 
# rows: i, alpha, beta, sigma
samplesmat = mcmcChain
nsamp = nrow(samplesmat)

xComb = seq(xLim[1],xLim[2],length=201)
# next, determine number of samples to take into account [here 1000] for 
# calculation of y at each x; here we take 1000 samples into account
for ( i in round(seq(from=1,to=nsamp,length=1000)) ) {
  lines( xComb , 
         samplesmat[i,"alpha"] + samplesmat[i,"beta"]*xComb, 
         col=adjustcolor("black", alpha=0.01))
}

# PLOT CREDIBLE REGION [red solid lines]
# the credible region is obtained by computing regression lines using
# all sampled slopes and intercepts; since we use all samples, regardless
# of the scatter intrscat, we are effectively marginalizing over this
# parameter; from the sampled regression lines, we can compute y values
# at each point of our x-grid; we can thus define suitable quantiles
# at each x value; this defines the credible region for the regression;
#
# at any given x value, the corresponding true y value most likely is
# located between the two red solid lines

# select first two colums [alpha and beta] from matrix
samplesmat2 <- as.matrix(samplesmat[,1:2])

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
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.025), lty=1, lwd=2, col='red')
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.975), lty=1, lwd=2, col='red')
# interpretation: we are 95% sure that the true regression line is within the
# credible interval

# PLOT MEDIAN REGRESSION LINE [blue solid line]
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.50), lty=1, lwd=2, col='blue')

# PLOT MEAN INTRINSIC SCATTER LINES [dashed lines]
# if we aquire new data under the same conditions, they will most likely be 
# located between the two blue dashed lines

lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.50) + mean(samplesmat[,3]), 
      lty=2, lwd=2, col='black')
lines(newdat$xx, apply(fitmat, 1, quantile, prob=0.50) - mean(samplesmat[,3]), 
      lty=2, lwd=2, col='black')

dev.off()

######################################################################
# PLOT CORRELATION FOR GIVEN PAIR WITH CREDIBLE CONTOURS
######################################################################
pdf("test_c.pdf")

plot(samplesmat2, main="", 
     cex.lab=1.1, cex.axis=1.0, cex.main=1.0, 
     xlab="alpha", ylab="beta", pch=".", 
     col=adjustcolor("black", alpha=0.1), 
     xlim=c(4,16), ylim=c(0.8,1.2))

# now add credible regions in 2D scatter plot of parameters;
# finding the highest density ellipse corresponds to finding the smallest 
# ellipse containing p% of the probability a.k.a. the minimum volume ellipse; 

library(emdbook)
HPDregionplot(samplesmat, prob = c(0.95, 0.75, 0.50), 
              col=c("blue", "green", "red"), lwd=3, add=TRUE)
legend("topleft", legend = c("95%", "75%", "50%"), 
       col = c("blue", "green", "red"), lty=c(1,1,1), lwd=c(3,3,3))

#points(sigma.clus.sim, vcentsim, pch=3, cex=3)
dev.off()

######################################################################
# PREDICTION 
######################################################################
# plot posterior predictive distribution at a given x by marginalization 
# over all parameters; we will use all sampled regression lines [i.e., 
# slope and intercept samples], calculate y-values for all slope/intercept 
# samples at given x; this set of y-values represents our predicted 
# posterior at x
pdf("test_d.pdf")

# define x value for which we would like to predict y
xchoice <- 0.02

# declare vector with y values
fitvec <- vector(mode = "numeric", length = nsamp)   

# calculate y values from all regression samples at given x value
for(i in 1:nsamp) fitvec[i] <- samplesmat2[i,1] + xchoice * samplesmat2[i,2]

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
