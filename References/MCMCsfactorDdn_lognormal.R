######################################################################
# MCMCsfactorDdn_lognormal.R
#
# in this code we have a lognormal likelihood,
# where the lognormal parameters mu and sigma are given 
# by Eq. (27) of Longland et al., NPA 841, 1 (2010):
#
# obsy1[i] ~ dlnorm(yl1[i],pow(corrl.err1[i], -2))                
# yl1[i] <- log(y1[i])-0.5*log(1+pow(corr.err1[i],2)/pow(y1[i],2)) ## lognormal mu
# corrl.err1[i] <- sqrt(log(1+pow(corr.err1[i],2)/pow(y1[i],2)))   ## lognormal sigma
#
#
# conditions:
# - relationship between x and y given by theory [table] 
# - theory relation from av8p3b-ddn.dat
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
#set.seed(123)
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i]

##We will work with energy is in MeV and S-factor in MeVb 
##(need to convert keVb to MeVb of data from the tables)

########################
### d(d,n)3He
########################
## DATA SET 1: Leo06
tab1 <- read.table("dd_n_master_Leo06.dat", header=FALSE)

obsx1 <- c(tab1[,1])

obsy1 <- c(tab1[,3]) * 1e-3

errobsy1 <- c(tab1[,4]) * 1e-3
########################
## DATA SET 2: Gre95
tab2 <- read.table("dd_n_master_Gre95.dat", header=FALSE)

obsx2 <- c(tab2[,1])

obsy2 <- c(tab2[,3]) * 1e-3

errobsy2 <- c(tab2[,4]) * 1e-3
########################
## DATA SET 3: Bro90
tab3 <- read.table("dd_n_master_Bro90.dat", header=FALSE)

obsx3 <- c(tab3[,1])

obsy3 <- c(tab3[,3]) * 1e-3

errobsy3 <- c(tab3[,4]) * 1e-3
########################
## DATA SET 4: Kra87 (B)
tab4 <- read.table("dd_n_master_Kra87B.dat", header=FALSE)

obsx4 <- c(tab4[,1])

obsy4 <- c(tab4[,3]) * 1e-3

errobsy4 <- c(tab4[,4]) * 1e-3
########################
## DATA SET 5: Kra87 (M)
tab5 <- read.table("dd_n_master_Kra87M.dat", header=FALSE)

obsx5 <- c(tab5[,1])

obsy5 <- c(tab5[,3]) * 1e-3

errobsy5 <- c(tab5[,4]) * 1e-3
######################################################################
# input of theoretical S-factor from Arai et al.(2011); E is in MeV, S is in keVb;
# convert latter to MeV b
theory <- read.table("av8p3b-ddn.dat", header=FALSE)

interp.x <- theory[,1]
interp.y <- theory[,2] * 1e-3

# we will be using JAGS interp.lin function to use this theoretical S-factor:
#
# - the columns of this table define vectors x and y
# - a single point is given by x_i, y_i
# - interp.lin gives the y value for the x value provided as argument e,
#   interp.lin(e, x, y)

######################################################################
# set up Bayesian model [with parameters] 
#
# import jags package 
library('rjags')
## for block updating [we do not need to center predictor variables]
load.module("glm")  
# 
cat('model {
# modeling of y variable:
# - first, we assume a [theoretical] model dependence between obsy and z, where
#   z is the true y-value; the scaling, a.scale, is a model parameter [it is
#   the same for all data sets since we are seeking the best fit to all data]
# - second, a new variable, y, is calculated by multiplying z with "n.norm",
#   i.e., the systematic normalization factor estimate affecting all points of 
#   a given data set equally
# - third, observed y value is drawn from a Gaussian with mean=y[i] and
#   sd=individual [i.e., statistical] errors
#
# 1. LIKELIHOOD
    # !!! in a for loop, make sure **all** variables on the LEFT of an 
    # !!! expression have the index [i]
    # yl1[] and corrl.err1[] are lognormal parameters (mu and sigma)
  for (i in 1:length(obsx1)) {
   obsy1[i] ~ dlnorm(yl1[i],pow(corrl.err1[i], -2))                 
   yl1[i] <- log(y1[i])-0.5*log(1+(pow(corr.err1[i],2)/pow(y1[i],2))) ## lognormal mu
   corrl.err1[i] <- sqrt(log(1+(pow(corr.err1[i],2)/pow(y1[i],2))))   ## lognormal sigma
   p.alt1[i] ~ dcat(p1[])
   corr.err1[i] <- errobsy1[i]*phi[p.alt1[i]]
   y1[i] <- n.norm1 * z1[i]    # systematic error as normalization factor
  ### y1[i] <- n.norm1 + z1[i]    # systematic error as offset  
# 2. REALTIONSHIP BETWEEN TRUE VARIABLES
   z1[i] <- a.scale * interp.lin(obsx1[i], interp.x, interp.y)     
  }

# 1. LIKELIHOOD
  for (i in 1:length(obsx2)) {
   obsy2[i] ~ dlnorm(yl2[i],pow(corrl.err2[i], -2)) 
   yl2[i] <- log(y2[i])-0.5*log(1+(pow(corr.err2[i],2)/pow(y2[i],2)))
   corrl.err2[i] <- sqrt(log(1+(pow(corr.err2[i],2)/pow(y2[i],2))))
   p.alt2[i] ~ dcat(p2[])
   corr.err2[i] <- errobsy2[i]*phi[p.alt2[i]]
   y2[i] <- n.norm2 * z2[i]    # systematic error as normalization factor
  ### y2[i] <- n.norm2 + z2[i]    # systematic error as offset  
# 2. REALTIONSHIP BETWEEN TRUE VARIABLES
   z2[i] <- a.scale * interp.lin(obsx2[i], interp.x, interp.y)
  }

# 1. LIKELIHOOD
  for (i in 1:length(obsx3)) {
   obsy3[i] ~ dlnorm(yl3[i],pow(corrl.err3[i], -2)) 
   yl3[i] <- log(y3[i])-0.5*log(1+(pow(corr.err3[i],2)/pow(y3[i],2)))
   corrl.err3[i] <- sqrt(log(1+(pow(corr.err3[i],2)/pow(y3[i],2))))
   p.alt3[i] ~ dcat(p3[])
   corr.err3[i] <- errobsy3[i]*phi[p.alt3[i]]
   y3[i] <- n.norm3 * z3[i]    # systematic error as normalization factor
  ### y3[i] <- n.norm3 + z3[i]    # systematic error as offset  
# 2. REALTIONSHIP BETWEEN TRUE VARIABLES
   z3[i] <- a.scale * interp.lin(obsx3[i], interp.x, interp.y)
  }

# 1. LIKELIHOOD
  for (i in 1:length(obsx4)) {
   obsy4[i] ~ dlnorm(yl4[i],pow(corrl.err4[i], -2)) 
   yl4[i] <- log(y4[i])-0.5*log(1+(pow(corr.err4[i],2)/pow(y4[i],2)))
   corrl.err4[i] <- sqrt(log(1+(pow(corr.err4[i],2)/pow(y4[i],2))))
   p.alt4[i] ~ dcat(p4[])
   corr.err4[i] <- errobsy4[i]*phi[p.alt4[i]]
   y4[i] <- n.norm4 * z4[i]    # systematic error as normalization factor
  ### y4[i] <- n.norm4 + z4[i]    # systematic error as offset  
# 2. REALTIONSHIP BETWEEN TRUE VARIABLES
   z4[i] <- a.scale * interp.lin(obsx4[i], interp.x, interp.y)
  }

# 1. LIKELIHOOD
  for (i in 1:length(obsx5)) {
    obsy5[i] ~ dlnorm(yl5[i],pow(corrl.err5[i], -2)) 
    yl5[i] <- log(y5[i])-0.5*log(1+(pow(corr.err5[i],2)/pow(y5[i],2)))
    corrl.err5[i] <- sqrt(log(1+(pow(corr.err5[i],2)/pow(y5[i],2))))
    p.alt5[i] ~ dcat(p5[])
    corr.err5[i] <- errobsy5[i]*phi[p.alt5[i]]
    y5[i] <- n.norm5 * z5[i]    # systematic error as normalization factor
    ### y5[i] <- n.norm5 + z5[i]    # systematic error as offset  
    # 2. REALTIONSHIP BETWEEN TRUE VARIABLES
    z5[i] <- a.scale * interp.lin(obsx5[i], interp.x, interp.y)
    }

# 3. PRIORS
  ### a.scale ~ dnorm(0.0, pow(2, -2))T(0,)    #  "a.scale" cannot become negative!
  a.scale ~ dnorm(0.0, pow(100, -2))T(0,)    #  "a.scale" cannot become negative!
  
  # systematic errors as normalization factors
  # [see Table II in Coc et al., PRD 92, 123526 (2015)] 
  ### lognormal density:
  n.norm1 ~ dlnorm(logmu1, pow(logsigma1, -2))
  logmu1 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma1 <- log(1.02)  # factor uncertainty is 1.02, i.e. 2% for Leo06
  
  n.norm2 ~ dlnorm(logmu2, pow(logsigma2, -2))
  logmu2 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma2 <- log(1.03)  # factor uncertainty is 1.03, i.e., 3% for Gre95

  n.norm3 ~ dlnorm(logmu3, pow(logsigma3, -2))
  logmu3 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma3 <- log(1.013)  # factor uncertainty is 1.013, i.e., 1.3% for Bro90
  
  n.norm4 ~ dlnorm(logmu4, pow(logsigma4, -2))
  logmu4 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma4 <- log(1.064) # factor uncertainty is 1.064, i.e., 6.4% for Kra87 (B) 

  n.norm5 ~ dlnorm(logmu5, pow(logsigma5, -2))
  logmu5 <- log(1.0)      # median of factor uncertainty is 1.0
  logsigma5 <- log(1.082) # factor uncertainty is 1.082, i.e., 8.2% for Kra87 (M) 
  
  # for case that measured errors are correct:
  phi[1] <- 1
  # for case that measured errors are too optimistic:
  phi[2] ~ dunif(1,50)
  
  p1[1] ~ dunif(0,1)
  p1[2] <- 1-p1[1]
  
  p2[1] ~ dunif(0,1)
  p2[2] <- 1-p2[1]
  
  p3[1] ~ dunif(0,1)
  p3[2] <- 1-p3[1]
  
  p4[1] ~ dunif(0,1)
  p4[2] <- 1-p4[1]

  p5[1] ~ dunif(0,1)
  p5[2] <- 1-p5[1]

}', file={f <- tempfile()})

######################################################################
# n.adapt:  number of iterations in the chain for adaptation (n.adapt) 
#           [JAGS will use to choose the sampler and to assure optimum 
#           mixing of the MCMC chain; will be discarded] 
# n.udpate: number of iterations for burnin; these will be discarded to 
#           allow the chain to converge before iterations are stored
# n.iter:   number of iterations to store in the final chain as samples 
#           from the posterior distribution 
# n.chains: number of mcmc chains
# n.thin:   store every n.thin element [=1 keeps all samples]

n.adapt  <- 1000  
n.update <- 2000  
n.iter   <- 7500  
n.chains <- 3
n.thin   <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model;
#
ourmodel <- jags.model(f,
                data = list('obsx1' = obsx1, ## jags wants all data in a list
                            'obsy1' = obsy1,
                            'errobsy1' = errobsy1,
                            'obsx2' = obsx2, 
                            'obsy2' = obsy2,
                            'errobsy2' = errobsy2,
                            'obsx3' = obsx3, 
                            'obsy3' = obsy3,
                            'errobsy3' = errobsy3,
                            'obsx4' = obsx4, 
                            'obsy4' = obsy4,
                            'errobsy4' = errobsy4,
                            'obsx5' = obsx5, 
                            'obsy5' = obsy5,
                            'errobsy5' = errobsy5,
                            'interp.x' = interp.x,
                            'interp.y' = interp.y),
                   n.chains = n.chains,
                   n.adapt = n.adapt)
# burnin 
update(ourmodel, n.update, progress.bar="none") 
    
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                        variable.names=c('a.scale',
                        'n.norm1', 'n.norm2', 'n.norm3', 'n.norm4', 'n.norm5',
                        'p.alt1', 'p.alt2', 'p.alt3', 'p.alt4', 'p.alt5'), 
                        n.iter=n.iter, n.thin=n.thin)
######################################################################
# output results on screen
cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n") 
# output
cat("REJECTION RATES:", "\n")  
cat("", "\n")
show(rejectionRate(mcmcChain[,1:6]))
cat("-------------------------------------------", "\n")
cat("GELMAN DIAGNOSTICS:", "\n")
cat("", "\n")
show(gelman.diag(mcmcChain[,1:6]))
cat("-------------------------------------------", "\n")
cat("HEIDELBERGER DIAGNOSTIC:", "\n")
cat("", "\n")
show(heidel.diag(mcmcChain[,1:6]))
cat("-------------------------------------------", "\n")
######################################################################
######################################################################
# PLOT TRACES, DENSITIES, AND GELMAN
pdf("MCMCsfactorDdn_a.pdf")
plot(mcmcChain[,1:6])
dev.off()

pdf("MCMCsfactorDdn_b.pdf")
## traceplot(mcmcChain)
## densplot(mcmcChain)
gelman.plot(mcmcChain[,1:6])
dev.off()
######################################################################
# PLOT S-FACTOR GRAPH
## make single-panel graph;
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

# first determine plot ranges
pdf("MCMCsfactorDdn_c.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(5e-3,0.5)
yLim = c(4e-2,2e-1)
#####################################
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
legend(0.18, 0.08, legend=c("Leo06", "Gre95", "Bro90", "Kra87(B)", "Kra87(M)"), pch=c(1, 5, 0, 6, 2), col=c(68,91,258,"darkorchid","darkorange"))

text(0.01, 0.17, labels=expression(paste("d(d,n)",NULL^"3","He")), cex=1.3)
######################################################################
# PLOT BELIEVABLE S-FACTORS

# matrix samplesmat contains the samples from all n.chains with rows:
# i, a.scale,...
# !!! make sure to check the order of the parameters in the MCMC output !!!
samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)

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
######################################################################
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
capture.output( print(samplesmat2, print.gap=3), file="MCMCsamplesDdn")

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

######################################################################
# PLOT MEDIAN REGRESSION LINE [blue solid line]
lines(xComb, apply(fitmat, 1, quantile, prob=0.50), lty=1, lwd=2, 
          col=adjustcolor("blue", alpha=0.5))
######################################################################
# function for error bars;
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

# add data Leo06 - circles
points( obsx1, obsy1, col=68, pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="black" )  

# add data Gre95 - diamonds
points( obsx2, obsy2, col=91, pch=5, cex=1.2 )
add.error.bars( obsx2, obsy2, 0.0, errobsy2, 0.0, col="black" )

# add data Bro90 - squares
points( obsx3, obsy3, col=258, pch=0, cex=1.2 )
add.error.bars( obsx3, obsy3, 0.0, errobsy3, 0.0, col="black" )

# add data Kra87(B) - triangles down
points( obsx4, obsy4, col="darkorchid", pch=6, cex=1.2 )
add.error.bars( obsx4, obsy4, 0.0, errobsy4, 0.0, col="black" )

# add data Kra87(M) - triangles up
points( obsx5, obsy5, col="darkorange", pch=2, cex=1.2 )
add.error.bars( obsx5, obsy5, 0.0, errobsy5, 0.0, col="black" )


######################################################################
# MARK DATA POINTS THAT ARE OUTLIERS

samplesmat3 <- samplesmat[,-c(1:6)]

# subtract 1 from each element (so that 1,2 is transformed to 0,1)
samplesmat4 <- samplesmat3 - 1

n.length <- length(obsx1)+length(obsx2)+length(obsx3)+length(obsx4)+length(obsx5)

prob <- vector(, n.length)
# sum across rows and divide by total number of rows (samples)
for (i in 1:n.length){
   prob[i] <- sum(samplesmat4[,i])/nrow(samplesmat4)
}

# output
cat(" PROBABILITY OF OUTLIERS:     ", "\n")
print(prob)

dev.off()
######################################################################
# PREDICTION: plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

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
pdf("MCMCsfactorDdn_e.pdf")
plot(density(fitvec))
dev.off()
######################################################################
######################################################################


