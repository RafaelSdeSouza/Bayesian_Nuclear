# He3dp analysis
#
# purpose: ARTIFICIAL DATA
#
# - 3 parameters are assumed: Er, gamma_d^2, gamma_n^2 [e1, gin, gout]
#
# - uses the function sfactorTdn_fast(obsx1[i], e1, gin, gout), which
#   is a C++ version of a Fortran code that includes Coulomb wave
#   function calculations; JAGS has been recompiled with this C++ function 
#
# Fortran code tdn_plot.f is needed in this R script for plotting the
# S-factor only
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
#set.seed(123)
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical 
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb; 
######################################################################
## ARTIFICIAL DATA GENERATION 

N <- 30

obsx1 <- exp(runif(N,log(0.001),log(1)))

res <- vector()
obsy1 <- vector()
errobsy1 <- vector()

# Barker values:
# Er  = 0.0912 MeV
# g^2_in = 2.93 MeV         ! reduced width of deuteron
# g^2_out = 0.0794 MeV      ! reduced width of neutron

#res[2] <-  0.35779  # resonance energy
res[2] <-  0.35  # resonance energy
res[3] <-  1.0085    # reduced width incoming
res[4] <- 0.025425   # reduced width outgoing

for (i in 1:length(obsx1)){
   res[1] <- obsx1[i]
   write.table(res, file="He3dp_AD.in", quote=TRUE, 
                        row.names=FALSE, col.names=FALSE)

   # Load the fortran code needed to calculate S-factor curve
   if(!is.loaded("He3dp_AD_Sub"))
   dyn.load("He3dp_AD.so") 
   .Fortran("He3dp_AD_Sub")

   tab1 <- read.table("He3dp_AD.out", header=FALSE)

   errobsy1[i] <- 0.1
   obsy1[i] <- rnorm( 1, tab1[1,2], errobsy1[i] )
}

######################################################################                 
# import jags package 
library(rjags)
library(runjags)
library(R2jags)
library(mcmcplots)
## for block updating [we do not need to center predictor variables]
#load.module("glm")  
load.module("nuclear")  
# 
######################################################################                 
cat('model {

# LIKELIHOOD
for (i in 1:length(obsx1)) {
  obsy1[i] ~ dnorm(y1[i], pow(tau, -2))
  y1[i] <- sfactor3Hedp(obsx1[i], e1, gin, gout)    
}    
# PRIORS
# e1, gin, gout are defined as in tdn.f (by Alain Coc):
# resonance energy, initial reduced width, final reduced 
# width;

## Physical priors:
  # Hyperpriors

##
   
   tau ~dgamma(1e-3,1e-3)
   e1 ~ dunif(0.0, 10.0)
   gin ~ dnorm(0, pow(10, -2))T(0,)
    gout ~ dnorm(0, pow(10, -2))T(0,)



#  gin ~ dunif(0, 100)
#  gout ~ dunif(0, 100)
#  gin ~ dunif(0, 10)
#  gout ~ dunif(gin, 10)

#  gin ~ dgamma(0.5,100)
#  gout ~ dgamma(0.5,100)
  
#  gin ~ dgamma(0.5,rl1)
#  gout ~ dgamma(0.5,rl2)
#  rl1 ~ dunif(0,500)
#  rl2 ~ dunif(0,500)
  
#  gin ~ dchisqr(1)T(gout,)
#  gin ~ dchisqr(1)
#  gout ~ dchisqr(1)

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

n.adapt  <- 5000   
n.update <- 1000  
n.iter   <- 10000  
n.chains <- 3
n.thin   <- 10
inits <- function () { list(e1 = runif(1,0,0.3),gin=runif(1,2,10),gout=runif(1,0.01,1)) }
# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 



# JAGS model with R2Jags;
out <- jags(data = list('obsx1' = obsx1, ## jags wants all data in a list
                        'obsy1' = obsy1),
              inits = inits,
              parameters = c("e1", "gin", "gout"),
              model.file = f,
              n.thin = 50,
              n.chains = 3,
              n.burnin = 7500,
              n.iter = 15000)
denplot(out)
mcmcplot(out)
traplot(out)

# JAGS model with runjags;
results <- run.jags(model=f,  data = list('obsx1' = obsx1, ## jags wants all data in a list
                                          'obsy1' = obsy1,
                                          'errobsy1' = errobsy1),
                    monitor=c("e1", "gin", "gout"), n.chains=4,
                    thin=10,adapt=10000,
                    inits=list(inits1=inits(),inits2=inits(),inits3=inits(),inits4=inits()))

mcmcChain <- as.mcmc(out)[,-1]
# JAGS model;
#
ourmodel <- jags.model(f,
                data = list('obsx1' = obsx1, ## jags wants all data in a list
                            'obsy1' = obsy1,
                            'errobsy1' = errobsy1),
                    inits = list(e1 = 1, gin = 6, gout = 0.1),
                    n.chains = n.chains,
                    n.adapt = n.adapt)
# burnin 
update(ourmodel, n.update, progress.bar="none") 

# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                    variable.names=c('e1', 'gin', 'gout' 
#                    ,'rl1','rl2'
#                    ,'s1', 's2', 'r1', 'r2'
                    ), 
                    n.iter=n.iter, n.thin=n.thin)
######################################################################
# output results on screen
cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n") 
# PLOT TRACES AND DENSITIES
pdf("..//figures/MCMCsfactorTdn_a.pdf")
plot(mcmcChain)
dev.off()
pdf("MCMCsfactorTdn_b.pdf")
## traceplot(mcmcChain)
## densplot(mcmcChain)
gelman.plot(mcmcChain)
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
# first determine plot ranges
pdf("MCMCHe3dp.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
    las=1)

xLim = c(2e-3,1)
yLim = c(0,30)
######################################################################
# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
      axes=FALSE, main="", xlab = "", ylab = "",
      cex=1.5 , cex.lab=1.3, cex.axis=1.0,log="x",
      cex.main=1.0 )
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
##legend(0.12, 3e-7, legend=c("MA97", "SCH97", "CAS02", "BYS08"), pch=c(1, 0, 6, 5))

text(0.04, 15, labels=expression(paste(NULL^"3","He(d,p)",NULL^"4","He")), cex=1.3)
######################################################################
# PLOT BELIEVABLE S-FACTORS

# matrix samplesmat contains the samples from all n.chains with rows:
# i, a.scale,...
# !!! make sure to check the order of the parameters in the MCMC output !!!
samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)

# define grid of 201 x values on a log scale for plotting of credible solutions;
# lseq is appropriate for log scale
# xComb = lseq(xLim[1],xLim[2],length=201)

# calculate for these energies and the set of Bayesian samples for
# Er, g^2_i, g^2_f the S_factor curve using Fortran code

for ( i in round(seq(from=1,to=nsamp,length=1000)) ) {
  # output samples to file for S-factor calculation; for each set of samples,
  # the values for the 3 parameters are written to file on separate lines 
  cat(samplesmat[i,], fill=1, file="He3dp.in")
  
  # Load the fortran code needed to calculate S-factor curve
  if(!is.loaded("He3dp_pSub"))
    dyn.load("He3dp_plot.so") 
  .Fortran("He3dp_pSub")
  
  # read data from file
  mydat <- read.table("He3dp.out", header=FALSE)
  # One line for each scale (approximate, just using mean value)
  lines( mydat[,1], mydat[,2], col=adjustcolor("orange", alpha=0.02))
}
######################################################################
# function for error bars;
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

# add data - circles
points( obsx1, obsy1, col="gray", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="gray" )  

dev.off()
######################################################################
# output samples to file for rate calculation
write.table(samplesmat, file="MCMCsamplesTdn", quote=TRUE, 
                        row.names=FALSE, col.names=FALSE, sep="   ")




