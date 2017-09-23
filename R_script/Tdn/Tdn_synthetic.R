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
######################################################################
## ARTIFICIAL DATA GENERATION 

N <- 35

#obsx1 <- runif(N,0,0.7)
obsx1 <- exp(runif(N,log(1e-3),log(0.7)))


res <- vector()
obsy1 <- vector()
errobsy1 <- vector()

# Barker values:
# Er  = 0.0912 MeV
# g^2_in = 2.93 MeV         ! reduced width of deuteron
# g^2_out = 0.0794 MeV      ! reduced width of neutron

res[2] <- 0.0912   # resonance energy
res[3] <-  2.93    # reduced width incoming
res[4] <- 0.0794   # reduced width outgoing

for (i in 1:length(obsx1)){
  res[1] <- obsx1[i]
  write.table(res, file="tdn_AD.in", quote=TRUE, 
              row.names=FALSE, col.names=FALSE)
  
  # Load the fortran code needed to calculate S-factor curve
  if(!is.loaded("tdn_AD_Sub"))
    dyn.load("tdn_AD.so") 
  .Fortran("tdn_AD_Sub")
  
  tab1 <- read.table("tdn_AD.out", header=FALSE)
  
  errobsy1[i] <- 1
  obsy1[i] <- rnorm( 1, tab1[1,2], errobsy1[i] )
}

M <- 150
xx <- seq(min(obsx1),max(obsx1),length.out = M)
model.data <- list(obsy = obsy1,    # Response variable
                   obsx =  obsx1,   # Predictors
                   erry = errobsy1,
                   N = N, # Sample size
                   M = M, 
                   xx = xx
)





######################################################################                 
# import jags package 
library(rjags)
library(runjags)
library(R2jags)
library(mcmcplots)
require(RcppGSL)
source("https://raw.githubusercontent.com/johnbaums/jagstools/master/R/jagsresults.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")  
load.module("nuclear")  
# 
######################################################################                 
cat('model {

# LIKELIHOOD
for (i in 1:N) {
  obsy[i] ~ dnorm(y1[i], pow(erry[i], -2))
#  y1[i] ~ dnorm(sfactorTdn(obsx[i], e1, gin, gout),pow(tau,-2))  
  y1[i] <- sfactorTdn(obsx[i], e1, gin, gout)  
}

# Predicted values

for (j in 1:M){
mux[j] <- sfactorTdn(xx[j], e1, gin, gout)
#yx[j] ~ dnorm(mux[j],pow(tau,-2))
}

# PRIORS
# e1, gin, gout are defined as in tdn.f (by Alain Coc):
# resonance energy, initial reduced width, final reduced 
# width;

## Physical priors:
  # Hyperpriors

##
   

#  tau ~dgamma(0.01,0.01)

e1 ~ dnorm(0,1)T(0,)
gin ~ dunif(0,10) 
gout ~ dunif(0,10)

# e1 ~   dgamma(sh0,ra0)
# sh0 <- pow(m0,2) / pow(sd0,2)
# ra0 <- m0/pow(sd0,2)
# m0  <- 0.1
#sd0 ~  dunif(0.1,1) 

# gin ~ dgamma(0.01,0.01)
#   gin ~ dgamma(sh2,ra2)
#    sh2 <- pow(m2,2) / pow(sd2,2)
#    ra2 <- m2/pow(sd2,2)
#   m2 <- 3 
#   sd2 <- 1

#gout ~ dgamma(0.01,0.01)
#   gout ~ dnorm(sh,ra)
#   sh <- pow(m,2) / pow(sd,2)
#    ra <- m/pow(sd,2)
#   m  <- 0.1
#   sd <- 1 

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

n.burnin  <- 4500   
n.iter   <- 9000  
n.chains <- 4
n.thin   <- 5
inits <- function () { list(e1 = runif(1,0,0.5),gin=runif(1,2,10),gout=runif(1,0.01,1)) }
# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 



# JAGS model with R2Jags;
out <- jags(data = model.data,
              inits = inits,
              parameters = c("e1", "gin", "gout","mux"),
              model.file = f,
              n.thin = n.thin,
              n.chains = n.chains,
              n.burnin = n.burnin,
              n.iter = n.iter)


traplot(out ,c("e1", "gin", "gout"),style="plain")
denplot(out ,c("e1", "gin", "gout"),style="plain")


# Plot
y <- jagsresults(x=out, params=c('mux'),probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
x <- xx
gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],upr1=y[,"75%"],upr2=y[,"97.5%"])
gobs <- data.frame(obsx1,obsy1,errobsy1)



ggplot(gobs,aes(x=obsx1,y=obsy1))+
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.95, fill = c("#fac901"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.95, fill=c("#225095"),show.legend=FALSE) +
  geom_point()+
  geom_errorbar(data=gobs,mapping=aes(x=obsx1,y=obsy1,ymin=obsy1-errobsy1,ymax=obsy1+errobsy1),alpha=0.85,
                colour="#dd0100",width=0.005)+
  geom_line(data=gdata,aes(x=xx,y=mean),colour="#ffffff",linetype="dashed",size=1,show.legend=FALSE)+
  theme_bw() + xlab("E") + ylab("S") + scale_x_log10() +coord_cartesian(xlim=c(2e-3,1),ylim=c(0,35))






# Plot posterior predicted y at x probe:
# Convert coda object to matrix: 
mcmcChain <- as.mcmc(out)[,2:4]
mcmcMat = as.matrix(mcmcChain)
# Find the xProbe and predicted y columns:
xPcols = grep( "xP" , colnames(mcmcMat) , value=FALSE )
yPcols = grep( "yP" , colnames(mcmcMat) , value=FALSE )
# Find the extreme predicted values for graph axis limits:
xLim = quantile( mcmcMat[,yPcols] , probs=c(0.005,0.995) )
# Make the plots of the posterior predicted values:
for ( i in 1:length(xPcols) ) {
  openGraph(width=4,height=3)
  plotPost( mcmcMat[,yPcols[i]] , xlab="y" , xlim=xLim , cenTend="mean" ,
            main=bquote( "Posterior Predicted y for x = "
                         * .(mcmcMat[1,xPcols[i]]) )  )
}











mcmcplot(out)
traplot(out)




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
yLim = c(0,35)
######################################################################
# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
      axes=FALSE, main="", xlab = "", ylab = "",
      cex=1.5 , cex.lab=1.3, cex.axis=1.0,
      cex.main=1.0, log="x")
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

text(0.4, 30, labels=expression(paste("t(d,n)",NULL^"4","He")), cex=1.3)
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

for ( i in round(seq(from=1,to=nsamp,length=1200)) ) {
  # output samples to file for S-factor calculation; for each set of samples,
  # the values for the 3 parameters are written to file on separate lines 
  cat(samplesmat[i,], fill=1, file="tdn.in")
  
  # Load the fortran code needed to calculate S-factor curve
  if(!is.loaded("tdn_pSub"))
    dyn.load("tdn_plot.so") 
  .Fortran("tdn_pSub")
  
  # read data from file
  mydat <- read.table("tdn.out", header=FALSE)
  
  lines( mydat[,1], mydat[,2], col=adjustcolor("red",alpha=1)) 
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
points( obsx1, obsy1, col="black", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="black" )  
dev.off()
######################################################################
# output samples to file for rate calculation
write.table(samplesmat, file="MCMCsamplesTdn", quote=TRUE, 
                        row.names=FALSE, col.names=FALSE, sep="   ")




