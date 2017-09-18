# tdn analysis
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

# Er=0.35779 MeV
# g_i^2=1.0085 MeV
# g_f^2=0.025425 MeV

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
# import jags package 
library(rjags)
library(runjags)
library(R2jags)
library(mcmcplots)
library(magrittr)
library(dplyr)
## for block updating [we do not need to center predictor variables]
load.module("glm")  
load.module("nuclear")  
# 
    
######################################################################
## Read DATA GENERATION 
ensamble <- read.csv("ensamble.csv",header = T)  %>%
  mutate(Stat=replace(Stat,Stat==0,0.1)) %>%
  filter(.,dat!="Lac05")  %>% droplevels(.)
#%>%
#  filter(.,dat!="gei99b") %>%
#  filter(.,dat!="gei99d") %>%
#  filter(.,dat!="Mol80") %>%
#  filter(.,dat!="Kra87m") %>%
#  filter(.,dat!="zhi77b")

re <- as.numeric(ensamble$dat)
Nre <- length(unique(ensamble$dat))   


#index <- sample(seq(1:125),50,replace = F)
#ensamble <- ensamble[index,]

obsy = ensamble$S    # Response variable
obsx =  ensamble$E   # Predictors
erry = ensamble$Stat
syst = ensamble$Syst

syst <- unique(syst)
syst <- round(c(syst,syst[5]),3)

model.data <- list(obsy = obsy,    # Response variable
                   obsx =  obsx,   # Predictors
                   erry = erry,
                   re = re,  # random effect
                   N = nrow(ensamble),    # Sample size
#                   a0 = rep(0,Nre), # priors for scale parameters
#                   A0 = diag(Nre), 
                   syst = syst
)


######################################################################
Model <- "model{
# LIKELIHOOD
for (i in 1:N) {
  obsy[i] ~ dnorm(y[i], pow(erry[i], -2))
  y[i] <-  scale[re[i]]*sfactor3Hedp(obsx[i], e1, gin, gout) 
#+ a[re[i]]
  }

# PRIORS
# e1, gin, gout are defined as in tdn.f (by Alain Coc):
# resonance energy, initial reduced width, final reduced 
# width;
#    tau ~  dgamma(0.01,0.01)
#    e1 ~   dgamma(0.01,0.01)

e1 ~   dgamma(sh0,ra0)
sh0 <- pow(m0,2) / pow(sd0,2)
ra0 <- m0/pow(sd0,2)
m0  <- 0.35
sd0 ~  dunif(0.025,0.05) 

#gin ~ dgamma(0.01,0.01)
gin ~ dgamma(sh2,ra2)
sh2 <- pow(m2,2) / pow(sd2,2)
ra2 <- m2/pow(sd2,2)
m2 ~ dunif(0.75,1.25) 
sd2 ~ dunif(0.1,0.5)

#gout ~ dgamma(0.01,0.01)
gout ~ dnorm(sh,ra)
sh <- pow(m,2) / pow(sd,2)
ra <- m/pow(sd,2)
m ~ dunif(0.01,0.03) 
sd ~ dunif(0.1,0.5) 


# Priors for random intercept groups
#a ~ dmnorm(a0, tau.plot * A0[,])
# Priors for the two sigmas and taus
#tau.plot <- 1 / (sigma.plot * sigma.plot)
#sigma.plot ~ dunif(0.001, 10)

# This is not performing so well, as it still alows things like 1.8, which is undesirable

#for (k in 1:6){
#scale[k] ~ dlnorm(log(1.0),pow(log(1+syst[k]),-2))
#}

for (k in 1:6){
scale[k] ~ dunif(1-3*syst[k],1+3*syst[k])
}

}"


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


inits <- function () { list(e1 = runif(1,0.3,1),gin=runif(1,0.5,1.5),gout=runif(1,0.001,0.05),
                            a = rnorm(6, 0, 1),scale=runif(6, 1, 1.01)) }
# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 



# JAGS model with R2Jags;
Normfit <- jags(data = model.data,
              inits = inits,
              parameters = c("e1", "gin", "gout","scale"),
              model = textConnection(Model),
              n.thin = 10,
              n.chains = 5,
              n.burnin = 5000,
              n.iter = 10000)
mcmcChain <- as.mcmc(Normfit)[,2:4]
tmp = Reduce('+', mcmcChain)
result = tmp/length(mcmcChain)
#result[,1] <- mcmcChain[[1]][,1]
mcmcChain <- result


# Plots
denplot(Normfit,c("e1", "gin", "gout"),style="plain")
mcmcplot(Normfit)
traplot(Normfit,c("e1", "gin", "gout"),style="plain")


# variable.names are variables to be recorded in output file of samples

write.matrix(as.matrix(mcmcChain[,-1]),"samples.dat")


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
pdf("MCMCHe3dp.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

xLim = c(2e-3,1)
yLim = c(0,20)
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
   
   lines( mydat[,1], Normfit$BUGSoutput$mean$scale[1]*mydat[,2], col=adjustcolor("red", alpha=0.02))
   lines( mydat[,1], Normfit$BUGSoutput$mean$scale[2]*mydat[,2], col=adjustcolor("blue", alpha=0.02)) 
   lines( mydat[,1], Normfit$BUGSoutput$mean$scale[3]*mydat[,2], col=adjustcolor("green", alpha=0.02)) 
   lines( mydat[,1], Normfit$BUGSoutput$mean$scale[4]*mydat[,2], col=adjustcolor("orange", alpha=0.02)) 
   lines( mydat[,1], Normfit$BUGSoutput$mean$scale[5]*mydat[,2], col=adjustcolor("brown", alpha=0.02)) 
   lines( mydat[,1], Normfit$BUGSoutput$mean$scale[6]*mydat[,2], col=adjustcolor("gray", alpha=0.02)) 
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
points( obsx, obsy, col=ensamble$dat, pch=1, cex=1.2 )
add.error.bars(obsx, obsy, 0.0, erry, 0.0, col=ensamble$dat )  

dev.off()
######################################################################
# output samples to file for rate calculation
write.table(samplesmat, file="MCMCsamplesTdn", quote=TRUE, 
                        row.names=FALSE, col.names=FALSE, sep="   ")



#syt <- data.frame(E=obsx1,S=obsy1,Stat =0, Syst=0,dat="synthetic")
#syt2 <- rbind(syt,ensamble)

#ggplot(data=syt2,aes(x=E,y=S,color=dat))+geom_point(size=2.85)+
 # scale_x_log10() +
 # coord_cartesian(xlim=c(5e-3,1))
  
