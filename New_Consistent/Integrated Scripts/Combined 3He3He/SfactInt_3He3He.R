#######################################
# MCMC Chains & Integration
# Integrated into 1 Script
#######################################
# 3He3He.R
#######################################
# Preparing Workspace
rm(list=ls())
# Import jags package and other relevant packages
library(rjags)
library(magicaxis)
library(MASS)
library(ggplot2)
require(magrittr)
## For block updating [we do not need to center predictor variables]
load.module("glm")  
# Random number seed
set.seed(123)

######################################################################
## Functions
######################################################################
# Error bars
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
## Reading functions that will be used for integrating out the reaction rates
# Change the relevant lines on setwd to the folder with the relevant files
setwd(paste("C:/Users/Hong Kiat/Desktop/UROPS Materials/UROPS 1 - Thermonuclear"
, "Reaction Rates/Working Folder/auxiliar_functions",sep = " "))
source("jagsresults.R")
source("pair_wise_plot.R")
setwd(paste("C:/Users/Hong Kiat/Desktop/UROPS Materials/UROPS 1 - Thermonuclear",
            "Reaction Rates/Working Folder/Integrated 3He3He/Combined 3He3He",sep = " "))
source("numRates_3He3He.R")
source("numRates_3He3He_Table.R")

######################################################################
# Data Input
######################################################################
# Data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, S-factor in MeVb

## Set this to the working directory where you can find ensamble_3He3He.csv. This 
# is also the directory where the pdf that are outputted will be found in
setwd(paste("C:/Users/Hong Kiat/Desktop/UROPS Materials/UROPS 1 - Thermonuclear",
      "Reaction Rates/Working Folder/Integrated 3He3He/Combined 3He3He",sep = " "))

## Change the directory accordingly
ensamble <- read.csv("ensamble_3He3He.csv",header = T,fileEncoding="UTF-8-BOM") 

re <- as.numeric(ensamble$lab) # Change the label to a numeric vector
# Note that the numbers are assigned by "sorting" the labels in alphabetical order
Nre <- length(unique(ensamble$lab))
# Unique removes duplicated vector, we want to know how many groups of 
# data are there
N <- nrow(ensamble) # Total No of data sets
obsy <- ensamble$S    # Response variable in MeV
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat # Error in MeV
set <- ensamble$lab # Get the labels as a vector

syst = c(log(1.057),log(1.082),log(1.037),log(1.045),log(1.038))
# In accordance to 
# 1 Bon99 5.7%
# 2 Daw71 8.2%
# 3 Jun98 3.7%
# 4 Kra87 4.5%
# 5 Kud04 3.8%
M <- 500
xx <- seq(0.9*min(obsx),1.1*max(obsx),length.out = M) 
# Sequence of interpolating x-values

####################################
# Initialising the Rjags model
####################################
cat("model{
## LIKELIHOOD

for (i in 1:N) {
  obsy[i] ~ dnorm(ya[i], pow(erry[i], -2))
# Propagating the errors in observations of y
  
  ya[i] ~ dnorm(yt[i],pow(y.scat[re[i]], -2))
  yt[i] <- y.norm[re[i]]*mut[i]
# Temp mean (mut) multiplied by normalising constant, 
# re[i] refers to the index that we would use to refer to the category
  
  mut[i] <- (exp(2.429819509 * i.screen * (obsx[i]^(-1.5)))) * 
            (alpha + beta * obsx[i] + gamma * (obsx[i]^2))
}

## PRIORS  
### polynomial parameters

alpha ~ dnorm(0.0,pow(100, -2))T(0,)
beta ~ dnorm(0.0,pow(100, -2))
gamma ~ dnorm(0.0,pow(100, -2))
i.screen ~ dnorm(0.0,pow(100, -2))T(0,) 

for (k in 1:Nre){
# Systematic Uncertainty as a highly informative prior
y.norm[k] ~ dlnorm(log(1.0),pow(syst[k],-2))
y.scat[k] ~ dnorm(mt, pow(5,-2))T(0,)
}

mt ~ dnorm(0, pow(5,-2))T(0,) }", file={f <- tempfile()})

model.data <- list(obsy = obsy,    # Response variable
                   obsx =  obsx,   # Predictors
                   erry = erry,
                   N = N, # Sample size
                   syst = syst,
                   Nre = Nre, 
                   re = re # This is used to "iterate"
)

######################################################################
# Running the MCMC Chains
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
n.adapt  <- 1000   
n.burn   <- 1000 
n.iter   <- 7500  
thin     <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model
#
ourmodel <- jags.model(f, data = model.data, n.chains = n.chains, 
                       n.adapt = n.adapt, 
                       inits = list(alpha = 1, beta = -1, gamma = 1, i.screen = 1e-6))

update(ourmodel, n.burn) 
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                          variable.names=c("y.norm", "y.scat", "alpha", "beta",
                                           "gamma", "i.screen"), 
                          n.iter=n.iter, thin=thin)

## Outputting the Results to the Console

cat("", "\n")    # output empty line
# sample size adjusted for autocorrelation
effectiveChainLength = effectiveSize(mcmcChain) 
show(effectiveChainLength)

cat("", "\n")    # output empty line
cat("SUMMARY:", "\n")
show(summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975)))
cat("-------------------------------------------", "\n")
# Or 
jagsresults(x=mcmcChain , params=c('y.scat',"y.norm","alpha", "beta", "gamma",
                                   "i.screen"),
            probs=c(0.005,0.025, 0.16, 0.5, 0.84, 0.975,0.995))


samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)

####################################################
# Performing Integration to Obtain Reaction Rates
####################################################
a <- unlist(samplesmat[,"alpha"])
b <- unlist(samplesmat[,"beta"])
g <- unlist(samplesmat[,"gamma"])
i.s <- unlist(samplesmat[,"i.screen"])
# Histogram for T9 = 0.01, to get a sense of the distribution of the rates
hist(numRates_3He3He(a,b,g,i.s,1))

# Loops over the entire temperature grid. The default value is at 
# 1000. Here, we try with N = 75000. To reconfirm the data points, we will
# attempt to integrate through ALL of the chains

Nrate  <- numRates_3He3He_table(samplesmat, N = 5000)

require(nuclear)

sfact <- matrix(ncol=length(a) ,nrow=length(xx) )
for (i in 1:length(xx)){
  sfact[i,] <- sfactorHe3He3(xx[i],a,b,g,i.s) 
}
# This allow us to obtain the S-fact at each value of xx to plot the graphs out.
# sfactorHe3He3 is a function in the nuclear package 
# Each column refers to an interpolation point xx, while each row is an iteration

y <- matrix(nrow = nrow(sfact),ncol = 7)
colnames(y) <- c("0.5%","2.5%","16%","mean","84%","97.5%","99.5%")
for (i in 1:nrow(sfact)){
  y[i,] <- quantile(sfact[i,],probs = c(0.005,0.025, 0.16, 0.5, 0.84, 0.975,0.995))
}
# Each row represents the posterior of S(E) at xx[i]



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
library(dplyr)
require(ggmcmc)
require(ggthemes)
require(scales)
require(ggridges)
require(viridis)

obsx1 <- obsx[c(sapply(seq(1:length(obsx)), function(x){re[x] == 1}))]
obsy1 <- obsy[c(sapply(seq(1:length(obsx)), function(x){re[x] == 1}))]
errobsy1 <- erry[c(sapply(seq(1:length(obsx)), function(x){re[x] == 1}))]
obsx2 <- obsx[c(sapply(seq(1:length(obsx)), function(x){re[x] == 2}))]
obsy2 <- obsy[c(sapply(seq(1:length(obsx)), function(x){re[x] == 2}))]
errobsy2 <- erry[c(sapply(seq(1:length(obsx)), function(x){re[x] == 2}))]
obsx3 <- obsx[c(sapply(seq(1:length(obsx)), function(x){re[x] == 3}))]
obsy3 <- obsy[c(sapply(seq(1:length(obsx)), function(x){re[x] == 3}))]
errobsy3 <- erry[c(sapply(seq(1:length(obsx)), function(x){re[x] == 3}))]
obsx4 <- obsx[c(sapply(seq(1:length(obsx)), function(x){re[x] == 4}))]
obsy4 <- obsy[c(sapply(seq(1:length(obsx)), function(x){re[x] == 4}))]
errobsy4 <- erry[c(sapply(seq(1:length(obsx)), function(x){re[x] == 4}))]
obsx5 <- obsx[c(sapply(seq(1:length(obsx)), function(x){re[x] == 5}))]
obsy5 <- obsy[c(sapply(seq(1:length(obsx)), function(x){re[x] == 5}))]
errobsy5 <- erry[c(sapply(seq(1:length(obsx)), function(x){re[x] == 5}))]

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
## Grammar of Plots, ggplot2


gobs <- data.frame(obsx,obsy,erry,set)
gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"16%"],lwr2=y[,"2.5%"],lwr3=y[,"0.5%"],upr1=y[,"84%"],
                    upr2=y[,"97.5%"],upr3=y[,"99.5%"])
ggplot(gobs,aes(x=obsx,y=obsy))+
  #  Bare
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("gray80"),show.legend=FALSE)+
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("gray60"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("gray40"),show.legend=FALSE) +
  #
  #  
  
  geom_point(data=gobs,aes(x=obsx,y=obsy,group=set,color=set,shape=set),size=2.75)+
  geom_errorbar(show.legend=FALSE,data=gobs,aes(x=obsx,y=obsy,ymin=obsy-erry,ymax=obsy+erry,group=set,color=set),
                width=0.01,alpha=0.4)+
  
  scale_colour_stata(name="")+
  scale_shape_manual(values=c(0,19,8,10,4,17,3),name="") +
  theme_bw() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  +
  scale_y_log10() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.8,0.675),
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size=14,colour = "black"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size=18.5),
        axis.text  = element_text(size=13),
        axis.ticks = element_line(size = 0.75),
        axis.line = element_line(size = 0.5, linetype = "solid")) 
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
legend(0.5, 27, legend=c("Bon99", "Daw71", "Jun98", "Kra87", "Kud04"), 
       pch=c(1, 5, 0, 6, 2), col=c(68,91,258,"darkorchid","darkorange"))

text(6e-2, 20, labels=expression(paste(NULL^"3","He(",NULL^"3","He,", "2p",")",NULL^"4","He")), 
     cex=2.0)

# add data BON99 - circles
points( obsx1, obsy1, col=68, pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col=68 )  

# add data DAW71 - diamonds
points( obsx2, obsy2, col=91, pch=5, cex=1.2 )
add.error.bars( obsx2, obsy2, 0.0, errobsy2, 0.0, col=91 )

# add data JUN98 - squares
points( obsx3, obsy3, col=258, pch=0, cex=1.2 )
add.error.bars( obsx3, obsy3, 0.0, errobsy3, 0.0, col=258 )

# add data KRA87 - triangles
points( obsx4, obsy4, col="darkorchid", pch=6, cex=1.2 )
add.error.bars( obsx4, obsy4, 0.0, errobsy4, 0.0, col="darkorchid" )

# add data KUD04 - triangles
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
plot(density(a), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
title(xlab=expression(paste(S(0), "  (MeV b)")), line=3.0, cex.lab=2.0)
title(ylab="Probability density", line=3.5, cex.lab=2.0)
polygon(density(samplesmat[,1]), col=adjustcolor("blue", alpha=0.5))

# plot S'(0) -> beta
plot(density(b), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
title(xlab=expression(paste(S,"'",(0), "  (b)")), line=3.0, cex.lab=2.0)
polygon(density(samplesmat[,2]), col=adjustcolor("blue", alpha=0.5))

# plot S''(0) -> 2 * gamma
plot(density(2*g), main="", xlab="", ylab="",
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

dens <- density(1e6*i.s)

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
       c("Bon99", "Daw71", "Jun98", "Kra87", "Kud04"), 
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
plot(density(samplesmat[,10]), main="", xlab="", ylab="", xlim=c(0, 7),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
#title(ylab="density", line=2.5, cex.lab=2.3)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3)
legend("topright", legend="Bon99", pch=NA, cex=1.5)
polygon(density(samplesmat[,10]), col=adjustcolor("blue", alpha=0.5))

# plot #2                  
plot(density(samplesmat[,11]), main="", xlab="", ylab="", xlim=c(0, 0.5),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
      yaxs='i', xaxs='i')
#title(ylab="Probability", line=2.5, cex.lab=2.3)
legend("topright", legend="Daw71", pch=NA, cex=1.5)
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
plot(density(samplesmat[,14]), main="", xlab="", ylab="", xlim=c(0, 1.0),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
      yaxs='i', xaxs='i')
legend("topright", legend="Kud04", pch=NA, cex=1.5)
polygon(density(samplesmat[,14]), col=adjustcolor("blue", alpha=0.5))

dev.off()

######################################################################
# CORRELATION PLOT
######################################################################
pdf("MCMC_3He3He_h.pdf", width=6, height=6, onefile=F)
samplesmat2 <- samplesmat[,1:3]
pairs(~alpha+beta+gamma, col=adjustcolor("red", alpha=0.5),  
      data=samplesmat2[sample(nrow(samplesmat2), size=1000, replace=FALSE),], 
      main="Simple Scatterplot Matrix")

dev.off()


###############################################
# For Ease of Reading the Outputs
###############################################

sto <- summary(mcmcChain, quantiles = c(0.025, 0.16, 0.5, 0.84, 0.975))
quan <- sto$"quantiles"
res <- matrix(nrow = length(quan[,1]), ncol = 3)
res[,1] <- quan[,3] - quan[,2]
res[,2] <- quan[,3]
res[,3] <- quan[,4] - quan[,3]
rownames(res) <- rownames(quan)
colnames(res) <- c('-','value','+')
format(res, scientific = F)
format(sto$statistics, scientific = F)

###############################################
# S(E0) Estimates
##############################################

pdf("MCMC_3He3He_i.pdf")

ref.E <- 0.02194

fitvec <- (exp(0.98951013*2*2/2*sqrt(3/2*3.0160293)*samplesmat[,4]*(ref.E)^(1.5)))*
  (samplesmat[,1] + samplesmat[,2]*ref.E + samplesmat[,3]*ref.E^2)

# define quantiles of y at xchoice
a16 <- quantile(fitvec, prob = 0.16)
a50 <- quantile(fitvec, prob = 0.50)
a84 <- quantile(fitvec, prob = 0.84)

# output
cat("", "\n") 
cat("PREDICTION FOR ref.E=",ref.E, "\n") 
cat("  16%        50%        84%", "\n")
cat(quantile(fitvec, probs = 0.16),quantile(fitvec, probs = 0.50),
    quantile(fitvec, probs = 0.84), "\n")

# plot density at xchoice
## plot(density(fitvec))
par(mfcol=c(1,1), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
    las=1)

#plot(density(fitvec), main="", xlab="", ylab="",
#     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, yaxs='i', xaxs='i')
plot(1, type="n", lwd=2, col="black",  ylim=c(0, 6), xlim=c(4.5, 5.8),
     main="", xlab="", ylab="", axes=FALSE,
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, log="", yaxs='i', xaxs='i')

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.3,0), minorn=0, cex=2.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

polygon(density(fitvec), col=adjustcolor("green", alpha=0.5))

title(xlab=expression(paste(S(E['o']), "  (MeV b)")), line=4.0, cex.lab=2.0)
title(ylab="Probability density", line=4.5, cex.lab=2.0)


dev.off()


###############################################
# S(E0) Estimates
##############################################
pdf("MCMC_3He3He_i.pdf")

ref.E <- 0.02194

# Gamow Factor ~ 1

fitvec <- (exp(0.98951013*2*2/2*sqrt(3/2*3.0160293*1.66*10^(-27))*samplesmat[,4]*1.602*10^(-19)/(ref.E*1.602*10^(-13))^(1.5)))*(samplesmat[,1] + samplesmat[,2]*ref.E + samplesmat[,3]*ref.E^2)

# define quantiles of y at xchoice
a16 <- quantile(fitvec, prob = 0.16)
a50 <- quantile(fitvec, prob = 0.50)
a84 <- quantile(fitvec, prob = 0.84)

# output
cat("", "\n") 
cat("PREDICTION FOR ref.E=",ref.E, "\n") 
cat("  16%        50%        84%", "\n")
cat(quantile(fitvec, probs = 0.16),quantile(fitvec, probs = 0.50),
    quantile(fitvec, probs = 0.84), "\n")

# plot density at xchoice
## plot(density(fitvec))
par(mfcol=c(1,1), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
    las=1)

#plot(density(fitvec), main="", xlab="", ylab="",
#     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, yaxs='i', xaxs='i')
plot(1, type="n", lwd=2, col="black",  ylim=c(0, 6), xlim=c(4.5, 5.8),
     main="", xlab="", ylab="", axes=FALSE,
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, log="", yaxs='i', xaxs='i')

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.3,0), minorn=0, cex=2.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

polygon(density(fitvec), col=adjustcolor("green", alpha=0.5))

title(xlab=expression(paste(S(E['o']), "  (MeV b)")), line=4.0, cex.lab=2.0)
title(ylab="Probability density", line=4.5, cex.lab=2.0)


dev.off()


###############################################
#  Reaction Rates Posteriors
###############################################

require(pracma)
# For 6 different temperature points
tpts <- c(0.001,0.05,1,0.01,0.1,10)
data <- matrix(nrow = n.iter*n.chains,ncol = 6)
params <- matrix(ncol=6,nrow = 2)
for (i in 1:6){
  data[,i] <- numRates_3He3He(a,b,g,i.s,tpts[i])
  params[1,] <- apply(log(data),2,mean) # log mu
  params[2,] <- apply(log(data),2,sd) # f.u = exp(log var)
}

# plot histogram of rate samples
pdf("plotpdf.pdf")
par(mfcol=c(3,2), mar=c(3.5,3.5,0.5,0.5), oma=c(0.2,0.2,0.2,0.2), tck=0.05, 
    las=1, mgp=c(1.8,0.2,0))

# for linear x-axis scale, take out log="x"
for (j in 1:6){
  histRes <- hist(data[,j], plot=FALSE, breaks=490)
  xvals <- histRes$breaks
  print(xvals)
  yvals <- histRes$counts
  xvals <- c(xvals,xvals[length(xvals)])
  yvals <- c(0,yvals,0)
  yvals <- yvals/trapz(xvals,yvals)       # normalizes by area
  plot(xvals, yvals, type="S", col="red", xlab="Reaction rate", ylab="", log="x")
  # don't touch the 'x'
  curve(dlnorm(x, params[1,j], params[2,j]), add=TRUE, col="blue", lwd=2)
  legend("topright", c(paste("T9=", as.character(tpts[j]), sep=""), 
                       paste("mu=", as.character(round(params[1,j],4)), sep=""), 
                       paste("sig=", as.character(round(params[2,j],4)), sep=""))) 
}
dev.off()





###############################################
# For Ease of Reading the Outputs
###############################################

sto <- summary(mcmcChain, quantiles = c(0.005 ,0.025, 0.16, 0.5, 0.84, 0.975, 0.995))
quan <- sto$"quantiles"
res <- matrix(nrow = length(quan[,1]), ncol = 3)
res[,1] <- quan[,3] - quan[,2]
res[,2] <- quan[,3]
res[,3] <- quan[,4] - quan[,3]
rownames(res) <- c(rownames(quan))
colnames(res) <- c('-','value','+')
format(res, scientific = F)
format(sto$statistics, scientific = F)

#### Optional: Obtain the code to be used to output the table in LaTeX
# Print in latex format
require(xtable)
#write.csv(Nrate,"rates_3He3He.csv",row.names = F)
print(xtable(Nrate, type = "latex",display=c("e","g","E","E","E",
                                             "g"),
             digits=4), include.rownames = FALSE)
# g E E E g means real, scientific * 3 then real
# since include.rownames = FALSE then the first argument is somewhat irrelvant?
####
