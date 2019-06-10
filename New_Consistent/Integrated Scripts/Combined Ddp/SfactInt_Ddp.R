#######################################
# MCMC Chains & Integration
# Integrated into 1 Script
#######################################
# Ddp.R
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
            "Reaction Rates/Working Folder/Integrated Ddp/Combined Ddp",sep = " "))
source("numRates_ddp.R")
source("numRates_ddp_Table.R")

######################################################################
# Data Input
######################################################################
# Data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, S-factor in MeVb

## Set this to the working directory where you can find ensamble_Ddp.csv. This 
# is also the directory where the pdf that are outputted will be found in
setwd(paste("C:/Users/Hong Kiat/Desktop/UROPS Materials/UROPS 1 - Thermonuclear",
            "Reaction Rates/Working Folder/Integrated Ddp/Combined Ddp",sep = " "))
## Change the directory accordingly

ensamble <- read.csv("ensamble_Ddp.csv",header = T,fileEncoding="UTF-8-BOM") 

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

syst = c(log(1.013),log(1.03),log(1.064),log(1.082),log(1.020))
# In accordance to 
# 1 = 1.3% for Bro90
# 2 = 3% for Gre95
# 3 = 6.4% for Kra87 (B)
# 4 = 8.2% for Kra87 (M)
# 5 = 2% for Leo06
M <- 500
xx <- seq(0.9*min(obsx),1.1*max(obsx),length.out = M) 
# Sequence of interpolating x-values

######################################################################
# Input of Theoretical Model [MARCUCCI et al.]
######################################################################
# E is in MeV, S is in eV b; convert to MeV and MeV b, respectively
theory <- read.table("Arai_ddp_2011.dat", header=FALSE)

interp.x <- theory[,1] 
interp.y <- theory[,2]*1e-3

# we will use JAGS interp.lin function to use this theoretical S-factor:
#
# - the columns of this table define vectors x and y
# - a single point is given by x_i, y_i
# - interp.lin gives the y value for the x value provided as argument e,
#   interp.lin(e, x, y)

####################################
# Initialising the Rjags model
####################################

cat("model{
## LIKELIHOOD
for (i in 1:N) {
  obsy[i] ~ dnorm(ya[i], pow(erry[i], -2))
# Propagating the errors in observations of y
  ya[i] ~ dnorm(y.norm[re[i]]*mut[i],pow(y.scat[re[i]], -2))
# Temp mean (mut) multiplied by normalising constant, 
# re[i] refers to the index that we would use to refer to the category
  mut[i] <- a.scale*interp.lin(obsx[i], interp.x, interp.y)
# Interpolate and scale
}

## PRIORS
# Scaling factor of theory 
a.scale ~ dnorm(0, pow(5,-2))T(0,)    

for (k in 1:Nre){
# Systematic Uncertainty as a highly informative prior
y.norm[k] ~ dlnorm(log(1.0),pow(syst[k],-2))
y.scat[k] ~  dnorm(mt, pow(5,-2))T(0,)
}
mt ~  dnorm(0, pow(5,-2))T(0,) # Hyperprior
}", file={f <- tempfile()})

model.data <- list(obsy = obsy,    # Response variable
                   obsx =  obsx,   # Predictors
                   erry = erry,
                   N = nrow(ensamble), # Sample size
                   syst = syst,
                   Nre = Nre, 
                   re = re, # This is used to "iterate"
                   interp.x = interp.x,
                   interp.y = interp.y
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
n.iter   <- 5000  
thin     <- 1

# "f": is the model specification from above; 
# data = list(...): define all data elements that are referenced in the 
# JAGS model
#
ourmodel <- jags.model(f, data = model.data, n.chains = n.chains, 
                       n.adapt = n.adapt)

update(ourmodel, n.burn) 
# variable.names are variables to be recorded in output file of samples
mcmcChain <- coda.samples(ourmodel, 
                          variable.names=c("a.scale", "y.norm","y.scat"), 
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
jagsresults(x=mcmcChain , params=c('y.scat',"a.scale", "y.norm"),
            probs=c(0.005,0.025, 0.16, 0.5, 0.84, 0.975,0.995))


samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)

####################################################
# Performing Integration to Obtain Reaction Rates
####################################################
apost <- unlist(samplesmat[,"a.scale"])

# Histogram for T9 = 0.01, to get a sense of the distribution of the rates
hist(numRates_ddp(apost,0.01))

# Loops over the entire temperature grid. The default value is at 
# 1000. Here, we try with N = 75000. To reconfirm the data points, we will
# attempt to integrate through ALL of the chains

Nrate  <- numRates_ddp_table(samplesmat,vars = "a.scale", N = 5000)

mux <- outer(apost,(approx(interp.x, interp.y, xx, method="linear")$y))
# Each column refers to an interpolation point xx, while each row is an iteration

y <- matrix(nrow = ncol(mux),ncol = 7)
colnames(y) <- c("0.5%","2.5%","16%","mean","84%","97.5%","99.5%")
for (i in 1:ncol(mux)){
  y[i,] <- quantile(mux[,i],probs = c(0.005,0.025, 0.16, 0.5, 0.84, 0.975,0.995))
}

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
pdf("MCMC_Ddp_a.pdf")
plot(mcmcChain)
dev.off()

######################################################################
# S-FACTOR FIT + DATA 
######################################################################
pdf("MCMC_Ddp_b.pdf",width=10,height=5,onefile=F)
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
        legend.position = c(0.1,0.675),
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
pdf("MCMC_Ddp_c.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
    las=1)

# determine plot ranges
xLim = c(1e-2,0.5)
yLim = c(0.03,0.2)

# plot axes only...add lines...then data
# cex.lab: controls axis labels
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim, 
      axes=FALSE, main="", xlab = "", ylab = "", log="xy")
# line: controls distance between axis and label 
title(xlab="Energy (MeV)", line=2.5, cex.lab=2.0)
title(ylab="S-Factor (MeV b)", line=2.5, cex.lab=2.0)
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
# cex.axis: size of tick mark labels
magaxis(mgp=c(0,0.2,0), cex.axis=1.3)
box()
# plot legend
legend(0.2, 0.08, legend=c( "Bro90", "Gre95" , "Kra87B" , "Kra87M", "Leo06" ), 
       pch=c(1, 0, 6, 5), col=c("blue", "red", "green4", "black", "orange"))

text(4e-3, 2.0e-6, labels=expression(paste("D(d,n)",NULL^"3","He")), 
     cex=2.0)

# add data Bro90 - circles
points( obsx1, obsy1, col="blue", pch=1, cex=1.2 )
add.error.bars( obsx1, obsy1, 0.0, errobsy1, 0.0, col="blue" )  

# add data Gre95 - diamonds
points( obsx2, obsy2, col="red", pch=0, cex=1.2 )
add.error.bars( obsx2, obsy2, 0.0, errobsy2, 0.0, col="red" )

# add data Kra87B - squares
points( obsx3, obsy3, col="green4", pch=6, cex=1.2 )
add.error.bars( obsx3, obsy3, 0.0, errobsy3, 0.0, col="green4" )

# add data Kra87M - triangles
points( obsx4, obsy4, col="black", pch=5, cex=1.2 )
add.error.bars( obsx4, obsy4, 0.0, errobsy4, 0.0, col="black" )

# add data Leo06 - pluses
points( obsx5, obsy5, col="black", pch=3, cex=1.2 )
add.error.bars( obsx5, obsy5, 0.0, errobsy5, 0.0, col="orange" )

dev.off()

######################################################################
# POSTERIOR OF S-FACTOR AT GIVEN ENERGY 
######################################################################
# plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

pdf("MCMC_Ddp_d.pdf")

# define x value for which we would like to predict y
xchoice <- 0.00

# declare vector with y values
fitvec <- vector(mode = "numeric", length = nsamp)   

# calculate y values from all regression samples at given x value
predmat <- approx(interp.x, interp.y, xchoice, method="linear")

for(i in 1:nsamp) fitvec[i] <- samplesmat[i,1] * predmat$y

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
## plot(density(fitvec))
par(mfcol=c(1,1), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
    las=1)

#plot(density(fitvec), main="", xlab="", ylab="",
#     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, yaxs='i', xaxs='i')
plot(1, type="n", lwd=2, col="black",  ylim=c(0, 10), xlim=c(5.0, 5.7),
     main="", xlab="", ylab="", axes=FALSE,
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, log="", yaxs='i', xaxs='i')

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.3,0), minorn=0, cex=2.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

polygon(density(fitvec*1e2), col=adjustcolor("blue", alpha=0.5))

title(xlab=expression(paste(S(0), "  (MeV b)")), line=4.0, cex.lab=2.0)
title(ylab="Probability density", line=4.5, cex.lab=2.0)

mtext(expression(paste("x10"^{-2})), side=1, line=2.0, adj=1.0, cex=1.5)

dev.off()

######################################################################
# THEORY SCALE FACTOR 
######################################################################
# plot posterior predictive distribution at a given x by 
# marginalization over all parameters; we will use all credible
# lines, calculate y-values for all samples at given x; this set of 
# y-values represents our predicted posterior at x

pdf("MCMC_Ddp_e.pdf")

# plot density at xchoice
## plot(density(fitvec))
par(mfcol=c(1,1), mar=c(5.5,6.5,4.5,0), oma=c(2.0,5.0,0.5,2.0), tck=0.02, 
    las=1)

plot(1, type="n", lwd=2, col="black",  
     main="", xlab="", ylab="", axes=FALSE, xlim=c(0.87, 1.1), ylim=c(1,50),
     cex=10.5, cex.main=10.0, cex.lab=2.0, cex.axis=1.3, log="", yaxs='i', xaxs='i')

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.3,0), minorn=0, cex=2.0, cex.lab=1.3, cex.main=1.0, cex.axis=1.5)
box()

polygon(density(samplesmat[,1]), col=adjustcolor("blue", alpha=0.5))

title(xlab="Theory scale factor", line=3.0, cex.lab=2.0)
title(ylab="Probability density", line=3.0, cex.lab=2.0)

dev.off()

######################################################################
# DENSITIES OF S-FACTOR NORMALIZATION FACTORS
######################################################################
pdf("MCMC_Ddp_f.pdf",width=10, height=6, onefile=F)
par(mfcol=c(1,1), mar=c(5.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
    las=1)

xLim = c(0.8, 1.15)

# plot density in first panel           
plot(density(samplesmat[,5]), main="", xlab="", ylab="",
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, 
     xlim=xLim, yaxs='i', xaxs='i',ylim=c(0,50))
title(ylab="Probability density", line=3.5, cex.lab=2.3)
title(xlab="S-factor normalization", line=4.0, cex.lab=2.3)

polygon(density(samplesmat[,2]),  
        col=adjustcolor("blue", alpha=0.5))
polygon(density(samplesmat[,3]), 
        col=adjustcolor("red", alpha=0.5))
polygon(density(samplesmat[,4]),  
        col=adjustcolor("green4", alpha=0.5))
polygon(density(samplesmat[,5]),  
        col=adjustcolor("black", alpha=0.5))
polygon(density(samplesmat[,6]),  
        col=adjustcolor("orange", alpha=0.5))

legend("topleft", inset=.01, 
       c( "Bro90", "Gre95" , "Kra87B" , "Kra87M", "Leo06" ), 
       fill=adjustcolor(c("blue", "red", "green4", "black", "orange"), alpha=0.5), 
       horiz=FALSE, cex=1, box.lty=0)

dev.off()

######################################################################
# POSTERIOR EXTRINSIC S-FACTOR SCATTER
######################################################################
pdf("MCMC_Ddp_g.pdf",width=10, height=6, onefile=F)
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
par(mfcol=c(2,3), mar=c(5.5,7.0,1.0,2.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
    las=1)

# plot #1          
plot(density(samplesmat[,7]*1e2), main="", xlab="", ylab="", xlim=c(0, 0.4),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
#title(ylab="density", line=2.5, cex.lab=2.3)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3)
legend("topright", legend="Bro90", pch=NA, cex=1.5)
polygon(density(samplesmat[,7]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)

# plot #2                  
plot(density(samplesmat[,8]*1e2), main="", xlab="", ylab="", xlim=c(0, 2),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
      yaxs='i', xaxs='i')
#title(ylab="Probability", line=2.5, cex.lab=2.3)
legend("topright", legend="Gre95", pch=NA, cex=1.5)
polygon(density(samplesmat[,8]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)

# plot #3                  
plot(density(samplesmat[,9]*1e2), main="", xlab="", ylab="", xlim=c(0, 1.5),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
      yaxs='i', xaxs='i')
legend("topright", legend="Kra87B", pch=NA, cex=1.5)
polygon(density(samplesmat[,9]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)

# plot #4                  
plot(density(samplesmat[,10]*1e2), main="", xlab="", ylab="", xlim=c(0, 1.5),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
      yaxs='i', xaxs='i')
legend("topright", legend="Kra87M", pch=NA, cex=1.5)
polygon(density(samplesmat[,10]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)

# plot #5                  
plot(density(samplesmat[,11]*1e2), main="", xlab="", ylab="", xlim=c(0, 2),
     cex=10.5, cex.main=10.0, cex.lab=2.3, cex.axis=1.5, yaxs='i', xaxs='i'
)
title(xlab=expression(paste(sigma [S], " (MeVb)")), line=4.0, cex.lab=2.3, 
      yaxs='i', xaxs='i')
legend("topright", legend="Leo06", pch=NA, cex=1.5)
polygon(density(samplesmat[,11]*1e2), col=adjustcolor("blue", alpha=0.5))

mtext(expression(paste("x10"^{-2})), side=1, line=2.7, adj=1.05, cex=1.2)

dev.off()

#################################################
# Plotting the Autocorrelation Plots
##################################################
## Import MCMC samples into a ggs object that can be used by all ggs_ graphical functions
# %>% pipes to tibble
pdf("MCMC_Ddp_h.pdf",width=10, height=6, onefile=F)
Sp <- ggs(mcmcChain[,c("a.scale", "y.norm[1]","y.norm[2]","y.norm[3]","y.norm[4]",
                       "y.scat[1]","y.scat[2]","y.scat[3]","y.scat[4]")]) %>% as_tibble()

pair_wise_plot(Sp)

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
  data[,i] <- numRates_ddp(apost,tpts[i])
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
#write.csv(Nrate,"rates_Ddp.csv",row.names = F)
print(xtable(Nrate, type = "latex",display=c("e","g","E","E","E",
                                             "g"),
             digits=4), include.rownames = FALSE)
# g E E E g means real, scientific * 3 then real
# since include.rownames = FALSE then the first argument is somewhat irrelvant?
####
