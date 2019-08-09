######################################################################
# Author: Christian Iliadis (05/13/2019) and HK
######################################################################
#
# this script splits each data sets into two parts: a "relative cross
# section" part, and an "absolute cross section" part; the relative
# part is modeled with a normalization factor [e.g., y.norm1] for which
# the prior is sampled from a very broad density; the absolute part
# consists of just one data point. It does have a statistical uncertainty.
# It is modeled with anormalization factor [e.g., y.norm11] for which
# the prior is sampled from a highly informative lognormal density, where
# the lognormal parameters mu and sigma a chosen according to the
# systematic uncertainty reported for a given experiment;
#
# In the plot of the reduced cross section, the red data points show
# the "absolute" cross section of each data set
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package
library(nimble)
library(magicaxis)
require(gsl)
require(RcppGSL)
library(emdbook)
require(nuclear)


sigma7Benp7mod <- function(ecm,
                           e0_1, ga_1, gb_1, ra_1, rb_1,
                           e0_2, ga_2, gb_2, ra_2, rb_2,
                           e0_3, ga_3, gb_3, ra_3, rb_3,
                           e0_4, ga_4, gb_4, ra_4, rb_4,
                           e0_5, ga_5, gb_5, ra_5, rb_5,
                           e0_6, ga_6, gb_6, ra_6, rb_6,
                           e0_7, ga_7, gb_7, ra_7, rb_7){

  SF1 <-  sigma7Benp(ecm, e0_1, ga_1, gb_1, ra_1, rb_1, jr = 2, la = 0, lb = 0)
  SF2 <-  sigma7Benp(ecm, e0_2, ga_2, gb_2, ra_2, rb_2, jr = 3, la = 1, lb = 1)
  SF3 <-  sigma7Benp(ecm, e0_3, ga_3, gb_3, ra_3, rb_3, jr = 3, la = 1, lb = 1)
  SF4 <-  sigma7Benp(ecm, e0_4, ga_4, gb_4, ra_4, rb_4, jr = 1, la = 0, lb = 0)
  SF5 <-  sigma7Benp(ecm, e0_5, ga_5, gb_5, ra_5, rb_5, jr = 4, la = 3, lb = 3)
  SF6 <-  sigma7Benp(ecm, e0_6, ga_6, gb_6, ra_6, rb_6, jr = 2, la = 1, lb = 1)
  SF7 <-   sigma7Benp(ecm, e0_7, ga_7, gb_7, ra_7, rb_7, jr = 0, la = 1, lb = 1)
  SF <- SF1 + SF2 + SF3 + SF4 + SF5 + SF6 + SF7
  return(SF = SF)
}

# Vectorised sigma function
sigmaBe7mod <- nimbleRcall(function(ecm = double(1),
                                    e0_1 = double(0), ga_1 = double(0), gb_1 = double(0),
                                    ra_1 = double(0), rb_1 = double(0),
                                    e0_2 = double(0), ga_2 = double(0), gb_2 = double(0),
                                    ra_2 = double(0), rb_2 = double(0),
                                    e0_3 = double(0), ga_3 = double(0), gb_3 = double(0),
                                    ra_3 = double(0), rb_3 = double(0),
                                    e0_4 = double(0), ga_4 = double(0), gb_4 = double(0),
                                    ra_4 = double(0), rb_4 = double(0),
                                    e0_5 = double(0), ga_5 = double(0), gb_5 = double(0),
                                    ra_5 = double(0), rb_5 = double(0),
                                    e0_6 = double(0), ga_6 = double(0), gb_6 = double(0),
                                    ra_6 = double(0), rb_6 = double(0),
                                    e0_7 = double(0), ga_7 = double(0), gb_7 = double(0),
                                    ra_7 = double(0), rb_7 = double(0)
){},
Rfun = "sigma7Benp7mod", returnType = double(1))

######################################################################
## DATA SETS
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the
# latter is the individual statistical error of each datum [i];
# energy is in MeV, sqrt(Ecm)*sigma is in sqrt(MeV)b

data(Be7np)

re <- as.numeric(as.factor(Be7np$dat))
Nre <- length(unique(Be7np$dat))
# Unique removes duplicated vector, we want to know how many groups of
# data are there
N <- nrow(Be7np) # Total No of data sets
obsy <- Be7np$S    # Response variable in MeV
obsx <-  Be7np$E   # Predictors
erry <- Be7np$Stat # Error in MeV
set <- Be7np$dat # Get the labels as a vector
fu <- log(c(1.085,2.25,1.10,2.25,1.050,1.051,2.25,1.020,1.051,1.032))
scat <- c(1e-6,2,1e-6,2,1e-6,1e-6,2,1e-6,2,1e-6)

samplerCode <- nimbleCode({

  ###################
  # LIKELIHOODS
  ###################
  # - careful: dnorm is differently defined in R and JAGS!
  # - precision=sigma^(-2)
  # - in a for loop, make sure **all** variables on the LEFT of an
  #   expression has the index [i]
  # - systematic error as normalization factor y.norm...

  ## Calling sigmaBe7mod once!
  sigmaBe7modT[1:N] <- sigmaBe7mod(obsx[1:N],
                                 e0_1, ga_1, gb_1, ra, rb,
                                 e0_2, ga_2, gb_2, ra, rb,
                                 e0_3, ga_3, gb_3, ra, rb,
                                 e0_4, ga_4, gb_4, ra, rb,
                                 e0_5, ga_5, gb_5, ra, rb,
                                 e0_6, ga_6, gb_6, ra, rb,
                                 e0_7, ga_7, gb_7, ra, rb)

  for (i in 1:N){
    obsy[i] ~ dnorm(ym[i], sd = erry[i])
     ym[i] ~ dnorm(yt[i],sd = y.scat[re[i]])
     yt[i] <- y.norm[re[i]]*sqrt(obsx[i])*sigmaBe7modT[i]
     }




  ###################
  # PRIORS
  ###################
  # parameters: ecm, e0, ga, gb, ra, rb, xj, xla, xlb

  # channel radii
  ra ~ T(dnorm(4.0, pow(0.5, -2)),0,Inf)
  rb ~ T(dnorm(4.0, pow(0.5, -2)),0,Inf)

  # Wigner limit: wl = hbar^2/(m_red a_c^2) = 41.80159/(M_red a_c^2)
  #
  # neutron channel: wl_n = 41.80159/(0.88186 a_c^2) = 47.40160/a_c^2
  # proton channel:  wl_p = 41.80159/(0.88079 a_c^2) = 47.45920/a_c^2
  #
  wl_n <- 47.40160*pow(ra, -2)
  wl_p <- 47.45920*pow(rb, -2)


  ###################################################################
  # RESONANCE 1: e0=0 MeV
  # resonance spin, orbital angular momenta   # 2-
  # energy eigenvalue
  e0_1 ~ T(dnorm(0.0, pow(0.1, -2)),0,Inf)

  # reduced widths
  ga_1 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_1 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 2: e0=0.150 MeV
  # resonance spin, orbital angular momenta # 3+


  # energy eigenvalue
  e0_2 ~ T(dnorm(0.15, pow(0.025, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_2 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_2 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)


  ##############
  # RESONANCE 3: e0=0.336 MeV
  # resonance spin, orbital angular momenta  # 3+


  # energy eigenvalue
  e0_3 ~ T(dnorm(0.336, pow(0.010, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_3 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_3 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 4: e0=0.510 MeV
  # resonance spin, orbital angular momenta # 1-

  # energy eigenvalue
  e0_4 ~ T(dnorm(0.51, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_4 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_4 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 5: e0=0.960 MeV
  # resonance spin, orbital angular momenta # 4+

  # energy eigenvalue
  e0_5 ~ T(dnorm(0.96, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_5 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_5 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 6: e0=1.23 MeV
  # resonance spin, orbital angular momenta # 2+


  # energy eigenvalue
  e0_6 ~ T(dnorm(1.23, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak

  # reduced widths
  ga_6 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_6 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)

  ##############
  # RESONANCE 7: e0 = 1.32 MeV
  # resonance spin, orbital angular momenta # 0+


  # energy eigenvalue
  e0_7 ~ T(dnorm(1.32, pow(0.1, -2)),0,Inf)         # positive since we see sigma peak
  # reduced widths
  ga_7 ~ T(dnorm(0.0, pow(wl_n/2, -2)),0,Inf)
  gb_7 ~ T(dnorm(0.0, pow(wl_p/2, -2)),0,Inf)
  ##################################################################


  for (k in 1:Nre) {
    # Systematic Uncertainty as a highly informative prior
    y.norm[k] ~ dlnorm(log(1.0), pow(fu[k], -2))
    y.scat[k] ~ T(dnorm(0,sd=scat[k]),0,Inf)
  }

#  k ~ T(dnorm(0,sd=0.15),0,Inf)
  r_1 <- ga_1/gb_1
  r_4 <- ga_4/gb_4
})

samplerData <- list(obsy = obsy,
                    re = re, # Sample size
                    erry = erry,
                    obsx = obsx,
                    fu = fu,
                    scat = scat)


samplerConst <- list(N = N,
                     Nre = Nre)


samplerInits <- list(y.norm = rep(1,Nre),
                     y.scat = rep(0.001,Nre),
                     e0_1 = 0.05, gb_1 = 0.6, ga_1 = 0.6,
                     e0_2 = 0.15, ga_2 = 0.6, gb_2 = 0.6,
                     e0_3 = 0.336, ga_3 = 0.6, gb_3 = 0.6,
                     e0_4 = 0.51, ga_4 = 0.6, gb_4 = 0.6,
                     e0_5 = 0.96, ga_5 = 0.6, gb_5 = 0.6,
                     e0_6 = 1.23, ga_6 = 0.6, gb_6 = 0.6,
                     e0_7 = 1.32, ga_7 = 0.6, gb_7 = 0.6,
                     ra = 4, rb = 4,mt=0
 #                    ,k = 0.1
                     )


###############################################################
# (Alternative) If invoking Nimble MCMC stepwise (but more customisable)
##############################################################
ourmodel <- nimbleModel(code = samplerCode, constants = samplerConst,
                        data = samplerData, inits = samplerInits, check = FALSE
)
compileNimble(ourmodel)
# Always compile the model after you are done setting up with it

ourmodel$calculate()
# Calculate the log likelihood (logProb). If the value is not NA,
# we have successfully initialised the model (FULLY)
# One iteration: simulate -> new values -> calculate

conf <- configureMCMC(ourmodel,print=TRUE)
# print = TRUE tells you what kind of samplers are being used for each stochastic node

conf$addMonitors(c('e0_1','ga_1','gb_1','e0_2','ga_2','gb_2',
                   'e0_3', 'ga_3', 'gb_3', 'e0_4', 'ga_4', 'gb_4',
                   'e0_5', 'ga_5', 'gb_5', 'e0_6', 'ga_6', 'gb_6',
                   'e0_7', 'ga_7', 'gb_7',
                   'r_1', 'r_4', 'ra', 'rb',
                   'y.norm',
                   'y.scat'))
# Add the parameters to monitor throughout the MCMC process

samplerMCMC <- buildMCMC(conf)
# Note that if no configuration is required, we can directly jump towards buildMCMC.
# But you know, then theres no point in using this entire chunk of code
# where you can just run nimbleMCMC

compiledMCMC <- compileNimble(samplerMCMC,project = ourmodel,showCompilerOutput = TRUE)
# Complie the configured model again once we are done with configuring it;
# Just before we perform the MCMC runs
# Can include project = ourmodel; to include all the nimbleFunctions, nimbleModels that
# you would want to inherit the functions and models from
# resetFunctions = TRUE; if you would want to reset all the previously created functions
# in order to addd the new MCMC

n.chains = 1
n.iter = 6000
n.burnin = 2500

system.time(
  mcmcChain <- runMCMC(compiledMCMC,niter = n.iter, nchain = n.chains, nburnin = n.burnin,
                       setSeed=9,samplesAsCodaMCMC = TRUE)
)

save(mcmcChain, file = "MCMCBe7.RData")

pdf("Final11.pdf")
plot(mcmcChain)
dev.off()





samp <- as.matrix(mcmcChain)
xx <- exp(seq(log(1e-8),log(2),length.out = 250))
# resonance 1
y1 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_1"],samp [,"ga_1"],samp [,"gb_1"],
                                          samp [,"ra"],
                                          samp [,"rb"],jr = 2, la = 0, lb = 0),probs=c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))
  y1 <- rbind(y1,mux)
}
gr1 <-  data.frame(x =xx, mean = y1[,"50%"],lwr1=y1[,"16%"],lwr2=y1[,"2.5%"],lwr3=y1[,"0.15%"],
                   upr1=y1[,"84%"],
                   upr2=y1[,"97.5%"],upr3=y1[,"99.85%"])
# resonance 2
y2 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_2"],samp[,"ga_2"],samp [,"gb_2"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 3, la = 1, lb = 1),
                                probs=c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))
y2 <- rbind(y2,mux)
}

gr2 <-  data.frame(x =xx, mean = y2[,"50%"],lwr1=y2[,"16%"],lwr2=y2[,"2.5%"],lwr3=y2[,"0.15%"],
                   upr1=y2[,"84%"],
                   upr2=y2[,"97.5%"],upr3=y2[,"99.85%"])

# resonance 3
y3 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_3"],samp[,"ga_3"],samp [,"gb_3"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 3, la = 1, lb = 1),
                   probs=c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))
  y3 <- rbind(y3,mux)
}

gr3 <-  data.frame(x =xx, mean = y3[,"50%"],lwr1=y3[,"16%"],lwr2=y3[,"2.5%"],lwr3=y3[,"0.15%"],
                   upr1=y3[,"84%"],
                   upr2=y3[,"97.5%"],upr3=y3[,"99.85%"])


# resonance 4
y4 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_4"],samp[,"ga_4"],samp [,"gb_4"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 1, la = 0, lb = 0),
                   probs=c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))
  y4 <- rbind(y4,mux)
}

gr4 <-  data.frame(x =xx, mean = y4[,"50%"],lwr1=y4[,"16%"],lwr2=y4[,"2.5%"],lwr3=y4[,"0.15%"],
                   upr1=y4[,"84%"],
                   upr2=y4[,"97.5%"],upr3=y4[,"99.85%"])

# resonance 5
y5 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_5"],samp[,"ga_5"],samp [,"gb_5"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 4, la = 3, lb = 3),
                   probs=c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))
  y5 <- rbind(y5,mux)
}

gr5 <-  data.frame(x =xx, mean = y5[,"50%"],lwr1=y5[,"16%"],lwr2=y5[,"2.5%"],lwr3=y5[,"0.15%"],
                   upr1=y5[,"84%"],
                   upr2=y5[,"97.5%"],upr3=y5[,"99.85%"])


# resonance 6
y6 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_6"],samp[,"ga_6"],samp [,"gb_6"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 2, la = 1, lb = 1),
                   probs=c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))
  y6  <- rbind(y6,mux)
}

gr6 <-  data.frame(x =xx, mean = y6[,"50%"],lwr1=y6[,"16%"],lwr2=y6[,"2.5%"],lwr3=y6[,"0.15%"],
                   upr1=y6[,"84%"],
                   upr2=y6[,"97.5%"],upr3=y6[,"99.85%"])



yall <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp7mod(xx[j],samp[,"e0_1"],samp[,"ga_1"],samp [,"gb_1"],
                                          samp[,"ra"],samp[,"rb"],
                                          samp[,"e0_2"],samp[,"ga_2"],samp [,"gb_2"],
                                          samp[,"ra"],samp[,"rb"],
                                          samp[,"e0_3"],samp[,"ga_3"],samp [,"gb_3"],
                                          samp[,"ra"],samp[,"rb"],
                                          samp[,"e0_4"],samp[,"ga_4"],samp [,"gb_4"],
                                          samp[,"ra"],samp[,"rb"],
                                          samp[,"e0_5"],samp[,"ga_5"],samp [,"gb_5"],
                                          samp[,"ra"],samp[,"rb"],
                                          samp[,"e0_6"],samp[,"ga_6"],samp [,"gb_6"],
                                          samp[,"ra"],samp[,"rb"],
                                          samp[,"e0_7"],samp[,"ga_7"],samp [,"gb_7"],
                                          samp[,"ra"],samp[,"rb"]),
                   probs=c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))
  yall <- rbind(yall,mux)
}

grall <-  data.frame(x =xx, mean = yall[,"50%"],lwr1=yall[,"16%"],lwr2=yall[,"2.5%"],lwr3=yall[,"0.15%"],
                   upr1=yall[,"84%"],
                   upr2=yall[,"97.5%"],upr3=yall[,"99.85%"])



### plot Maxwell-Boltzmann "factor" at T=1 GK [arbitrary scale]
#for ( i in round(runif(500, min=1, max=nsamp)) ) {
#  lines(x1, 3e2*(x1*(2.718^(-x1/(0.086173*0.5)))))
#}


  
# Plot all

ggplot(Be7np,aes(x=E,y=S))+
  
  geom_ribbon(data=gr1,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("darkolivegreen4"),alpha=0.4,show.legend=FALSE) +
  geom_ribbon(data=gr1,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("darkolivegreen4"),alpha=0.8,show.legend=FALSE) +
  
  geom_ribbon(data=gr2,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#9e9ac8"),alpha=0.4,show.legend=FALSE) +
  geom_ribbon(data=gr2,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#984ea3"),alpha=0.8,show.legend=FALSE) +
  
  geom_ribbon(data=gr3,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("darkgoldenrod2"),alpha=0.4,show.legend=FALSE) +
  geom_ribbon(data=gr3,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("darkgoldenrod2"),alpha=0.8,show.legend=FALSE) +
  
  
  geom_ribbon(data=gr4,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("cyan3"),alpha=0.4,show.legend=FALSE) +
  geom_ribbon(data=gr4,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("cyan3"),alpha=0.8,show.legend=FALSE) +
  
  geom_ribbon(data=gr5,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("chocolate1"),alpha=0.4,show.legend=FALSE) +
  geom_ribbon(data=gr5,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("chocolate1"),alpha=0.8,show.legend=FALSE) +
  
  geom_ribbon(data=gr6,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("orange"),alpha=0.4,show.legend=FALSE) +
  geom_ribbon(data=gr6,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("orange"),alpha=0.8,show.legend=FALSE) +
  
  geom_ribbon(data=grall,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("red3"),alpha=0.4,show.legend=FALSE) +
  geom_ribbon(data=grall,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("red3"),alpha=0.8,show.legend=FALSE) +
  
  
  
  geom_point(data=Be7np,aes(x=obsx,y=obsy,group=dat,color=dat,shape=type),size=2.75,alpha=0.75)+
  geom_errorbar(show.legend = FALSE,aes(x=E,y=S,ymin=S-Stat,ymax=S+Stat,group=dat,color=dat),width=0.025)+
  scale_shape_stata(name="Dataset")+
  scale_color_stata(name="Dataset")+
  scale_x_log10() +
  #scale_y_log10() +
  theme_bw() + ylab(expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)"))) + 
  xlab("Energy (MeV)") + 
  theme(legend.position = "top",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  coord_cartesian(xlim=c(1e-8,1),ylim=c(0,10))























samplesmat = as.matrix(mcmcChain)
nsamp = nrow(samplesmat)
samplesmat2 <- samplesmat[,1:3]
library(sfsmisc)
library(plotrix)
library(emdbook)
library(magicaxis)

pdf("Final11b.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02,
    las=1)

xLim = c(4e-9, 10.0)
yLim = c(0,10)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xLim, ylim=yLim,
      axes=FALSE, main="", xlab = "", ylab = "", log="x", yaxs='i', xaxs='i' )
# control distance between axis and label [line=...]
title(xlab="Energy (MeV)", line=2, cex.lab=2.0)
title(ylab=expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)")), line=2, cex.lab=2.0)

# shade temperature region [low left corner to high right corner]
#rect(0.004, -2, 0.12, 32, col="snow2", border = NA)

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

magaxis(mgp=c(0,0.6,0), cex.axis=1.3)
box()

# plot reaction label
text(1e-7, 2, labels=expression(paste(NULL^"7","Be(n,p)",NULL^"7","Li")), cex=2.0)

# plot credible S-factors:
# calculate for these energies and the set of Bayesian samples for
# Er, g^2_i, g^2_f the S_factor curve using R code above

# define grid of energy values in MeV on a log scale for plotting of credible
# solutions; lseq is appropriate for log scale
x1 = lseq(1e-9, 5.0, length=500)

# ecm, E0, ga, gb, ra, rb, jr, la, lb

# 500 samples are randomly selected for plotting...

### ...individual resonances
##for ( i in round(runif(500, min=1, max=nsamp)) ) {
##   lines(x1, sqrt(x1)*(
##      Sigma7Benp(x1,samplesmat[i,1],samplesmat[i,2],samplesmat[i,3], 5, 5, 2, 0, 0)
##                      ),
##    col=adjustcolor("blue", alpha=0.02), lw=0.1)
##}


# ecm, E0, ga, gb, ra, rb, jr, la, lb
# 500 samples are randomly selected for plotting...
######### plot individual resonance curves
# resonance 1
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    samplesmat[i,22]*sigma7Benp(x1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15], 5, 5, 2, 0, 0)
  ),
  col=adjustcolor("blue", alpha=0.02), lw=0.1)
}

# resonance 2
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    samplesmat[i,23]*sigma7Benp(x1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16], 5, 5, 3, 1, 1)
  ),
  col=adjustcolor("black", alpha=0.02), lw=0.1)
}

# resonance 3
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    samplesmat[i,24]*sigma7Benp(x1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17], 5, 5, 3, 1, 1)
  ),
  col=adjustcolor("green4", alpha=0.02), lw=0.1)
}

# resonance 4
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    samplesmat[i,25]*sigma7Benp(x1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18], 5, 5, 1, 0, 0)
  ),
  col=adjustcolor("chocolate1", alpha=0.02), lw=0.1)
}

# resonance 5
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    samplesmat[i,26]*sigma7Benp(x1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19], 5, 5, 4, 3, 3)
  ),
  col=adjustcolor("violet", alpha=0.02), lw=0.1)
}

# resonance 6
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    samplesmat[i,27]*sigma7Benp(x1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20], 5, 5, 2, 1, 1)
  ),
  col=adjustcolor("orange", alpha=0.02), lw=0.1)
}


# resonance 7
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    samplesmat[i,28]*sigma7Benp(x1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21], 5, 5, 0, 1, 1)
  ),
  col=adjustcolor("skyblue", alpha=0.02), lw=0.1)
}


### sum of all resonances
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, sqrt(x1)*(
    samplesmat[i,22]*sigma7Benp(x1,samplesmat[i,1],samplesmat[i,8],samplesmat[i,15], 5, 5, 2, 0, 0)
    + samplesmat[i,23]*sigma7Benp(x1,samplesmat[i,2],samplesmat[i,9],samplesmat[i,16], 5, 5, 3, 1, 1)
    + samplesmat[i,24]*sigma7Benp(x1,samplesmat[i,3],samplesmat[i,10],samplesmat[i,17], 5, 5, 3, 1, 1)
    + samplesmat[i,25]*sigma7Benp(x1,samplesmat[i,4],samplesmat[i,11],samplesmat[i,18], 5, 5, 1, 0, 0)
    + samplesmat[i,26]*sigma7Benp(x1,samplesmat[i,5],samplesmat[i,12],samplesmat[i,19], 5, 5, 4, 3, 3)
    + samplesmat[i,27]*sigma7Benp(x1,samplesmat[i,6],samplesmat[i,13],samplesmat[i,20], 5, 5, 2, 1, 1)
    + samplesmat[i,28]*sigma7Benp(x1,samplesmat[i,7],samplesmat[i,14],samplesmat[i,21], 5, 5, 0, 1, 1)
  ) + samplesmat[i,29],
  col=adjustcolor("red", alpha=0.02), lw=0.1)
}


### plot legend
legend(0.5, 9.8, legend=c("Dam18", "Gib59" , "Mar19", "Koe88",
                          "Dam18", "Koe88", "Mar19", "Gib59", "Cer89", "Tom19"),
       pch=c(1, 2, 0, 4, 18, 19, 15, 17, 21, 22),
       pt.cex=c(1, 1, 1, 1, 1, 1, 1.1, 1.4))

### add data

## DATA SET 1: damone 2018;
points( obsx1, obsy1, col="gray40", pch=1, cex=1.2 )
add.error.bars(obsx1, obsy1, 0.0, errobsy1, 0.0, col="gray40" )

## DATA SET 2: gibbons & macklin 1959;
points( obsx2, obsy2, col="gray40", pch=2, cex=1.2 )
add.error.bars(obsx2, obsy2, 0.0, errobsy2, 0.0, col="gray40" )

## DATA SET 3: martin-hernandez 2018;
points( obsx3, obsy3, col="gray40", pch=0, cex=1.2 )
add.error.bars(obsx3, obsy3, 0.0, errobsy3, 0.0, col="gray40" )

## DATA SET 4: koehler 1988;
points( obsx4, obsy4, col="gray40", pch=4, cex=1.2 )
add.error.bars(obsx4, obsy4, 0.0, errobsy4, 0.0, col="gray40" )


## DATA SET 10: Damone 2018;
points( obsx10, obsy10, col="gray40", pch=18, cex=1.2 )
add.error.bars(obsx10, obsy10, 0.0, errobsy10, 0.0, col="gray40" )

## DATA SET 11: Koehler;
points( obsx11, obsy11, col="gray40", pch=19, cex=1.2 )
add.error.bars(obsx11, obsy11, 0.0, errobsy11, 0.0, col="gray40" )

## DATA SET 12: Hernandexs 2019;
points( obsx12, obsy12, col="gray40", pch=15, cex=1.2 )
add.error.bars(obsx12, obsy12, 0.0, errobsy12, 0.0, col="gray40" )

## DATA SET 13: Gibbons 1959;
points( obsx13, obsy13, col="gray40", pch=17, cex=1.2 )
add.error.bars(obsx13, obsy13, 0.0, errobsy13, 0.0, col="gray40" )

## DATA SET 14: cervena 1989;
points( obsx14, obsy14, col="gray40", pch=21, cex=1.2 )
add.error.bars(obsx14, obsy14, 0.0, errobsy14, 0.0, col="gray40" )

## DATA SET 15: tomandl 2019;
points( obsx15, obsy15, col="gray40", pch=22, cex=1.5 )
add.error.bars(obsx15, obsy15, 0.0, errobsy15, 0.0, col="gray40" )


### plot Maxwell-Boltzmann "factor" at T=1 GK [arbitrary scale]
for ( i in round(runif(500, min=1, max=nsamp)) ) {
  lines(x1, 3e2*(x1*(2.718^(-x1/(0.086173*0.5)))))
}

dev.off()
