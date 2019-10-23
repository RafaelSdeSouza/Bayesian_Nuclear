# preparation: remove all variables from the work space
rm(list=ls())
# import jags package
library(nimble)
require(gsl)
require(RcppGSL)
require(nuclear)
require(forcats)
require(utils)
require(stats)
require(methods)
require(ggplot2)
require(ggthemes)
path <- getwd()
setwd(path)
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


Be7np <- read.csv("Be7np.csv")


Be7np$dat <- as.factor(Be7np$dat)
Be7np$dat <- fct_relevel(Be7np$dat, "Dam18","Gib59","Mar19","Koe88","Koe88b","Dam18b","Gib59b",
                                    "Her19","Cer89","Tom19")  

Be7np$type <- as.factor(Be7np$type)

re <- as.numeric(Be7np$dat) 
Nre <- length(unique(Be7np$dat))
# Unique removes duplicated vector, we want to know how many groups of
# data are there
N <- nrow(Be7np) # Total No of data sets
obsy <- Be7np$S    # Response variable in MeV
obsx <-  Be7np$E   # Predictors
erry <- Be7np$Stat # Error in MeV
set <- Be7np$dat # Get the labels as a vector
fu <- log(c(1.020,1.10,1.050,1.051,1.085,1.032))

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
       obsy[i] ~ dnorm(yt[i], pow(sqrt(erry[i]^2+y.scat[re[i]]^2),-2))
        yt[i] <- y.norm[re[i]]*(sqrt(obsx[i])*sigmaBe7modT[i] + hbg)
        }
  
 # for (i in 1:N){
#    obsy[i] ~ dnorm(yt[i], sd = sqrt(erry[i]^2 + y.scat[re[i]]^2))
#     obsy[i] ~ dnorm(yt[i], sd = sqrt(erry[i]^2 + y.scat^2))
#     yt[i] <- y.norm[re[i]]*sqrt(obsx[i])*sigmaBe7modT[i]
#  }
  


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


 
  hbg ~ dgamma(0.1,0.1)

  for (k in 1:4) {
    y.scat[k]  <- 0
  }
  for (k in 5:10) {
   y.scat[k] ~ dgamma(1,20)
  }

  
  
   for (k in 1:4) {
    # Systematic Uncertainty as a weakly informative prior
     nf[k] ~ dunif(-1,1)
     y.norm[k] <- 10^nf[k]
   }
  
  for (k in 5:10) {
    # Systematic Uncertainty as a highly informative prior
    y.norm[k] ~ dlnorm(log(1.0), pow(fu[k-4], -2))
  }
  


#  k ~ T(dnorm(0,sd=0.15),0,Inf)
  r_1 <- ga_1/gb_1
  r_4 <- ga_4/gb_4
})

samplerData <- list(obsy = obsy,
                    re = re, # Sample size
                    erry = erry,
                    obsx = obsx,
                    fu = fu)


samplerConst <- list(N = N,
                     Nre = Nre)


samplerInits <- list(y.norm = rep(1,Nre),
                     y.scat = rep(0.01,Nre),
                     nf = rep(0,4), 
                     y.norm = rep(1,Nre),
                     e0_1 = 0.05, gb_1 = 0.1, ga_1 = 0.1,
                     e0_2 = 0.15, ga_2 = 0.1, gb_2 = 0.1,
                     e0_3 = 0.336, ga_3 = 0.1, gb_3 = 0.1,
                     e0_4 = 0.51, ga_4 = 0.1, gb_4 = 0.1,
                     e0_5 = 0.96, ga_5 = 0.1,  gb_5 = 0.6,
                     e0_6 = 1.23, ga_6 = 0.1, gb_6 = 0.1,
                     e0_7 = 1.32, ga_7 = 0.1, gb_7 = 0.1,
                     ra = 4, rb = 4,
                     yt = rep(1,N),
                     hbg = 0.01
                     )



#nimbleOptions(oldConjugacyChecking = FALSE) 
#nimbleOptions(useNewConfigureMCMC = TRUE)

###############################################################
# (Alternative) If invoking Nimble MCMC stepwise (but more customisable)
##############################################################
ourmodel <- nimbleModel(code = samplerCode, constants = samplerConst,
                        data = samplerData, inits = samplerInits, check = FALSE
)

nimbleOptions(oldConjugacyChecking = FALSE)
nimbleOptions(useNewConfigureMCMC = TRUE)

compileNimble(ourmodel)
# Always compile the model after you are done setting up with it

ourmodel$calculate()
# Calculate the log likelihood (logProb). If the value is not NA,
# we have successfully initialised the model (FULLY)
# One iteration: simulate -> new values -> calculate

conf <- configureMCMC(ourmodel,print=TRUE,thin = 5)
# print = TRUE tells you what kind of samplers are being used for each stochastic node

conf$addMonitors(c('e0_1','ga_1','gb_1','e0_2','ga_2','gb_2',
                   'e0_3', 'ga_3', 'gb_3', 'e0_4', 'ga_4', 'gb_4',
                   'e0_5', 'ga_5', 'gb_5', 'e0_6', 'ga_6', 'gb_6',
                   'e0_7', 'ga_7', 'gb_7',
                   'r_1', 'r_4', 'ra', 'rb',
                   'y.norm','y.scat','hbg'))


conf$removeSampler(c('e0_1','ga_1','ga_2','ga_3','gb_1','gb_2','ga_4','gb_3', 'gb_4',
                     'ra', 'rb'))

conf$addSampler(target = c('e0_1','ga_1','ga_2','ga_3','gb_1','gb_2','ga_4','gb_3', 'gb_4',
                           'ra', 'rb'),
               type = "AF_slice")


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
n.iter =   1500
n.burnin = 1000

system.time(
  mcmcChain <- runMCMC(compiledMCMC,niter = n.iter, nchains = n.chains, nburnin = n.burnin,
                       samplesAsCodaMCMC = TRUE)
)

samp <- as.matrix(mcmcChain)


write.csv(samp,"MCMC_ApJ_ultimaterun.csv")

pdf("Chains_slice_ApJ.pdf")
plot(mcmcChain)
dev.off()




xx <- exp(seq(log(1e-9),log(3),length.out = 400))
# resonance 1
y1 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_1"],samp [,"ga_1"],samp [,"gb_1"],
                                          samp [,"ra"],
                                          samp [,"rb"],jr = 2, la = 0, lb = 0),probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y1 <- rbind(y1,mux)
}
gr1 <-  data.frame(x =xx, mean = y1[,"50%"],lwr1=y1[,"25%"],lwr2=y1[,"2.5%"],
                   upr1=y1[,"75%"],
                   upr2=y1[,"97.5%"])
# resonance 2
y2 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_2"],samp[,"ga_2"],samp [,"gb_2"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 3, la = 1, lb = 1),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y2 <- rbind(y2,mux)
}

gr2 <-  data.frame(x =xx, mean = y2[,"50%"],lwr1=y2[,"25%"],lwr2=y2[,"2.5%"],
                   upr1=y2[,"75%"],
                   upr2=y2[,"97.5%"])

# resonance 3
y3 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_3"],samp[,"ga_3"],samp [,"gb_3"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 3, la = 1, lb = 1),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y3 <- rbind(y3,mux)
}

gr3 <-  data.frame(x =xx, mean = y3[,"50%"],lwr1=y3[,"25%"],lwr2=y3[,"2.5%"],
                   upr1=y3[,"75%"],
                   upr2=y3[,"97.5%"])


# resonance 4
y4 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_4"],samp[,"ga_4"],samp [,"gb_4"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 1, la = 0, lb = 0),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y4 <- rbind(y4,mux)
}

gr4 <-  data.frame(x =xx, mean = y4[,"50%"],lwr1=y4[,"25%"],lwr2=y4[,"2.5%"],
                   upr1=y4[,"75%"],
                   upr2=y4[,"97.5%"])

# resonance 5
y5 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_5"],samp[,"ga_5"],samp [,"gb_5"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 4, la = 3, lb = 3),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y5 <- rbind(y5,mux)
}

gr5 <-  data.frame(x =xx, mean = y5[,"50%"],lwr1=y5[,"25%"],lwr2=y5[,"2.5%"],
                   upr1=y5[,"75%"],
                   upr2=y5[,"97.5%"])


# resonance 6
y6 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_6"],samp[,"ga_6"],samp [,"gb_6"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 2, la = 1, lb = 1),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y6  <- rbind(y6,mux)
}

gr6 <-  data.frame(x =xx, mean = y6[,"50%"],lwr1=y6[,"25%"],lwr2=y6[,"2.5%"],
                   upr1=y6[,"75%"],
                   upr2=y6[,"97.5%"])


# resonance 7
y7 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_7"],samp[,"ga_7"],samp [,"gb_7"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 0, la = 1, lb = 1),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y7  <- rbind(y7,mux)
}

gr7 <-  data.frame(x =xx, mean = y7[,"50%"],lwr1=y7[,"25%"],lwr2=y7[,"2.5%"],
                   upr1=y7[,"75%"],
                   upr2=y7[,"97.5%"])



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
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  yall <- rbind(yall,mux)
}

grall <-  data.frame(x =xx, mean = yall[,"50%"],lwr1=yall[,"25%"],lwr2=yall[,"2.5%"],
                     upr1=yall[,"75%"],
                     upr2=yall[,"97.5%"])



### plot Maxwell-Boltzmann "factor" at T=1 GK [arbitrary scale]
#for ( i in round(runif(500, min=1, max=nsamp)) ) {
#  lines(x1, 3e2*(x1*(2.718^(-x1/(0.086173*0.5)))))
#}

MB <- function(x1){3e2*(x1*(2.718^(-x1/(0.086173*0.5))))}
MBD <- data.frame(x=xx,y=MB(xx))

Be7npG <-  Be7np 
library(plyr)
Be7npG$dat <- revalue(Be7npG$dat, c("Koe88b"="Koe88","Dam18b"="Dam18",Gib59b="Gib59"))

# Ressonances

dr=data.frame(x=c(0,0.15,0.34,0.51,0.96,1.23,1.32), 
              y=rep(10,7), 
              vy=rep(9.5,7))

# Plot all
pdf("Be7_ApJ_ultimate.pdf", width=7.5, height=5)
ggplot(Be7npG,aes(x=E,y=S)) +
  
  
  geom_area(data=MBD,aes(x=x,y=y),color="#a6cee3",fill="#a6cee3",
            size=0,alpha=0.4) +
  
  #
  geom_ribbon(data=gr1,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#739a43"),alpha=0.6,show.legend=FALSE) +
  geom_ribbon(data=gr1,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#739a43"),alpha=0.95,show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 5.75+2, yend = 5.75+2,size=1.5,
           colour = "#739a43") +
  annotate(geom="text", x = 1.25,  y = 5.75+2,
           label=expression(2^"-"*", 0.00"),size=4) +
  
  geom_ribbon(data=gr2,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#6780d8"),alpha=0.6,show.legend=FALSE) +
  geom_ribbon(data=gr2,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#6780d8"),alpha=0.95,show.legend=FALSE) +
  annotate("segment",  x = 0.2, xend = 0.35, y = 5.25+2, yend = 5.25+2,size=1.5,
           colour = "#6780d8") +
  annotate(geom="text", x = 1.25,  y = 5.25+2,
           label=expression(3^"+"*", 0.15"),size=4) +
  
  
  geom_ribbon(data=gr3,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#b84c7d"),alpha=0.6,show.legend=FALSE) +
  geom_ribbon(data=gr3,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#b84c7d"),alpha=0.95,show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 4.75+2, yend = 4.75+2,size=1.5,
           colour = "#b84c7d") +
  annotate(geom="text", x = 1.25,  y = 4.75+2,
           label=expression(3^"+"*", 0.34"),size=4) +
  
  
  geom_ribbon(data=gr4,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#df7357"),alpha=0.6,show.legend=FALSE) +
  geom_ribbon(data=gr4,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#df7357"),alpha=0.95,show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 4.25+2, yend = 4.25+2,size=1.5,
           colour = "#df7357") +
  annotate(geom="text", x = 1.25,  y = 4.25+2,
           label=expression(1^"-"*", 0.51"),size=4) +
  
  
  
  geom_ribbon(data=gr6,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#c1bd5c"),alpha=0.6,show.legend=FALSE) +
  geom_ribbon(data=gr6,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#c1bd5c"),alpha=0.95,show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 3.25+2, yend = 3.25+2,size=1.5,
           colour = "#c1bd5c") +
  annotate(geom="text",  x = 1.25,  y = 3.25+2,
           label=expression(2^"+"*", 1.23"),size=4) +
  
  
  geom_ribbon(data=gr5,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#5d398d"),alpha=0.6,show.legend=FALSE) +
  geom_ribbon(data=gr5,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#5d398d"),alpha=0.95,show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 3.75+2, yend = 3.75+2,size=1.5,
           colour = "#5d398d") +
  annotate(geom="text", x = 1.25,  y = 3.75+2,
           label=expression(4^"+"*", 0.96"),size=4) +
  
  
  
  geom_ribbon(data=gr7,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#c19c3c"),alpha=0.6,show.legend=FALSE) +
  geom_ribbon(data=gr7,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#c19c3c"),alpha=0.95,show.legend=FALSE) +
  annotate("segment",x = 0.2, xend = 0.35, y = 2.75+2, yend = 2.75+2,size=1.5,
           colour = "#c19c3c") +
  annotate(geom="text",  x = 1.25,  y = 2.75+2,
           label=expression(0^"+"*", 1.32"),size=4) +
  
  geom_ribbon(data=grall,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("black"),alpha=0.5,show.legend=FALSE) +
  geom_ribbon(data=grall,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("black"),alpha=0.85,show.legend=FALSE) +
  
  
  
  geom_errorbar(show.legend = FALSE,aes(x=E,y=S,ymin=S-Stat,ymax=S+Stat,group=dat,color=type),alpha=0.75,width=0.025)+
  geom_point(data=Be7npG,aes(x=obsx,y=obsy,group=dat,shape=dat,color=type,fill=type),size=2)+
  
  scale_shape_manual(values=c(22,4,3,24,25,23,21),name="",
                     guide = guide_legend(nrow = 2))+
  scale_color_manual(name="",values=c("red","cyan3"),guide="none")+
  scale_alpha_manual(name="",values=c(1,0.5),guide="none")+
  scale_fill_manual(values=c("red","cyan3"),name="",guide="none") +
  scale_x_log10(breaks = c(1e-6,1e-3,1),
                labels=c(expression(10^-6),expression(10^-3),"1")) +
  theme_economist_white() + ylab(expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)"))) + 
  xlab("Energy (MeV)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.45,0.95),
        legend.text = element_text(size=13),
        legend.text.align = 1,
        
        legend.background = element_rect(colour = "white",fill=NA),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="black", size=17.5),
        axis.text  = element_text(size=16),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  coord_cartesian(xlim=c(1e-8,2.25),ylim=c(0.175,9.75)) +
  annotate(geom="text",1e-7, 1,
           label=expression(paste(NULL^"7","Be(n,p)",NULL^"7","Li")),
           size=6) 
dev.off()




