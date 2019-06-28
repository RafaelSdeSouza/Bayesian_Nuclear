# 3Hedp analysis
#
# purpose: Real  DATA
#
# - 5 parameters are assumed: Er, gamma_d^2, gamma_n^2 [e1, gin, gout]
#
# - uses the function sfactorHe3dp(obsx1[i], e1, gin, gout), which
#   is a C++ version of a Fortran code that includes Coulomb wave
#   function calculations; JAGS has been recompiled with this C++ function
#
######################################################################
# preparation: remove all variables from the work space
#rm(list=ls())
set.seed(123)
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb;

######################################################################
# import packages

require(RcppGSL);require(ggplot2);require(ggthemes)
require(nuclear);library(magrittr);library(wesanderson)
library(dplyr);require(gsl);library(latex2exp)
require(nimble)
######################################################################
## ARTIFICIAL DATA GENERATION

data("He3dp")

re <- as.numeric(as.factor(He3dp$dat))
Nre <- length(unique(He3dp$dat))
ik <- as.numeric(as.factor(He3dp$Kinematics))
Nik <- length(unique(He3dp$Kinematics))

N <- nrow(He3dp)
obsy <- He3dp$S    # Response variable
obsx <-  He3dp$E   # Predictors
erry <- He3dp$Stat
set <- as.factor(He3dp$dat)
lab <- He3dp$Kinematics
syst <- as.numeric(unlist(exp(aggregate(He3dp$Syst, by=list(He3dp$dat), FUN=mean)[2])))

M <- 500
xx <- seq(min(obsx),max(obsx),length.out = M)


model.data <- list(obsy = obsy,    # Response variable
                   obsx =  obsx,   # Predictors
                   erry = erry,
                   syst = syst,
                   xx = xx
)

samplerConst <- list(N = N, # Sample size
                     M = M,
                     Nre = Nre, 
                     re = re, # This is used to "iterate"
                     Nik = Nik,
                     ik  = ik
)


library(compiler)
sfactorHe3dpc <- cmpfun(sfactorHe3dp)

  
sfactorHe3dpNimble <- nimbleRcall(function(ecm = double(0),
                              e0 = double(0),gi = double(0),
                              gf = double(0),ri = double(0),
                              rf = double(0),ue = double(0)){},
                              Rfun = 'sfactorHe3dpc',
                              returnType = double(0))



library(microbenchmark)
res <- microbenchmark(sfactorHe3dpNimble(obsx,0.35,1,0.02,6,5,0),
               sfactorHe3dpc(obsx,0.35,1,0.02,6,5,0),
               sfactorHe3dp(obsx,0.35,1,0.02,6,5,0),times=1000L)



model <- nimbleCode({
   for (i in 1:N) {
    obsy[i] ~ dnorm(y[i],pow(erry[i], -2))
     y[i] ~ dnorm(scale[re[i]]*sfactorHe3dpNimble(obsx[i], E0,  gd2, gp2, ad, ap, ue[ik[i]]), pow(tau, -2))
}

  # Predicted values
  for (j in 1:M){

    # Bare...

    mux0[j] <- sfactorHe3dpNimble(xx[j], E0,  gd2, gp2, ad, ap,0)

    # No inverse Kinematics

    mux1[j] <- sfactorHe3dpNimble(xx[j], E0,  gd2, gp2, ad, ap,ue[1])
    yx1[j] ~ dnorm(mux1[j],pow(tau,-2))

    # With inverse Kinematics
    mux2[j] <- sfactorHe3dpNimble(xx[j], E0,  gd2, gp2, ad, ap,ue[2])
    yx2[j] ~ dnorm(mux1[j],pow(tau,-2))

  }

  for (k in 1:Nre){
    scale[k] ~ dlnorm(log(1.0),pow(syst[k],-2))
  }

  for (z in 1:Nik){
    ue[z] ~ T(dnorm(0,pow(100,-2)),0,Inf)
  }

  # PRIORS

  # Case I
  tau ~  dgamma(0.1,0.1)
  E0  ~  dgamma(0.1,0.1)
  gd2 ~  dgamma(0.1,0.1)
  gp2 ~  dgamma(0.1,0.1)
  ad  ~ T(dnorm(3.5,pow(0.5,-2)),0,Inf)
  ap  ~ T(dnorm(5.5,pow(1,-2)),0,Inf)
#  ue_ev[1] <-1e6*ue[1]
#  ue_ev[2] <-1e6*ue[2]
})
inits <- list(E0  = runif(1,0.01,1),gd2=0.01,gp2=runif(1,0.01,1),
              ad = 5, ap = 5, ue = c(100,100),scale = runif(7,0.9,1.1),
              tau  = runif(1,0.01,1),
              y =  sfactorHe3dp(obsx,0.35,1,0.02,6,5,0),
              yx1 = sfactorHe3dp(xx,1,1,2,6,5,0.001),
              yx2 = sfactorHe3dp(xx,2,6,5,0.001))


Rmodel <- nimbleModel(code = model,data = model.data,constants = samplerConst,
                      inits = inits,check = FALSE)
compileNimble(Rmodel)

mcmcConf <- configureMCMC(Rmodel,
                          monitors = c("E0","gd2", "gp2","ue","tau", "ad","ap",
                                               "mux0","mux1","mux2","scale"))
mcmc_CL <- buildMCMC(mcmcConf)
CRmodel <- compileNimble(mcmc_CL,project = Rmodel)


mcmcChain <- runMCMC(CRmodel,niter = 30000, nchains = 1, nburnin = 10000,samplesAsCodaMCMC = TRUE)




S <- ggs(mcmcChain[,c("E0","gd2","gp2","ue[1]","ue[2]")])


ggs_traceplot(S)






mcmc.output <- nimbleMCMC(Rmodel, data = model.data, inits = inits,
                          monitors = c("e1", "gin", "gout","sd"), thin = 10,
                          niter = 20000, nburnin = 1000, nchains = 3,
                          summary = TRUE, WAIC = TRUE)





mcmcConf <- configureMCMC(Rmodel, monitors = c("e1", "gin", "gout","sd"))
mcmc_CL <- buildMCMC(mcmcConf)
samplesList <- runMCMC(mcmc_CL, niter = 50, nchains = 3, inits = inits)


Rmcmc <- buildMCMC(model)
Rmodel$setData(model.data)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)






