N <- 500

#obsx1 <- runif(N,0,0.7)
obsx1 <- runif(N,0,10)
sd <- 2
obsy1 <- rnorm(N, 3 + 2*obsx1,sd=sd)

plot(obsx1,obsy1)


model.data <- list(obsy = obsy1,    # Response variable
                   obsx =  obsx1
)
#


model <- nimbleCode({
  for (i in 1:500) {
    obsy[i] ~ dnorm(mu[i],tau)
    mu[i] <- a + b*obsx[i]
    
  }
  
  tau ~  dunif(1e-3,10)
  a ~  dunif(1e-3,10)
  b ~  dunif(1e-3,10)
  sd <- 1/sqrt(tau)  
  
})
  inits <- list(a = runif(1,0.01,1),b=1,
              sd = runif(1,0.01,1)) 

Rmodel <- nimbleModel(code = model,data = model.data, inits = inits,check = FALSE)
compileNimble(Rmodel)

mcmcConf <- configureMCMC(Rmodel, monitors = c("a", "b","sd"))
mcmc_CL <- buildMCMC(mcmcConf)
CRmodel <- compileNimble(mcmc_CL,project = Rmodel)

mcmcChain <- runMCMC(CRmodel ,niter = 10000, nchain = 3, nburnin = 5000,
                     setSeed=15,samplesAsCodaMCMC = TRUE)


S <- ggs(mcmcChain)

ggs_histogram(S)
ggs_traceplot(S)