library(nimble)

######################################################################
# GENERATE ARTIFICIAL DATA
######################################################################
# sample data input/generation; the input is of the form: obsx, obsy, 
# errobsy, where the latter is the individual error of each datum [i]

N <- 500
## DATA SET
# x values are sampled uniformly
obsx <- runif(N, 0, 10)

# y values:
# - mean=linear model of slope=0.5 and intercept=0
# - y errors are drawn from a uniform distribution with sd=individual 
#   errors
# - finally, option to shift all data by a normalization factor [for
#   offset, change '*' into '+'

alpha = 0
beta  = 0.5
ynorm = 1.1
z = alpha + beta*obsx  
y = ynorm*(z)
errobsy <-  runif(N, 0.2, 0.5)
obsy =  rnorm(N, y, errobsy)


Modeldata <- list(errobsy = errobsy,
                  obsy = obsy,
                  obsx = obsx
)

# Nimble Code
NormalFit <-  nimbleCode({
    for (i in 1:500) {
      z[i] <- alpha + beta*obsx[i];
      y[i] <- ynorm*z[i];
      obsy[i] ~ dnorm(y[i], sd = pow(errobsy[i],-2));
      }
  alpha ~ dnorm(0,1e-2)
  beta ~  dnorm(0,1e-2)
  ynorm ~ dlnorm(0,pow(log(1.2),-2))
})


inits <- list(alpha = 1,
              beta = 1,
              ynorm = 1)

model <- nimbleModel(code = NormalFit, name = "NormalFit",
                         data = Modeldata,inits = inits )

mcmc.out <- nimbleMCMC(code = NormalFit, 
                       data = Modeldata,
                       nchains = 10, niter = 50000,nburnin = 20000,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c('alpha','beta','ynorm'))

mcmc.out$summary


