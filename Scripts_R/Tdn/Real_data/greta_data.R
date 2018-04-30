library(greta)
require(gsl)
source("Sfactor3.R")
N <- 150

## Read DATA 
ensamble <- read.csv("ensamble_Tdn.csv",header = T)

# data
N <- nrow(ensamble)
obsy <- ensamble$S    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat
set <- ensamble$dat
lab <- ensamble$invK
syst = c(unique(ensamble$Syst))



# variables and priors
Er = normal(0, 1,truncation = c(0, Inf))
gi = normal(0, 1,truncation = c(0, Inf))
gf = normal(0, 1,truncation = c(0, Inf))
sd = normal(0, 1,truncation = c(0, Inf))
# operations

Model <- 
# LIKELIHOOD
obsy[i] ~ dnorm(y[i], pow(erry[i], -2))
y[i] ~ dnorm(scale[re[i]]*sfactorTdn(obsx[i], e1,ex, gin, gout,ri,rf,0),pow(tau, -2))
res[i] <- obsy[i]-sfactorTdn(obsx[i], e1,ex, gin, gout,ri,rf,0)

# LIKELIHOOD
mean <- Sfactor3(obsx1,Er,0.0912,gi,gf,6,5,0)
distribution(obsy) = normal(y, erry)

# defining the model
m <- model(Er, gi, gf,sd)

# plotting
plot(m)

# sampling
draws <- mcmc(m, warmup = 5000,n_samples = 10000,chains = 3)


library (bayesplot)
mcmc_trace(draws)
mcmc_intervals(draws)