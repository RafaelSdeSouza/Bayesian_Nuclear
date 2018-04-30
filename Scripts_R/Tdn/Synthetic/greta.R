library(greta)
require(gsl)
source("Sfactor3.R")
N <- 150


# data
obsx1 <- exp(runif(N,log(1e-3),log(1)))
sd <- 1
# Artificial values, here we will just fit 3 parameters Er, gi, gf
y <- rnorm(N, Sfactor3(obsx1,0.0912,0.0912,2.93,0.0794,6,5,0),sd = sd)


# variables and priors
Er = normal(0, 1,truncation = c(0, Inf))
gi = normal(0, 1,truncation = c(0, Inf))
gf = normal(0, 1,truncation = c(0, Inf))
sd = normal(0, 1,truncation = c(0, Inf))
# operations
mean <- Sfactor3(obsx1,Er,0.0912,gi,gf,6,5,0)

# likelihood
distribution(y) = normal(mean, sd)

# defining the model
m <- model(Er, gi, gf,sd)

# plotting
plot(m)

# sampling
draws <- mcmc(m, warmup = 5000,n_samples = 10000,chains = 3)


library (bayesplot)
mcmc_trace(draws)
mcmc_intervals(draws)