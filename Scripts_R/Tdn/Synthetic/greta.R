library(greta)
require(gsl)
source("Sfactor3.R")
N <- 150


# data
obsx1 <- exp(runif(N,log(1e-3),log(1)))

# Artificial values, here we will just fit 3 parameters Er, gi, gf
y <- rnorm(N, sfactorTdn_5p(obsx1,0.0912,2.93,0.0794,6,5),sd = 1)


# variables and priors
Er = uniform(min = 0, max = 50)
gi = uniform(min = 0, max = 50)
gf = uniform(min = 0, max = 50)
sd =  uniform(min = 0, max = 50)
# operations

if (is.na(as.matrix(Er)[[1]])){
  xx <- 0.0912  
  mean = sfactorTdn_5p(obsx1,xx,gi,gf,6,5)
}else{
  mean = sfactorTdn_5p(obsx1,Er,gi,gf,6,5) 
}



# likelihood
distribution(y) = normal(mean, sd)

# defining the model
m <- model(Er, gi, gf,sd)

# plotting
plot(m)

# sampling
draws <- mcmc(m, warmup = 500,n_samples = 1000,chains = 1)


library (bayesplot)
mcmc_trace(draws)
mcmc_intervals(draws)