library(greta)
require(gsl)
source("Sfactor3.R")

## Read DATA 
ensamble <- read.csv("ensamble_Tdn.csv",header = T)

# data
N <- nrow(ensamble)
obsy <- as_data(ensamble$S)    # Response variable
obsx <- ensamble$E[]   # Predictors
erry <- as_data(ensamble$Stat)
set <- ensamble$dat
lab <- ensamble$invK
syst = c(unique(ensamble$Syst))



# variables and priors
Er = normal(0, 5,truncation = c(0, 5))
Ex = normal(0, 5,truncation = c(0, 5))
gi = normal(0, 10,truncation = c(0, Inf))
gf = normal(0, 10,truncation = c(0, Inf))
sd = normal(0, 10,truncation = c(0, Inf))
epsilon = normal(0, 10)

# operations

mean <- sfactorTdn_5p(obsx,Er,gi,gf,6,5)
#distribution(y) = normal(mean,sd)
distribution(obsy) = normal(mean,sd)


# defining the model
m <- model(Er, gi, gf,sd)

# plotting
plot(m)

# sampling
draws <- mcmc(m, warmup = 1000,n_samples = 5000,chains = 3)

summary(draws)
library (bayesplot)
mcmc_trace(draws)
mcmc_intervals(draws)