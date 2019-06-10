Tgrid = c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,
          0.013,0.014,0.015,0.016,0.018,0.020,0.025,0.030,0.040,0.050,0.060,0.070,
          0.080,0.090,0.100,0.110,0.120,0.130,0.140,0.150,0.160,0.180,0.200,0.250,0.300,
          0.350,0.400,0.450,0.500,0.600,0.700,0.800,0.900,1.000,1.250,1.5,1.75,2,
          2.5,3,3.5,4,5,6,7,8,9,10)

numRates_ddp_table <- function(mcmcobj, vars, N = 1000){
  mcdat_I <- as.data.frame(unlist(as.matrix(mcmcobj)[,vars])) 
  # In context of DPN-jags, we are extracting a.scale
  # 1) Convert to mcmc "matrix", extract "vars" column
  # 2) Unlist to combine them all
  # 3) Convert these into a data.frame
  
  index <- sample(1:nrow(mcdat_I),size=N,replace=FALSE)
  # Random Sample of 1000 data sets to do integration
  a  <- mcdat_I[index,]

  gdat <-  sapply(Tgrid,function(Tgrid){numRates_ddp(a,T9 = Tgrid)})
  # Multiple Integration
  # Different temperature at different columns 
  # Recall that numRates_ddp does integration and return a list of the different reaction rates
  # for EACH chain
  
  gg <-  as.data.frame(gdat)
  
  gg2 <- apply(gg, 2, quantile, probs=c(0.16, 0.5, 0.84), na.rm=TRUE)
  # 2nd argument of apply refers to margin -> i.e 2 = Column
  
  fu <- function(x){exp(sqrt(log(1+var(x)/mean(x)^2)))}
  ## f.u = exp(sigma), which the formula is obtained from equation (3) of Bayesian estimation of 
  # Thermonuclear Reaction Rates
  
  fu_I <- apply(gg, 2, fu)
  
  gg2data <- data.frame(T9 = Tgrid, mean = gg2["50%",],lower = gg2["16%",], 
                        upper = gg2["84%",])
  gg2data$fu <- fu_I
  rownames( gg2data) <- NULL
  return(gg2data)
}

