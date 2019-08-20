## Integration Script

## Load your mcmcChain here
load(file = "Final3.RData")

## Functions to be used to calculate rate at a given temperature

sigma7Benp <- function(ecm, e0, ga, gb, ra, rb, jr, la, lb){
  # input masses, charges, angular momenta
  m1_i = 7.0147344 
  m2_i = 1.0086649158   # masses (amu) of 7Be and n
  m1_f = 7.014357697 
  m2_f = 1.007276452    # masses (amu) of 7Li and p
  z1_i = 4 
  z2_i = 0              # charges of 7Be and n
  z1_f = 3 
  z2_f = 1              # charges of 7Li and p
  jt = 1.5              # spins of target, projectile
  jp = 0.5 
  Q = 1.6442402         # reaction Q-value (MeV)
  
  # reduced masses
  mue_i <- (m1_i*m2_i)/(m1_i+m2_i)
  mue_f <- (m1_f*m2_f)/(m1_f+m2_f)
  
  # constants
  pek <- 6.56618216e-1/mue_i
  omega <- (2*jr+1)/((2*jt+1)*(2*jp+1))
  
  ### CALCULATE S-FACTOR
  ## incoming channel   
  eta_a=0.1574854*z2_i*z1_i*sqrt(mue_i)
  rho_a=0.218735*ra*sqrt(mue_i)      
  eta_i=eta_a/(sqrt(ecm))
  rho_i=rho_a*(sqrt(ecm))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k=0)
  # penetration and shift factor 
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  s_i <- rho_i*(P3$val_F*P3$val_Fp + P3$val_G*P3$val_Gp)/(P3$val_F^2 + P3$val_G^2)
  # shift factor at energy E0 [eigenvalue]
  xeta_i=eta_a/(sqrt(e0))
  xrho_i=rho_a*(sqrt(e0))
  PX1 <- coulomb_wave_FG(xeta_i, xrho_i, la, k=0)
  b_i <- xrho_i*(PX1$val_F*PX1$val_Fp + PX1$val_G*PX1$val_Gp)/(PX1$val_F^2 + PX1$val_G^2)
  # partial width
  Ga <- 2*ga*p_i
  
  ## outgoing channel
  eta_b=0.1574854*z2_f*z1_f*sqrt(mue_f)
  rho_b=0.218735*rb*sqrt(mue_f)      
  eta_f=eta_b/(sqrt(ecm+Q))
  rho_f=rho_b*(sqrt(ecm+Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k=0)
  # penetration and shift factor
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  s_f <- rho_f*(P4$val_F*P4$val_Fp + P4$val_G*P4$val_Gp)/(P4$val_F^2 + P4$val_G^2)
  # shift factor at energy E0+Q
  xeta_f=eta_b/(sqrt(e0+Q))
  xrho_f=rho_b*(sqrt(e0+Q)) 
  PX2 <- coulomb_wave_FG(xeta_f, xrho_f, lb, k=0)
  b_f <- xrho_f*(PX2$val_F*PX2$val_Fp + PX2$val_G*PX2$val_Gp)/(PX2$val_F^2 + PX2$val_G^2)
  # partial width
  Gb <- 2*gb*p_f
  
  tapp <- (s_i-b_i)*ga+(s_f-b_f)*gb
  
  s1=pek*omega*Ga*Gb
  s2=((e0-ecm-tapp)^2)+0.25*((Ga+Gb)^2)
  SF <- (s1/s2)*(1/ecm)
  
  return(SF = SF)
}

sigma7Benp7mod <- function(ecm, 
                           e0_1, ga_1, gb_1, ra_1, rb_1, jr_1, la_1, lb_1, i_1,
                           e0_2, ga_2, gb_2, ra_2, rb_2, jr_2, la_2, lb_2, i_2,
                           e0_3, ga_3, gb_3, ra_3, rb_3, jr_3, la_3, lb_3, i_3,
                           e0_4, ga_4, gb_4, ra_4, rb_4, jr_4, la_4, lb_4, i_4,
                           e0_5, ga_5, gb_5, ra_5, rb_5, jr_5, la_5, lb_5, i_5,
                           e0_6, ga_6, gb_6, ra_6, rb_6, jr_6, la_6, lb_6, i_6,
                           e0_7, ga_7, gb_7, ra_7, rb_7, jr_7, la_7, lb_7, i_7){
  
  SF1 <- if(i_1 == 0) 0 else sigma7Benp(ecm, e0_1, ga_1, gb_1, ra_1, rb_1, jr_1, la_1, lb_1)
  SF2 <- if(i_2 == 0) 0 else sigma7Benp(ecm, e0_2, ga_2, gb_2, ra_2, rb_2, jr_2, la_2, lb_2)
  SF3 <- if(i_3 == 0) 0 else sigma7Benp(ecm, e0_3, ga_3, gb_3, ra_3, rb_3, jr_3, la_3, lb_3)
  SF4 <- if(i_4 == 0) 0 else sigma7Benp(ecm, e0_4, ga_4, gb_4, ra_4, rb_4, jr_4, la_4, lb_4)
  SF5 <- if(i_5 == 0) 0 else sigma7Benp(ecm, e0_5, ga_5, gb_5, ra_5, rb_5, jr_5, la_5, lb_5)
  SF6 <- if(i_6 == 0) 0 else sigma7Benp(ecm, e0_6, ga_6, gb_6, ra_6, rb_6, jr_6, la_6, lb_6)
  SF7 <- if(i_7 == 0) 0 else sigma7Benp(ecm, e0_7, ga_7, gb_7, ra_7, rb_7, jr_7, la_7, lb_7)
  SF <- SF1 + SF2 + SF3 + SF4 + SF5 + SF6 + SF7
  return(SF = SF)
}



numRates_7Benp <- Vectorize(function(e0_1, ga_1, gb_1, ra_1, rb_1, jr_1, la_1, lb_1, i_1,
                                     e0_2, ga_2, gb_2, ra_2, rb_2, jr_2, la_2, lb_2, i_2,
                                     e0_3, ga_3, gb_3, ra_3, rb_3, jr_3, la_3, lb_3, i_3,
                                     e0_4, ga_4, gb_4, ra_4, rb_4, jr_4, la_4, lb_4, i_4,
                                     e0_5, ga_5, gb_5, ra_5, rb_5, jr_5, la_5, lb_5, i_5,
                                     e0_6, ga_6, gb_6, ra_6, rb_6, jr_6, la_6, lb_6, i_6,
                                     e0_7, ga_7, gb_7, ra_7, rb_7, jr_7, la_7, lb_7, i_7,
                                     K, T9){
  
  
  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------
  
  integrand <- function(E,T9) {(sqrt(E) * sigma7Benp7mod(E,e0_1, ga_1, gb_1, ra_1, rb_1, jr_1, la_1, lb_1, i_1,
                                                        e0_2, ga_2, gb_2, ra_2, rb_2, jr_2, la_2, lb_2, i_2,
                                                        e0_3, ga_3, gb_3, ra_3, rb_3, jr_3, la_3, lb_3, i_3,
                                                        e0_4, ga_4, gb_4, ra_4, rb_4, jr_4, la_4, lb_4, i_4,
                                                        e0_5, ga_5, gb_5, ra_5, rb_5, jr_5, la_5, lb_5, i_5,
                                                        e0_6, ga_6, gb_6, ra_6, rb_6, jr_6, la_6, lb_6, i_6,
                                                        e0_7, ga_7, gb_7, ra_7, rb_7, jr_7, la_7, lb_7, i_7) + K)*
      exp(-E/(0.086173324*T9))}
  
  # CALCULATE Nuclear rate
  
  m1 = 7.0147344 
  m2 = 1.0086649158   # masses (amu) of 7Be and n
  mue = (m1*m2)/(m1+m2)
  
  Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 1e-10, upper = 2,
                                                                      abs.tol = 0L,
                                                                      T9 = Temp)$value}
  
  # Note to self, the limits of integration, in some sense, the scale should be appropriate.
  # From HELP, the first argument MUST BE integrated. The optional argument T9 is used to be substituted in
  # Nasv <-> N A <sigma v >
  
  out <- Nasv(T9)
  return(Nasv=out)
}
)

## Single Temperature value:
# Important Temperature 0.5 G K

# Chain reduction to 2000 iterations (randomly sampled) if needed
index <- sample(1:nrow(mcmcChain),size=2000,replace=FALSE)
mcmcChaint <- mcmcChain[index,]

NRate <- numRates_7Benp(mcmcChaint[,"e0_1"], mcmcChaint[,"ga_1"], mcmcChaint[,"gb_1"], 5, 5, 2, 0, 0, mcmcChaint[,"i_1"],
               mcmcChaint[,"e0_2"], mcmcChaint[,"ga_2"], mcmcChaint[,"gb_2"], 5, 5, 3, 1, 1, mcmcChaint[,"i_2"],
               mcmcChaint[,"e0_3"], mcmcChaint[,"ga_3"], mcmcChaint[,"gb_3"], 5, 5, 3, 1, 1, mcmcChaint[,"i_3"],
               mcmcChaint[,"e0_4"], mcmcChaint[,"ga_4"], mcmcChaint[,"gb_4"], 5, 5, 1, 0, 0, mcmcChaint[,"i_4"],
               mcmcChaint[,"e0_5"], mcmcChaint[,"ga_5"], mcmcChaint[,"gb_5"], 5, 5, 4, 3, 3, mcmcChaint[,"i_5"],
               mcmcChaint[,"e0_6"], mcmcChaint[,"ga_6"], mcmcChaint[,"gb_6"], 5, 5, 2, 1, 1, mcmcChaint[,"i_6"],
               mcmcChaint[,"e0_7"], mcmcChaint[,"ga_7"], mcmcChaint[,"gb_7"], 5, 5, 0, 1, 1, mcmcChaint[,"i_7"],
               mcmcChaint[,"k"],0.5)

quan <- quantile(NRate, probs = c(0.16,0.5,0.84))



## Latex Table Output

require(xtable)
#write.csv(Nrate,"rates_dpg.csv",row.names = F)
print(xtable(Nrate, type = "latex",display=c("e","g","E","E","E",
                                             "g"),
             digits=4), include.rownames = FALSE)


## If we are obtaining the rates at different temperatures of interest
Tgrid = c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,
          0.013,0.014,0.015,0.016,0.018,0.020,0.025,0.030,0.040,0.050,0.060,0.070,
          0.080,0.090,0.100,0.110,0.120,0.130,0.140,0.150,0.160,0.180,0.200,0.250,0.300,
          0.350,0.400,0.450,0.500,0.600,0.700,0.800,0.900,1.000,1.250,1.5,1.75,2,
          2.5,3,3.5,4,5,6,7,8,9,10)

numRates_7Benp_table <- function(mcmcobj, N = 1000){
  # Condition: The variables are named as such:
  # e0_1, e0_2, ..., ga_1, gb_1, ..., i_1, ..., i_7, k, 
  mcdat <- as.matrix(mcmcobj) 
  # In context of DPG-jags, we are extracting a.scale
  # 1) Convert to mcmc "matrix", extract "vars" column
  # 2) Unlist to combine them all
  # 3) Convert these into a data.frame
  
  index <- sample(1:nrow(mcdat),size=N,replace=FALSE)
  # Random Sample of 1000 data sets to do integration
  mcdatless  <- mcdat[index,]
  
  gdat <-  sapply(Tgrid,function(Tgrid){numRates_7Benp(mcdatless[,"e0_1"], mcdatless[,"ga_1"], mcdatless[,"gb_1"], 5, 5, 2, 0, 0, mcdatless[,"i_1"],
                                                       mcdatless[,"e0_2"], mcdatless[,"ga_2"], mcdatless[,"gb_2"], 5, 5, 3, 1, 1, mcdatless[,"i_2"],
                                                       mcdatless[,"e0_3"], mcdatless[,"ga_3"], mcdatless[,"gb_3"], 5, 5, 3, 1, 1, mcdatless[,"i_3"],
                                                       mcdatless[,"e0_4"], mcdatless[,"ga_4"], mcdatless[,"gb_4"], 5, 5, 1, 0, 0, mcdatless[,"i_4"],
                                                       mcdatless[,"e0_5"], mcdatless[,"ga_5"], mcdatless[,"gb_5"], 5, 5, 4, 3, 3, mcdatless[,"i_5"],
                                                       mcdatless[,"e0_6"], mcdatless[,"ga_6"], mcdatless[,"gb_6"], 5, 5, 2, 1, 1, mcdatless[,"i_6"],
                                                       mcdatless[,"e0_7"], mcdatless[,"ga_7"], mcdatless[,"gb_7"], 5, 5, 0, 1, 1, mcdatless[,"i_7"],
                                                       mcdatless[,"k"],T9 = Tgrid)})
  # Multiple Integration
  # Different temperature at different columns 
  # Recall that numRates_dpg does integration and return a list of the different reaction rates
  # for EACH chain
  
  gg <-  as.data.frame(gdat)
  
  gg2 <- apply(gg, 2, quantile, probs=c(0.16, 0.5, 0.84), na.rm=TRUE)
  
  gg2data <- data.frame(T9 = Tgrid, median = gg2["50%",],lower = gg2["16%",], 
                        upper = gg2["84%",])
  return(gg2data)
}












