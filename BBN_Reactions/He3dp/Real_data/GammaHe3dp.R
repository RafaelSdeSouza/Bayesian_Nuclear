require(gsl)
GammaHe3dp <- function(Normfit){
  mcdat_I <- as.data.frame(do.call(rbind, as.mcmc(Normfit)[,c("E0_b","gd2_b", "gp2_b","ad_b","ap_b")]))  
  e0 <- mcdat_I[,1]
  gi <- mcdat_I[,2]
  gf <- mcdat_I[,3]
  ri <- mcdat_I[,4]
  rf <- mcdat_I[,5]
  # input masses, charges, angular momenta
  m1_i = 3.01550 
  m2_i = 2.01355        # masses (amu) of t and d
  m1_f = 4.00151 
  m2_f = 1.008664       # masses (amu) of n and 4He
  z1_i = 2 
  z2_i = 1              # charges of t and d
  z1_f = 2 
  z2_f = 1				# charges of n and 4He
  jt = 0.5                # spins of target, projectile, resonance
  jp = 1.0 
  jr = 1.5                
  Q = 18.353053		# reaction Q-value (MeV)
  la = 0 
  lb = 2				# orbital angular momenta of d and n	
  
  # reduced masses
  mue_i <- (m1_i*m2_i)/(m1_i+m2_i)
  mue_f <- (m1_f*m2_f)/(m1_f+m2_f)
  
  # constants
  omega <- (2*jr+1)/((2*jt+1)*(2*jp+1))
  
  ## incoming channel   
  eta_a=0.1574854*z2_i*z1_i*sqrt(mue_i)
  rho_a=0.218735*ri*sqrt(mue_i)      
  eta_i=eta_a/(sqrt(e0))
  rho_i=rho_a*(sqrt(e0))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k = 0)
  # penetration and shift factor 
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  # partial width
  Ga <- 2*gi*p_i
  
  ## outgoing channel
  eta_b=0.1574854*z2_f*z1_f*sqrt(mue_f)
  rho_b=0.218735*rf*sqrt(mue_f)      
  eta_f=eta_b/(sqrt(e0+Q))
  rho_f=rho_b*(sqrt(e0+Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k = 0)
  # penetration and shift factor
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  # partial width
  Gb <- 2*gf*p_f
  
  return(list(Gd = Ga, Gp = Gb))
}