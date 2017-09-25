Sfactor_He3dp <- function(ECM,ER,gi,gf){
  require(gsl)  
  # Constants
  m1_i = 3.01493; m2_i = 2.01355;		# masses (amu) of t and d
  m1_f = 4.00151; m2_f = 1.007277;	# masses (amu) of n and 4He
  z1_i = 2; z2_i = 1;			# charges of t and d
  z1_f = 2; z2_f = 1;				#charges of n and 4He
  r_i = 6.0; r_f = 5.0;			# channel radii (fm)
  la = 0; lb = 2;					#orbital angular momenta of d and n	
  Q = 18.353053;						#reaction Q-value (MeV)
  hm = 20.9          # "hm" is Pierre's constant (hbar.c)^2/(2u) i.e. (197.3)^2/(2*931) = 20.9
  jt = 0.5; jp=1.0; jr=1.5;			#spins of target, projectile, resonance

  #   DEFINITIONS
  
  mue_i <- (m1_i*m2_i)/(m1_i+m2_i);
  mue_f <- (m1_f*m2_f)/(m1_f+m2_f);
  pek <- 6.56618216e-1/mue_i;
  omega <- (2*jr+1)/((2*jt+1)*(2*jp+1));
  
  #     ----------------------------------------------------
  #     PENETRABILITY AND SHIFT FUNCTION AT ER
  #     ----------------------------------------------------
  
  
  eta_a=.15748927*z2_i*z1_i*sqrt(mue_i)
  rho_a=.218735097*r_i*sqrt(mue_i)      
  reta_i=eta_a/(sqrt(ER))
  rrho_i=rho_a*(sqrt(ER))
  
  P1 <- coulomb_wave_FG(reta_i, rrho_i, la, k=0)
  pr_i <- rrho_i/(P1$val_F^2 + P1$val_G^2)
  sr_i <- rrho_i*(P1$val_F*P1$val_Fp + P1$val_G*P1$val_Gp)/(P1$val_F^2 + P1$val_G^2)
  ga <- 2*gi*pr_i
  
  
  eta_b=.15748927*z2_f*z1_f*sqrt(mue_f)
  rho_b=.218735097*r_f*sqrt(mue_f)      
  reta_f=eta_b/(sqrt(ER+Q))
  rrho_f=rho_b*(sqrt(ER+Q))
  
  P2 <- coulomb_wave_FG(reta_f, rrho_f, lb, k=0)
  pr_f <- rrho_f/(P2$val_F^2 + P2$val_G^2)
  sr_f <- rrho_f*(P2$val_F*P2$val_Fp + P2$val_G*P2$val_Gp)/(P2$val_F^2 + P2$val_G^2)
  gb <- 2*gf*pr_f
  
  # CALCULATE S-FACTOR   
  
  etpe_i=exp(0.989534267*z1_i*z2_i*sqrt(mue_i/ECM))
  eta_i=eta_a/(sqrt(ECM))
  rho_i=rho_a*(sqrt(ECM))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k=0)
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  s_i <- rho_i*(P3$val_F*P3$val_Fp + P3$val_G*P3$val_Gp)/(P3$val_F^2 + P3$val_G^2)
  prat_i=p_i/pr_i
  
  eta_f=eta_b/(sqrt(ECM+Q))
  rho_f=rho_b*(sqrt(ECM+Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k=0)
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  s_f <- rho_f*(P4$val_F*P4$val_Fp + P4$val_G*P4$val_Gp)/(P4$val_F^2 + P4$val_G^2)
  prat_f=p_f/pr_f
  
  tapp <- (s_i-sr_i)*gi+(s_f-sr_f)*gf
  
  s1=pek*etpe_i*omega*prat_i*prat_f*ga*gb
  s2=((ER-ECM-tapp)^2)+0.25*((ga*prat_i+gb*prat_f)^2)
  SF <- s1/s2
  return(SF = SF)
}