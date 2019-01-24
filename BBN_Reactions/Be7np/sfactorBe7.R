# Rafael de Souza, UNC
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
sigmaE2Be7 <- function(ecm,e0,er,gi,gf, ri = 6, rf =5,ue=0){
  # Constants
  m1_i = 3.01493; m2_i = 2.01355;		# masses (amu) of t and d
  m1_f = 4.00151; m2_f = 1.007277;	# masses (amu) of n and 4He
  z1_i = 2; z2_i = 1;			# charges of t and d
  z1_f = 2; z2_f = 1;				#charges of n and 4He
#  rd = 6.0; rp = 5.0;			# channel radii (fm)
  la = 0; lb = 2;					#orbital angular momenta of d and n
  Q = 18.353053;						#reaction Q-value (MeV)
  jt = 0.5; jp=1.0; jr=1.5;			#spins of target, projectile, resonance

  #   DEFINITIONS

  mue_i <- (m1_i*m2_i)/(m1_i+m2_i);
  mue_f <- (m1_f*m2_f)/(m1_f+m2_f);
  pek <- 6.56618216e-1/mue_i;
  omega <- (2*jr+1)/((2*jt+1)*(2*jp+1));

  ### CALCULATE S-FACTOR
  ## incoming channel
  etpe_i=exp(0.98951013*z1_i*z2_i*sqrt(mue_i/ecm))
  eta_a=0.1574854*z2_i*z1_i*sqrt(mue_i)
  rho_a=0.218735*ri*sqrt(mue_i)
  eta_i=eta_a/(sqrt(ecm))
  rho_i=rho_a*(sqrt(ecm))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k=0)
  # penetration and shift factor
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  s_i <- rho_i*(P3$val_F*P3$val_Fp + P3$val_G*P3$val_Gp)/(P3$val_F^2 + P3$val_G^2)
  # shift factor at energy Er
  xeta_i=eta_a/(sqrt(er))
  xrho_i=rho_a*(sqrt(er))
  PX1 <- coulomb_wave_FG(xeta_i, xrho_i, la, k=0)
  b_i <- xrho_i*(PX1$val_F*PX1$val_Fp + PX1$val_G*PX1$val_Gp)/(PX1$val_F^2 + PX1$val_G^2)
  # partial width
  Ga <- 2*gi*p_i

  ## outgoing channel
  eta_b=0.1574854*z2_f*z1_f*sqrt(mue_f)
  rho_b=0.218735*rf*sqrt(mue_f)
  eta_f=eta_b/(sqrt(ecm+Q))
  rho_f=rho_b*(sqrt(ecm+Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k=0)
  # penetration and shift factor
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  s_f <- rho_f*(P4$val_F*P4$val_Fp + P4$val_G*P4$val_Gp)/(P4$val_F^2 + P4$val_G^2)
  # shift factor at energy Er+Q
  xeta_f=eta_b/(sqrt(er+Q))
  xrho_f=rho_b*(sqrt(er+Q))
  PX2 <- coulomb_wave_FG(xeta_f, xrho_f, lb, k=0)
  b_f <- xrho_f*(PX2$val_F*PX2$val_Fp + PX2$val_G*PX2$val_Gp)/(PX2$val_F^2 + PX2$val_G^2)
  # partial width
  Gb <- 2*gf*p_f

  tapp <- (s_i-b_i)*gi+(s_f-b_f)*gf

  s1=pek*etpe_i*omega*Ga*Gb
  s2=((e0-ecm-tapp)^2)+0.25*((Ga+Gb)^2)
  SF <- exp( 0.5*0.98951013e0*z1_i*z2_i*sqrt(mue_i)*(1e-6*ue)*ecm^(-1.5) )*s1/s2

  return(SF = SF)
}
