Sfactor <- function(ECM,ER,gi,gf){
require(gsl)  
# Constants
m1_i = 3.016; m2_i = 2.014;		# masses (amu) of t and d
m1_f = 4.0026; m2_f = 1.0087;	# masses (amu) of n and 4He
z1_i = 1; z2_i = 1;			# charges of t and d
z1_f = 2; z2_f = 0;				#charges of n and 4He
r_i = 6.0; r_f = 5.0;			# channel radii (fm)
jt = 0.5; jp=1.0; jr=1.5;			#spins of target, projectile, resonance
Q = 17.589;						#reaction Q-value (MeV)
la = 0; lb = 2;					#orbital angular momenta of d and n	

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


N <- 30
x1 <- exp(runif(N, log(1e-3), log(1)))
res <- vector()
y1 <- vector()
# Barker values:
# Er  = 0.0912 MeV
# g^2_in = 2.93 MeV         ! reduced width of deuteron
# g^2_out = 0.0794 MeV      ! reduced width of neutron

res[2] <- 0.0912   # resonance energy
res[3] <-  2.93    # reduced width incoming
res[4] <- 0.0794   # reduced width outgoing

for (i in 1:length(x1)){
  res[1] <- x1[i]
  write.table(res, file="tdn_AD.in", quote=TRUE, 
              row.names=FALSE, col.names=FALSE)
  
  # Load the fortran code needed to calculate S-factor curve
  if(!is.loaded("tdn_AD_Sub"))
    dyn.load("tdn_AD.so") 
  .Fortran("tdn_AD_Sub")
  
  tab1 <- read.table("tdn_AD.out", header=FALSE)
  

  y1[i] <- tab1[1,2]
}


# Test
plot(x1,Sfactor(x1,0.0912,2.93,0.0794),log="x",ylim=c(0,35),
     xlim=c(1e-3,1),col="red",cex=0.5)
par(new=TRUE)
plot(x1,y1,col="green",log="x",ylim=c(0,35),xlim=c(1e-3,1),cex=0.5)


hist((Sfactor(x1,0.0912,2.93,0.0794)-y1))
