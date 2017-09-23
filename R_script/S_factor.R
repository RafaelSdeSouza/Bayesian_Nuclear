require(gsl)
# Constants
m1_i = 3.016; m2_i = 2.014;		# masses (amu) of t and d
m1_f = 4.0026; m2_f = 1.0087;	# masses (amu) of n and 4He
z1_i = 1; z2_i = 1;			# charges of t and d
z1_f = 0; z2_f = 2;				#charges of n and 4He
r_i = 6.0; r_f = 5.0;			# channel radii (fm)
jt=0.5; jp=1.0; jr=1.5;			#spins of target, projectile, resonance
Q = 17.589;						#reaction Q-value (MeV)
la = 0; lb = 2;					#orbital angular momenta of d and n	
     
#   DEFINITIONS

mue_i <- (m1_i*m2_i)/(m1_i+m2_i);
mue_f <- (m1_f*m2_f)/(m1_f+m2_f);
pek <- 6.56618216e-1/mue_i;
omega <- (2*jr+1)/((2*jt+1)*(2*jp+1));



# Function to calcuate the penetration factor
# INPUT: E, L, Mass0, Mass1, Charge0, Charge1
# OUTPUT: P, S
PenFactor <- function(E, L, R, mue, qQ) {
rho = 0.218735097*R*sqrt(mue*E)  
eta = 0.15748927*qQ*sqrt(mue/E)  

jj <- coulomb_wave_FG(eta, rho, L, k=0)

Fp.val <- jj$val_Fp
Gp.val <- jj$val_Gp

F_l <- jj$val_F*exp(jj$exp_F)
G_l <- jj$val_G*exp(jj$exp_G)
P <- rho/(F_l^2 + G_l^2)
S <- rho*( F_l^2*Fp.val^2 + G_l^2*Gp.val^2)/
  (F_l^2 + G_l^2);  
return(list(P=P,S=S))  
}


Sfactor <- function(E,e1,gi,gf){
Pen1 <- PenFactor(e1, la, r_i, mue_i, z1_i*z2_i)
pr_i <- Pen1$P
sr_i <- Pen1$S
ga <- 2*gi*pr_i

Pen2 <- PenFactor(e1+Q, lb, r_f, mue_f, z1_f*z2_f);
pr_f <- Pen2$P
sr_f <- Pen2$S
gb <- 2.e0*gf*pr_f

#   CALCULATE S-FACTOR FOR EACH ENERGY
etpe_i <- exp( 0.989534267e0*z1_i*z2_i*sqrt(mue_i/E))

Pen3 <- PenFactor(E, la, r_i, mue_i, z1_i*z2_i);
p_i <- Pen3$P
s_i <- Pen3$S
prat_i <- p_i/pr_i
  

Pen4 <- PenFactor(E+Q, lb, r_f, mue_f, z1_f*z2_f);
p_f <- Pen4$P
s_f <- Pen4$S
prat_f <- p_f/pr_f

tapp <- (s_i-sr_i)*gi+(s_f-sr_f)*gf

s1=pek*etpe_i*omega*prat_i*prat_f*ga*gb
s2=( (e1-E-tapp)^2 )+0.25e0*( (ga*prat_i+gb*prat_f)^2 )

return(S = s1/s2)
}
