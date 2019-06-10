numRates_3He3He <- Vectorize(function(alpha,beta,gamma,i.screen,T9){

  # Constants
  M0 = 3.0160293; M1 = 3.0160293;		# masses (amu) of p and d
  Z0 = 2; Z1 = 2 ;			# charges of t and d
  
  #   DEFINITIONS
  mue <- (M0*M1)/(M0 + M1) # Reduced Center of Mass
  dpieta <- function(E){0.98951013*Z0*Z1*sqrt(mue/E)} # 2 pi eta
  
  sfactor3He3He <- function(E,alpha,beta,gamma,i.screen){
    (exp(0.98951013*Z0*Z1*sqrt(mue)/(2*E^(3/2))*i.screen))*
      (alpha + beta*E + 0.5*gamma*E^2)}
  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------
  
  integrand <- function(E,T9) {exp(-dpieta(E))*sfactor3He3He(E,alpha,beta,gamma,i.screen)*
      exp(-E/(0.086173324*T9))}
  # Integrand, this means that we are looking at Gamow*S-factor*Boltzmann Factor
  
  # CALCULATE Nuclear rate
  
  Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 0, upper = Inf,
                                                                      T9 = Temp)$value}
  
  # Note to self, the limits of integration, in some sense, the scale should be appropriate.
  # From HELP, the first argument MUST BE integrated. The optional argument T9 is used to be substituted in
  # Nasv <-> N A <sigma v >
  
  out <- Nasv(T9)
  return(Nasv=out)
}
)










