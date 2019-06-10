theory <- read.table("Marcucci2005.dat", header = FALSE)
th <- approxfun(theory[,1], theory[,2])
# Piecewise Linear Interpolation? Can we do better?
sfactordpg <- function(E,a.scale){a.scale*th(E)*1e-6}

## APPLY AFTER YOU RUN THE FUNCTION!!
# Okay now I'm not sure what this function does
# vnumRates_dpg <- Vectorize(numRates_dpg)


numRates_dpg <- Vectorize(function(a.scale,T9){

  # Constants
  M0 = 1.007276466812; M1 = 2.01355318262;		# masses (amu) of p and d
  Z0 = 1; Z1 = 1 ;			# charges of t and d
  
  #   DEFINITIONS
  mue <- (M0*M1)/(M0 + M1) # Reduced Center of Mass
  dpieta <- function(E){0.98951013*Z0*Z1*sqrt(mue/E)} # 2 pi eta
  
  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------
  
  integrand <- function(E,T9) {exp(-dpieta(E))*sfactordpg(E,a.scale)*exp(-E/(0.086173324*T9))}
  # Integrand, this means that we are looking at Gamow*S-factor*Boltzmann Factor
  
  # CALCULATE Nuclear rate
  
  Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 0, upper = Inf,
                                                                      abs.tol = 0L,T9 = Temp)$value}
  
  # Note to self, the limits of integration, in some sense, the scale should be appropriate.
  # From HELP, the first argument MUST BE integrated. The optional argument T9 is used to be substituted in
  # Nasv <-> N A <sigma v >
  
  out <- Nasv(T9)
  return(Nasv=out)
}
)










