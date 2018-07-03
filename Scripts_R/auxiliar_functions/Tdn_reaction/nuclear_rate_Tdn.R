nuclear_rate_Tdn <- function(e0, er, gi, gf, ri, rf, ue, T9){
  # Constants
  M0 = 3.01550; M1 = 2.01355332;		# masses (amu) of t and d
  Z0 = 1; Z1 = 1 ;			# charges of t and d

  #   DEFINITIONS
  mue <- (M0*M1)/(M0 + M1)
  dpieta <- function(E){0.98951013*Z0*Z1*sqrt(mue/E)}

  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------

  integrand <- function(E,T9) {exp(-dpieta(E))*SfacTdn(E, e0, er, gi, gf, ri, rf, ue)*exp(-E/(0.086173324*T9))}

  # CALCULATE Nuclear rate

  Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 1e-5, upper = Inf,
          abs.tol = 0L,T9 = Temp)$value}

  out <- Nasv(T9)
  return(Nasv=out)
}
