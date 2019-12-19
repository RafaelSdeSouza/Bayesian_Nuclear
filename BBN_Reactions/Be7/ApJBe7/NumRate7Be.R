## Integration Script
require(xtable)
require(gsl)
require(RcppGSL)
## Auxiliar Function



fu <- function(x){exp(sqrt(log(1+var(x)/mean(x)^2)))}

##



## Load your mcmcChain here
mat <- read.csv("MCMC_ApJ_ultimaterun.csv",header = T) 
#mat$hbg <- 0

## Functions to be used to calculate rate at a given temperature

sigma7Benp  <-  function(ecm, e0, ga, gb, ra, rb, jr, la, lb){
  # input masses, charges, angular momenta
  m1_i = 7.01473482886
  m2_i = 1.00866491582  # masses (amu) of 7Be and n
  m1_f = 7.01435791572
  m2_f = 1.00727646658  # masses (amu) of 7Li and p
  z1_i = 4
  z2_i = 0              # charges of 7Be and n
  z1_f = 3
  z2_f = 1              # charges of 7Li and p
  jt = 1.5              # spins of target, projectile
  jp = 0.5
  Q = 1.644425          # reaction Q-value (MeV) [from nuclear masses]

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

  tapp <- (s_i-b_i)*ga + (s_f-b_f)*gb

  s1=pek*omega*Ga*Gb
  s2=((e0-ecm-tapp)^2)+0.25*((Ga+Gb)^2)
  SG <- (s1/s2)*(1/ecm)

  return(SG = SG)
}


sigma7Benp7mod <- function(ecm,M){

  SF1 <-  sigma7Benp(ecm, e0 = M[["e0_1"]], ga = M[["ga_1"]], gb = M[["gb_1"]], ra = M[["ra"]], rb = M[["rb"]], jr = 2, la = 0, lb = 0)
  SF2 <-  sigma7Benp(ecm, e0 = M[["e0_2"]], ga = M[["ga_2"]], gb = M[["gb_2"]], ra = M[["ra"]], rb = M[["rb"]], jr = 3, la = 1, lb = 1)
  SF3 <-  sigma7Benp(ecm, e0 = M[["e0_3"]], ga = M[["ga_3"]], gb = M[["gb_3"]], ra = M[["ra"]], rb = M[["rb"]], jr = 3, la = 1, lb = 1)
  SF4 <-  sigma7Benp(ecm, e0 = M[["e0_4"]], ga = M[["ga_4"]], gb = M[["gb_4"]], ra = M[["ra"]], rb = M[["rb"]], jr = 1, la = 0, lb = 0)
  SF5 <-  sigma7Benp(ecm, e0 = M[["e0_5"]], ga = M[["ga_5"]], gb = M[["gb_5"]], ra = M[["ra"]], rb = M[["rb"]], jr = 4, la = 3, lb = 3)
  SF6 <-  sigma7Benp(ecm, e0 = M[["e0_6"]], ga = M[["ga_6"]], gb = M[["gb_6"]], ra = M[["ra"]], rb = M[["rb"]], jr = 2, la = 1, lb = 1)
  SF7 <-  sigma7Benp(ecm, e0 = M[["e0_7"]], ga = M[["ga_7"]], gb = M[["gb_7"]], ra = M[["ra"]], rb = M[["rb"]], jr = 0, la = 1, lb = 1)
  SF <- SF1 + SF2 + SF3 + SF4 + SF5 + SF6 + SF7
  return(SF = SF)
}

NumRate7Benp   <- function(x, T9){


  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------

  integrand <- function(E,T9) {
    E * (sigma7Benp7mod(E,x) + x["hbg"]/sqrt(E)) * exp(-11.605*E/T9)
#    (sqrt(E) * sigma7Benp7mod(E,x) + x["hbg"])*
#      exp(-E/(0.086173324*T9))
    }

  # CALCULATE Nuclear rate

  m1 = 7.01473482886
  m2 = 1.00866491582   # masses (amu) of 7Be and n
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



NumRate7BenpTable <- function(mat, N = 1000,T9){
  
  index <- sample(1:nrow(mat),size=N,replace=FALSE)
  Mat  <- mat[index,]
  

  Mat2 <-  as.data.frame(sapply(T9,function(Tgrid){apply(Mat,1,NumRate7Benp,T9=Tgrid)}))


  gg2 <- apply(Mat2, 2, quantile, probs=c(0.16, 0.5, 0.84), na.rm=TRUE)
  fu_I <- apply(Mat2, 2, fu)

  gg2data <- data.frame(T9 =T9, lower = gg2["16%",], median = gg2["50%",], upper = gg2["84%",],
                        fu = fu_I)
#  gg2data$fu <- fu_I
  rownames(gg2data) <- c()
  return(gg2data)
}

# Chain reduction to 2000 iterations (randomly sampled) if needed
Tgrid = c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,
          0.013,0.014,0.015,0.016,0.018,0.020,0.025,0.030,0.040,0.050,0.060,0.070,
          0.080,0.090,0.100,0.110,0.120,0.130,0.140,0.150,0.160,0.180,0.200,0.250,0.300,
          0.350,0.400,0.450,0.500,0.600,0.700,0.800,0.900,1.000,1.250,1.5,1.75,2,
         2.5,3,3.5,4,5,6,7,8,9,10)


NRate <- NumRate7BenpTable(mat,N = 1e3,T9=Tgrid)

df1 = data.frame(NRate[1:30,c(1,3,2,4,5)])
df2 = data.frame(NRate[31:60,c(1,3,2,4,5)])
dfFull = data.frame(df1,df2)

## Latex Table Output


write.csv(NRate,"MCMCrates_Be7.csv",row.names = F)
print(xtable(dfFull[,c(1,2,5,6,7,10)],
             type = "latex",display=c("e","g","E","g","g","E","g"),
             digits=4),include.rownames = FALSE)



rata_tab <- read.csv("MCMCrates_Be7.csv",header = T)

write.matrix(rata_tab[c("T9","median","fu")],"rateBe7.dat")




