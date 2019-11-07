
Gamma7Benp <- function(e0, ga, gb, ra, rb, jr, la, lb){
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
  
  ## incoming channel   
  eta_a=0.1574854*z2_i*z1_i*sqrt(mue_i)
  rho_a=0.218735*ra*sqrt(mue_i)      
  eta_i=eta_a/(sqrt(e0))
  rho_i=rho_a*(sqrt(e0))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k=0)
  # penetration and shift factor 
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  # partial width
  Ga <- 2*ga*p_i
  
  ## outgoing channel
  eta_b=0.1574854*z2_f*z1_f*sqrt(mue_f)
  rho_b=0.218735*rb*sqrt(mue_f)      
  eta_f=eta_b/(sqrt(e0+Q))
  rho_f=rho_b*(sqrt(e0+Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k=0)
  # penetration and shift factor
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  # partial width
  Gb <- 2*gb*p_f
  
  return(list(Ga = Ga, Gb = Gb))
}


M <- read.csv("MCMC_ApJ_ultimaterun.csv",header = T) 

probBe7 <- function(x){
  quantile(x,probs=c(0.16, 0.5, 0.84))
}
formtab <- function(tab){
  low <-  tab[2,] - tab[1,]
  up <- tab[3,] - tab[2,]
  mean <- tab[2,]
  out <- data.frame(low,mean,up)
  return(out)
}


G1 <- Gamma7Benp(e0 = M[["e0_1"]], ga = M[["ga_1"]], gb = M[["gb_1"]],
                 ra = M[["ra"]], rb = M[["rb"]], jr = 2, la = 0, lb = 0) %>%
                  as.data.frame()
G2 <- Gamma7Benp(e0 = M[["e0_2"]], ga = M[["ga_2"]], gb = M[["gb_2"]],
                 ra = M[["ra"]], rb = M[["rb"]], jr = 3, la = 1, lb = 1) %>%
  as.data.frame()

G3 <- Gamma7Benp(e0 = M[["e0_3"]], ga = M[["ga_3"]], gb = M[["gb_3"]],
                 ra = M[["ra"]], rb = M[["rb"]], jr = 3, la = 1, lb = 1) %>%
  as.data.frame()

G4 <- Gamma7Benp(e0 = M[["e0_4"]], ga = M[["ga_4"]], gb = M[["gb_4"]],
                 ra = M[["ra"]], rb = M[["rb"]], jr = 1, la = 0, lb = 0) %>%
  as.data.frame()
G5 <- Gamma7Benp(e0 = M[["e0_5"]], ga = M[["ga_5"]], gb = M[["gb_5"]],
                 ra = M[["ra"]], rb = M[["rb"]], jr = 4, la = 3, lb = 3) %>%
  as.data.frame()

G6 <- Gamma7Benp(e0 = M[["e0_6"]], ga = M[["ga_6"]], gb = M[["gb_6"]],
                 ra = M[["ra"]], rb = M[["rb"]], jr = 2, la = 1, lb = 1) %>%
  as.data.frame()

G7 <- Gamma7Benp(e0 = M[["e0_7"]], ga = M[["ga_7"]], gb = M[["gb_7"]],
                 ra = M[["ra"]], rb = M[["rb"]], jr = 0, la = 1, lb = 1) %>%
  as.data.frame()

tabG1 <- apply(G1,2,probBe7)
formtab(tabG1)

tabG2 <- apply(G2,2,probBe7)
formtab(tabG2)

tabG3 <- apply(G3,2,probBe7)
formtab(tabG3)

tabG4 <- apply(G4,2,probBe7)
formtab(tabG4)

tabG5 <- apply(G5,2,probBe7)
formtab(tabG5)

tabG6 <- apply(G6,2,probBe7)
formtab(tabG6)

tabG7 <- apply(G7,2,probBe7)
formtab(tabG7)

Ex1 = M[,c("e0_1","e0_2",
           "e0_3","e0_4",
           "e0_5","e0_6",
           "e0_7")] + rnorm(40000,18898.64,0.08)