require(dplyr)
require(nuclear)
require(gsl)
require(purrr)
require(ggplot2)
require(reshape)
source('theme_GAMMA.R')
source('asinh.R')
colpal <- c("#e31a1c","#008000","#1f78b4","#b15928",
            "#ff7f00","#fdbf6f","#6a3d9a")

M <- read.csv("MCMC_ApJ_ultimaterun.csv",header = T) 



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



probBe7 <- function(x){
  quantile(x,probs=c(0.16, 0.5, 0.84))
}

uplimit <- function(x){
  quantile(x,probs=c(0.95))
}

formtab <- function(tab){
  low <-  tab[2,] - tab[1,]
  up <- tab[3,] - tab[2,]
  mean <- tab[2,]
  out <- data.frame(low,mean,up)
  return(out)
}

formtabs <- function(tab){
  low <-  tab[2] - tab[1]
  up <- tab[3] - tab[2]
  mean <- tab[2]
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

G1s <- (G1$Ga + G1$Gb) %>% as.data.frame()
tabG1s <- apply(G1s,2,probBe7)
formtabs(tabG1s)
apply(G1s,2,uplimit)



tabG2 <- apply(G2,2,probBe7)
formtab(tabG2)

G2s <- (G2$Ga + G2$Gb) %>% as.data.frame()
tabG2s <- apply(G2s,2,probBe7)
formtabs(tabG2s)



tabG3 <- apply(G3,2,probBe7)
formtab(tabG3)
G3s <- (G3$Ga + G3$Gb) %>% as.data.frame()
tabG3s <- apply(G3s,2,probBe7)
formtabs(tabG3s)




tabG4 <- apply(G4,2,probBe7)
formtab(tabG4)
G4s <- (G4$Ga + G4$Gb) %>% as.data.frame()
tabG4s <- apply(G4s,2,probBe7)
formtabs(tabG4s)



tabG5 <- apply(G5,2,probBe7)
formtab(tabG5)
G5s <- (G5$Ga + G5$Gb) %>% as.data.frame()
tabG5s <- apply(G5s,2,probBe7)
formtabs(tabG5s)

tabG6 <- apply(G6,2,probBe7)
formtab(tabG6)
G6s <- (G6$Ga + G6$Gb) %>% as.data.frame()
tabG6s <- apply(G6s,2,probBe7)
formtabs(tabG6s)

tabG7 <- apply(G7,2,probBe7)
formtab(tabG7)
G7s <- (G7$Ga + G7$Gb) %>% as.data.frame()
tabG7s <- apply(G7s,2,probBe7)
formtabs(tabG7s)

<<<<<<< HEAD
GAMMASUM <- data.frame(G1s,G2s, G3s, G4s, G5s, G6s, G7s) %>% 
  set_names(c("G1","G2","G3","G4","G5","G6","G7")) %>%
  melt() %>%
  mutate(variable  = factor(variable, labels = c("Gamma[1]","Gamma[2]","Gamma[3]",
                                                   "Gamma[4]","Gamma[5]",
                                                   "Gamma[6]","Gamma[7]")))  


pdf("GAMMA.pdf", width=9, height=5)
ggplot(GAMMASUM,aes(x= value,color=variable)) +
  geom_vline(xintercept = 1e-1, color="gray45",linetype="dashed",size=0.75) +
  stat_ecdf(geom = "step",size=1.25) +
  scale_x_log10(breaks = c(1e-4,1e-3,1e-2,1e-1,1),
                labels=c(expression(10^-4),
                         expression(10^-3),expression(10^-2),
                         expression(10^-1),"1")) + 
  scale_color_manual(values = colpal) +
  facet_wrap(.~variable,nrow=2,labeller = "label_parsed") +
  scale_y_continuous(breaks = c(1e-4,1e-3,1e-2,1e-1,1),
                labels=c(expression(10^-4),
                         expression(10^-3),expression(10^-2),
                         expression(10^-1),"1"),trans="log") + 
  theme_GAMMA() +
  xlab(expression(paste(E["c.m."]," (MeV)"))) +
  ylab("ECDF")
dev.off()
=======
GAMMASUM <- data.frame(G1s,G2s, G3s, G4s, G5s, G6s, G7s) %>% set_names(c("G1","G2","G3","G4","G5",
                                                                         "G6","G7")) %>%
  melt()


ggplot(GAMMASUM,aes(x= value,fill=variable)) +
  geom_density() +
  facet_wrap(.~variable,scales="free")
>>>>>>> 0c3d9e1f3365f33ca4ff94939b76bb912cc10efa





Ex = M[,c("e0_1","e0_2",
           "e0_3","e0_4",
           "e0_5","e0_6",
           "e0_7")] + rnorm(40000,18898.64,0.08)*1e-3

tabEx <- apply(Ex ,2,probBe7)
formtab(tabEx)






gan <- M[,c('ga_1','ga_2', 'ga_3', 'ga_4', 'ga_5', 'ga_6', 'ga_7')]
gbn <- M[,c('gb_1','gb_2', 'gb_3', 'gb_4', 'gb_5', 'gb_6', 'gb_7')]

apply(gan,2,uplimit)
