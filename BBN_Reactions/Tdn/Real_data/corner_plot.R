# Corner plot Tdn
require(dplyr)
require(ggmcmc)
require(msm);require(mcmcplots);
require(ggridges);require(plyr);require(MASS)
source("..//..//auxiliar_functions/pair_wise_plot.R")

ssDat <- read.table("Tdn_DREAM.dat",header = T)

Corr_chain <- ggs(as.mcmc(ssDat[,c("e0","er","gd2","gn2","ad","an","ue")])) %>% 
  as_tibble() %>%
  mutate(Parameter = factor(Parameter, levels = c("e0","er","gd2","gn2","ad","an","ue"))) %>%
  mutate(Parameter  = factor(Parameter, labels = c("E[0]~(MeV)","E[r]~(MeV)","gamma[d]^2~(MeV)", 
                                                   "gamma[n]^2~(MeV)","a[d]~(fm)",
                                                   "a[n]~(fm)","U[e]~(eV)")))

pdf("Tdn_corner.pdf",height = 9,width = 9)
pair_wise_plot(Corr_chain)
dev.off()
