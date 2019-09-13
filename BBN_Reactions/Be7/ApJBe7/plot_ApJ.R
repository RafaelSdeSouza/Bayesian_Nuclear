require(ggmcmc)
require(mcmcplots)
require(coda)
require(ggridges)
source("pair_wise_plot.R")
source("external_functions/plot_normfactors_DREAM.R")

# Plotting routines
mat <- read.csv("MCMCBe7.csv",header = T)

# Corner plot
Corr_chain <- ggs(as.mcmc(mat[,c("e0_1","e0_2","e0_3","e0_4","e0_5","e0_6","e0_7")])) %>%
  mutate(Parameter = factor(Parameter, levels = c("e0_1","e0_2","e0_3","e0_4","e0_5","e0_6","e0_7"))) %>%
  mutate(Parameter  = factor(Parameter, labels = c("E[1]~(MeV)","E[2]~(MeV)","E[3]~(MeV)","E[4]~(MeV)",
                                                   "E[5]~(MeV)","E[6]~(MeV)","E[7]~(MeV)")))

pdf("Be7_corr.pdf",height = 9,width = 9)
pair_wise_plot(Corr_chain)
dev.off()

ggs_pairs(Corr_chain )

ggs_pairs(Corr_chain, 
          labeller = "label_parsed",
          diag=list(continuous = my_hist),
          upper = "blank",
          lower = list(continuous = my_bin),
          switch="both",
          showStrips=TRUE
) 