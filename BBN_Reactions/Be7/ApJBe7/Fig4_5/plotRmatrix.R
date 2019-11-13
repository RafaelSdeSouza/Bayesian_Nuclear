#plot R-matrix
require(truncnorm)
source('plot_Er.R')
source('plot_gb.R')
source('plot_ga.R')
source('plot_normfactor.R')
source('theme_rmatrix.R')
colpal <- c("#e31a1c","#008000","#1f78b4","#b15928",
                               "#ff7f00","#fdbf6f","#6a3d9a")
samp <- read.csv("..//Chains_ApJ/MCMC_ApJ_ultimaterun.csv",header = T) 

theta1 <- plot_Er(samp)
theta2 <- plot_ga(samp)
theta3 <- plot_gb(samp)

pdf("Be7_Rmatrix.pdf", width=21, height=4.25*3)
plot_grid(
  theta1,theta2,theta3 ,
  align = "hv", axis = "tb",
  nrow = 3)  
dev.off()



pdf("Be7_ac.pdf", width=7, height=6)
plot_ac(samp)
dev.off()


source('plot_normfactor.R')
pdf("Be7_norm.pdf", width=3.55, height=4.25)
plot_normfactors(samp)
dev.off()