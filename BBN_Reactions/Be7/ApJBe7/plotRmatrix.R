#plot R-matrix

samp <- read.csv("MCMC_ApJ_ultimaterun.csv",header = T) 
source('plot_Er.R')
source('plot_gb.R')
source('plot_ga.R')
theta1 <- plot_Er(samp)
theta2 <- plot_ga(samp)
theta3 <- plot_gb(samp)

pdf("Be7_Rmatrix.pdf", width=23.5, height=4.25*3)
plot_grid(
  theta1,theta2,theta3 ,
  align = "hv", axis = "tb",
  nrow = 3)  
dev.off()

pdf("Be7_ac.pdf", width=7, height=6)
plot_ac(samp)
dev.off()



pdf("Be7_norm.pdf", width=3.5, height=4)
plot_normfactors(samp)
dev.off()