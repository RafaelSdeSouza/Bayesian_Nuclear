plot_normfactors_DREAM <- function(Normfit){

  Sa <- ggs(as.mcmc(Normfit),family="scale")
  
  Sa$Parameter <- revalue(Sa$Parameter,
                          c("scale1" = "Ali01a","scale2" = "Ali01b","scale3" = "Cos00",
                            "scale4" = "Gei99","scale5" = "Kra87","scale6" = "Mol80",
                            "scale7" = "Zhi77"))


require(ggridges)

gg <- ggplot(Sa, aes(x = value, y = Parameter, fill=factor(..quantile..),alpha=factor(..quantile..))) +

  #   geom_density_ridges(scale = 2.5,panel_scaling=F) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, alpha=0.5, 
                      quantile_lines = TRUE,quantiles = c(0.025, 0.16,  0.84, 0.975)) +
  theme_economist_white() +

  scale_fill_manual(name = "Probability", values = c("#A9A9A9D9", "#6E6E6ED9", "#191919D9",
                                                     "#6E6E6ED9","#A9A9A9D9"))+
  geom_vline(xintercept = 1,linetype="dashed",color="red") +
  #  scale_fill_manual(values=c(rep("gray75",7))) +
  # geom_point(size=1,color="red") +
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="black", size=15),
        axis.text  = element_text(size=10),
        strip.background = element_rect("white")) +
  ylab("") +
  #  xlab("Highest Probability Interval")
  xlab(expression(paste("Normalization factors ")))


return(gg)
}
