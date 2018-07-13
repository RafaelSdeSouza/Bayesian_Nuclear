plot_normfactors <- function(Normfit){

Sa <- ggs(as.mcmc(Normfit),family="scale")

Sa$Parameter <- revalue(Sa$Parameter,
                        c("scale[1]" = "Ali01a","scale[2]" = "Ali01b","scale[3]" = "Cos00",
                          "scale[4]" = "Gei99","scale[5]" = "Kra87","scale[6]" = "Mol80",
                          "scale[7]" = "Zhi77"))

jagsresults(x = Normfit, params = c("scale"),probs = c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))


require(ggridges)

gg <- ggplot(Sa, aes(x = value, y = Parameter, fill=factor(..quantile..),alpha=factor(..quantile..))) +

  #   geom_density_ridges(scale = 2.5,panel_scaling=F) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, alpha=0.3,
                      quantile_lines = TRUE,quantiles = c(0.025, 0.16,  0.84, 0.975)) +
  theme_economist_white() +

  scale_fill_manual(name = "Probability", values = c("#440154BF",  "#21908CBF","#FDE725BF",
                                                     "#21908CBF","#440154BF"))+
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
  xlab("Normalization factors")


return(gg)
}
