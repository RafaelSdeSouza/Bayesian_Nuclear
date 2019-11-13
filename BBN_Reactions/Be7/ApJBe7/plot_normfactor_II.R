require(plyr)
require(ggridges)
plot_normfactors <- function(Normfit){

  prior <- ggs(as.mcmc(data.frame(
    ga1 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
    ga2 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
    ga3 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
    ga4 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
    ga5 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
    ga6 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
    ga7 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5)
  ))) %>% mutate(type="Prior")  

#  Sa <- ggs(as.mcmc(Normfit[,c("y.norm.1.","y.norm.2.",
#                               "y.norm.3.","y.norm.4.",
#                               "y.norm.5.","y.norm.6.",
#                               "y.norm.7.",
#                               "y.norm.8.",
#                               "y.norm.9.","y.norm.10.")]))
  
 Sa <- ggs(as.mcmc(Normfit[,c("y.norm.5.","y.norm.6.",
                                 "y.norm.7.",
                                 "y.norm.8.",
                                 "y.norm.9.","y.norm.10.")]))
  
  
 Sa$Parameter <- revalue(Sa$Parameter,
                          c("y.norm.5." = "Koe88",
                            "y.norm.6." = "Dam18",
                            "y.norm.7." = "Gib59",
                            "y.norm.8." = "Her19",
                            "y.norm.9." = "Cer89",
                            "y.norm.10." = "Tom19"))

 gg <- ggplot(Sa, aes(x = value, y = Parameter)) +
   geom_vline(xintercept = 1,linetype="dashed",color="gray30") +
  stat_halfeyeh(slab_type = "pdf",adjust=3,normalize="panels",fill="#e41a1c",alpha=0.9) + 
  
#  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, 
#                     quantile_lines = FALSE,bandwidth = 0.01,
 #                    quantiles = c(0.025,0.975)) +
  
#  geom_density_ridges(stat = "binline", binwidth = 0.005,
#                      scale = 1.5,panel_scaling=F, color="black",
#                      alpha = 1) +
 
  theme_economist_white() +

#  scale_fill_manual(name = "Probability", values = c("#d95f02","#e41a1c","#377eb8","#984ea3",
#                                                     "#4daf4a",
#                                                     "#e41a1c","#e41a1c","#e41a1c",
#                                                     "#e41a1c","#e41a1c")) +
  scale_alpha_manual(values=c(0,1,0)) +
  
  #  scale_fill_manual(values=c(rep("gray75",7))) +
  # geom_point(size=1,color="red") +
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="black", size=12),
        axis.text  = element_text(size=10),
        strip.background = element_rect("white")) +
  ylab("") +
  #  xlab("Highest Probability Interval")
  xlab(expression(paste("Normalization factors ")))


return(gg)
}
