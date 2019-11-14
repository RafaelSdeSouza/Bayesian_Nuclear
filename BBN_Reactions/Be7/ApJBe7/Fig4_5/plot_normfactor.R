require(plyr)
require(ggridges)
require(ggmcmc)
require(tidybayes)
require(coda)
require(ggthemes)
plot_normfactors <- function(Normfit){

   prior <- ggs(as.mcmc(data.frame(
    Koe88 = rlnorm(1e4, meanlog   = 0, sdlog  = log(1.020)),
    Dam18 = rlnorm(1e4, meanlog   = 0, sdlog  = log(1.10)),
    Gib59 = rlnorm(1e4, meanlog   = 0, sdlog  = log(1.050)),
    Mar19 = rlnorm(1e4, meanlog   = 0, sdlog  = log(1.051)),
    Cer89 = rlnorm(1e4, meanlog   = 0, sdlog  = log(1.085)),
    Tom19 = rlnorm(1e4, meanlog   = 0, sdlog  = log(1.032))
  ))) %>% mutate(type="Prior")  

#  Sa <- ggs(as.mcmc(Normfit[,c("y.norm.1.","y.norm.2.",
#                               "y.norm.3.","y.norm.4.",
#                               "y.norm.5.","y.norm.6.",
#                               "y.norm.7.",
#                               "y.norm.8.",
#                               "y.norm.9.","y.norm.10.")]))
  
   post <- ggs(as.mcmc(Normfit[,c("y.norm.5.","y.norm.6.",
                                 "y.norm.7.",
                                 "y.norm.8.",
                                 "y.norm.9.","y.norm.10.")])) %>%
                                  mutate(type="Posterior")  
  
  
   post$Parameter <- revalue(post$Parameter,
                          c("y.norm.5." = "Koe88",
                            "y.norm.6." = "Dam18",
                            "y.norm.7." = "Gib59",
                            "y.norm.8." = "Mar19",
                            "y.norm.9." = "Cer89",
                            "y.norm.10." = "Tom19"))

 
joind <- rbind(prior,post) %>% mutate(type=as.factor(type))  %>%
  mutate(Parameter=as.factor(Parameter))   
 
 gg <- ggplot(joind, aes(x = value, y=Parameter, fill=type,color=type,alpha = type)) +
   geom_vline(xintercept = 1,linetype="dashed",color="black") +
   stat_halfeyeh(slab_type = "pdf",adjust=3,normalize="panels",alpha=0.9) + 
  
  theme_economist_white() +

  scale_fill_manual(name = "Probability", values = c("#fb6a4a","gray80")) +
  scale_color_manual(name = "Probability", values = c("#67000d","gray40")) +
  scale_alpha_manual(values=c(1,1)) +
  
  #  scale_fill_manual(values=c(rep("gray75",7))) +
  # geom_point(size=1,color="red") +
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="black", size = 12),
        axis.text.x   = element_text(size=11.5),
        axis.text.y   = element_text(size=12),
        strip.background = element_rect("white")) +
  ylab("") +
  #  xlab("Highest Probability Interval")
  xlab(expression(paste("Normalization factor ",f[k])))


return(gg)
}
