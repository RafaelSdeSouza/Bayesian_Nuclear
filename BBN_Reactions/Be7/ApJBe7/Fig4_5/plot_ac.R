#prior posterior comparison
require(truncnorm)
plot_ac <- function(M){prior <- ggs(as.mcmc(data.frame(
  an = rtruncnorm(1e5, a = 0, b = Inf, mean = 4, sd = 0.5),
  ap = rtruncnorm(1e5, a = 0, b = Inf, mean = 4, sd = 0.5)
))) %>% mutate(type="Prior")

post <- ggs(as.mcmc(M[,c("ra","rb")])) %>% mutate(type="Posterior")

post$Parameter <- revalue(post$Parameter,
                          c("ra" = "an","rb" = "ap"))

joind <- rbind(prior,post) %>% mutate(type=as.factor(type)) %>%
  mutate(Parameter = factor(Parameter, levels = c("an","ap"))) %>%
  mutate(Parameter  = factor(Parameter, labels = c("a[n]","a[p]")))  



gg <- ggplot(joind, aes(x = value, y=type, fill=Parameter,alpha=type)) +
  # stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, 
  #                      quantile_lines = FALSE,bandwidth = 0.01,
  #                      quantiles = c(0.025,0.975)) +
  
  #  geom_density(aes(y=..scaled..),size=1.25) +
  stat_halfeyeh(slab_type = "pdf",normalize="groups") + 
  # stat_ecdf() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4,min.n = 4)) +
  theme_economist_white() +
  
  scale_fill_manual(name = "Probability",
                    values = c("#e31a1c","#e31a1c")) +
  scale_color_manual(name = "Probability",
                     values = c("#e31a1c","#e31a1c")) +
  scale_alpha_manual(values=c(1,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="black", size=23),
        axis.text  = element_text(size=18),
        strip.text = element_text(size=25),
        strip.background = element_rect("gray85")) +
  ylab("") +
  #  xlab("Highest Probability Interval")
  xlab("Chanel radii (fm)") + 
  facet_wrap(.~Parameter,scale="free_x",nrow=1,labeller = "label_parsed")
return(gg)
}