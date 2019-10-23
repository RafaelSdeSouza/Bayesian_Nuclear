require(truncnorm)
plot_gb <- function(M){     prior <- ggs(as.mcmc(data.frame(
  gb1 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  gb2 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  gb3 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  gb4 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  gb5 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  gb6 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  gb7 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5)
))) %>% mutate(type="Prior")

post <- ggs(as.mcmc(M[,c("gb_1","gb_2",
                         "gb_3","gb_4",
                         "gb_5","gb_6",
                         "gb_7")])) %>% mutate(type="Posterior")

post$Parameter <- revalue(post$Parameter,
                          c("gb_1" = "gb1","gb_2" = "gb2","gb_3" = "gb3",
                            "gb_4" = "gb4","gb_5" = "gb5","gb_6" = "gb6",
                            "gb_7" = "gb7"))

joind <- rbind(prior,post) %>% mutate(type=as.factor(type)) %>% mutate(type=as.factor(type)) %>%
  mutate(Parameter = factor(Parameter, levels = c("gb1","gb2","gb3",
                                                  "gb4","gb5","gb6","gb7"))) %>%
  mutate(Parameter  = factor(Parameter, labels = c("gamma[p[1]]^2","gamma[p[2]]^2","gamma[p[3]]^2",
                                                   "gamma[p[4]]^2","gamma[p[5]]^2",
                                                   "gamma[p[6]]^2","gamma[p[7]]^2")))  


gg <- ggplot(joind, aes(x = value,y=type, fill=Parameter,alpha=type)) +
  # stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, 
  #                      quantile_lines = FALSE,bandwidth = 0.01,
  #                      quantiles = c(0.025,0.975)) +
  
  stat_halfeyeh(slab_type = "pdf",normalize="groups") + 
  # stat_ecdf() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  
  theme_economist_white() +
  
  scale_fill_manual(name = "Probability",
                    values = c("#411BEA","#008000","#BF40BF","#A6761D",
                               "#ff7f00","#f781bf","cyan3")) +
  scale_color_manual(name = "Probability",
                     values = c("#411BEA","#008000","#BF40BF","#A6761D",
                                "#ff7f00","#f781bf","cyan3")) +
  scale_alpha_manual(values=c(0.9,0.7)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="black", size=15),
        axis.text  = element_text(size=11),
        strip.text = element_text(size=15),
        strip.background = element_rect("gray85")) +
  ylab("") +
  #  xlab("Highest Probability Interval")
  xlab(expression(paste("Energy (MeV)"))) +
  facet_wrap(.~Parameter,scale="free_x",nrow=1,labeller = "label_parsed")
return(gg)
}