require(truncnorm)
plot_ga <- function(M){     prior <- ggs(as.mcmc(data.frame(
  ga1 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  ga2 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  ga3 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  ga4 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  ga5 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  ga6 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5),
  ga7 = rtruncnorm(1e4, a = 0, b = Inf, mean = 0, sd = 2.5)
))) %>% mutate(type="Prior")

post <- ggs(as.mcmc(M[,c("ga_1","ga_2",
                         "ga_3","ga_4",
                         "ga_5","ga_6",
                         "ga_7")])) %>% mutate(type="Posterior")

post$Parameter <- revalue(post$Parameter,
                          c("ga_1" = "ga1","ga_2" = "ga2","ga_3" = "ga3",
                            "ga_4" = "ga4","ga_5" = "ga5","ga_6" = "ga6",
                            "ga_7" = "ga7"))

joind <- rbind(prior,post) %>% mutate(type=as.factor(type)) %>% mutate(type=as.factor(type)) %>%
  mutate(Parameter = factor(Parameter, levels = c("ga1","ga2","ga3",
                                                  "ga4","ga5","ga6","ga7"))) %>%
  mutate(Parameter  = factor(Parameter, labels = c("gamma[n[1]]^2","gamma[n[2]]^2","gamma[n[3]]^2",
                                                   "gamma[n[4]]^2","gamma[n[5]]^2",
                                                   "gamma[n[6]]^2","gamma[n[7]]^2")))  


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