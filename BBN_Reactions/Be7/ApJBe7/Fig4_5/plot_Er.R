#prior posterior comparison
require(truncnorm)
plot_Er <- function(M){     prior <- ggs(as.mcmc(data.frame(
                            E1 = rtruncnorm(1e5, a = 0, b = Inf, mean = 0, sd = 0.1),
                            E2 = rtruncnorm(1e5, a = 0, b = Inf, mean = 0.15, sd = 0.025),
                            E3 = rtruncnorm(1e5, a = 0, b = Inf, mean = 0.336, sd = 0.010),
                            E4 = rtruncnorm(1e5, a = 0, b = Inf, mean = 0.51, sd = 0.1),
                            E5 = rtruncnorm(1e5, a = 0, b = Inf, mean = 0.96, sd = 0.1),
                            E6 = rtruncnorm(1e5, a = 0, b = Inf, mean = 1.23, sd = 0.1),
                            E7 = rtruncnorm(1e5, a = 0, b = Inf, mean = 1.32, sd = 0.1)
                            ))) %>% mutate(type="Prior")

post <- ggs(as.mcmc(M[,c("e0_1","e0_2",
                         "e0_3","e0_4",
                          "e0_5","e0_6",
                          "e0_7")])) %>% mutate(type="Posterior")

post$Parameter <- revalue(post$Parameter,
                        c("e0_1" = "E1","e0_2" = "E2","e0_3" = "E3",
                          "e0_4" = "E4","e0_5" = "E5","e0_6" = "E6",
                          "e0_7" = "E7"))

joind <- rbind(prior,post) %>% mutate(type=as.factor(type)) %>%
mutate(Parameter = factor(Parameter, levels = c("E1","E2","E3",
                                                "E4","E5","E6","E7"))) %>%
mutate(Parameter  = factor(Parameter, labels = c("E[1]","E[2]","E[3]",
                                                 "E[4]","E[5]",
                                                 "E[6]","E[7]")))  



gg <- ggplot(joind, aes(x = value, y=type, fill=Parameter,alpha=type)) +
  stat_halfeyeh(slab_type = "pdf",normalize="groups") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4,min.n = 4)) +
  scale_fill_manual(name = "Probability",values = colpal) +
  scale_color_manual(name = "Probability",values = colpal) +
  scale_alpha_manual(values=c(1,1)) +
  theme_rmatrix() +
     ylab("") +
  xlab(expression(paste(E["c.m."]," (MeV)"))) + 
  facet_wrap(.~Parameter,scale="free_x",nrow=1,labeller = "label_parsed")
return(gg)
}