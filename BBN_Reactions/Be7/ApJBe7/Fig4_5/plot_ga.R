
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