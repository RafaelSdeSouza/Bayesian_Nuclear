# Plot data
require(nuclear)
require(dplyr)
require(ggplot2)
require(ggthemes)
require(forcats)
require(tidybayes)
require(gsl)
require(cowplot)
library(plyr)
require(tidyr)
require(emojifont)

Be7np <- read.csv("Be7np.csv")


Be7np$dat <- as.factor(Be7np$dat)
Be7np$dat <- fct_relevel(Be7np$dat, "Dam18","Gib59","Mar19","Koe88","Koe88b","Dam18b","Gib59b",
                         "Mar19b","Cer89","Tom19") 

Be7np$type <- as.factor(Be7np$type)

re <- as.numeric(Be7np$dat) 
Nre <- length(unique(Be7np$dat))
# Unique removes duplicated vector, we want to know how many groups of
# data are there
N <- nrow(Be7np) # Total No of data sets
obsy <- Be7np$S    # Response variable in MeV
obsx <-  Be7np$E   # Predictors
erry <- Be7np$Stat # Error in MeV
set <- Be7np$dat # Get the labels as a vector
fu <- log(c(1.020,1.10,1.050,1.051,1.085,1.032))

xx <- exp(seq(log(1e-9),log(3),length.out = 400))

MB <- function(x1){3e2*(x1*(2.718^(-x1/(0.086173*0.5))))}
MBD <- data.frame(x=xx,y=MB(xx))

Be7npG <-  Be7np 

Be7npG$dat <- revalue(Be7npG$dat, c("Koe88b"="Koe88","Dam18b"="Dam18",Gib59b="Gib59",Mar19b = "Mar19"))



Be7npG2 <- Be7npG %>% unite(comb, c(dat,type),remove=F) %>%
  mutate(comb =as.factor(comb))



#Be7npG2 <- Be7npG2 %>% mutate_cond(type == "abs", E = E + runif(6,-1e-5,1e-5))

absdat <- filter(Be7npG2, type=="abs") %>% mutate(syst = exp(fu)*S-S)


Be7npG2$E[155:160] <- c(1.672947e-08, 2.569022e-08, 1.247648e-02, 1.790060e-03, 0.85*2.416324e-08, 1.5*1.833538e-08)
absdat$E <- Be7npG2$E[155:160] 


# Plot all
pdf("data.pdf", width=7.5, height=5)
ggplot(Be7npG2,aes(x=E,y=S)) +
  
  
  geom_area(data=MBD,aes(x=x,y=y),color="#BbC2C2",fill="#BbC2C2",
            size=0,alpha=0.4) +
  
  
  geom_errorbar(data=filter(Be7npG2,type=="rel"),show.legend = FALSE,aes(x=E,ymin=S-Stat,ymax=S+Stat,group=dat,color=comb),size=0.5,width=0.01)+
  
  geom_errorbar(data=absdat,show.legend = FALSE,aes(x=E,ymin=S-syst,ymax=S+syst),color="#e41a1c",size=0.5,width=0.01)+
  
  geom_point(data=Be7npG2,aes(x=E,y=S,group=dat,shape=dat,color=comb,fill=type,size=type))+
 
  scale_shape_manual(values=c(22,24,21,25,8,23),name="",
                     guide = guide_legend(ncol = 1,
                                          override.aes = list(size = 2.5))) +
  scale_size_manual(values=c(2.85,2.6),name="",guide="none")+
  #  scale_color_discrete_qualitative(name="",guide="none") +
  scale_color_manual(name="",guide="none",
                     values=c("#e41a1c","black","#d95f02","black","#377eb8","black",
                              "#4daf4a","black","#984ea3","black")) +
  scale_alpha_manual(name="",values=c(1,0.5),guide="none")+
  scale_fill_manual(values=c("#e41a1c","#ffffff"),name="",guide="none") +
  scale_x_log10(breaks = c(1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1),
                labels=c(expression(10^-8),expression(10^-7),expression(10^-6),
                         expression(10^-5),expression(10^-4),
                         expression(10^-3),expression(10^-2),
                         expression(10^-1),"1")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  ylab(expression(paste(sqrt(E["c.m."]), sigma["np"], " (", sqrt(MeV), "b)"))) + 
  xlab(expression(paste(E["c.m."]," (MeV)"))) + 
 theme_cowplot() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.85,0.75),
        legend.text = element_text(size=13),
        legend.text.align = 1,
        legend.background = element_rect(colour = "white",fill=NA),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="black", size=17.5),
        axis.text  = element_text(size=16),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  coord_cartesian(xlim=c(2e-8,1.75),ylim=c(0.175,9.75)) +
  annotate(geom="text",1e-7, 1,
           label=expression(paste(NULL^"7","Be(n,p)",NULL^"7","Li")),
           size=6) 
dev.off()



