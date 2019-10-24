require(ggmcmc)
require(coda)
require(ggridges)
require(nuclear)
require(dplyr)
require(ggplot2)
require(ggthemes)
require(forcats)
require(tidybayes)
require(truncnorm)
require(nuclear)
require(gsl)

source("pair_wise_plot.R")
source("plot_normfactor_II.R")
source("sigma7Benp7mod.R")
source('plot_Er.R')
source('plot_gb.R')
source('plot_ga.R')
Be7np <- read.csv("Be7np.csv")


Be7np$dat <- as.factor(Be7np$dat)
Be7np$dat <- fct_relevel(Be7np$dat, "Dam18","Gib59","Mar19","Koe88","Koe88b","Dam18b","Gib59b",
                         "Her19","Cer89","Tom19")  

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


# Plotting routines
samp <- read.csv("MCMC_ApJ_ultimaterun.csv",header = T) 
en <- samp[,c('e0_1','e0_2',
              'e0_3', 'e0_4',
              'e0_5', 'e0_6', 
              'e0_7')]

pdf("Be7_norm.pdf", width=4, height=5.5)
plot_normfactors(samp)
dev.off()

pdf("Be7_Er.pdf", width=18.25, height=3.5)
plot_Er(samp)
dev.off()

pdf("Be7_ga.pdf", width=18.25, height=3.5)
plot_ga(samp)
dev.off()

pdf("Be7_gb.pdf", width=18.25, height=3.5)
plot_gb(samp)
dev.off()

SE <- ggs(as.mcmc(en))
ggs_traceplot(SE)

gan <- samp[,c('ga_1','ga_2', 'ga_3', 'ga_4', 'ga_5', 'ga_6', 'ga_7')]
Sa <- ggs(as.mcmc(gan))
ggs_traceplot(Sa)

gbn <- samp[,c('gb_1','gb_2', 'gb_3', 'gb_4', 'gb_5', 'gb_6', 'gb_7')]
Sb <- ggs(as.mcmc(gbn))
ggs_traceplot(Sb)



xx <- exp(seq(log(1e-9),log(3),length.out = 400))
# resonance 1
y1 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_1"],samp [,"ga_1"],samp [,"gb_1"],
                                          samp [,"ra"],
                                          samp [,"rb"],jr = 2, la = 0, lb = 0),probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y1 <- rbind(y1,mux)
}
gr1 <-  data.frame(x =xx, mean = y1[,"50%"],lwr1=y1[,"25%"],lwr2=y1[,"2.5%"],
                   upr1=y1[,"75%"],
                   upr2=y1[,"97.5%"])
# resonance 2
y2 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_2"],samp[,"ga_2"],samp [,"gb_2"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 3, la = 1, lb = 1),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y2 <- rbind(y2,mux)
}

gr2 <-  data.frame(x =xx, mean = y2[,"50%"],lwr1=y2[,"25%"],lwr2=y2[,"2.5%"],
                   upr1=y2[,"75%"],
                   upr2=y2[,"97.5%"])

# resonance 3
y3 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_3"],samp[,"ga_3"],samp [,"gb_3"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 3, la = 1, lb = 1),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y3 <- rbind(y3,mux)
}

gr3 <-  data.frame(x =xx, mean = y3[,"50%"],lwr1=y3[,"25%"],lwr2=y3[,"2.5%"],
                   upr1=y3[,"75%"],
                   upr2=y3[,"97.5%"])


# resonance 4
y4 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_4"],samp[,"ga_4"],samp [,"gb_4"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 1, la = 0, lb = 0),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y4 <- rbind(y4,mux)
}

gr4 <-  data.frame(x =xx, mean = y4[,"50%"],lwr1=y4[,"25%"],lwr2=y4[,"2.5%"],
                   upr1=y4[,"75%"],
                   upr2=y4[,"97.5%"])

# resonance 5
y5 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_5"],samp[,"ga_5"],samp [,"gb_5"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 4, la = 3, lb = 3),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y5 <- rbind(y5,mux)
}

gr5 <-  data.frame(x =xx, mean = y5[,"50%"],lwr1=y5[,"25%"],lwr2=y5[,"2.5%"],
                   upr1=y5[,"75%"],
                   upr2=y5[,"97.5%"])


# resonance 6
y6 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_6"],samp[,"ga_6"],samp [,"gb_6"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 2, la = 1, lb = 1),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y6  <- rbind(y6,mux)
}

gr6 <-  data.frame(x =xx, mean = y6[,"50%"],lwr1=y6[,"25%"],lwr2=y6[,"2.5%"],
                   upr1=y6[,"75%"],
                   upr2=y6[,"97.5%"])


# resonance 7
y7 <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp(xx[j],samp[,"e0_7"],samp[,"ga_7"],samp [,"gb_7"],
                                          samp[,"ra"],
                                          samp[,"rb"],jr = 0, la = 1, lb = 1),
                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  y7  <- rbind(y7,mux)
}

gr7 <-  data.frame(x =xx, mean = y7[,"50%"],lwr1=y7[,"25%"],lwr2=y7[,"2.5%"],
                   upr1=y7[,"75%"],
                   upr2=y7[,"97.5%"])



yall <- NULL
for(j in 1:length(xx)){
  mux <-  quantile(sqrt(xx[j])*sigma7Benp7mod(xx[j],
                                              samp[,"e0_1"],samp[,"ga_1"],samp[,"gb_1"],samp[,"ra"],samp[,"rb"],
                                              samp[,"e0_2"],samp[,"ga_2"],samp[,"gb_2"],samp[,"ra"],samp[,"rb"],
                                              samp[,"e0_3"],samp[,"ga_3"],samp[,"gb_3"],samp[,"ra"],samp[,"rb"],
                                              samp[,"e0_4"],samp[,"ga_4"],samp[,"gb_4"],samp[,"ra"],samp[,"rb"],
                                              samp[,"e0_5"],samp[,"ga_5"],samp[,"gb_5"],samp[,"ra"],samp[,"rb"],
                                              samp[,"e0_6"],samp[,"ga_6"],samp[,"gb_6"],samp[,"ra"],samp[,"rb"],
                                              samp[,"e0_7"],samp[,"ga_7"],samp[,"gb_7"],samp[,"ra"],samp[,"rb"]),
                                              probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
     yall <- rbind(yall,mux)
     }

grall <-  data.frame(x = xx, mean = yall[,"50%"],lwr1=yall[,"25%"],lwr2=yall[,"2.5%"],
                     upr1=yall[,"75%"],
                     upr2=yall[,"97.5%"])

MB <- function(x1){3e2*(x1*(2.718^(-x1/(0.086173*0.5))))}
MBD <- data.frame(x=xx,y=MB(xx))

Be7npG <-  Be7np 
library(plyr)
Be7npG$dat <- revalue(Be7npG$dat, c("Koe88b"="Koe88","Dam18b"="Dam18",Gib59b="Gib59"))

# Ressonances

dr=data.frame(x=c(0,0.15,0.34,0.51,0.96,1.23,1.32), 
              y=rep(10,7), 
              vy=rep(9.5,7))

# Plot all
pdf("Be7_slice.pdf", width=7.5, height=5)
ggplot(Be7npG,aes(x=E,y=S)) +
  
  
#  geom_area(data=MBD,aes(x=x,y=y),color="#a6cee3",fill="#a6cee3",
#            size=0,alpha=0.4) +
  
#
#  geom_ribbon(data=gr1,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#6baed6"),show.legend=FALSE) +
  geom_ribbon(data=gr1,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#411BEA"),show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 5.75+2, yend = 5.75+2,size=1.5,
           colour = "#411BEA") +
  annotate(geom="text", x = 1.25,  y = 5.75+2,
           label=expression(2^"-"*", 0.00"),size=4) +
  
#  geom_ribbon(data=gr2,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#74c476"),show.legend=FALSE) +
  geom_ribbon(data=gr2,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#008000"),show.legend=FALSE) +
  annotate("segment",  x = 0.2, xend = 0.35, y = 5.25+2, yend = 5.25+2,size=1.5,
           colour = "#008000") +
  annotate(geom="text", x = 1.25,  y = 5.25+2,
           label=expression(3^"+"*", 0.15"),size=4) +
  
  
#  geom_ribbon(data=gr3,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#BF40BF"),alpha=0.7,show.legend=FALSE) +
  geom_ribbon(data=gr3,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#BF40BF"),alpha=0.95,show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 4.75+2, yend = 4.75+2,size=1.5,
           colour = "#BF40BF") +
  annotate(geom="text", x = 1.25,  y = 4.75+2,
           label=expression(3^"+"*", 0.34"),size=4) +
  
  
#  geom_ribbon(data=gr4,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#E6AB02"),show.legend=FALSE) +
  geom_ribbon(data=gr4,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#A6761D"),show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 4.25+2, yend = 4.25+2,size=1.5,
           colour = "#A6761D") +
  annotate(geom="text", x = 1.25,  y = 4.25+2,
           label=expression(1^"-"*", 0.51"),size=4) +
  
  
  
#  geom_ribbon(data=gr6,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#fccde5"),show.legend=FALSE) +
  geom_ribbon(data=gr6,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#f781bf"),show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 3.25+2, yend = 3.25+2,size=1.5,
           colour = "#f781bf") +
  annotate(geom="text",  x = 1.25,  y = 3.25+2,
           label=expression(2^"+"*", 1.23"),size=4) +
  
  

  
  
#  geom_ribbon(data=gr7,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("cyan"),show.legend=FALSE) +
  geom_ribbon(data=gr7,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("cyan3"),show.legend=FALSE) +
  annotate("segment",x = 0.2, xend = 0.35, y = 2.75+2, yend = 2.75+2,size=1.5,
           colour = "cyan3") +
  annotate(geom="text",  x = 1.25,  y = 2.75+2,
           label=expression(0^"+"*", 1.32"),size=4) +
  
#  geom_ribbon(data=gr5,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#fdb462"),show.legend=FALSE) +
  geom_ribbon(data=gr5,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#ff7f00"),show.legend=FALSE) +
  annotate("segment", x = 0.2, xend = 0.35, y = 3.75+2, yend = 3.75+2,size=1.5,
           colour = "#ff7f00") +
  annotate(geom="text", x = 1.25,  y = 3.75+2,
           label=expression(4^"+"*", 0.96"),size=4) +
  
  
 # geom_ribbon(data=grall,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#fb6a4a"),show.legend=FALSE) +
  geom_ribbon(data=grall,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#99000d"),show.legend=FALSE) +
  
  
  
  geom_errorbar(show.legend = FALSE,aes(x=E,ymin=S-Stat,ymax=S+Stat,group=dat,color=type),alpha=0.75,width=0.025)+
  geom_point(data=Be7npG,aes(x=obsx,y=obsy,group=dat,shape=dat,color=type,fill=type,size=type,alpha=type)
             )+
  
  scale_shape_manual(guide = guide_legend(nrow = 2,
                    override.aes = list(size = 2.5)),
                    values=c(22,24,21,4,25,8,23),name="")+
  scale_size_manual(values=c(2.7,2.5),name="",guide="none")+
  scale_color_manual(name="",values=c("black","gray45"),guide="none")+
  scale_alpha_manual(name="",values=c(1,0.75),guide="none")+
  scale_fill_manual(values=c("black","white"),name="",guide="none") +
  scale_x_log10(breaks = c(1e-6,1e-3,1),
                labels=c(expression(10^-6),expression(10^-3),"1")) +
  theme_economist_white() + ylab(expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)"))) + 
  xlab("Energy (MeV)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.45,0.95),
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
  coord_cartesian(xlim=c(1e-8,2.25),ylim=c(0.175,9.75)) +
  annotate(geom="text",1e-7, 1,
           label=expression(paste(NULL^"7","Be(n,p)",NULL^"7","Li")),
           size=6) 
dev.off()



Be7npG2 <- Be7npG %>% unite(comb, c(dat,type),remove=F) %>%
  mutate(comb =as.factor(comb))

# Plot all
pdf("data.pdf", width=7.5, height=5)
ggplot(Be7npG,aes(x=E,y=S)) +
  
  
  geom_area(data=MBD,aes(x=x,y=y),color="#a6cee3",fill="#a6cee3",
            size=0,alpha=0.4) +
  
 
  geom_errorbar(data=Be7npG2,show.legend = FALSE,aes(x=E,y=S,ymin=S-Stat,ymax=S+Stat,group=dat,color=comb),width=0.025)+
  geom_point(data=Be7npG2,aes(x=obsx,y=obsy,group=dat,shape=dat,color=comb,fill=type,size=type))+
  
  scale_shape_manual(values=c(22,24,21,4,25,8,23),name="",
                     guide = guide_legend(ncol = 1,
                             override.aes = list(size = 2.5))) +
  scale_size_manual(values=c(2.75,2.5),name="",guide="none")+
#  scale_color_discrete_qualitative(name="",guide="none") +
  scale_color_manual(name="",guide="none",
  values=c("#e41a1c","black","#d95f02","black","#377eb8","black",
           "#e41a1c","#4daf4a","#984ea3","black")) +
  scale_alpha_manual(name="",values=c(1,0.5),guide="none")+
  scale_fill_manual(values=c("#e41a1c","white"),name="",guide="none") +
  scale_x_log10(breaks = c(1e-6,1e-3,1),
                labels=c(expression(10^-6),expression(10^-3),"1")) +
  theme_economist_white() + ylab(expression(paste(sqrt(E), sigma, " (", sqrt(MeV), "b)"))) + 
  xlab("Energy (MeV)") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.875,0.75),
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
  coord_cartesian(xlim=c(1e-8,2.25),ylim=c(0.175,9.75)) +
  annotate(geom="text",1e-7, 1,
           label=expression(paste(NULL^"7","Be(n,p)",NULL^"7","Li")),
           size=6) 
dev.off()




# Corner plot
Corr_chain <- ggs(as.mcmc(mat[,c("e0_1","e0_2","e0_3","e0_4","e0_5","e0_6","e0_7")])) %>%
  mutate(Parameter = factor(Parameter, levels = c("e0_1","e0_2","e0_3","e0_4","e0_5","e0_6","e0_7"))) %>%
  mutate(Parameter  = factor(Parameter, labels = c("E[1]~(MeV)","E[2]~(MeV)","E[3]~(MeV)","E[4]~(MeV)",
                                                   "E[5]~(MeV)","E[6]~(MeV)","E[7]~(MeV)")))

pdf("Be7_corr.pdf",height = 9,width = 9)
pair_wise_plot(Corr_chain)
dev.off()

ggs_pairs(Corr_chain )

ggs_pairs(Corr_chain, 
          labeller = "label_parsed",
          diag=list(continuous = my_hist),
          upper = "blank",
          lower = list(continuous = my_bin),
          switch="both",
          showStrips=TRUE
) 