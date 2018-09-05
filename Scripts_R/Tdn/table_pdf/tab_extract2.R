# import packages
library(magrittr);
library(dplyr);require(plyr)
require(ggplot2)
D <- read.csv("tabula93_Bosch.csv",header = F)

ff <- c(1e-26,1e-25,1e-24,1e-23,rep(1e-22,2),rep(1e-21,2),rep(1e-20,3),
        rep(1e-19,4),rep(1e-18,2),rep(1e-17,3),rep(1e-16,7))
D$V2 <- D$V2*ff


colnames(D) <- c("T9","Adopted")
D$Lower <- 0.975*D$Adopted
D$Upper <- 1.025*D$Adopted

# Temporary

D <- D[-c(9,12),]


Norm <- read.table("numRates2.dat",header = TRUE)$median

new <- read.table("numRates2.dat",header = TRUE) %>%
  select(c("x","median","lower","upper")) %>%
  set_colnames(c("T9","Adopted","Lower","Upper")) %>%
  mutate(data = "present") %>% 
  mutate(Adopted = Adopted/Norm) %>% 
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)



old <- D %>%  
  select(c("T9","Adopted","Lower","Upper"))  %>%
  mutate(data="previous") %>% 
  mutate(Adopted = Adopted/Norm) %>% 
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)
  
joint <- rbind(old,new)


pdf("rate_ratioii.pdf",height = 7,width = 10)
ggplot(joint,aes(x=T9,y=Adopted, group=data,fill=data,linetype=data,alpha=0.3)) +
  geom_ribbon(aes(x=T9,ymin=Lower, ymax=Upper),show.legend=F) +
  geom_line() +
  coord_cartesian(xlim=c(0.2,50),ylim=c(0.95,1.08)) +
  theme_bw() + xlab("Thermal energy kT (keV)") + ylab("Reactivity ratio") +
  scale_fill_manual(values=c("#606060","green3"))+
  scale_x_log10(breaks = c(0.001,0.01,0.1,1))  +
  annotation_logticks(short = unit(0.2, "cm"), mid = unit(0.3, "cm"), long = unit(0.4, "cm"),
                      sides = "b",size = 1) +
  scale_linetype_manual(guide=F,values=c("dashed","solid")) +
  theme(panel.border = element_rect(size=1.5),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(size = 3, linetype='dashed',colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size = 22),
        axis.ticks = element_line(size = 1),
 #       axis.line = element_line(size = 0.75, linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3, "mm")) 
        
dev.off()

