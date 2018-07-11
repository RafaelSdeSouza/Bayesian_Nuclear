# import packages
library(magrittr);
library(dplyr);require(plyr)
require(ggplot2)

Norm <- read.table("rates.dat",header = TRUE)$median

new <- read.table("rates.dat",header = TRUE) %>%
  select(c("x","median","lower","upper")) %>%
  set_colnames(c("T9","Adopted","Lower","Upper")) %>%
  mutate(data = "present") %>% 
  mutate(Adopted = Adopted/Norm) %>% 
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)

old <- read.csv("tabula-table_tdn.csv",header = TRUE) %>%  
  select(c("T9","Adopted","Lower","Upper"))  %>%
  mutate(data="previous") %>% 
  mutate(Adopted = Adopted/Norm) %>% 
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)
  
joint <- rbind(old,new)


pdf("rate_ratio.pdf",height = 7,width = 10)
ggplot(joint,aes(x=T9,y=Adopted, group=data,fill=data,linetype=data,alpha=0.3)) +
  geom_ribbon(aes(x=T9,ymin=Lower, ymax=Upper),show.legend=FALSE) +
  geom_line() +
  coord_cartesian(xlim=c(0.00125,1),ylim=c(0.95,1.05)) +
  theme_bw() + xlab("Temperature (GK)") + ylab("Reaction rate ratio") +
  scale_fill_manual(values=c("#606060","#6600CC"))+
  scale_x_log10(breaks = c(0.001,0.01,0.1,1))  +
  annotation_logticks(short = unit(0.2, "cm"), mid = unit(0.3, "cm"), long = unit(0.4, "cm"),
                      sides = "b") +
  annotation_logticks(base=2.875,sides = "l") +
  scale_linetype_manual(guide=F,values=c("dashed","solid")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(size = 3, linetype='dashed',colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size = 22),
        axis.ticks = element_line(size = 0.75),
 #       axis.line = element_line(size = 0.75, linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3, "mm")) 
        
dev.off()

