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

ggplot(joint,aes(x=T9,y=Adopted, group=data,fill=data)) +
  geom_ribbon(aes(x=T9,ymin=Lower, ymax=Upper),show.legend=FALSE) +
  geom_line() +
  coord_cartesian(xlim=c(5e-3,1),ylim=c(0.95,1.05)) +
  theme_bw() + xlab("Temperature (GK)") + ylab("Reaction") +
  scale_fill_brewer()+
  scale_x_log10()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.925,0.7),
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size=22),
        axis.text  = element_text(size=18),
        axis.ticks = element_line(size = 0.75),
        axis.line = element_line(size = 0.75, linetype = "solid"))


