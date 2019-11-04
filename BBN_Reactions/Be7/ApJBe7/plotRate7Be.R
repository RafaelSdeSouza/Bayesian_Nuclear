require(purrr)
require(forcats)
require(dplyr)
require(tidyverse)
require(ggplot2)
require(ggthemes)
require(colorspace)
NABe7 <- read.csv("MCMCrates_Be7.csv",header = TRUE)
Norm <- NABe7$median

NAI_new <- NABe7[,c("T9","median","lower","upper")]  %>%
  set_names(c("T9","Adopted","Lower","Upper")) %>%
  mutate(Adopted = Adopted/Norm) %>%
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)  %>%
  mutate(data="presentI")


old <- read.csv("Be7Li.csv",header = TRUE) 
old <- old[,c("T.9","Adopted","Lower","Upper")]  %>%
  set_names(c("T9","Adopted","Lower","Upper"))  %>%
#  filter(T9 <=1) %>%
  mutate(data="previous") %>%
  mutate(Adopted = Adopted/Norm) %>%
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)


joint <- rbind(old,NAI_new)
joint$data <- as.factor(joint$data)
joint$data <- factor(joint$data, levels = c("previous","presentI"))


pdf("Be7_NAratio.pdf", width=7.5, height=5)
ggplot(joint,aes(x=T9,y=Adopted, group=data,fill=data,linetype=data)) +
  geom_ribbon(aes(x=T9,ymin=Lower, ymax=Upper,alpha=data),show.legend=FALSE) +
  scale_alpha_manual(values=c(1,0.375))+
  geom_line(size=0.5) +
  coord_cartesian(ylim=c(0.85,1.05),xlim=c(0.00125,10)) +
  theme_cowplot() + xlab("Temperature (GK)") + ylab("Reaction rate ratio") +
  scale_fill_discrete_qualitative(name="") +
  scale_x_log10(breaks = c(0.001,0.01,0.1,1,10),labels=c(expression(10^-3),expression(10^-2),
                                                         expression(10^-1),"1","10"))  +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
   #  annotation_logticks(base=2.875,sides = "l") +
  scale_linetype_manual(guide=F,values=c("dashed","solid"),name="") +
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
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85"),
        axis.text.y = element_text(size = 16, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 16, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")))
dev.off()
