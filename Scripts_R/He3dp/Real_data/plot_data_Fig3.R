# import packages
require(ggplot2);require(ggthemes)
library(magrittr);library(wesanderson)
library(dplyr);require(ggsci);library(latex2exp)

######################################################################
## Read DATA 
ensamble <- read.csv("ensamble.csv",header = T) %>%
  mutate(Syst=replace(Syst,Syst==0.06,0.078))


obsy <- ensamble$S    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat
set <- ensamble$dat
lab <- ensamble$invK



gobs <- data.frame(obsx,obsy,erry,set,lab)
gobs$set <- as.factor(gobs$set)


pdf("plot/He3dp_data.pdf",height = 7,width = 10)
ggplot(gobs,aes(x=obsx,y=obsy))+
  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="gray90",alpha=0.6) +
  geom_point(data=gobs,aes(x=obsx,y=obsy,group=lab,color=set,shape=set),size=3)+
  geom_errorbar(show.legend=FALSE,data=gobs,aes(x=obsx,y=obsy,ymin=obsy-erry,ymax=obsy+erry,group=set,color=set),
                width=0.01)+
  scale_color_manual(name="",values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
          "#a65628","gray30"))+
  scale_shape_manual(values=c(0,19,8,10,4,17,3),name="")+
  coord_cartesian(xlim=c(4.25e-3,0.85),ylim=c(0,20)) +
  theme_bw() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + 
  scale_x_log10()  +
  annotation_logticks(sides = "b") +
  annotation_logticks(base=2.875,sides = "l") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.675),
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size=16,colour = set),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(size=0.75, fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size=18.5),
        axis.text  = element_text(size=13),
        axis.ticks = element_line(size = 0.75),
        axis.line = element_line(size = 0.5, linetype = "solid"))
dev.off()






dodge <- position_dodge(width=0.25)


mtcars$cyl2 <- factor(mtcars$cyl, labels = c("alpha", "beta", "gamma"))

labels <- c("4He" = "4He", "D" = "D/H", "3He" = "3He/H","7Li" = "7Li/H")

 ggplot(df1, aes(x=iso,y=resp,group=group, colour = group)) +
   geom_errorbar(aes(ymin = resp-lower, ymax = resp+upper), width = 0.15,
                  position = dodge) +
   geom_point(position = dodge)+
   
  facet_wrap(~iso,scales = "free",labeller=labeller(iso = labels)) +
   xlab("") + scale_color_tableau(name="") + theme_bw() +
   ylab("") +  coord_flip() +
   theme(axis.text.y = element_blank())
