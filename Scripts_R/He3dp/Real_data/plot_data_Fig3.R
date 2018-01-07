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
  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="gray90",alpha=0.3) +
  geom_point(data=gobs,aes(x=obsx,y=obsy,group=lab,color=set,shape=set),size=2.85)+
  geom_errorbar(show.legend=FALSE,data=gobs,aes(x=obsx,y=obsy,ymin=obsy-erry,ymax=obsy+erry,group=set,color=set),
                width=0.01)+
  scale_color_manual(name="",values=c("cyan2","red","brown","violet","green",
          "orange2","blue"))+
  scale_shape_manual(values=c(0,1,2,5,25,11,13),name="")+
  coord_cartesian(xlim=c(4.25e-3,0.85),ylim=c(0,20)) +
  theme_bw() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + 
  scale_x_log10()  +
#  annotation_logticks(sides = "b",size=1) +
#  annotation_logticks(base=2.875,sides = "l",size=1) +
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




# Bar plot

#X(4He):               0.2484+-0.0002               0.2449+-0.0040
#D/H(x10^-5):       2.45+-0.05                       2.53+-0.04
#3He/H(x10^-5):   1.07+-0.03                       <1.3
#7Li/H(x10^-10):   5.61+-0.26                      1.58+0.35-0.28
pred <-  c(0.2484,2.45,1.07,5.61)
obs <- c(0.2449,2.53,1.3,1.58)
id <- c(1,2,3,4)

d1 <- data.frame(pred, obs,id)
d1 <- melt(d1,id.vars ="id")

ggplot(d1,aes(x=id,y=value,group=variable,fill=variable)) + 
  geom_bar(stat="identity", position="dodge",aes(group=variable))


df1 <- data.frame(
  iso = factor(c("4He","D","3He","7Li")),
  resp = c(0.2484,2.45,1.07,5.61,0.2449,2.53,1.3,1.58),
  group = factor(c("pred","pred", "pred", "pred","obs","obs","obs","obs")),
  upper = c(0.0002, 0.05, 0.03, 0.26,0.0040,0.04,0,0.35),
  lower = c(0.0002, 0.05, 0.03, 0.26,0.0040,0.04,1.3,0.28)
)


p + geom_linerange(aes(ymin = lower, ymax = upper))
p + geom_pointrange(aes(ymin = lower, ymax = upper))
p + geom_crossbar(aes(ymin = lower, ymax = upper), width = 0.2)

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
