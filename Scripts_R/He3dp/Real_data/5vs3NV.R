
NV <- read.csv("NV.csv")
NV5 <- read.csv("5p_NV.csv")
NV$case <- "3p"
NV5$case <- "5p"

NV_all <- rbind(NV,NV5)
NV_all$case <- as.factor(NV_all$case)


cmb3<- read.csv("df.csv")
cmb5 <- read.csv("df_5.csv")
cmb5W <- read.csv("df_5W.csv")
cmb3$case <- "3p"
cmb5$case <- "5p"
cmb5W$case <- "5pW"

cmb_all <- rbind(cmb3,cmb5,cmb5W)

ggplot(cmb_all, aes(x,y,group=case,alpha=0.75)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=case)) +
  scale_fill_tableau()+
  geom_line() + theme_wsj() + xlab(expression(N[A]~sigma*v)) + ylab("Density") +
  theme(legend.position = "top",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He")))


cmb <- as.data.frame(filter(NV_all, x <= 1.1 & x >= 0.99))



g1 <- ggplot(NV_all,aes(x=x,y=mean,group=case,fill=case,alpha=0.8,color=case))+
  theme_bw()  +
  
#  geom_ribbon(data=NV_all,aes(x=x,ymin=lwr3, ymax=upr3,y=NULL), show.legend=FALSE) +
 # geom_ribbon(data=NV_all,aes(x=x,ymin=lwr2, ymax=upr2,y=NULL), show.legend=FALSE) +
 geom_ribbon(data=NV_all,aes(x=x,ymin=lwr1, ymax=upr1,y=NULL),show.legend=FALSE) +
  geom_line(size=1,colour="#ffffff",linetype="dashed",size=1,show.legend=FALSE) +
  theme_wsj() + xlab("Temperature (GK)") + ylab(TeX('$N_A\\sigma v$')) +
  scale_colour_npg() +
  scale_fill_npg() +
  #  scale_y_log10() + scale_x_log10()+
  theme(legend.position = "none",
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He"))) +
coord_cartesian(xlim=c(0,1.5))

box <- data.frame(x=NV$x,y=(NV$mean/NV5$mean))

bb <- exp(seq(log(1e-3),log(10),length.out =  9))
box$x <- cut(box$x,breaks=bb)

ggplot(box[-1,],aes(x=x,y=y)) +
  geom_boxplot(fill="blue") +xlab("T (GK)") +
  ylab("Ratio") + theme_bw()

hist((NV$mean-NV5$mean)/NV$mean)
g1