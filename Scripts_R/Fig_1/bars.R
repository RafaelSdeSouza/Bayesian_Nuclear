# Bar plot
require(reshape2)
require(dplyr)
library(htmltab)

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
  lower = c(0.0002, 0.05, 0.03, 0.26,0.0040,0.04,0,0.28)
)

labels <- c("4He" = "4He", "D" = "D/H", "3He" = "3He/H","7Li" = "7Li/H")


ann_line <-data.frame(xmid="3He",xmin="3He",xmax="3He",y0=2,y2=0,y=1,
                     group=factor("obs",levels=c("pred","obs")))

df1$seg <- c(0,0,0,0,0,0,1,0)
dodge <- position_dodge(width=0.35)

pdf("plot/abunda.pdf",height = 5,width = 7)
 ggplot(df1, aes(x=iso,y=resp,group=group, colour = group)) +
  geom_errorbar(aes(ymin = resp-lower, ymax = resp+upper), width = 0.35,
                position = dodge,size=1.25) +
  geom_point(position = dodge,size=3) +
  facet_wrap(~iso,scales = "free",labeller=labeller(iso = labels,label_parsed)) +
  xlab("") + scale_color_tableau(name="") + theme_rafa() +
  ylab("") +  
  theme(axis.text.x = element_blank(),legend.position = "none")
dev.off()


ukLang <- htmltab(doc = "/Users/Rafael/Desktop/nn.html")
colnames(ukLang)

plot(ukLang[,1],ukLang[,3])

ukLang[,3]