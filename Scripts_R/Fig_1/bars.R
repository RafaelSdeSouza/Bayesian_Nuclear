# Bar plot
require(reshape2)
require(dplyr)
library(htmltab)
require(ggplot2)
require(ggthemes)
source("/Users/Rafael/Documents/GitHub/JAGS_UNC/Scripts_R/auxiliar_functions/theme_rafa.R")

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

#labels <- c("4He" = "4He", "D" = "D/H", "3He" = "3He/H","7Li" = "7Li/H")
labels <- c("4He" = "", "D" = "", "3He" = "","7Li" = "")



df1$seg <- as.factor(c(0,0,0,0,0,0,1,0))
dodge <- position_dodge(width=0.35)

df1$group <- factor(df1$group, levels = c("pred", "obs"))

require(ggthemr)
ggthemr_reset()
ggthemr('fresh', type = 'inner')
darken_swatch(amount = 0.3)



pdf("abunda.pdf",height = 5,width = 7)
 ggplot(df1, aes(x=iso,y=resp,group=group, colour = group,fill=group)) +
#   geom_bar(stat="identity", position="dodge") +  
   geom_errorbar(aes(ymin = resp-lower, ymax = resp+upper), width = 0.35,
                 position = dodge,size=1.25) +
  geom_point(position = dodge,aes(size=seg)) +
   scale_size_manual(values=c(3,0))+
   scale_color_manual(name="",values=c("#143054","#805715")) +
   scale_fill_manual(name="",values=c("#143054","#805715")) +
  facet_wrap(~iso,scales = "free",labeller = labeller(iso = labels,label_parsed)) +
  xlab("")  + 
 # theme_rafa() +
  ylab("") +  
  theme(axis.text.x = element_blank(),legend.position = "none") 
dev.off()


ukLang <- htmltab(doc = "/Users/Rafael/Desktop/nn.html")
colnames(ukLang)

plot(ukLang[,1],ukLang[,3])

ukLang[,3]