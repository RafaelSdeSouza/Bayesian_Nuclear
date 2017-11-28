# Read 3HeDp 
require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)

#--Get File names
path <- "/Users/Rafael/Documents/GitHub/JAGS_UNC/data/3Hedp"
files <- list.files(path = path)

# Data of interest


Ali01 <-  read.table(paste(path,"Ali01.dat",sep="/"),header = F) %>%
  set_colnames(.,c("E","S","Stat","E2","S2","Stat2")) %>%
  mutate(.,Syst=0.03) %>%  
  mutate(.,dat="Ali01") %>% 
  select(.,c("E","S","Stat","Syst","dat"))



Cos00 <-  read.table(paste(path,"data_cos00a.dat",sep="/"),header = F) %>%
  set_colnames(.,c("Ebeam","E","S","Stat","Syst")) %>%
  mutate(.,Syst=0.055) %>%  
  mutate(.,dat="Cos00") %>% 
  select(.,c(-Ebeam))




#Lac05 <-  read.table(paste(path,"hdp_thm05b.dat",sep="/"),header = F) %>%
 # set_colnames(.,c("E","E_rr","S","Stat")) %>%
#  mutate(.,Syst=0.15) %>%  
#  mutate(.,dat="Lac05") %>% 
#  select(.,c(-E_rr))


Mol80 <- read.table(paste(path,"data_moeller.dat",sep="/"),header = F) %>%
         set_colnames(.,c("E","S","Stat")) %>%
         mutate(.,Syst= mean(c(0.28,0.55,0.65,0.72,0.67,0.69,0.58,
                          0.52,0.44,0.4,0.29,0.2,0.15,0.11,
                          0.09,0.08,0.07,0.06)/S)) %>%
        select(.,c("E","S","Stat","Syst")) %>%
        mutate(.,dat="Mol80")

Kra87 <- read.table(paste(path,"data_krauss.dat",sep="/"),header = F,sep="") %>% 
          set_colnames(.,c("E","S","Stat")) %>%
          mutate(.,Syst= 0.078) %>% 
          mutate(.,dat="Kra87")

#Kra87b <- read.csv(paste(path,"hdp_kra87b.csv",sep="/"),header = T,sep="")  %>% 
#             mutate(.,Syst= 0.078) %>%  
#             mutate(.,dat="Kra87b")

zhi77b <- read.table(paste(path,"data_zhichang.dat",sep="/"),header = F) %>%
          set_colnames(c("E","S","Stat"))  %>%  
          mutate(.,Syst=0.034) %>%  
          mutate(.,dat="Zhi77") 

#gei99b <- read.table(paste(path,"hdp_gei99b.dat",sep="/"),header = F) %>%
#  set_colnames(.,c("E","Syst","S","Stat")) %>%
#  mutate(.,Syst= 0.043) %>%
 # select(.,c("E","S","Stat","Syst")) %>%
#  mutate(.,dat="gei99")


gei99 <- read.table(paste(path,"data_geist.dat",sep="/"),header = F) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst= 0.043) %>%
  select(.,c("E","S","Stat","Syst")) %>%
  mutate(.,dat="Gei99")


ensamble <- rbind(Ali01,Cos00,Mol80,Kra87,zhi77b,gei99 )

write.csv(ensamble,"ensamble.csv",row.names = F)

ggplot(ensamble,aes(x=E,y=S,color=dat))+geom_point(size=2.85)+
  theme_bw()+scale_x_log10()+scale_color_hc() +
  coord_cartesian(xlim=c(5e-3,1))+geom_vline(xintercept = 0.2, linetype="dotted", 
  color = "blue", size=1.5)


pdf("plot/3Hedp_data.pdf",height = 7,width = 8)
ggplot(ensamble,aes(x=E,y=S,color=dat,group=dat,shape=dat))+
  geom_point(size=2)+
  geom_errorbar(show.legend = FALSE,aes(x=E,y=S,ymin=S-Stat,ymax=S+Stat,group=dat,color=dat),width=0.025)+
  scale_colour_stata(name="Dataset")+
  scale_shape_stata(name="Dataset")+
  theme_wsj() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + scale_x_log10()  +
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
dev.off()

