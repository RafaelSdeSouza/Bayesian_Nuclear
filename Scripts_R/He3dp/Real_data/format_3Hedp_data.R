# Read 3HeDp 
require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)

#--Get File names
path <- "lab_data"


# Data of interest


Ali01a <-  read.table(paste(path,"data_aliotta_a.dat",sep="/"),header = F) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.03) %>% 
  mutate(.,invK = "No") %>%
  mutate(.,dat="Ali01a") 

Ali01b <-  read.table(paste(path,"data_aliotta_b.dat",sep="/"),header = F) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.03) %>% 
  mutate(.,invK = "Yes") %>%
  mutate(.,dat="Ali01b")

Cos00 <-  read.table(paste(path,"data_costantini.dat",sep="/"),header = F) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.055) %>%  
  mutate(.,invK = "Yes") %>%
  mutate(.,dat="Cos00") 



Mol80 <- read.table(paste(path,"data_moeller.dat",sep="/"),header = F) %>%
         set_colnames(.,c("E","S","Stat")) %>%
         mutate(.,Syst = 0.039) %>%
         mutate(.,invK = "No") %>%
         mutate(.,dat="Mol80")

Kra87 <- read.table(paste(path,"data_krauss.dat",sep="/"),header = F,sep="") %>% 
          set_colnames(.,c("E","S","Stat")) %>%
          mutate(.,Syst= 0.078) %>% 
          mutate(.,invK = "No") %>%
          mutate(.,dat="Kra87")


zhi77b <- read.table(paste(path,"data_zhichang.dat",sep="/"),header = F) %>%
          set_colnames(c("E","S","Stat"))  %>%  
          mutate(.,Syst=0.034) %>%  
          mutate(.,invK = "No") %>%
          mutate(.,dat="Zhi77") 


gei99 <- read.table(paste(path,"data_geist.dat",sep="/"),header = F) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst= 0.043) %>%
  mutate(.,invK = "No") %>%
  mutate(.,dat="Gei99")


ensamble <- rbind(Ali01a,Ali01b,Cos00,Mol80,Kra87,zhi77b,gei99 )

write.csv(ensamble,"ensamble.csv",row.names = F)


pdf("plot/3Hedp_data.pdf",height = 7,width = 8)
ggplot(ensamble,aes(x=E,y=S,color=invK,group=dat,shape=dat))+
  geom_point(size=2)+
  geom_errorbar(show.legend = FALSE,aes(x=E,y=S,ymin=S-Stat,ymax=S+Stat,group=invK,color=invK),width=0.025)+
  scale_colour_stata(name="Inverse Kinematics")+
  scale_shape_stata(name="Dataset")+
  theme_economist_white() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + scale_x_log10()  +
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

