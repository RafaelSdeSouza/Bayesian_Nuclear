# Read 3HeDp 
require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)

#--Get File names


# Data of interest

Bame57 <-  read.table("data_bame.dat") %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=mean(0.02*S)) %>%  
  mutate(.,dat="Bam57") 


Brown87 <- read.table("data_brown.dat") %>%
         set_colnames(.,c("E","S","Stat")) %>%
         mutate(.,Syst= mean(0.014*S)) %>%
        mutate(.,dat="Bro87")

Jarmie86 <- read.table("data_jarmie.dat") %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst= mean(0.0126*S)) %>%
  mutate(.,dat="Jar86 ")


ensamble <- rbind(Bame57, Brown87,Jarmie86)

write.csv(ensamble,"ensamble_Tdn.csv",row.names = F)

ggplot(ensamble,aes(x=E,y=S,color=dat))+geom_point(size=2.85)+
  theme_bw()+scale_x_log10()+scale_color_hc() 


