# Read Tdn
require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)

#--Get File names


# Data of interest

# extra

Mag75  <-  read.table("tdn_mag75b.dat") %>%
  set_colnames(.,c("E","E_stat","S","Stat")) %>%
  select(.,c("E","S","Stat")) %>%
   mutate(.,Syst=mean(0.5*S)) %>%  
  mutate(.,dat="Mag75") 




Bam57  <-  read.table("tdn_bam57b.dat") %>%
  set_colnames(.,c("E","E_stat","S","Stat")) %>%
  select(.,c("E","E_stat","S","Stat")) %>%
  mutate(.,Syst=mean(0.5*S)) %>%  
  mutate(.,dat="Bam57") 


# Current
Arn53  <-  read.table("data_arnold.dat") %>%
  set_colnames(.,c("E","E_stat","S","Stat")) %>%
  mutate(.,Syst=mean(0.02*S)) %>%  
  mutate(.,dat="Arn53") 


Bro87 <- read.table("data_brown.dat") %>%
  set_colnames(.,c("E","E_stat","S","Stat")) %>%
         mutate(.,Syst= mean(0.013*S)) %>%
        mutate(.,dat="Bro87")

Jar84 <- read.table("data_jarmie.dat") %>%
  set_colnames(.,c("E","E_stat","S","Stat")) %>%
  mutate(.,Syst= mean(0.0126*S)) %>%
  mutate(.,dat="Jar84")


Kob66 <- read.table("data_kobzev.dat") %>%
  set_colnames(.,c("E","E_stat","S","Stat")) %>%
  mutate(.,Syst = mean(0.025*S)) %>%
  mutate(.,dat="Kob66")


Con52 <- read.table("data_conner.dat") %>%
  set_colnames(.,c("E","E_stat","S","Stat")) %>%
  mutate(.,Syst = mean(0.018*S)) %>%
  mutate(.,dat="Con52")


# Remove bad points

ensamble <- rbind(Arn53, Bro87,Jar84,Kob66,Con52)

#ensamble <- rbind(Bame57, Brown87,Jarmie86)


write.csv(ensamble,"ensamble_Tdn_extra.csv",row.names = F)

ggplot(ensamble,aes(x=E,y=S,color=dat))+geom_point(size=2.85)+
  theme_bw()+scale_x_log10()+scale_color_hc() 


