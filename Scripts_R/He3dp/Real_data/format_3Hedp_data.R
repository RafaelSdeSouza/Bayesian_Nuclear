# Read 3HeDp 
require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)

#--Get File names
path <- "../data/3Hedp/hdp_data-2017"
files <- list.files(path = path)

# Data of interest

Lac05 <-  read.table(paste(path,"hdp_thm05b.dat",sep="/"),header = F) %>%
  set_colnames(.,c("E","E_rr","S","Stat")) %>%
  mutate(.,Syst=0.15) %>%  
  mutate(.,dat="Lac05") %>% 
  select(.,c(-E_rr))


Mol80 <- read.table(paste(path,"hdp_mol80b.dat",sep="/"),header = F) %>%
         set_colnames(.,c("E","Syst","S","Stat")) %>%
         mutate(.,Syst= mean(c(0.28,0.55,0.65,0.72,0.67,0.69,0.58,
                          0.52,0.44,0.4,0.29,0.2,0.15,0.11,
                          0.09,0.08,0.07,0.06)/S)) %>%
        select(.,c("E","S","Stat","Syst")) %>%
        mutate(.,dat="Mol80")

Kra87m <- read.csv(paste(path,"hdp_kra87m.csv",sep="/"),header = T,sep="") %>% 
          mutate(.,Syst= 0.078) %>% 
          mutate(.,dat="Kra87m")

Kra87b <- read.csv(paste(path,"hdp_kra87b.csv",sep="/"),header = T,sep="")  %>% 
              mutate(.,Syst= 0.078) %>%  
              mutate(.,dat="Kra87b")

zhi77b <- read.table(paste(path,"hdp_zhi77b.dat",sep="/"),header = F) %>%
          set_colnames(c("E","E_rr","S","Stat"))  %>%  
          mutate(.,Syst=0.034) %>%  
          mutate(.,dat="zhi77b") %>% 
          select(.,c(-E_rr))

gei99b <- read.table(paste(path,"hdp_gei99b.dat",sep="/"),header = F) %>%
  set_colnames(.,c("E","Syst","S","Stat")) %>%
  mutate(.,Syst= 0.043) %>%
  select(.,c("E","S","Stat","Syst")) %>%
  mutate(.,dat="gei99b")


gei99d <- read.table(paste(path,"hdp_gei99d.dat",sep="/"),header = F) %>%
  set_colnames(.,c("E","Syst","S","Stat")) %>%
  mutate(.,Syst= 0.043) %>%
  select(.,c("E","S","Stat","Syst")) %>%
  mutate(.,dat="gei99d")


ensamble <- rbind(Lac05,Mol80,Kra87m,Kra87b,zhi77b,gei99b,gei99d)

write.csv(ensamble,"ensamble.csv",row.names = F)

ggplot(ensamble,aes(x=E,y=S,color=dat))+geom_point(size=2.85)+
  theme_bw()+scale_x_log10()+scale_color_hc() +
  coord_cartesian(xlim=c(5e-3,1))+geom_vline(xintercept = 0.2, linetype="dotted", 
  color = "blue", size=1.5)



