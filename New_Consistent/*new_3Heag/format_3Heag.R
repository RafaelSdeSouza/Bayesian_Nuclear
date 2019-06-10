# Read Dpg
require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)

# E is in MeV, S is in keVb; convert latter to MeV b
He3ag_th <- read.table("TheoryNeff.dat", header = FALSE) %>%
  set_colnames(c("E","S")) 

save(He3ag_th,file="He3ag_th.RData")


# Data of interest

## DATA SET 1: di Leva 2009 [recoil only]
obsx1    <- c(0.701, 0.802, 0.902, 1.002, 1.002, 1.102, 
              1.102, 1.103, 1.203, 1.203, 1.353, 1.403, 
              1.403, 1.403, 1.504)

obsy1    <- c(0.393e-3, 0.385e-3, 0.339e-3, 0.351e-3, 0.333e-3, 0.334e-3, 
              0.339e-3, 0.334e-3, 0.333e-3, 0.333e-3, 0.327e-3, 0.343e-3, 
              0.340e-3, 0.343e-3, 0.339e-3)

errobsy1 <- c(0.069e-3, 0.021e-3, 0.015e-3, 0.013e-3, 0.011e-3, 0.003e-3, 
              0.006e-3, 0.009e-3, 0.007e-3, 0.012e-3, 0.008e-3, 0.004e-3, 
              0.009e-3, 0.011e-3, 0.010e-3)


Dil09 <- data.frame(obsx1,obsy1,errobsy1) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.05) %>%
  mutate(.,Experiment= "Recoil") %>%
  mutate(.,dat="Dil09")


## DATA SET 2: nara singh 2004 [activation]
obsx2    <- c(0.4200, 0.5060, 0.6145, 0.9500)
obsy2    <- c(0.420e-3, 0.379e-3, 0.362e-3, 0.316e-3)
errobsy2 <- c(0.014e-3, 0.015e-3, 0.010e-3, 0.006e-3)


Nar04 <- data.frame(obsx2,obsy2,errobsy2) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.051) %>%
  mutate(.,Experiment= "Activation") %>%
  mutate(.,dat="Nar04")



## DATA SET 3: brown 2007 [activation only]
obsx3    <- c(0.3274, 0.4260, 0.5180, 0.5815, 0.7024, 0.7968, 
              1.2337, 1.2347)
obsy3    <- c(0.495e-3, 0.458e-3, 0.440e-3, 0.400e-3, 0.375e-3, 0.363e-3, 
              0.330e-3, 0.324e-3)
errobsy3 <- c(0.015e-3, 0.010e-3, 0.010e-3, 0.011e-3, 0.010e-3, 0.007e-3, 
              0.006e-3, 0.006e-3)

Bro07 <- data.frame(obsx3,obsy3,errobsy3) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.03) %>%
  mutate(.,Experiment= "Activation") %>%
  mutate(.,dat="Bro07")


## DATA SET 4: costantini 2008 [activation only]
obsx4    <- c(0.0929, 0.1057, 0.1265, 0.1477, 0.1689, 0.1695)
obsy4    <- c(0.534e-3, 0.493e-3, 0.514e-3, 0.499e-3, 0.482e-3, 0.507e-3)
errobsy4 <- c(0.016e-3, 0.015e-3, 0.020e-3, 0.017e-3, 0.020e-3, 0.010e-3)


Cos08 <- data.frame(obsx4,obsy4,errobsy4) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.031) %>%
  mutate(.,Experiment= "Activation") %>%
  mutate(.,dat="Cos08")

He3ag  <- rbind(Dil09,Nar04,Bro07,Cos08)

save(He3ag, file = "He3ag.RData")


#pdf("plot/3Hedp_data.pdf",height = 7,width = 8)
ggplot(He3ag,aes(x=E,y=S,group=dat,shape=dat,color=dat))+
  geom_point(size=2)+
  geom_errorbar(show.legend = FALSE,aes(x=E,y=S,ymin=S-Stat,ymax=S+Stat),width=0.025)+
  scale_shape_stata(name="Dataset")+
  scale_color_stata(name="Dataset")+
  theme_bw() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + 
  theme(legend.position = "top",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) 
#dev.off()

