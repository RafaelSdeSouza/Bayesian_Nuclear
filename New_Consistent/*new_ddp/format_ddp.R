# Read 3HeDp
require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)

# E is in MeV, S is in keVb; convert latter to MeV b
ddp_th <- read.table("Arai_ddp_2011.dat", header = FALSE) %>%
  set_colnames(c("E","S"))  %>% mutate(S = S*1e-3)

save(ddp_th,file="ddp_th.RData")


# Data of interest

## DATA SET 1: Leo06
obsx1    <- c(5.6100E-02, 8.6400E-02, 1.1640E-01, 1.5730E-01, 1.8960E-01,
              2.3510E-01, 2.7865E-01, 3.2305E-01)
obsy1    <- c(7.3361E+01, 7.4820E+01, 8.0457E+01, 8.8004E+01, 9.8428E+01,
              1.0653E+02, 1.2432E+02, 1.3815E+02) * 1e-3
errobsy1 <- c(3.1541E+00, 2.2281E+00, 1.7523E+00, 2.0381E+00, 1.7427E+00,
              1.4576E+00, 1.7546E+00, 1.9829E+00) * 1e-3

Leo06 <- data.frame(obsx1,obsy1,errobsy1) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.02) %>%
  mutate(.,dat="Leo06")


obsx2    <- c(2.4660E-02, 3.9620E-02, 4.1000E-02, 5.0000E-02, 5.2000E-02,
              7.6000E-02, 1.0200E-01, 1.2800E-01)
obsy2    <- c(6.4500E+01, 6.8500E+01, 7.0100E+01, 6.9000E+01, 6.7500E+01,
              7.4000E+01, 7.8500E+01, 8.6600E+01) * 1e-3
errobsy2 <- c(2.7000E+00, 2.9000E+00, 3.3000E+00, 2.5000E+00, 2.4000E+00,
              2.5000E+00, 2.7000E+00, 3.3000E+00) * 1e-3


Gre95 <- data.frame(obsx2,obsy2,errobsy2) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.03) %>%
  mutate(.,dat="Gre95")



## DATA SET 3: Bro90
obsx3    <- c(1.9964E-02, 2.4961E-02, 2.9959E-02, 3.4957E-02, 3.9956E-02,
              4.4956E-02, 4.9955E-02, 5.4955E-02, 5.8455E-02)
obsy3    <- c(5.9390E+01, 5.9960E+01, 6.1130E+01, 6.1980E+01, 6.2610E+01,
              6.4210E+01, 6.4990E+01, 6.5850E+01, 6.6770E+01) * 1e-3
errobsy3 <- c(5.2000E-01, 4.7000E-01, 4.5000E-01, 4.2000E-01, 4.0000E-01,
              2.7000E-01, 3.8000E-01, 3.7000E-01, 5.0000E-01) * 1e-3


Bro90 <- data.frame(obsx3,obsy3,errobsy3) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.013) %>%
  mutate(.,dat="Bro90")


## DATA SET 4: Kra87 (B)
obsx4    <- c(1.9880E-02, 2.2360E-02, 2.9820E-02, 3.4780E-02, 3.9790E-02,
              4.4720E-02, 4.9670E-02)
obsy4    <- c(5.2900E+01, 5.3700E+01, 5.5100E+01, 5.6200E+01, 5.6200E+01,
              5.8200E+01, 5.6400E+01) * 1e-3
errobsy4 <- c(3.4000E+00, 3.4000E+00, 3.5000E+00, 3.6000E+00, 3.6000E+00,
              3.7000E+00, 3.6000E+00) * 1e-3

Kra87B <- data.frame(obsx4,obsy4,errobsy4) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.064) %>%
  mutate(.,dat="Kra87B")


## DATA SET 5: Kra87 (M)
obsx5    <- c(1.9600E-02, 2.4000E-02, 2.9000E-02, 3.4800E-02, 3.9700E-02,
              4.4600E-02, 4.9500E-02, 5.4300E-02, 5.9300E-02, 6.4300E-02,
              7.4000E-02, 8.3800E-02, 9.3600E-02, 1.0340E-01, 1.1320E-01,
              1.2300E-01, 1.3300E-01, 1.4250E-01, 1.5250E-01, 1.6250E-01)
obsy5    <- c(5.6800E+01, 5.5700E+01, 5.5300E+01, 5.4800E+01, 5.6200E+01,
              5.6800E+01, 5.6400E+01, 5.9200E+01, 6.3000E+01, 6.3100E+01,
              6.6800E+01, 7.0000E+01, 7.2600E+01, 7.6100E+01, 7.7800E+01,
              8.0800E+01, 8.4400E+01, 8.8300E+01, 8.5100E+01, 9.1400E+01) * 1e-3
errobsy5 <- c(8.1000E+00, 7.4000E+00, 5.8000E+00, 5.7000E+00, 4.7000E+00,
              4.6000E+00, 4.6000E+00, 4.8000E+00, 5.1000E+00, 5.1000E+00,
              5.4000E+00, 5.7000E+00, 5.9000E+00, 6.2000E+00, 6.3000E+00,
              6.6000E+00, 6.9000E+00, 7.2000E+00, 6.9000E+00, 7.4000E+00) * 1e-3

Kra87M <- data.frame(obsx5,obsy5,errobsy5) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.082) %>%
  mutate(.,dat="Kra87M")


ddp  <- rbind(Leo06,Gre95,Bro90,Kra87B,Kra87M)

save(ddp, file = "ddp.RData")


#pdf("plot/3Hedp_data.pdf",height = 7,width = 8)
ggplot(ddp,aes(x=E,y=S,group=dat,shape=dat,color=dat))+
  geom_point(size=2)+
  geom_errorbar(show.legend = FALSE,aes(x=E,y=S,ymin=S-Stat,ymax=S+Stat),width=0.025)+
  scale_shape_stata(name="Dataset")+
  scale_color_stata(name="Dataset")+
  theme_bw() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + scale_x_log10()  +
  theme(legend.position = "top",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  ggtitle(expression(paste("d(d,p)",NULL^"3","H")))
#dev.off()

