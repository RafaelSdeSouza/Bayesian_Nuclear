# Read d(d,n)

# E is in MeV, S is in keVb; convert latter to MeV b
ddn_th <- read.table("Arai_ddn_2011.dat", header = FALSE) %>%
  set_colnames(c("E","S"))  %>% mutate(S = S*1e-3)

save(ddn_th,file="ddn_th.RData")


require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)


# Data of interest

## DATA SET 1: Leo06
obsx1    <- c(5.6100E-02, 8.6400E-02, 1.1640E-01, 1.5730E-01, 1.8960E-01,
              2.3510E-01, 2.7865E-01, 3.2305E-01)
obsy1    <- c(7.8185E+01, 8.3631E+01, 9.2168E+01, 1.0769E+02, 1.1606E+02,
              1.2672E+02, 1.4662E+02, 1.5629E+02) * 1e-3
errobsy1 <- c(3.4139E+00, 2.5573E+00, 2.0515E+00, 2.5573E+00, 2.0765E+00,
              1.8949E+00, 2.1567E+00, 2.9465E+00) * 1e-3

Leo06 <- data.frame(obsx1,obsy1,errobsy1) %>%
         set_colnames(.,c("E","S","Stat")) %>%
         mutate(.,Syst=0.02) %>%
         mutate(.,dat="Leo06")


## DATA SET 2: Gre95
obsx2    <- c(2.4660E-02, 3.9620E-02, 4.1000E-02, 5.0000E-02, 5.2000E-02,
              7.6000E-02, 1.0200E-01, 1.2800E-01)
obsy2    <- c(6.8700E+01, 7.3300E+01, 7.3000E+01, 7.5900E+01, 7.3300E+01,
              8.3000E+01, 9.0000E+01, 1.0060E+02) * 1e-3
errobsy2 <- c(2.7000E+00, 2.9000E+00, 2.8000E+00, 3.0000E+00, 2.8000E+00,
              2.8000E+00, 3.2000E+00, 3.5000E+00) * 1e-3

Gre95 <- data.frame(obsx2,obsy2,errobsy2) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.03) %>%
  mutate(.,dat="Gre95")



## DATA SET 3: Bro90
obsx3    <- c(1.9964E-02, 2.4961E-02, 2.9959E-02, 3.4957E-02, 3.9956E-02,
              4.4956E-02, 4.9955E-02, 5.4955E-02, 5.8495E-02)
obsy3    <- c(5.9990E+01, 6.2240E+01, 6.4310E+01, 6.5360E+01, 6.5320E+01,
              6.8410E+01, 6.9790E+01, 7.2500E+01, 7.2270E+01) * 1e-3
errobsy3 <- c(7.2000E-01, 6.2000E-01, 6.4000E-01, 5.6000E-01, 5.4000E-01,
              3.1000E-01, 5.5000E-01, 5.4000E-01, 6.4000E-01) * 1e-3

Bro90 <- data.frame(obsx3,obsy3,errobsy3) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.013) %>%
  mutate(.,dat="Bro90")


## DATA SET 4: Kra87 (B)
obsx4    <- c(1.9880E-02, 2.2360E-02, 2.9820E-02, 3.4780E-02, 3.9790E-02,
              4.4720E-02, 4.9670E-02)
obsy4    <- c(5.4600E+01, 5.6500E+01, 5.8200E+01, 5.9400E+01, 5.8900E+01,
              6.0900E+01, 5.9300E+01) * 1e-3
errobsy4 <- c(3.8000E+00, 3.8000E+00, 3.9000E+00, 3.9000E+00, 3.9000E+00,
              4.1000E+00, 3.9000E+00) * 1e-3

Kra87B <- data.frame(obsx4,obsy4,errobsy4) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.064) %>%
  mutate(.,dat="Kra87B")


## DATA SET 5: Kra87 (M)
obsx5    <- c(1.9600E-02, 2.4000E-02, 2.9000E-02, 3.4800E-02, 3.9700E-02,
              4.4600E-02, 4.9500E-02, 5.4300E-02, 5.9300E-02, 6.4300E-02,
              7.4000E-02, 8.3800E-02, 9.3600E-02, 1.0340E-01, 1.1320E-01,
              1.2300E-01, 1.3300E-01, 1.4250E-01, 1.5250E-01, 1.6250E-01)
obsy5    <- c(5.8300E+01, 5.6400E+01, 5.4400E+01, 5.4200E+01, 5.8900E+01,
              5.9600E+01, 6.1900E+01, 6.4600E+01, 6.8500E+01, 7.0900E+01,
              7.8600E+01, 8.3000E+01, 8.6700E+01, 8.9600E+01, 9.3200E+01,
              9.8200E+01, 1.0570E+02, 1.0930E+02, 1.0790E+02, 1.1420E+02) * 1e-3
errobsy5 <- c(8.1000E+00, 7.5000E+00, 5.7000E+00, 5.6000E+00, 4.9000E+00,
              4.8000E+00, 5.0000E+00, 5.2000E+00, 5.6000E+00, 5.8000E+00,
              6.4000E+00, 6.7000E+00, 7.0000E+00, 7.3000E+00, 7.6000E+00,
              8.0000E+00, 8.6000E+00, 8.9000E+00, 8.8000E+00, 9.3000E+00) * 1e-3

Kra87M <- data.frame(obsx5,obsy5,errobsy5) %>%
  set_colnames(.,c("E","S","Stat")) %>%
  mutate(.,Syst=0.082) %>%
  mutate(.,dat="Kra87M")


ddn  <- rbind(Leo06,Gre95,Bro90,Kra87B,Kra87M)

save(ddn, file = "ddn.RData")


#pdf("plot/3Hedp_data.pdf",height = 7,width = 8)
ggplot(ddndata,aes(x=E,y=S,group=dat,shape=dat,color=dat))+
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
  ggtitle(expression(paste("d(d,n)",NULL^"3","H")))
#dev.off()

