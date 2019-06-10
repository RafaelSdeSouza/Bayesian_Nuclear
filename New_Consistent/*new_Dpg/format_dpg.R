# Read Dpg
require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)

# E is in MeV, S is in keVb; convert latter to MeV b
dpg_th <- read.table("Marcucci2005.dat", header = FALSE) %>%
  set_colnames(c("E","S"))  %>% mutate(S = S*1e-6)

save(dpg_th,file="dpg_th.RData")


# Data of interest

## DATA SET 1: ma97a
obsx1    <- c(7.490E-02, 1.070E-01, 1.330E-01, 1.730E-01)
obsy1    <- c(6.850E-07, 7.080E-07, 9.560E-07, 1.260E-06)
errobsy1 <- c(7.020E-08, 6.840E-08, 8.400E-08, 9.820E-08)


ma97a <- data.frame(obsx1,obsy1,errobsy1) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.09) %>%
  mutate(.,dat="ma97a")


## DATA SET 2: bys08a
obsx2    <- c(8.280E-03, 9.490E-03, 10.10E-03)
obsy2    <- c(2.370E-07, 2.770E-07, 2.980E-07)
errobsy2 <- c(7.100E-08, 6.400E-08, 6.500E-08)


bys08a <- data.frame(obsx2,obsy2,errobsy2) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.08) %>%
  mutate(.,dat="bys08a")



## DATA SET 3: sch97a
obsx3    <- c(1.000E-02, 1.670E-02, 2.330E-02, 3.000E-02, 3.670E-02, 4.330E-02,
              5.000E-02)
obsy3    <- c(2.425E-07, 2.740E-07, 3.452E-07, 3.974E-07, 4.452E-07, 4.738E-07,
              4.744E-07)
errobsy3 <- c(1.250E-08, 7.500E-09, 6.500E-09, 6.100E-09, 5.700E-09, 7.200E-09,
              6.400E-09)

sch97a <- data.frame(obsx3,obsy3,errobsy3) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.09) %>%
  mutate(.,dat="sch97a")


## DATA SET 4: cas02a
obsx4    <- c( 2.600E-03, 3.200E-03, 3.200E-03, 3.800E-03, 4.200E-03,
               4.600E-03, 4.600E-03, 5.200E-03, 5.800E-03, 6.100E-03, 6.500E-03,
               6.700E-03, 7.200E-03, 7.800E-03, 8.400E-03, 8.500E-03, 9.100E-03,
               9.100E-03, 9.700E-03, 9.800E-03, 1.040E-02, 1.040E-02, 1.050E-02,
               1.060E-02, 1.120E-02, 1.120E-02, 1.190E-02, 1.190E-02, 1.230E-02,
               1.240E-02, 1.320E-02, 1.320E-02, 1.380E-02, 1.440E-02, 1.460E-02,
               1.500E-02, 1.580E-02, 1.590E-02, 1.630E-02, 1.720E-02, 1.720E-02,
               1.760E-02, 1.850E-02, 1.860E-02, 1.910E-02, 1.970E-02, 1.990E-02,
               1.980E-02, 2.040E-02, 2.120E-02)
obsy4    <- c( 1.590E-07, 2.370E-07, 2.560E-07, 2.420E-07, 2.330E-07,
               2.290E-07, 2.510E-07, 2.490E-07, 2.430E-07, 2.520E-07, 2.620E-07,
               2.610E-07, 2.690E-07, 2.730E-07, 2.730E-07, 2.640E-07, 2.620E-07,
               2.710E-07, 2.700E-07, 2.790E-07, 3.010E-07, 2.850E-07, 2.900E-07,
               2.860E-07, 2.880E-07, 2.800E-07, 2.900E-07, 2.890E-07, 2.770E-07,
               2.890E-07, 3.000E-07, 2.680E-07, 3.110E-07, 2.900E-07, 3.460E-07,
               2.990E-07, 2.820E-07, 3.250E-07, 2.960E-07, 3.140E-07, 3.390E-07,
               3.020E-07, 3.240E-07, 3.550E-07, 3.250E-07, 3.280E-07, 3.280E-07,
               3.320E-07, 3.090E-07, 3.280E-07)
errobsy4 <- c( 4.900E-08, 4.000E-08, 3.800E-08, 2.900E-08, 1.300E-08,
               1.500E-08, 1.500E-08, 1.800E-08, 1.500E-08, 1.300E-08, 1.100E-08,
               1.300E-08, 1.100E-08, 9.700E-09, 9.170E-08, 1.700E-08, 4.200E-08,
               1.900E-08, 1.800E-08, 1.600E-08, 6.200E-08, 1.900E-08, 5.800E-08,
               1.500E-08, 3.000E-08, 1.400E-08, 3.100E-08, 2.800E-08, 2.400E-08,
               2.500E-08, 2.500E-08, 2.700E-08, 2.200E-08, 1.900E-08, 2.400E-08,
               2.400E-08, 3.400E-08, 1.600E-08, 2.000E-08, 2.500E-08, 1.400E-08,
               1.800E-08, 2.200E-08, 1.400E-08, 1.900E-08, 1.800E-08, 2.200E-08,
               3.100E-08, 1.700E-08, 1.200E-08)


cas02a <- data.frame(obsx4,obsy4,errobsy4) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.045) %>%
  mutate(.,dat="cas02a")

dpg  <- rbind(ma97a,bys08a,sch97a,cas02a)

save(dpg, file = "dpg.RData")


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

