# Read Dpg
require(ggplot2);require(ggthemes);require(dplyr);require(magrittr)

# Data of interest

## DATA SET 1: kudomi 2004
obsx1    <- c(0.0312, 0.0331, 0.0352, 0.0373,    
              0.0393, 0.0413, 0.0433, 0.0453 )
obsy1    <- c(6.40, 5.48, 5.62, 5.46,
              5.69, 5.51, 5.43, 5.39 )
errobsy1 <- c(0.39, 0.22, 0.21, 0.20,
              0.25, 0.18, 0.14, 0.09 )

Kud04 <- data.frame(obsx1,obsy1,errobsy1) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.038) %>%
  mutate(.,dat="Kud04")



## DATA SET 2: bonetti 1999
obsx2    <- c(0.01699, 0.01846, 0.01898, 0.01946, 0.01993, 0.02143,
              0.02337, 0.02436 )
obsy2    <- c(13.15, 7.86, 8.25, 7.67,  5.10, 4.72,
              7.31, 5.44 )
errobsy2 <- c(4.98, 2.97, 2.29, 2.22, 1.70, 0.65,
              0.63, 0.34 )


Bon99 <- data.frame(obsx2,obsy2,errobsy2) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.057) %>%
  mutate(.,dat="Bon99")



## DATA SET 3: junker 1998
obsx3    <- c(0.02076, 0.02123, 0.02175, 0.02228, 0.02233, 0.02278,    
              0.02282, 0.02315, 0.02321, 0.02370, 0.02425, 0.02430,    
              0.02452, 0.02470, 0.02480, 0.04582, 0.05064, 0.05594,    
              0.06106, 0.06606, 0.07122, 0.07629, 0.08150, 0.08651,    
              0.09170 )
obsy3    <- c(6.80, 7.15, 7.63, 5.85, 7.27, 5.97,
              7.21, 6.82, 7.50, 6.87, 6.66, 6.90,
              7.10, 6.23, 5.96, 6.14, 5.63, 5.50,
              5.41, 5.43, 5.43, 5.32, 5.33, 5.23,
              5.15 )
errobsy3 <- c(0.82, 1.06, 0.91, 0.89, 1.05, 0.64,
              0.84, 1.47, 1.02, 0.74, 0.74, 0.72,
              0.79, 0.37, 0.62, 0.23, 0.14, 0.16,
              0.14, 0.15, 0.14, 0.11, 0.12, 0.11,
              0.11 )


Jun98 <- data.frame(obsx3,obsy3,errobsy3) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.037) %>%
  mutate(.,dat="Jun98")


## DATA SET 4: krauss 1987
obsx4    <- c(0.02451, 0.02655, 0.02900, 0.03145, 0.03390, 0.03634, 0.03909,    
              0.04124, 0.04373, 0.04648, 0.04808, 0.04900, 0.04932, 0.05440,    
              0.05940, 0.06440, 0.06460, 0.06800, 0.06930, 0.07270, 0.07340,    
              0.07780, 0.07940, 0.08450, 0.08630, 0.08950, 0.09160, 0.09400,    
              0.09720, 0.10340, 0.10920, 0.11600, 0.12150, 0.13360, 0.14130,    
              0.14600, 0.15630, 0.15790, 0.16890, 0.17050, 0.19540, 0.21980,    
              0.24430, 0.26880, 0.29330, 0.31790, 0.34250 )
obsy4    <- c(5.07, 5.18, 5.23, 5.45, 5.26, 5.35, 5.77,
              5.03, 4.88, 4.98, 5.08, 5.06, 5.86, 5.71,
              5.10, 5.18, 5.56, 5.39, 5.93, 5.30, 5.55,
              5.27, 5.26, 5.12, 4.92, 5.31, 4.69, 4.86,
              4.97, 4.93, 4.77, 4.89, 4.67, 4.56, 4.62,
              4.97, 4.63, 4.56, 4.67, 4.73, 4.68, 4.35,
              4.57, 4.73, 5.09, 4.40, 4.41 )
errobsy4 <- c(1.34, 1.06, 0.58, 0.45, 0.52, 0.41, 0.35,
              0.43, 0.24, 0.26, 0.16, 0.19, 0.32, 0.32,
              0.36, 0.20, 0.23, 0.31, 0.17, 0.18, 0.25,
              0.20, 0.18, 0.18, 0.13, 0.30, 0.07, 0.08,
              0.08, 0.10, 0.16, 0.08, 0.08, 0.13, 0.09,
              0.10, 0.05, 0.08, 0.05, 0.05, 0.21, 0.22,
              0.19, 0.26, 0.28, 0.26, 0.24 )

Kra87 <- data.frame(obsx4,obsy4,errobsy4) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.045) %>%
  mutate(.,dat="Kra87")

## DATA SET 5: dawarakanath 1971
obsx5    <- c(0.088, 0.126,	0.155, 0.193, 0.234, 0.288, 0.338,	   
              0.379, 0.435,	0.488, 0.591, 0.691, 0.746, 0.792,	   
              0.895, 0.997, 1.081 )
obsy5    <- c(4.864, 4.882, 4.956, 4.513, 4.680, 4.449, 4.342,
              4.503, 4.214, 4.131, 3.901, 3.762, 3.699, 3.503,
              3.511, 3.661, 3.504 )
errobsy5 <- c(0.340, 0.340, 0.350, 0.316, 0.328, 0.311, 0.174,
              0.180, 0.169, 0.165, 0.156, 0.150, 0.148, 0.140,
              0.140, 0.146, 0.140 )

Daw71 <- data.frame(obsx5,obsy5,errobsy5) %>%
  set_colnames(c("E","S","Stat")) %>%
  mutate(.,Syst=0.082) %>%
  mutate(.,dat="Daw71")




He3He3  <- rbind(Kud04,Bon99,Jun98,Kra87,Daw71)

save(He3He3, file = "He3He3.RData")


#pdf("plot/3Hedp_data.pdf",height = 7,width = 8)
ggplot(He3He3,aes(x=E,y=S,group=dat,shape=dat,color=dat))+
  geom_point(size=2)+
  geom_errorbar(show.legend = FALSE,aes(x=E,y=S,ymin=S-Stat,ymax=S+Stat),width=0.025)+
  scale_shape_stata(name="Dataset")+
  scale_color_stata(name="Dataset")+
  scale_x_log10() + scale_y_log10() +
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

