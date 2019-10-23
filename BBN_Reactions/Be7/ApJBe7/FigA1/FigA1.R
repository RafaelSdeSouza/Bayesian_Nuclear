# Reduced cross section, reciprocity theorem
require(readr)
dtop <- read.table("data1.txt",header=F) %>% 
  mutate(type = rep( c("class. kinematics w/o B.E.",
                       "class. kinematics with B.E.",
                        "rel. kinematics with B.E."),each=38)) %>%
  mutate(reaction = "Gib59")

dtop2 <- read.table("data2.txt",header=F) %>% 
  mutate(type = rep( c("class. kinematics w/o B.E.",
                       "class. kinematics with B.E.",
                       "rel. kinematics with B.E."),each=46)) %>%
  mutate(reaction = "Her19")

dd <- rbind(dtop,dtop2)

anno <- data.frame(x1 = c(9.813e-3, 1.987e-3), x2 = c(9.813e-3, 1.987e-3), 
                   y1 = c(3, 3.5), y2 = c(4.75, 8.5),
                   reaction = c("Gib59", "Her19"))


pdf("fig_compData.pdf",height = 7.5,width = 6)
ggplot(dd,aes(x=V1,y=V2,group=type)) +
  geom_point(aes(shape=type,color=type),size=2) +
  scale_x_log10(breaks = c(1e-5,1e-4,1e-3,1e-2,1e-1,1),
                labels=c(expression(10^-5),
                         expression(10^-4),expression(10^-3),
                         expression(10^-2),expression(10^-1),"1")) +
  coord_cartesian(xlim=c(1e-5,1))+
  scale_color_discrete_qualitative(name="") + 
  scale_shape_manual(values=c(2,0,5),name="",
  guide = guide_legend(ncol = 1,override.aes = list(size = 2.5))) +
  theme_economist_white() +
  geom_segment(data = anno, aes(x = x1, xend = x2, 
                                y = y1, yend = y2,group=reaction),
               colour = "black",linetype="dashed") +
  facet_wrap(.~ reaction,ncol=1,scale="free_y") +
  ylab(expression(paste(sigma[np],sqrt(E[paste(NULL^7,"Be+n")]), " (", sqrt(MeV), "b)"))) +
  xlab(expression(paste(E[paste(NULL^7,"Be+n")], " (", "MeV)"))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.725,0.375),
        legend.text = element_text(size=13),
        legend.text.align = 1,
        legend.background = element_rect(colour = "white",fill=NA),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="black", size=18),
        axis.text  = element_text(size=16),
        strip.text = element_text(size=15),
        strip.background = element_rect("gray95"),
        axis.title.x = element_text(margin = margin(t = 10)))
dev.off()