# New Figure
require(ggplot2)
require(dplyr)
require(viridis)
require(cowplot)
require(reshape2)
require(nuclear)
require("plot3D")

dream_dat <- read.table("He3dp_DREAM.dat",header = T)
ssDat <- dream_dat[,c("e0","er","gd2","gp2","ad","ap","ue1","ue2")]


xx <- 1e-4
S0 <- sfactor3Hedp_5p(xx,ssDat[,1],ssDat[,2],ssDat[,3],ssDat[,4],ssDat[,5],ssDat[,6],0)
S0_prior <- sfactor3Hedp_5p(xx,ssDat[,1],ssDat[,2],ssDat[,3],ssDat[,4],ssDat[,5],seq(2,10,length.out = 1e4),0)


quantile(S0,prob=c(0.16, 0.5, 0.84))
#dd <- data.frame(S0,ssDat[,5:6]) %>% as.mcmc() %>% ggs() %>% as_tibble() 

dd <- data.frame(S0,ssDat[,5:6]) %>% melt(id.vars = "S0") %>% 
  mutate(case="posterior")
ddp <- data.frame(S0=S0_prior,ad=ssDat[,5],ap=seq(2,10,length.out = 1e4)) %>% 
  melt(id.vars = "S0") %>% 
  mutate(case="prior")


gdd <- rbind(dd,ddp)

levels(gdd$variable) <- c("a[d]","a[p]")


gg <- ggplot(data=gdd,aes(y=S0,x=value,group=variable)) +
#  geom_point() +
#  stat_density_2d(aes(fill = stat(nlevel)),n=500,
#                  geom="polygon",contour = TRUE) +
  stat_density_2d(
    aes(fill = ..level..,alpha = ..level..),
    geom = "polygon",
    n = 100,
    contour = TRUE
  ) +
#  stat_ellipse(type = "t",level = 0.95,color="black",linetype="dashed",alpha=0.5,size=1) +
#  stat_ellipse(type = "t",level = 0.68,color="black",linetype="dashed",alpha=0.5,size=1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_fill_viridis() +
  scale_alpha(range=c(0.2,1)) +
  theme_bw() +
  theme(strip.text = element_text(size=25),
        panel.border = element_rect(size=1),
        legend.position = "none",
        axis.ticks = element_line(size = 0.45),
        text = element_text(size=25),
        axis.text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("(fm)") +
  ylab(expression(S[0]~(MeV~b))) +
  facet_grid(case~variable,scale="free",labeller = label_parsed)


pdf("plot/S0_a.pdf",height = 0.75*7,width = 0.75*14)
gg
dev.off()


ad <- filter(dd,variable=="ad")

ap <- filter(dd,variable=="ap")

gad <- ggplot(data=ad,aes(y=S0,x=value)) +
  stat_density_2d(aes(fill =..level..,alpha =..level..),n=500,
                  geom="polygon") +
  stat_ellipse(type = "t",level = 0.95,color="black",linetype="dashed",alpha=0.5,size=1) +
  stat_ellipse(type = "t",level = 0.68,color="black",linetype="dashed",alpha=0.5,size=1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_fill_viridis() +
  scale_alpha(range=c(0.2,1)) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=22),
        axis.text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(expression(a[d]~(fm))) +
  ylab(expression(S[0]~(MeV~b)))


gap <- ggplot(data=ap,aes(y=S0,x=value)) +
  stat_density_2d(aes(fill =..level..,alpha =..level..),n=500,
                  geom="polygon") +
  stat_ellipse(type = "t",level = 0.95,color="black",linetype="dashed",alpha=0.5,size=1) +
  stat_ellipse(type = "t",level = 0.68,color="black",linetype="dashed",alpha=0.5,size=1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_fill_viridis() +
  scale_alpha(range=c(0.2,1)) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=22),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(expression(a[p]~(fm))) +
  ylab("")


dall <- data.frame(S0,ssDat) %>% melt(id.vars = "S0") 


ggplot(data=dall,aes(y=S0,x=value,group=variable)) +
  #  stat_density_2d(aes(fill = stat(nlevel)),n=500,
  #                  geom="polygon",contour = TRUE) +
  stat_density_2d(
    aes(fill = ..level..,alpha = ..level..),
    geom = "polygon",
    n = 100,
    contour = TRUE
  ) +
  stat_ellipse(type = "t",level = 0.95,color="black",linetype="dashed",alpha=0.5,size=1) +
  stat_ellipse(type = "t",level = 0.68,color="black",linetype="dashed",alpha=0.5,size=1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_fill_viridis() +
  scale_alpha(range=c(0.2,1)) +
  theme_bw() +
  theme(panel.border = element_rect(size=1),
        legend.position = "none",
        axis.ticks = element_line(size = 0.45),
        text = element_text(size=22),
        axis.text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab(expression(S[0]~(MeV~b))) +
  facet_wrap(.~variable,scale="free_x")



sfactor3Hedp_5p(xx,ssDat[,1],ssDat[,2],ssDat[,3],ssDat[,4],ssDat[,5],ssDat[,6],0)


grid.lines = 100
x <- seq(2.5, 5, length.out = grid.lines)
y <- seq(2.5,7, length.out = grid.lines)
f <- function(x,y){sfactor3Hedp_5p(xx,median(ssDat[,1]),median(ssDat[,2]),
                                   median(ssDat[,3]),
                                   median(ssDat[,4]),x,y,0)}
z <- outer(x,y)

persp(x,y,z,theta = 50, phi = 30,ticktype = "detailed",
      xlab="ad",
      ylab="ap",
      zlab="S0")

y.pred <- matrix(sfactor3Hedp_5p(xx,ssDat[,1],ssDat[,2],ssDat[,3],ssDat[,4],x1x2[,1],x1x2[,2],0)
,nrow = grid.lines, ncol = grid.lines)


# fitted points for droplines to surface

surf(x = x1.pred, y = x2.pred, z = y.pred,   col = viridis(200), shade = 0.35,
       lwd=1.25,lty=2,colkey = FALSE,
          theta = 80, phi = 30, ticktype = "detailed",bty = "b2",
          xlab="x1",
          ylab="x2",
          zlab="y", 
          expand =0.6)



