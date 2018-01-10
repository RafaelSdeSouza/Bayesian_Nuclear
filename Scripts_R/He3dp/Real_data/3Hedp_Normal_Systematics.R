# 3Hedp analysis
#
# purpose: Real  DATA
#
# - 5 parameters are assumed: Er, gamma_d^2, gamma_n^2 [e1, gin, gout]
#
# - uses the function sfactorHe3dp(obsx1[i], e1, gin, gout), which
#   is a C++ version of a Fortran code that includes Coulomb wave
#   function calculations; JAGS has been recompiled with this C++ function
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
set.seed(123)
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb;

######################################################################
# import packages
library(rjags);library(R2jags);library(mcmcplots)
require(RcppGSL);require(ggplot2);require(ggthemes)
require(nuclear);library(magrittr);library(wesanderson)
library(dplyr);require(ggsci);require(ggmcmc);require(plyr);library(latex2exp)
source("/Users/Rafael/Documents/GitHub/JAGS_UNC/Scripts_R/auxiliar_functions/jagsresults.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")
load.module("nuclear")


######################################################################
## Read DATA 
ensamble <- read.csv("ensamble.csv",header = T) %>%
  mutate(Syst=replace(Syst,Syst==0.06,0.078))


re <- as.numeric(ensamble$dat)
Nre <- length(unique(ensamble$dat))
ik <- as.numeric(ensamble$invK)
Nik <- length(unique(ensamble$invK))
# Radius
# r_i = 6
# r_f = 5

# Literature
#  0.35779   # resonance energy
#  1.0085    # reduced width incoming
#  0.025425   # reduced width outgoing


N <- nrow(ensamble)
obsy <- ensamble$S    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat
set <- ensamble$dat
lab <- ensamble$invK
syst = c(0.03,unique(ensamble$Syst))
#syst <- syst[-3]

M <- 500
xx <- seq(min(obsx),max(obsx),length.out = M)

model.data <- list(obsy = obsy,    # Response variable
                   obsx =  obsx,   # Predictors
                   erry = erry,
                   N = nrow(ensamble), # Sample size
                   syst = syst,
                   Nre = Nre,
                   re = re,
                   Nik = Nik,
                   ik  = ik,
                   M = M,
                   xx = xx
)


######################################################################
Model <- "model{
# LIKELIHOOD
for (i in 1:N) {
obsy[i] ~ dnorm(y[i], pow(erry[i], -2))
y[i] ~ dnorm(scale[re[i]]*sfactor3Hedp(obsx[i], e1, gin, gout,ri,rf,ue[ik[i]]),pow(tau, -2))
#y[i] <- scale[re[i]]*sfactor3Hedp(obsx[i], e1, gin, gout)
}

# Predicted values

for (j in 1:M){

# Bare...

mux0[j] <- sfactor3Hedp(xx[j], e1, gin, gout,ri,rf,0)
yx0[j] ~ dnorm(mux0[j],pow(tau,-2))

# No inverse Kinematics 

mux1[j] <- sfactor3Hedp(xx[j], e1, gin, gout,ri,rf,ue[1])
yx1[j] ~ dnorm(mux1[j],pow(tau,-2))

# With inverse Kinematics 
mux2[j] <- sfactor3Hedp(xx[j], e1, gin, gout,ri,rf,ue[2])
yx2[j] ~ dnorm(mux1[j],pow(tau,-2))

}



for (k in 1:Nre){
scale[k] ~ dlnorm(log(1.0),pow(log(1+syst[k]),-2))
}

for (z in 1:Nik){
ue[z] ~ dunif(0,1e-3)
}


# PRIORS
# e1, gin, gout are defined as in tdn.f (by Alain Coc):
# resonance energy, initial reduced width, final reduced
# width;

tau ~  dgamma(0.01,0.01)
e1 ~   dnorm(0,0.01)T(0,)
gin ~ dnorm(0,0.01)T(0,)
gout ~ dbeta(2,5)
#ri ~ ddexp(5,100)
#rf ~ ddexp(4,100)

rf <- 4
ri <- 5

#gb ~ dbeta(2,2)
#gin ~ dunif(0,20)

#gout ~ dbeta(2,5)

# Channel radius
#  rb ~ dbeta(2,2)
#  rb2 ~ dbeta(2,2)

#  ri <- 4*rb + 3
#  rf <- 4*rb2 + 3


}"


######################################################################
# n.adapt:  number of iterations in the chain for adaptation (n.adapt)
#           [JAGS will use to choose the sampler and to assure optimum
#           mixing of the MCMC chain; will be discarded]
# n.udpate: number of iterations for burnin; these will be discarded to
#           allow the chain to converge before iterations are stored
# n.iter:   number of iterations to store in the final chain as samples
#           from the posterior distribution
# n.chains: number of mcmc chains
# n.thin:   store every n.thin element [=1 keeps all samples]


inits <- function () { list(e1 = runif(1,0.15,1),gout=runif(1,0.01,1)) }
# "f": is the model specification from above;
# data = list(...): define all data elements that are referenced in the



# JAGS model with R2Jags;
Normfit <- jags(data = model.data,
                inits = inits,
                parameters = c("e1", "gin", "gout","ue","tau", "ri","rf","mux0","mux1","mux2","scale"),
                model = textConnection(Model),
                n.thin = 1,
                n.chains = 3,
                n.burnin = 5000,
                n.iter = 10000)

jagsresults(x=Normfit , params=c("e1", "gin", "gout","ue","tau","ri","rf"),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))

require(emdbook)

s <- as.mcmc(Normfit)
s1 <- lump.mcmc.list(s)
HPDinterval(s1, prob = 0.95)[c("e1"),]

traplot(Normfit  ,c("e1", "gin", "gout","ri","rf"),style="plain")
denplot(Normfit  ,c("e1", "gin", "gout","ri","rf","ue"),style="plain")
caterplot(Normfit,c("scale","tau"),style="plain")



# Plot

y <- jagsresults(x=Normfit , params=c('mux1'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))
x <- xx
gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],lwr3=y[,"0.5%"],upr1=y[,"75%"],
                    upr2=y[,"97.5%"],upr3=y[,"99.5%"])
gobs <- data.frame(obsx,obsy,erry,set,lab)
gobs$set <- as.factor(gobs$set)



y2 <- jagsresults(x=Normfit , params=c('mux2'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))
gdata2 <- data.frame(x =xx, mean = y2[,"mean"],lwr1=y2[,"25%"],lwr2=y2[,"2.5%"],lwr3=y2[,"0.5%"],upr1=y2[,"75%"],
                    upr2=y2[,"97.5%"],upr3=y2[,"99.5%"])


y0 <- jagsresults(x=Normfit , params=c('mux0'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))

gdata0 <- data.frame(x =xx, mean = y0[,"mean"],lwr1=y0[,"25%"],lwr2=y0[,"2.5%"],lwr3=y0[,"0.5%"],upr1=y0[,"75%"],
                    upr2=y0[,"97.5%"],upr3=y0[,"99.5%"])


pdf("plot/He3dp_syst.pdf",height = 7,width = 10)
ggplot(gobs,aes(x=obsx,y=obsy))+
  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="gray90",alpha=0.3) +

  
  
  geom_ribbon(data=gdata2,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("#e5f5e0"),show.legend=FALSE)+
  geom_ribbon(data=gdata2,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#a1d99b"),show.legend=FALSE) +
  geom_ribbon(data=gdata2,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#31a354"),show.legend=FALSE) +
  
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("#f0f0f0"),show.legend=FALSE)+
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#bdbdbd"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#636363"),show.legend=FALSE) +
  
  #  
  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("#ffeda0"),show.legend=FALSE)+
  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#feb24c"),show.legend=FALSE) +
  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#f03b20"),show.legend=FALSE) +
  #
  #  

  geom_point(data=gobs,aes(x=obsx,y=obsy,group=set,color=set,shape=set),size=2.75)+
  geom_errorbar(show.legend=FALSE,data=gobs,aes(x=obsx,y=obsy,ymin=obsy-erry,ymax=obsy+erry,group=set,color=set),
                width=0.01,alpha=0.4)+
#  geom_line(data=gdata,aes(x=xx,y=mean),linetype="dashed",size=0.75,show.legend=FALSE)+
#  geom_line(data=gdata2,aes(x=xx,y=mean),linetype="dashed",size=0.75,show.legend=FALSE)+
  scale_colour_manual(name="",values=c("#636363","#31a354","#31a354","#636363","#636363",
                                            "#636363","#636363"))+
  scale_shape_manual(values=c(0,1,2,5,25,11,13),name="")+
  coord_cartesian(xlim=c(5e-3,0.85),ylim=c(0,20)) +
  theme_bw() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + 
  scale_x_log10()  +
  annotation_logticks(sides = "b") +
  annotation_logticks(base=2.875,sides = "l") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.675),
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size=14,colour = set),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=17.5),
        axis.text  = element_text(size=13),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) 
dev.off()




S <- ggs(as.mcmc(Normfit)[,c("e1", "gin", "gout")])

S$Parameter <- revalue(S$Parameter, c("e1"="E[r]", "gin"="Gamma['d']",
                                      "gout"="Gamma['p']"))



pdf("plot/He3dp_posterior_syst.pdf",height = 7,width = 8)
ggplot(data=S,aes(x=value,group=Parameter,fill=Parameter)) +
  geom_histogram(aes(y=..count../sum(..count..)),bins = 10,color="black") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  geom_rug(sides = "b", aes(y = 0,colour = Parameter),size=0.1,alpha=0.1 )+
  scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb"))+
  scale_color_manual(values = c("#66c2a5","#fc8d62","#8da0cb"))+
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=15),
        strip.background = element_rect("#F0B27A"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Parameter,scales="free_x",ncol=2,nrow=2,labeller=label_parsed) +
  ylab("Posterior probability") + xlab("Parameter")
dev.off()






Sa <- ggs(as.mcmc(Normfit),family="scale")

Sa$Parameter <- revalue(Sa$Parameter, 
c("scale[1]" = "Ali01a","scale[2]" = "Ali01b","scale[3]" = "Cos00",
"scale[4]" = "Gei99","scale[5]" = "Kra87","scale[6]" = "Mol80",
"scale[7]" = "Zhi77"))

pdf("plot/He3dp_scale_syst.pdf",height = 5,width = 4.5)
ggs_caterpillar(Sa) + aes(size=5,color = Parameter,shape=Parameter) +
  theme_economist_white() +
  geom_vline(xintercept = 1,linetype="dashed",color="gray35") +
  scale_shape_manual(values=c(0,1,2,5,25,11,13),name="")+
  scale_color_manual(values=c(rep("gray25",7))) +
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=10),
        strip.background = element_rect("white")) +
    ylab("Data set") +
#  xlab("Highest Probability Interval")
  xlab("Systematic effects")
dev.off()




pdf("plot/He3dp_trace_syst.pdf",height = 7,width = 8)
ggplot(data=S,aes(x= Iteration,y=value,group=Parameter,color=factor(Chain))) +
  geom_line(alpha=0.5,size=0.25) +
  theme_wsj() +
  scale_color_economist()+
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=15),
        strip.background = element_rect("white")) +
  facet_wrap(~Parameter,scales="free",ncol=2,nrow=2,labeller=label_parsed) +
  ylab("Parameter value") + xlab("Iteration")+
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He")))
dev.off()
#


pdf("plot/He3dp_trace_scale_syst.pdf",height = 7,width = 8)
ggplot(data=Sa,aes(x= Iteration,y=value,group=Parameter,color=factor(Chain))) +
  geom_line(alpha=0.5,size=0.25) +
  theme_wsj() +
  scale_color_futurama()+
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=15),
        strip.background = element_rect("white")) +
  facet_wrap(~Parameter,scales="free",ncol=2,nrow=2,labeller=label_parsed) +
  ylab("Parameter value") + xlab("Iteration")+
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He")))
dev.off()


# Reaction rates

Nsamp <- 500
mcdat <- as.data.frame(do.call(rbind, as.mcmc(Normfit)[,c("e1","gin","gout")]))

index <- sample(1:nrow(mcdat),size=Nsamp,replace=FALSE)
mcdat <- mcdat[index,]
Tgrid <- exp(seq(log(1e-3),log(10),length.out =  100))


gdat <- list()
for(i in 1:Nsamp){
  y <- sapply(Tgrid,nuclear_rate3Hedp,ER = mcdat[i,1],gi = mcdat[i,2],gf = mcdat[i,3])
  dd <- data.frame(y)
  gdat[[i]] <- dd
}

gg <-  as.data.frame(gdat)
gg$x <- Tgrid

gg2<-apply(gg, 1, quantile, probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995), na.rm=TRUE)


gg2data <- data.frame(x =Tgrid, mean = gg2["50%",],lwr1=gg2["25%",],
                      lwr2 = gg2["2.5%",],lwr3=gg2["0.5%",],upr1=gg2["75%",],
                      upr2=gg2["97.5%",],upr3=gg2["99.5%",])

write.csv(gg2data,"NV.csv",row.names = F)

 


pdf("plot/He3dp_cross.pdf",height = 7,width = 10)
ggplot(gg2data,aes(x=x,y=mean))+
  theme_bw()  +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr3, ymax=upr3,y=NULL), fill=c("#ffeda0"),show.legend=FALSE) +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr2, ymax=upr2,y=NULL), fill=c("#feb24c"),show.legend=FALSE) +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr1, ymax=upr1,y=NULL), fill=c("#f03b20"),show.legend=FALSE) +
   xlab("Temperature (GK)") + ylab(expression(N[A]~symbol("\341")*sigma*nu*symbol("\361"))) +
 scale_y_continuous(breaks=c(0,5e7,1e8,1.5e8),labels=c(0,expression(5%*%10^7),
                              expression(10^8),expression(1.5%*%10^8))) +
  theme(legend.position = "none",
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "gray95"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=17.5),
        axis.text  = element_text(size=13),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
dev.off()




cmb <- as.data.frame(filter(gg, x <= 1.1 & x >= 0.9))

cmbhist <- data.frame(x=as.numeric(cmb)[1:Nsamp])
ggplot(cmbhist, aes(x)) +
  geom_histogram(aes(y=..count../sum(..count..)),bins = 10,fill="#f03b20",color="#feb24c") +
  theme_bw() +
#  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_x_continuous(breaks=c(1.15e8,1.2e8,1.25e8),labels=c(expression(1.15%*%10^8),
                              expression(1.20%*%10^8),expression(1.25%*%10^8))) +
  geom_rug(sides = "b", aes(y = 0),colour = "#feb24c") +
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=15),
        strip.background = element_rect("#F0B27A"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("Posterior probability") + xlab(expression(N[A]~symbol("\341")*sigma*nu*symbol("\361")))



dens <- density(cmbhist)
df <- data.frame(x=dens$x, y=dens$y)
probs=c(0.025, 0.25, 0.5, 0.75, 0.975)
quantiles <- quantile(df$x, prob=probs)
df$quant <- factor(findInterval(df$x,quantiles))

write.csv(df,"df.csv",row.names = F)

ggplot(df, aes(x,y)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
  scale_fill_manual(values=c("#deebf7","#9ecae1","#3182bd","#3182bd","#9ecae1","#deebf7")) +
  geom_line() + theme_wsj() + xlab(expression(N[A]~sigma*v)) + ylab("Density") +
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He")))




cmbquant <- apply(cmb, 1, quantile, probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995), na.rm=TRUE)
quant <- cut(cmbhist,breaks=as.numeric(cmbquant))
gcmb <- data.frame(cmbhist,quant)



