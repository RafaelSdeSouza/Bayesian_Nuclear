# 3Hedp analysis
#
# purpose: Real  DATA
#
# - 3 parameters are assumed: Er, gamma_d^2, gamma_n^2 [e1, gin, gout]
#
# - uses the function sfactorTdn_fast(obsx1[i], e1, gin, gout), which
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
require(nuclear);library(magrittr);require(runjags)
library(dplyr);require(ggsci);require(ggmcmc);require(plyr);require(latex2exp)
source("jagsresults.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")
load.module("nuclear")
load.runjagsmodule()



######################################################################
## Read DATA GENERATION
ensamble <- read.csv("ensamble_final.csv",header = T,stringsAsFactors=FALSE)  %>%
  mutate(Stat=replace(Stat,Stat==0,0.1)) %>%
  mutate(dat=replace(dat,dat %in% c("gei99b","gei99d"),"Gei99")) %>%
  mutate(dat=replace(dat,dat %in% c("Kra87m","Kra87b"),"Kra87")) %>%
  mutate(dat=replace(dat,dat == "zhi77b","Zhi77")) %>%
  filter(.,dat!="Lac05")  %>% droplevels(.) %>%
  mutate(dat=as.factor(dat))


re <- as.numeric(ensamble$dat)
Nre <- length(unique(ensamble$dat))

# Radius
# ri = 6
# rf = 5

# Literature
#  0.35779   # resonance energy
#  1.0085    # reduced width incoming
#  0.025425   # reduced width outgoing


N <- nrow(ensamble)
obsy <- ensamble$S    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat
set <- ensamble$dat
syst = unique(ensamble$Syst)
syst <- syst[-3]

# Literature radii
ri0 <- 6
rf0 <- 5

M <- 500
xx <- seq(min(obsx),max(obsx),length.out = M)

model.data <- list(obsy = obsy,    # Response variable
                   obsx =  obsx,   # Predictors
                   erry = erry,
                   N = nrow(ensamble), # Sample size
                   syst = syst,
                   Nre = Nre,
                   re = re,
                   M = M,
                   xx = xx
)


######################################################################
Model <- "model{
# LIKELIHOOD
for (i in 1:N) {
obsy[i] ~ dnorm(y[i], pow(erry[i], -2))
y[i] ~ dnorm(scale[re[i]]*sfactor3Hedp_5p(obsx[i], e1, gin, gout, ri, rf),pow(tau, -2))
}

# Predicted values

for (j in 1:M){
mux[j] <- sfactor3Hedp_5p(xx[j], e1, gin, gout,ri, rf)
yx[j] ~ dnorm(mux[j],pow(tau,-2))
}



for (k in 1:Nre){
scale[k] ~ dlnorm(log(1.0),pow(log(1+syst[k]),-2))
}


# PRIORS
# e1, gin, gout are defined as in tdn.f (by Alain Coc):
# resonance energy, initial reduced width, final reduced
# width;


#gout ~ dnorm(0.05,1)T(0,)
# gin ~  dnorm(6,1)T(0,)
#gin ~  dunif(1,10)

tau ~  dunif(0.01,10)
e1 ~   dunif(0,10)

#gout ~ dhalfcauchy(1)

#gin ~ dhalfcauchy(5)


gout ~ dbeta(2,2)
gb ~ dbeta(2,2)
gin <- 3*gb + 5

 
#gin ~ dnorm(6,0.25)T(0,)

#gin ~  dunif(0,20)



# Channel radius
 rb ~ dbeta(2,2)
 rb2 ~ dbeta(2,2)
 
ri <- 6*rb + 2
rf <- 6*rb + 2

#ri <- 3*rb + 3
#rf <- rb2 + 4.5

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


#inits <- function () { list(e1 = 0.38,gin= 6,gout=0.1,ri=4,rf=5) }

inits <- function () { list(e1 = 0.38) }

# "f": is the model specification from above;
# data = list(...): define all data elements that are referenced in the



# JAGS model with R2Jags;
Normfit <- jags(data = model.data,
                inits = inits,
                parameters.to.save  = c("e1", "gin", "gout","ri","rf","tau","mux","scale"),
                model.file  = textConnection(Model),
                n.thin = 10,
                n.chains = 4,
                n.burnin = 5000,
                n.iter = 10000)

jagsresults(x = Normfit , params = c("e1", "gin", "gout","ri","rf","tau"))

   

###################################  Nice plots with ggmcmc package

L.radon.intercepts <- data.frame(
  Parameter=paste("scale[", seq(1:4), "]", sep=""),
  Label=levels(ensamble$dat))
head(L.radon.intercepts)
S <- ggs(as.mcmc(Normfit),par_labels=L.radon.intercepts, family = "scale")
ggs_caterpillar(S) + theme_wsj() +
  theme(legend.position = "top",
  legend.background = element_rect(colour = "white", fill = "white"),
   plot.background = element_rect(colour = "white", fill = "white"),
   panel.background = element_rect(colour = "white", fill = "white"),
   legend.key = element_rect(colour = "white", fill = "white"),
   axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=15),
   axis.text  = element_text(size=12),
   strip.text = element_text(size=10),
  strip.background = element_rect("gray85")) + geom_vline(xintercept = 1,linetype="dashed",color="red",alpha=0.7) +
 ylab("Dataset")
####################################################


traplot(Normfit  ,c("e1", "gin", "gout","ri","rf"),style="plain")
denplot(Normfit  ,c("e1", "gin", "gout","ri","rf"),style="plain")
caterplot(Normfit,c("scale"),style="plain")


# Plot of predicted values
y <- jagsresults(x=Normfit , params=c('mux'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))
x <- xx
gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],lwr3=y[,"0.5%"],upr1=y[,"75%"],
                    upr2=y[,"97.5%"],upr3=y[,"99.5%"])
gobs <- data.frame(obsx,obsy,erry,set)
gobs$set <- as.factor(gobs$set)

pdf("plot/He3dp_syst_5p.pdf",height = 7,width = 8)
ggplot(gobs,aes(x=obsx,y=obsy))+
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),alpha=0.7,fill=c("#deebf7"),show.legend=FALSE)+
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),alpha=0.6,  fill = c("#9ecae1"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),alpha=0.4,fill=c("#3182bd"),show.legend=FALSE) +
  geom_point(data=gobs,aes(x=obsx,y=obsy,group=set,color=set,shape=set),size=2)+
  geom_errorbar(data=gobs,aes(x=obsx,y=obsy,ymin=obsy-erry,ymax=obsy+erry,group=set,color=set),width=0.025)+
  geom_line(data=gdata,aes(x=xx,y=mean),colour="white",linetype="dashed",size=1,show.legend=FALSE)+
  scale_colour_futurama(name="Dataset")+
  scale_shape_stata(name="Dataset")+
  theme_wsj() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + scale_x_log10()  +
  theme(legend.position = "top",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He")))
dev.off()


S <- ggs(as.mcmc(Normfit)[,c("e1", "gin", "gout","tau","ri","rf")])

S$Parameter <- revalue(S$Parameter, c("e1"="E[r]", "gin"="Gamma['in']",
                                      "gout"="Gamma['out']","tau"="tau","ri"="r[i]","rf"="r[f]"))



pdf("plot/He3dp_posterior_syst_5p.pdf",height = 7,width = 8)
ggplot(data=S,aes(x=value,group=Parameter,fill=Parameter)) +
  geom_density(alpha=0.75) + 
  theme_wsj() +
  scale_fill_futurama()+
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=15),
        strip.background = element_rect("white")) + 
  facet_wrap(~Parameter,scales="free",ncol=3,nrow=2,labeller=label_parsed) +
  ylab("Posterior probability") + xlab("Parameter value")+
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He"))) 
dev.off()


pdf("plot/He3dp_trace_syst_5.pdf",height = 7,width = 8)
ggplot(data=S,aes(x= Iteration,y=value,group=Parameter,color=factor(Chain))) +
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
  facet_wrap(~Parameter,scales="free",ncol=3,nrow=2,labeller=label_parsed) +
  ylab("Parameter value") + xlab("Iteration")+
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He"))) 
dev.off()
#


Sa <- ggs(as.mcmc(Normfit),family="scale")

Sa$Parameter <- revalue(Sa$Parameter, c("scale[1]" = "Gei99","scale[2]" = "Kra87","scale[3]" = "Mol80","scale[4]" = "Zhi77"))

pdf("plot/He3dp_scale_syst_5p.pdf",height = 7,width = 4)
ggs_caterpillar(Sa) + aes(color=Parameter) +
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
  #  facet_wrap(~Parameter,scales="free",ncol=2,nrow=2,labeller=label_parsed) +
  #  ylab("Parameter value") + 
  xlab("Highest Probability Density")+
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He"))) +
  geom_vline(xintercept = 1,linetype="dashed",color="gray")
dev.off()





# Reaction rates

Nsamp <- 500
mcdat <- as.data.frame(do.call(rbind, as.mcmc(Normfit)[,c("e1","gin","gout","ri","rf")]))

index <- sample(1:nrow(mcdat),size=Nsamp,replace=FALSE)
mcdat <- mcdat[index,]
Tgrid <- exp(seq(log(1e-3),log(10),length.out =  100))


gdat <- list()
for(i in 1:Nsamp){
  y <- sapply(Tgrid,nuclear_rate3Hedp_5p,ER = mcdat[i,1],gi = mcdat[i,2],gf = mcdat[i,3],r_i= mcdat[i,4],r_f=mcdat[i,5])
  dd <- data.frame(y)
  gdat[[i]] <- dd
}

gg <-  as.data.frame(gdat)
gg$x <- Tgrid

gg2<-apply(gg, 1, quantile, probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995), na.rm=TRUE)


gg2data <- data.frame(x =Tgrid, mean = gg2["50%",],lwr1=gg2["25%",],
                      lwr2 = gg2["2.5%",],lwr3=gg2["0.5%",],upr1=gg2["75%",],
                      upr2=gg2["97.5%",],upr3=gg2["99.5%",])

write.csv(gg2data,"NV_5.csv",row.names = F)

g1 <- ggplot(gg2data,aes(x=x,y=mean))+
  theme_bw()  +
  
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr3, ymax=upr3,y=NULL), fill=c("#deebf7"),show.legend=FALSE) +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr2, ymax=upr2,y=NULL), fill=c("#9ecae1"),show.legend=FALSE) +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr1, ymax=upr1,y=NULL), fill=c("#3182bd"),show.legend=FALSE) +
  geom_line(size=1,colour="#ffffff",linetype="dashed",size=1,show.legend=FALSE) +
  theme_wsj() + xlab("Temperature (GK)") + ylab(TeX('$N_A\\sigma v$')) +
  #  scale_y_log10() + scale_x_log10()+
  theme(legend.position = "none",
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He")))


g1


cmb <- as.data.frame(filter(gg, x <= 1.1 & x >= 0.99))

cmbhist <- as.numeric(cmb)[1:Nsamp]
dens <- density(cmbhist)
df <- data.frame(x=dens$x, y=dens$y)
probs=c(0.025, 0.25, 0.5, 0.75, 0.975)
quantiles <- quantile(df$x, prob=probs)
df$quant <- factor(findInterval(df$x,quantiles))

write.csv(df,"df_5W.csv",row.names = F)

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