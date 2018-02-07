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
source("..//..//auxiliar_functions/jagsresults.R")
source("..//..//auxiliar_functions/theme_rafa.R")
source("..//..//auxiliar_functions/pair_wise_plot.R")
source("..//..//auxiliar_functions/Gamma3Hedp.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")
load.module("nuclear")


######################################################################
## Read DATA 
ensamble <- read.csv("ensamble_Tdn.csv",header = T)


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
syst = c(unique(ensamble$Syst))
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
#y[i] ~ dnorm(scale[re[i]]*sfactorTdn3(obsx[i], e1, ex,gin, gout,ri,rf,ue[ik[i]]),pow(tau, -2))
y[i] ~ dnorm(scale[re[i]]*sfactorTdn3(obsx[i], e1,ex, gin, gout,ri,rf,0),pow(tau, -2))
res[i] <- obsy[i]-sfactorTdn3(obsx[i], e1,ex, gin, gout,ri,rf,0)
}

RSS <- sum(res^2)

# Predicted values

for (j in 1:M){

# Bare...

mux0[j] <- sfactorTdn3(xx[j], e1,ex, gin, gout,ri,rf,0)
yx0[j] ~ dnorm(mux0[j],pow(tau,-2))

# No inverse Kinematics 

#mux1[j] <- sfactorTdn(xx[j], e1, gin, gout,ri,rf,ue[1])
#yx1[j] ~ dnorm(mux1[j],pow(tau,-2))

# With inverse Kinematics 
#mux2[j] <- sfactorTdn(xx[j], e1, gin, gout,ri,rf,ue[2])
#yx2[j] ~ dnorm(mux1[j],pow(tau,-2))

}



for (k in 1:Nre){
scale[k] ~ dlnorm(log(1.0),1/log(1+pow(syst[k],2)))
}

#for (z in 1:Nik){
#ue[z] ~ dnorm(0,1e3)T(0,)
#}


# PRIORS
# Wigner limit: wl = hbar^2/(m_red a_c^2) = 41.80159/(M_red a_c^2)
#
# deuteron channel: wl_d = 41.80159/(1.207357 a_c^2) = 34.6224/a_c^2
# neutron channel:  wl_n = 41.80159/(0.805597 a_c^2) = 51.8889/a_c^2
# e1, gin, gout are defined as in tdn.f (by Alain Coc):
# resonance energy, initial reduced width, final reduced
# width;

tau ~  dt(0, pow(5,-2), 1)T(0,)
e1 ~  dnorm(0,0.1)T(0,)
ex ~  dnorm(0,0.1)T(0,)

rf ~  dnorm(5,pow(0.1,-2))T(0,)
ri ~  dnorm(6,pow(0.1,-2))T(0,)


gin ~  dnorm(0,pow(0.5,-2))T(0,)
gout ~ dnorm(0,pow(0.5,-2))T(0,)




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


inits <- function () { list(e1 = runif(1,0.01,1),gout=runif(1,0.01,10),gin=runif(1,0.01,10)) }
# "f": is the model specification from above;
# data = list(...): define all data elements that are referenced in the



# JAGS model with R2Jags;
Normfit <- jags(data = model.data,
                inits = inits,
                parameters.to.save  = c("e1","ex", "gin", "gout","ue","tau", "ri","rf","RSS","mux0","mux1","mux2","scale"),
                model.file  = textConnection(Model),
                n.thin = 4,
                n.chains = 10,
                n.burnin = 15000,
                n.iter = 30000)


jagsresults(x = Normfit , params = c("e1","ex", "gin", "gout","scale","tau","ri","rf"),probs = c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))

RSS <- as.matrix(as.mcmc(Normfit)[,c("RSS")])
rss0 <- function(x) crossprod(x-mean(x))[1]
1-mean(RSS)/rss0(obsy)


color_scheme_set("darkgray")
div_style <- parcoord_style_np(div_color = "green", div_size = 0.05, div_alpha = 0.4)
mcmc_parcoord(as.mcmc(Normfit),alpha = 0.05, regex_pars = c("e1", "gin", "gout","ri","rf"))


traplot(Normfit  ,c("e1","ex", "gin", "gout"),style="plain")
denplot(Normfit  ,c("e1", "gin", "gout","ue"),style="plain")
caterplot(Normfit,c("scale","tau"),style="plain")



# Plot

y <- jagsresults(x=Normfit , params=c('mux1'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))
x <- xx
gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],lwr3=y[,"0.5%"],upr1=y[,"75%"],
                    upr2=y[,"97.5%"],upr3=y[,"99.5%"])
gobs <- data.frame(obsx,obsy,erry,set)
gobs$set <- as.factor(gobs$set)



y2 <- jagsresults(x=Normfit , params=c('mux2'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))
gdata2 <- data.frame(x =xx, mean = y2[,"mean"],lwr1=y2[,"25%"],lwr2=y2[,"2.5%"],lwr3=y2[,"0.5%"],upr1=y2[,"75%"],
                    upr2=y2[,"97.5%"],upr3=y2[,"99.5%"])


y0 <- jagsresults(x=Normfit , params=c('mux0'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))

gdata0 <- data.frame(x =xx, mean = y0[,"mean"],lwr1=y0[,"25%"],lwr2=y0[,"2.5%"],lwr3=y0[,"0.5%"],upr1=y0[,"75%"],
                    upr2=y0[,"97.5%"],upr3=y0[,"99.5%"])


pdf("plot/He3dp_syst.pdf",height = 7,width = 10)
ggplot(gobs,aes(x=obsx,y=obsy))+
  geom_rect(aes(xmin=0.004, xmax=0.12, ymin=-2, ymax=32), fill="gray90",alpha=0.4) +

  # red
#  geom_ribbon(data=gdata2,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("#fdd0a2"),show.legend=FALSE)+
#  geom_ribbon(data=gdata2,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#fdae6b"),show.legend=FALSE) +
#  geom_ribbon(data=gdata2,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#d94801"),show.legend=FALSE) +
  
#  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("gray95"),show.legend=FALSE)+
#  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("gray75"),show.legend=FALSE) +
#  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("gray55"),show.legend=FALSE) +
  
  #  Bare
  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("#dadaeb"),show.legend=FALSE)+
  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#9e9ac8"),show.legend=FALSE) +
  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#984ea3"),show.legend=FALSE) +
  #
  #  

  geom_point(data=gobs,aes(x=obsx,y=obsy,group=set,color=set,shape=set),size=2.75,alpha=0.75)+
  geom_errorbar(show.legend=FALSE,data=gobs,aes(x=obsx,y=obsy,ymin=obsy-erry,ymax=obsy+erry,group=set,color=set),
                width=0.01,alpha=0.4)+
#  geom_line(data=gdata,aes(x=xx,y=mean),linetype="dashed",size=0.75,show.legend=FALSE)+
#  geom_line(data=gdata2,aes(x=xx,y=mean),linetype="dashed",size=0.75,show.legend=FALSE)+
#  scale_colour_manual(name="",values=c("gray20","#7f0000","#7f0000","gray20","gray20",
#                                       "gray20","gray20"))+
  scale_color_stata(name="")+
  scale_shape_manual(values=c(0,19,8,10,4,17,3),name="")+
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
        axis.title = element_text(size=18.5),
        axis.text  = element_text(size=13),
        axis.ticks = element_line(size = 0.75),
        axis.line = element_line(size = 0.5, linetype = "solid")) 
dev.off()




dummy <- as.mcmc(Normfit)[,c("e1", "gin", "gout","ri","rf")]

as.matrix(dummy)

#dummy  <- as.matrix(dummy)

for(i in 1:3){
dummy[[i]][,2] <- Gamma3Hedp(dummy[[i]][,1],dummy[[i]][,2],dummy[[i]][,3],dummy[[i]][,4],dummy[[i]][,5])$Ga
}
for(i in 1:3){
  dummy[[i]][,3] <- Gamma3Hedp(dummy[[i]][,1],dummy[[i]][,2],dummy[[i]][,3],dummy[[i]][,4],dummy[[i]][,5])$Gb
}


traplot(dummy)

ss <- ggs(dummy)
pair_wise_plot(ggs(dummy))

ggs_traceplot(ss)

Sp <- ggs(as.mcmc(Normfit)[,c("e1", "gin", "gout","ri","rf")]) 

DD <- as.matrix(as.mcmc(Normfit)[,c("e1", "gin", "gout","ri","rf")])




Sp0 <- Sp %>% as_tibble()  %>% mutate(value = ifelse(Parameter == 'e1', 10*value, value)) 







levels(Sp0$Parameter) <- as.factor(c("E[r]","Gamma[d]", "Gamma[p]","a[c]^i","a[c]^f"))

pdf("plot/He3dp_corr.pdf",height = 8,width =8)
pair_wise_plot(Sp0)
dev.off()





S <- ggs(as.mcmc(Normfit),family="ue")
S$value <- 1e6*S$value
levels(S$Parameter) <- c("U[e1]","U[e2]")

pdf("plot/He3dp_posterior_screen.pdf",height = 7,width = 8)
ggplot(data=S,aes(x=value,group=Parameter,fill=Parameter)) +
  geom_histogram(aes(y=..count../sum(..count..)),bins = 10,color="black") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  geom_rug(sides = "b", aes(y = 0,colour = Parameter),size=0.1,alpha=0.1 )+
  scale_fill_manual(values = c("gray70","#d7301f"))+
  scale_color_manual(values = c("gray70","#d7301f"))+
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=15),
        strip.background = element_rect("gray80"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Parameter,scales="free_x",ncol=2,nrow=2,labeller=label_parsed) +
  ylab("Posterior probability") + xlab("eV")
dev.off()






Sa <- ggs(as.mcmc(Normfit),family="scale")

Sa$Parameter <- revalue(Sa$Parameter, 
c("scale[1]" = "Ali01a","scale[2]" = "Ali01b","scale[3]" = "Cos00",
"scale[4]" = "Gei99","scale[5]" = "Kra87","scale[6]" = "Mol80",
"scale[7]" = "Zhi77"))

pdf("plot/He3dp_scale_syst.pdf",height = 5,width = 4.5)
ggs_caterpillar(Sa) + aes(size=0.1,color = Parameter,shape=Parameter) +
  theme_economist_white() +
  geom_vline(xintercept = 1,linetype="dashed",color="gray35") +
  scale_shape_manual(values=c(0,19,8,10,4,17,3),name="")+
  scale_color_manual(values=c(rep("gray55",7))) +
  geom_point(size=3,color="red") +
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



SS <- ggs(as.mcmc(Normfit)[,c("e1", "gin", "gout","ri","rf")])

ss$Chain <- as.factor(ss$Chain)

ggs_density(ss)
pdf("plot/He3dp_trace_syst.pdf",height = 7,width = 8)
ggplot(data=ss,aes(x= Iteration,y=value,group=Parameter,color=factor(Chain))) +
  geom_density(alpha=0.5,size=0.25) +
  theme_wsj() +
#  scale_color_economist()+
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=15),
        strip.background = element_rect("white")) +
  facet_wrap(~Parameter,scales="free",labeller=label_parsed) +
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



