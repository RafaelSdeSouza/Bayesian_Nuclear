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
#rm(list=ls())
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
require(MCMCvis);require(ggridges)
source("..//..//auxiliar_functions/jagsresults.R")
source("..//..//auxiliar_functions/theme_rafa.R")
source("..//..//auxiliar_functions/pair_wise_plot.R")
source("..//..//auxiliar_functions/Gamma3Hedp.R")
source("plot_Sfactor.R")
source("plot_normfactors.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")
load.module("nuclear")


######################################################################
## Read DATA
ensamble <- read.csv("ensamble.csv",header = T) %>%
  mutate(Syst=replace(Syst,Syst==0.06,0.078))  %>% filter(E <= 0.5)


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
                   obsy2 = obsy,    # Response variable
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
 #                  ap  = 5,
 #                  ad = 6

)


# Conservative case
######################################################################
Model <- "model{
# LIKELIHOOD informative
for (i in 1:N) {
obsy[i] ~ dnorm(y[i], pow(erry[i], -2))
y[i] ~ dnorm(scale[re[i]]*sfactor3Hedp(obsx[i], E0, Er, gd2, gp2, ad, ap, ue[ik[i]]), pow(tau, -2))
res[i] <- obsy[i]-sfactor3Hedp(obsx[i], E0, Er, gd2, gp2, ad, ap,0)
nres[i] <- res[i]/obsy[i]
}


# LIKELIHOOD broad
for (i in 1:N) {
obsy2[i] ~ dnorm(y_2[i], pow(erry[i], -2))
y_2[i] ~ dnorm(scale[re[i]]*sfactor3Hedp(obsx[i],  E0_b, Er_b, gd2_b, gp2_b, ad_b, ap_b, ue[ik[i]]),pow(tau_2, -2))
}
RSS <- sum(res^2)

# Predicted values
for (j in 1:M){

# Bare...

mux0[j] <- sfactor3Hedp(xx[j], E0, Er, gd2, gp2, ad, ap,0)


mux0_2[j] <- sfactor3Hedp(xx[j], E0_b, Er_b, gd2_b, gp2_b, ad_b, ap_b,0)

DeltaM[j] <- mux0[j]/mux0_2[j]



# No inverse Kinematics

mux1[j] <- sfactor3Hedp(xx[j], E0, Er, gd2, gp2, ad, ap,ue[1])
yx1[j] ~ dnorm(mux1[j],pow(tau,-2))

# With inverse Kinematics
mux2[j] <- sfactor3Hedp(xx[j], E0, Er, gd2, gp2, ad, ap,ue[2])
yx2[j] ~ dnorm(mux1[j],pow(tau,-2))

}

for (k in 1:Nre){
scale[k] ~ dlnorm(log(1.0),1/log(1+pow(syst[k],2)))
}

for (z in 1:Nik){
ue[z] ~ dnorm(0,pow(0.1,-2))T(0,)
}


# PRIORS

# Case I
tau ~  dnorm(0, pow(1,-2))T(0,)
E0  ~  dnorm(0, pow(1,-2))T(0,)
Er <-  E0
gd2 ~  dnorm(0, pow(3,-2))T(0,)
gp2 ~  dnorm(0, pow(3,-2))T(0,)

ad  ~  dnorm(6, pow(0.01,-2))T(0,)
ap  ~  dnorm(5, pow(0.01,-2))T(0,)



# Case II
tau_2  ~  dnorm(0, pow(1,-2))T(0,)
Er_b  ~   dbeta(0.5,0.5)
E0_b  ~   dbeta(2,5)

gp2_b ~  dnorm(0, pow(3,-2))T(0,)
gd2_b  ~ dnorm(0, pow(3,-2))T(0,)

ad_b  ~  dunif(2,10)
ap_b  ~  dunif(2,10)


ue_ev[1] <-1e6*ue[1]
ue_ev[2] <-1e6*ue[2]

S_0   <- sfactor3Hedp(1e-4, E0, Er, gd2, gp2, ad, ap,0)
S_0b  <- sfactor3Hedp(1e-4, E0_b, Er_b, gd2_b, gp2_b, ad_b, ap_b,0)


}"

inits <- function(){ list(E0 = runif(1,0.3,0.35),E0_b = 0.4,Er_b = 0.4,gd2 = 1,
                        gp2 = runif(1,0.01,0.1),gd2_b = 0.5) }


set.seed(24)
# JAGS model with R2Jags;
Normfit <- jags(data = model.data,
                inits = inits,
                parameters.to.save  = c("Er","E0","gd2", "gp2","ue_ev","tau", "ad","ap",
                                        "RSS","mux0","mux1","mux2","scale","DeltaM","S_0",
                                        "S_0b","E0_b","Er_b","gd2_b",
                                        "gp2_b","tau_2","ad_b","ap_b","res","nres"),
                model.file  = textConnection(Model),
                n.thin = 30,
                n.chains = 5,
<<<<<<< HEAD
                n.burnin = 7500,
                n.iter = 15000)
jagsresults(x = Normfit, params = c("E0","gd2", "gp2","ue","tau", "ad","ap","ue_ev","S_0"),probs = c(0.16, 0.5, 0.84))
temp <- Normfit
temp <- update(temp, n.thin = 50, n.iter = 50000)

getmcmc_var <- function(outjags=outjags,vars = vars){
as.data.frame(do.call(rbind, as.mcmc(outjags)[,vars]))
  }


sum(res[,"mean"]^2)
res <- jagsresults(x = Normfit, params = c("nres"),probs = c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))
res_data <- data.frame(x=obsx,sd=res[,"sd"], mean=res[,"50%"],lwr1=res[,"16%"],lwr2=res[,"2.5%"],lwr3=res[,"0.15%"],upr1=res[,"84%"],
                       upr2=res[,"97.5%"],upr3=res[,"99.85%"],set,lab)

ggplot(res_data,aes(x=obsx,y=mean,roup=set,color=set,shape=set)) + 
  geom_point(size=2.75,alpha=0.75) +
  geom_errorbar(show.legend=FALSE,aes(x=obsx,y=mean,ymin=lwr2,ymax=upr2),
                width=0.01,alpha=0.4)+
  coord_cartesian(xlim=c(0.03,0.8),ylim=c(-0.5,0.5)) +
  scale_colour_stata(name="") +
  scale_shape_manual(values=c(0,19,8,10,4,17,3),name="") + 
  theme_bw() + xlab("Energy (MeV)") + ylab("Residuals") +
  scale_x_log10(breaks = c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))  +
  annotation_logticks(short = unit(0.2, "cm"), mid = unit(0.25, "cm"), long = unit(0.3, "cm"),
                      sides = "b",size = 0.45) +
  #  annotation_logticks(base=2.875,
  #  short = unit(0.2, "cm"), mid = unit(0.25, "cm"), long = unit(0.3, "cm"),sides = "l",size = 0.45) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.925,0.65),
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size=14,colour = set),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size=22),
        axis.text  = element_text(size=18),
        axis.ticks = element_line(size = 0.45),
        #        axis.line = element_line(size = 0.45, linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3, "mm")) +
  geom_smooth(method = "lm")



dtc <- getmcmc_var(Normfit,c("E0","gd2","gp2","ad_b","ap_b","S_0","ue_ev[1]",
                             "ue_ev[2]"))
=======
                n.burnin = 5000,
                n.iter = 10000)


tab <- jagsresults(x = Normfit, params = c("E0","gd2", "gp2","ue","tau", "ad","ap","ue_ev","S_0"),probs = c(0.16, 0.5, 0.84))
tab <- as.data.frame(tab)
>>>>>>> e782103f99d3ed2da11c3ed8bed8204ba2d00510

tab$low <- tab[,4] - tab[,3]
tab$hi <-  tab[,5] - tab[,4]

jagsresults(x = Normfit, params = c("E0","gd2", "gp2","ue","tau", "ad","ap","ue_ev","S_0"),probs=c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))

tab2 <- jagsresults(x = Normfit , params = c("E0_b","Er_b","gd2_b", "gp2_b","tau_2","ad_b","ap_b","ue_ev","S_0b"),probs = c(0.16, 0.5, 0.84))
tab2 <- as.data.frame(tab2)
tab2$low <- tab2[,4] - tab2[,3]
tab2$hi <-  tab2[,5] - tab2[,4]


traplot(Normfit  ,c("E0","gd2", "gp2","ue"),style="plain")
denplot(Normfit  ,c("E0","gd2", "gp2","ue"),style="plain")


#temp <- Normfit
#temp <- update(temp, n.thin = 50, n.iter = 50000)

getmcmc_var <- function(outjags=outjags,vars = vars){
as.data.frame(do.call(rbind, as.mcmc(outjags)[,vars]))
  }
require(HDInterval)
dtc <- getmcmc_var(Normfit,c("E0","gd2","gp2","ad_b","ap_b","S_0","ue_ev[1]",
                             "ue_ev[2]"))





pdf("plot/ROPE_E0.pdf",height = 0.75*7,width = 0.75*10)
plotAreaInROPE(dtc$E0/0.35779, compVal = 1,maxROPEradius = 0.1)
dev.off()





E_cdf <- data.frame(E0=dtc$E0/0.35779)
dat <- with(density(E_cdf$E0,bw=0.00278,adjust = 1,n = 512), data.frame(x, y))
pdf("plot/NoROPE_E0.pdf",height = 0.75*7,width = 0.75*10)
ggplot(E_cdf,aes(x=E0)) +
  stat_density_ridges(aes(y=0,fill=factor(..quantile..),alpha=factor(..quantile..)),geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      quantiles = c(0.025, 0.16,  0.84, 0.975)) +
  scale_fill_manual(name = "Probability", values = c("#d9d9d9","#969696","#525252",
                                                     "#969696","#d9d9d9")) +
# geom_rect(aes( xmin=0.99,xmax=1.01,ymin=0,ymax=20),color="red",fill="white",alpha=0,linetype="dashed") +
#  geom_area(data = subset(E_cdf, E0 > 0.99 &  E0 < 1.01),fill="red")+
#  geom_density(fill="#21908CBF") +

  theme_bw() + coord_cartesian(xlim=c(0.925,1.01)) +
 
  xlab("Eigenenergy ratio") + 
  ylab("Probability density") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size=14,colour = set),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size=22),
        axis.text  = element_text(size=18),
        axis.ticks = element_line(size = 0.45),
        #        axis.line = element_line(size = 0.45, linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3, "mm"))  +
  annotate(geom = "text",label=paste("4.5% in ","ROPE"),x=1,y=20,size=7) +
  geom_errorbarh(aes(y=15, x=1, xmin=0.99, xmax=1.01),size=1,height=5,color="red2") +
  geom_segment(aes(x = 0.99, y = 0, xend = 0.99, yend = 12),linetype="dashed") +
  geom_segment(aes(x = 1.01, y = 0, xend = 1.01, yend = 12),linetype="dashed")
 # geom_ribbon(data=subset(dat, x  <= 1.01 & x >= 0.99),aes(x=x,ymax=y+0.5,ymin=0),fill = "red") 
                 


dev.off()
  

  



#hdi_jags <- function(mcmc=mcmc, par = par,credMass = 0.95,allowSplit=TRUE){
#  temp <- as.mcmc(mcmc)[,c(par)]
#  return(hdi(temp,allowSplit=TRUE))
#}
#hdi_jags(Normfit,par=c("Er","gd2", "gp2","ue_ev[1]","ue_ev[2]","S_0"),credMass = 0.95)



# Plot S-factor
pdf("plot/He3dp_syst_l.pdf",height = 7,width = 10)
plot_Sfactor(Normfit)
dev.off()

quantile((GammaHe3dp(Normfit)$Gd),probs = c(0.16, 0.5, 0.84))

quantile((GammaHe3dp(Normfit)$Gp),probs = c(0.16, 0.5, 0.84))

# Plot poir-wise plot Case-I


# Case I
Sp <- ggs(as.mcmc(Normfit)[,c("Er", "gd2", "gp2","ue_ev[1]","ue_ev[2]","S_0")])
Sp0 <- Sp %>% as_tibble()
Sp0$Parameter <- ordered(Sp0$Parameter, levels = c("Er", "gd2", "gp2","ue_ev[1]","ue_ev[2]","S_0"))
levels(Sp0$Parameter) <- as.factor(c("E[0]~(MeV)","gamma[d]^2~(MeV)", "gamma[p]^2~(MeV)","U[e1]~(eV)", "U[e2]~(eV)","S[0]~(MeV~b)"))
#
pdf("plot/He3dp_corr.pdf",height = 7,width = 7)
pair_wise_plot(Sp0)
dev.off()



pdf("plot/He3dp_scale_syst.pdf",height = 5,width = 4.5)
plot_normfactors(Normfit)
dev.off()


jagsresults(x = Normfit, params = c("scale"),probs = c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))






require(nuclear)

Tgrid = c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,
    0.013,0.014,0.015,0.016,0.018,0.020,0.025,0.030,0.040,0.050,0.060,0.070,
    0.080,0.090,0.100,0.110,0.120,0.130,0.140,0.150,0.160,0.180,0.200,0.250,0.300,
    0.350,0.400,0.450,0.500,0.600,0.700,0.800,0.900,1.000,1.250,1.500,1.750,2.000,2.500,3.000,3.500,4.000,5.000,
    6.000,7.000,8.000,9.000,10.000)


NAI <- table_reaction_He3dp(Normfit, vars=c("E0","Er","gd2", "gp2", "ad","ap"),N=800,T9=Tgrid )
NAII <- table_reaction_He3dp(Normfit, vars=c("E0_b","Er_b","gd2_b", "gp2_b", "ad_b","ap_b"),N=800,T9=Tgrid )



#write.csv(NAII ,"NA_II.csv",row.names = F)


Norm <- NAII$mean


NAI_new <- NAI  %>%  mutate(data = "present") %>% 
  select(c("T9","mean","lower","upper")) %>%
  set_colnames(c("T9","Adopted","Lower","Upper")) %>%
  mutate(Adopted = Adopted/Norm) %>% 
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)  %>%  
  mutate(data="presentI")  


NAII_new <- NAII  %>%  mutate(data = "present") %>% 
  select(c("T9","mean","lower","upper")) %>%
  set_colnames(c("T9","Adopted","Lower","Upper")) %>%
  mutate(Adopted = Adopted/Norm) %>% 
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)  %>%  
  mutate(data="presentII") 


old <- read.csv("tabula-tab_he3dp.csv",header = TRUE) %>%  
  select(c("T9","Adopted","Lower","Upper"))  %>%
  mutate(data="previous") %>% 
  mutate(Adopted = Adopted/Norm) %>% 
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)


joint <- rbind(old,NAI_new,NAII_new)
joint$data <- as.factor(joint$data)
joint$data <- factor(joint$data, levels = c("previous","presentI","presentII"))



#write.csv(joint ,"joint_rate.csv",row.names = F)

#joint <- read.csv("joint_rate.csv",header = T)

jointf <- filter(joint, data %in% c("previous","presentII"))
jointf$data <- factor(jointf$data, levels = c("presentII","previous"))

pdf("rate_ratio_he3dp.pdf",height = 0.75*7,width = 0.75*10)
ggplot(jointf,aes(x=T9,y=Adopted, group=data,fill=data,linetype=data)) +
#  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="#F0F8FF",alpha=0.4) +
  geom_ribbon(aes(x=T9,ymin=Lower, ymax=Upper),show.legend=FALSE,alpha=0.65) +
  geom_line(size=0.5) +
  coord_cartesian(ylim=c(0.9,1.06),xlim=c(0.00125,1)) +
  theme_bw() + xlab("Temperature (GK)") + ylab("Reaction rate ratio") +
  scale_fill_manual(values=c("#969696","#e41a1c"),name="") +
  scale_x_log10(breaks = c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))  +
  annotation_logticks(short = unit(0.2, "cm"), mid = unit(0.3, "cm"), long = unit(0.4, "cm"),
                      sides = "b") +
  #  annotation_logticks(base=2.875,sides = "l") +
  scale_linetype_manual(guide=F,values=c("dashed","solid"),name="") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
 #       plot.background = element_rect(size = 3, linetype='dashed',colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size = 22),
        axis.ticks = element_line(size = 0.75),
#        axis.line = element_line(linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3, "mm"),
       panel.border = element_rect(size = 0.5)
) 
dev.off()



jointI_II <- filter(joint, data %in% c("presentI","presentII"))
jointI_II$data <- factor(jointI_II$data, levels = c("presentII","presentI"))

pdf("rate_ratio_he3dp_I_II.pdf",height = 7,width = 10)
ggplot(jointI_II,aes(x=T9,y=Adopted, group=data,fill=data,linetype=data)) +
#  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="#F0F8FF",alpha=0.2) +
  geom_ribbon(aes(x=T9,ymin=Lower, ymax=Upper),show.legend=FALSE,alpha=0.65) +
  geom_line(size=0.75) +
  coord_cartesian(ylim=c(0.9,1.1),xlim=c(0.00125,1)) +
  theme_bw() + xlab("Temperature (GK)") + ylab("Reaction rate ratio") +
  scale_fill_manual(values=c("#969696","#e41a1c"),name="") +
  scale_x_log10(breaks = c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))  +
  annotation_logticks(short = unit(0.2, "cm"), mid = unit(0.3, "cm"), long = unit(0.4, "cm"),
                      sides = "b") +
#  annotation_logticks(base=2.875,sides = "l") +
  scale_linetype_manual(guide=F,values=c("dashed","solid"),name="",labels = c("Descouvemont", "Present","Case I")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(size = 3, linetype='dashed',colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size=22),
        axis.ticks = element_line(size = 0.75),
#        axis.line = element_line(size = 0.75, linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3, "mm"))
 #       panel.border = element_rect(size = 1.2)) 
dev.off()






Delta_y <- jagsresults(x=Normfit , params=c('DeltaM'),probs=c(0.025, 0.16, 0.5, 0.84, 0.975))
gdata02 <- data.frame(x =xx, mean = Delta_y[,"mean"],lwr1=Delta_y[,"16%"],lwr2=Delta_y[,"2.5%"],
                     upr1=Delta_y[,"84%"],
                     upr2=Delta_y[,"97.5%"])



pdf("plot/Delta_models.pdf",height = 0.75*7,width = 0.75*10)
ggplot(gdata02,aes(x=x,y=mean))+
  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="gray90",alpha=0.3) +


  # Delta  Bare
# geom_ribbon(data=gdata02,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("#dadaeb"),show.legend=FALSE)+
  geom_ribbon(data=gdata02,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.65, fill = c("#9e9ac8"),show.legend=FALSE) +
  geom_ribbon(data=gdata02,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),alpha=0.85,fill=c("#984ea3"),show.legend=FALSE) +
  geom_line(linetype="dashed")+
  #  Bare
#  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("#dadaeb"),show.legend=FALSE)+
#  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#9e9ac8"),show.legend=FALSE) +
#  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#984ea3"),alpha=0.5,show.legend=FALSE) +
  #
  #
 coord_cartesian(xlim=c(5e-3,0.325),ylim=c(0.9,1.055)) +
  theme_bw() + xlab("Energy (MeV)") + ylab("S-factor ratio") +
  scale_x_log10(breaks = c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))  +
  annotation_logticks(short = unit(0.2, "cm"), mid = unit(0.25, "cm"), long = unit(0.3, "cm"),
                      sides = "b",size = 0.45) +
  #  annotation_logticks(base=2.875,
  #  short = unit(0.2, "cm"), mid = unit(0.25, "cm"), long = unit(0.3, "cm"),sides = "l",size = 0.45) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.925,0.7),
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size=14,colour = set),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size=22),
        axis.text  = element_text(size=18),
        axis.ticks = element_line(size = 0.45),
        #        axis.line = element_line(size = 0.45, linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3, "mm")) 
dev.off()



####
#1-sum(dtc$ad_b/6 >= 0.95  &  dtc$ad_b/6 <= 1.05)/length(dtc$ad_b)
#sum(dtc$ad_b/6 >= 0.7  &  dtc$ad_b/6 <= 1.3)/length(dtc$ad_b)


#1-sum(dtc$ap_b/5 >= 0.9  &  dtc$ap_b/5 <= 1.1)/length(dtc$ap_b)
#sum(dtc$ap_b/5 >= 0.7  &  dtc$ap_b/5 <= 1.3)/length(dtc$ap_b)


#sum(dtc$S_0/5.9 >= 0.9  &  dtc$S_0/5.9 <= 1.1)/length(dtc$S_0/5.9)

#quantile(dtc$E0-0.35779,probs=c(0.0015,0.025, 0.16, 0.5, 0.84, 0.975,0.9985))
#sum(dtc$E0/0.35779 >= 0.951  &  dtc$E0/0.35779 <= 1.049)/length(dtc$E0)
#sum(dtc$gp2/0.025425 >= 0.95  &  dtc$gp2/0.025425 <= 1.05)/length(dtc$gp2)
#sum(dtc$gd2/1.0085 >= 0.95  &  dtc$gd2/1.0085 <= 1.05)/length(dtc$gd2)
#1-sum(dtc$`ue_ev[1]`/rnorm(835,219,7) >= 0.98  &  dtc$`ue_ev[1]`/rnorm(835,219,7) <= 1.02)/nrow(dtc)
#1 - sum(dtc$`ue_ev[2]`/rnorm(835,146,5) >= 0.85  &  dtc$`ue_ev[2]`/rnorm(835,146,5) <= 1.15)/nrow(dtc)
#1 - sum(dtc$`ue_ev[2]`/65 >= 0.6  &  dtc$`ue_ev[2]`/65 <= 1.4)/nrow(dtc)
#1-sum(dtc2$E0_b/0.35779 >= 0.975  &  dtc2$E0_b/0.35779 <= 1.025)/nrow(dtc2)
#1-sum(dtc2$Er_b/0.35779 >= 0.975  &  dtc2$Er_b/0.35779 <= 1.025)/nrow(dtc2)
#sum(dtc2$ad_b/5 >= 0.95  &  dtc2$ad_b/5 <= 1.05)/nrow(dtc2)
#1-sum(dtc2$ap_b/5 >= 0.95  &  dtc2$ap_b/5 <= 1.05)/nrow(dtc2)

ad_rope <- data.frame(ad=dtc2$ad_b)
pdf("plot/ad_rope.pdf",height = 0.75*7,width = 0.75*10)
ggplot(ad_rope,aes(x=ad)) +
  stat_density_ridges(aes(y=0,fill=factor(..quantile..),alpha=factor(..quantile..)),geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      quantiles = c(0.025, 0.16,  0.84, 0.975)) +
  scale_fill_manual(name = "Probability", values = c("#f2f0f7","#bcbddc","#756bb1",
                                                     "#bcbddc","#f2f0f7")) +
  geom_rect(aes( xmin=5.4,xmax=6.6,ymin=0,ymax=2),color="black",fill="skyblue",alpha=0.05,linetype="dashed") +
  #  geom_density(fill="#21908CBF") +
  
  theme_bw() + coord_cartesian(xlim=c(4,6.5)) +
  annotate(geom = "text",label=paste("2.5% in \n","ROPE"),x=6,y=1,size=7) +
  xlab(expression(paste(a[d]," (fm)"))) + 
  ylab("Probability density") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size=14,colour = set),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size=22),
        axis.text  = element_text(size=18),
        axis.ticks = element_line(size = 0.45),
        #        axis.line = element_line(size = 0.45, linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3, "mm")) 
dev.off()


