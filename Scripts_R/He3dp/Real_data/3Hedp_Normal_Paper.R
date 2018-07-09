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
tau_2  ~    dnorm(0, pow(1,-2))T(0,)
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

inits <- function () { list(E0 = runif(1,0.3,0.35),E0_b = 0.4,Er_b = 0.4,gd2 = 1,
                        gp2 = runif(1,0.01,0.1),gd2_b = 0.5) }


# JAGS model with R2Jags;
Normfit <- jags(data = model.data,
                inits = inits,
                parameters.to.save  = c("Er","E0","gd2", "gp2","ue_ev","tau", "ad","ap",
                                        "RSS","mux0","mux1","mux2","scale","DeltaM","S_0",
                                        "S_0b","E0_b","Er_b","gd2_b",
                                        "gp2_b","tau_2","ad_b","ap_b" ),
                model.file  = textConnection(Model),
                n.thin = 30,
                n.chains = 10,
                n.burnin = 7500,
                n.iter = 15000)

temp <- Normfit
temp <- update(temp, n.thin = 50, n.iter = 50000)

#hdi_jags <- function(mcmc=mcmc, par = par,credMass = 0.95,allowSplit=TRUE){
#  temp <- as.mcmc(mcmc)[,c(par)]
#  return(hdi(temp,allowSplit=TRUE))
#}
#hdi_jags(Normfit,par=c("Er","gd2", "gp2","ue_ev[1]","ue_ev[2]","S_0"),credMass = 0.95)


tab <- jagsresults(x = Normfit, params = c("E0","gd2", "gp2","ue","tau", "ad","ap","ue_ev","S_0"),probs = c(0.16, 0.5, 0.84))
tab <- as.data.frame(tab)

tab$low <- tab[,4] - tab[,3]
tab$hi <-  tab[,5] - tab[,4]

jagsresults(x = Normfit, params = c("E0","gd2", "gp2","ue","tau", "ad","ap","ue_ev","S_0"),probs = c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))

tab2 <- jagsresults(x = Normfit , params = c("E0_b","Er_b","gd2_b", "gp2_b","tau_2","ad_b","ap_b","ue_ev","S_0b"),probs = c(0.16, 0.5, 0.84))
tab2 <- as.data.frame(tab2)
tab2$low <- tab2[,4] - tab2[,3]
tab2$hi <-  tab2[,5] - tab2[,4]


traplot(Normfit  ,c("E0","gd2", "gp2","ue"),style="plain")
denplot(Normfit  ,c("E0","gd2", "gp2","ue"),style="plain")


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
levels(Sp0$Parameter) <- as.factor(c("E[0]~(MeV)","gamma[d]^2~(MeV)", "gamma[p]^2~(MeV)","U[e1]~(eV)", "U[e2]~(eV)","S[b](0)~(MeV~b)"))
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
#joint$data <- factor(joint$data, levels = c("previous","presentII","presentI"))



write.csv(joint ,"joint_rate.csv",row.names = F)

joint <- read.csv("joint_rate.csv",header = T)

jointf <- filter(joint, data %in% c("previous","presentII"))

pdf("rate_ratio_he3dp.pdf",height = 7,width = 10)
ggplot(jointf,aes(x=T9,y=Adopted, group=data,fill=data,linetype=data,alpha=0.3)) +
  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="gray90",alpha=0.4) +
  geom_ribbon(aes(x=T9,ymin=Lower, ymax=Upper),show.legend=FALSE) +
  geom_line() +
  coord_cartesian(ylim=c(0.9,1.1),xlim=c(0.00125,1)) +
  theme_bw() + xlab("Temperature (GK)") + ylab("Reaction rate ratio") +
  scale_fill_manual(values=c("#606060","#6600CC"),name="") +
  scale_x_log10(breaks = c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))  +
  annotation_logticks(sides = "b") +
  annotation_logticks(base=2.875,sides = "l") +
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
        axis.line = element_line(size = 0.75, linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3.75, "mm"),
        panel.border = element_rect(size = 1.2)) 
dev.off()



jointI_II <- filter(joint, data %in% c("presentI","presentII"))
jointI_II$data <- factor(jointI_II$data, levels = c("presentII","presentI"))

pdf("rate_ratio_he3dp_I_II.pdf",height = 7,width = 10)
ggplot(jointI_II,aes(x=T9,y=Adopted, group=data,fill=data,linetype=data,alpha=0.3)) +
  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="gray90",alpha=0.4) +
  geom_ribbon(aes(x=T9,ymin=Lower, ymax=Upper),show.legend=FALSE) +
  geom_line() +
  coord_cartesian(ylim=c(0.9,1.1),xlim=c(0.00125,1)) +
  theme_bw() + xlab("Temperature (GK)") + ylab("Reaction rate ratio") +
  scale_fill_manual(values=c("#606060","#ff7f00"),name="") +
  scale_x_log10(breaks = c(0.001,0.01,0.1,1),labels=c("0.001","0.01","0.1","1"))  +
  annotation_logticks(sides = "b") +
  annotation_logticks(base=2.875,sides = "l") +
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
        axis.line = element_line(size = 0.75, linetype = "solid"),
        axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(-3.75, "mm"),
        panel.border = element_rect(size = 1.2)) 
dev.off()













Delta_y <- jagsresults(x=Normfit , params=c('DeltaM'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))
gdata02 <- data.frame(x =xx, mean = Delta_y[,"mean"],lwr1=Delta_y[,"25%"],lwr2=Delta_y[,"2.5%"],
                      lwr3=Delta_y[,"0.5%"],upr1=Delta_y[,"75%"],
                     upr2=Delta_y[,"97.5%"],upr3=Delta_y[,"99.5%"])



pdf("plot/Delta_models.pdf",height = 7,width = 10)
ggplot(gdata02,aes(x=x,y=mean))+
  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="gray90",alpha=0.4) +


  # Delta  Bare
 geom_ribbon(data=gdata02,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("#dadaeb"),show.legend=FALSE)+
  geom_ribbon(data=gdata02,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#9e9ac8"),show.legend=FALSE) +
  geom_ribbon(data=gdata02,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#984ea3"),show.legend=FALSE) +

  #  Bare
#  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("#dadaeb"),show.legend=FALSE)+
#  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("#9e9ac8"),show.legend=FALSE) +
#  geom_ribbon(data=gdata0,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("#984ea3"),alpha=0.5,show.legend=FALSE) +
  #
  #
 coord_cartesian(xlim=c(5e-3,0.3425),ylim=c(0.9,1.1)) +
  theme_bw() + xlab("Energy (MeV)") + ylab(expression(delta["S"])) +
  scale_x_log10()  +
  annotation_logticks(sides = "b") +
  annotation_logticks(base=2.875,sides = "l") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.925,0.7),
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









# Diagnostics overleap between prior and posterior

PR_gin <- runif( 15000, 0, 5)
MCMCtrace(Normfit, params = 'gin', priors = PR_gin, pdf = FALSE)

PR_out <- runif( 15000, 0, 5)
MCMCtrace(Normfit, params = 'gout', priors = PR_out, pdf = FALSE,type = "density")

PR_er <- runif( 15000, 0, 5)
MCMCtrace(Normfit, params = 'e1', priors = PR_er, pdf = FALSE)


PR_ri <- rtruncnorm(6000, a=0, b=Inf, mean = 5, sd = 0.447)
MCMCtrace(Normfit, params = 'ri', priors = PR_ri, pdf = FALSE)

PR_rf <- rtruncnorm(6000, a=0, b=Inf, mean = 5, sd = 0.447)
MCMCtrace(Normfit, params = 'rf', priors = PR_rf, pdf = FALSE)

PR_ue <- rtruncnorm(6000, a=0, b=Inf, mean = 0, sd = 0.031)
MCMCtrace(Normfit, params = 'ue', priors = PR_ue, pdf = FALSE)





PR_gin <- rtruncnorm(6000, a=0, b=Inf, mean = 1, sd = 0.25)
MCMCtrace(Normfit, params = 'gin', priors = PR_gin, pdf = FALSE)

PR_out <- rtruncnorm(6000, a=0, b=Inf, mean = 0, sd = 0.05)
MCMCtrace(Normfit, params = 'gout', priors = PR_out, pdf = FALSE)

PR_ri <- rtruncnorm(6000, a=0, b=Inf, mean = 5, sd = 0.447)
MCMCtrace(Normfit, params = 'ri', priors = PR_ri, pdf = FALSE)

PR_rf <- rtruncnorm(6000, a=0, b=Inf, mean = 5, sd = 0.447)
MCMCtrace(Normfit, params = 'rf', priors = PR_rf, pdf = FALSE)

PR_ue <- rtruncnorm(6000, a=0, b=Inf, mean = 0, sd = 0.031)
MCMCtrace(Normfit, params = 'ue', priors = PR_ue, pdf = FALSE)


RSS <- as.matrix(as.mcmc(Normfit)[,c("RSS")])
rss0 <- function(x) crossprod(x-mean(x))[1]
1-mean(RSS)/rss0(obsy)


color_scheme_set("darkgray")
div_style <- parcoord_style_np(div_color = "green", div_size = 0.05, div_alpha = 0.4)
mcmc_parcoord(as.mcmc(Normfit),alpha = 0.05, regex_pars = c("e1", "gin", "gout","ri","rf"))



















ss <- ggs(as.mcmc(Normfit)[,c("e1", "gin", "gout")])

ss$Chain <- as.factor(ss$Chain)

ggs_density(ss)
pdf("plot/He3dp_trace_syst.pdf",height = 7,width = 8)
ggplot(data=ss,aes(x= Iteration,y=value,group=Parameter,color=factor(Chain))) +
  geom_line(alpha=0.5,size=0.25) +
  theme_rafa() +
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
  facet_wrap(~Parameter,scales="free",labeller=label_parsed,nrow=3) +
  ylab("Parameter value") + xlab("Iteration")
dev.off()
#





# Reaction rates

Nsamp <- 2500
Tgrid <- c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,
           0.013,0.014,0.015,0.016,0.018,0.020,0.025,0.030,0.040,0.050,0.060,0.070,
           0.080,0.090,0.100,0.110,0.120,0.130,0.140,0.150,0.160,0.180,0.200,0.250,0.300,
           0.350,0.400,0.450,0.500,0.600,0.700,0.800,0.900,1.000,1.250,1.500,1.750,2.000,2.500,3.000,3.500,4.000,5.000,
           6.000,7.000,8.000,9.000,10.000)








NA_I <-  table_reaction(Normfit,vars=c("Er", "gd", "gp","ad","ap"), N=1000)

NA_II <- table_reaction(Normfit,vars=c("Er_b","E0_b", "gd_b", "gp_b","ad_b","ap_b"),N=1000)

NA_I$Case <- "Case I"
NA_II$Case <- "Case II"
NA_total <-rbind(NA_I,NA_II)


pdf("plot/Delta_rates.pdf",height = 7,width = 10)
ggplot(Drate, aes(x=T,y=mean)) +
  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="gray90",alpha=0.4) +
  geom_line() +
#  geom_ribbon(aes(ymin= lower, ymax=upper)) +
  geom_ribbon(aes(ymin=lwr3, ymax=upr3),fill=c("#ece7f2"),show.legend=FALSE)+
  geom_ribbon(aes(ymin=lwr2, ymax=upr2),  fill = c("#a6bddb"),show.legend=FALSE) +
  geom_ribbon(aes(ymin=lwr1, ymax=upr1),fill=c("#2b8cbe"),show.legend=FALSE) +
 coord_cartesian(xlim=c(0.00125,7.7),ylim=c(0.9,1.05)) +
  scale_fill_tableau() +
  theme_bw() +
  scale_x_log10() +
   scale_alpha(guide = 'none') +
  xlab("Temperature (GK)") +
  ylab(expression(N[A]~I~symbol("\341")*sigma*nu*symbol("\361")/N[A]~II~symbol("\341")*sigma*nu*symbol("\361"))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.925,0.7),
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

write.csv(NA_I ,"NA_I.csv",row.names = F)
write.csv(NA_II ,"NA_II.csv",row.names = F)



#gg2 <- apply(gg, 1, quantile, probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995), na.rm=TRUE)





pdf("plot/He3dp_cross.pdf",height = 7,width = 10)
ggplot(gg2data,aes(x=x,y=mean))+
  theme_bw()  +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr3, ymax=upr3,y=NULL), fill=c("#ffeda0"),show.legend=FALSE) +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr2, ymax=upr2,y=NULL), fill=c("#feb24c"),show.legend=FALSE) +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr1, ymax=upr1,y=NULL), fill=c("#f03b20"),show.legend=FALSE) +
   xlab("Temperature (GK)") + ylab(expression(N[A]~symbol("\341")*sigma*nu*symbol("\361"))) +
 #scale_y_continuous(breaks=c(0,5e7,1e8,1.5e8),labels=c(0,expression(5%*%10^7),
 #                             expression(10^8),expression(1.5%*%10^8))) +
  scale_y_log10()+
  theme(legend.position = "none",
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "gray95"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=17.5),
        axis.text  = element_text(size=13),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()




cmb <- as.data.frame(filter(gg, x <= 0.105 & x >= 0.095))
cmbhist <- data.frame(x=as.numeric(cmb[1:Nsamp]))
pdf("plot/He3dp_hist_cmb.pdf",height = 7,width = 8)
ggplot(cmbhist, aes(x)) +
#  geom_histogram(aes(y=..count../sum(..count..)),bins = 15,fill="#4357a3",color="#d84951") +
  geom_density(fill="#4357a3",color="#d84951") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
#  scale_x_continuous(breaks=c(1.15e8,1.2e8,1.25e8),labels=c(expression(1.15%*%10^8),
 #                             expression(1.20%*%10^8),expression(1.25%*%10^8))) +
  geom_rug(sides = "b", aes(y = 0),colour = "#d84951") +
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=15),
        strip.text = element_text(size=15),
        strip.background = element_rect("#F0B27A"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("Posterior probability") + xlab(expression(N[A]~symbol("\341")*sigma*nu*symbol("\361")))
dev.off()


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



