# 3Hedp analysis
#
# purpose: Real  DATA
#
# - 7 R-matrix parameters are assumed: E_0, Er, gamma_d^2, gamma_p^2, ad, ap, ue  [e0,er,gd2, gp2,ad,ap, ue]
#
# - uses the function sfactor3Hedp_5p(obsx, e0,er,gd2, gp2,ad,ap,ue), which
#   is a R  version of a Fortran code that includes Coulomb wave
#   function calculations.
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
set.seed(27)
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb;

######################################################################
# import packages
require(RcppGSL);require(ggplot2);require(ggthemes)
require(nuclear);library(magrittr);
library(dplyr);require(lessR);library(BayesianTools)
require(msm);require(mcmcplots);require(ggmcmc);
require(ggridges);require(plyr);require(MASS)
source("..//..//auxiliar_functions/pair_wise_plot.R")
source("external_functions/plot_Sfactor_DREAM.R")
source("external_functions/plot_normfactors_DREAM.R")

######################################################################
## Read DATA
ensamble <- read.csv("data/ensamble.csv",header = T) %>%
  mutate(Syst=replace(Syst,Syst==0.06,0.078))  %>% filter(E <= 0.5)



He3dpdata <- read.csv("3Hedpdata.csv",header = T)


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
syst = 1 + c(0.03,unique(ensamble$Syst))

likelihood <- function(par){
  e0 = par[1]
  er = par[2]
  gd2 = par[3]
  gp2 = par[4]
  ad   = par[5]
  ap =  par[6]
  sigmax = par[7]
  scale = par[8:14]
  ue = par[15:16]
  y = par[17:(N + 16)]
  nu = par[231]

  llRandom = sum(dlnorm(scale,meanlog = log(1), sdlog = log(syst), log = T))
#  lly <- sum(dst(y,mu = scale[re]*sfactor3Hedp_5p(obsx, e0,er,gd2, gp2,ad,ap,ue = ue[ik]),sigma = sigmax,nu = nu, log = T))
  lly <- sum(dnorm(y,mean = scale[re]*sfactor3Hedp_5p(obsx, e0,er,gd2, gp2,ad,ap,ue = ue[ik]),sd = sigmax, log = T))
  llobs = sum(dnorm(obsy,mean = y,sd = erry,log = T))
  return(llRandom + llobs + lly)

}


low <- c(0.1,1e-3,rep(1e-4,2), 1,1,1e-4,rep(0.5,7),rep(0,2),obsy - 2*erry,1)
up <- c(0.4,1,rep(3,2),10,10,5,rep(1.5,7),rep(350,2),obsy + 2*erry,100)


createHedPrior <- function(lower, upper, best = NULL){
density = function(par){
  d1 = dunif(par[1], 0.1, 0.4,log = TRUE)
  d2 = dtnorm(par[2], 0, 1,log = TRUE)
  d3 = dtnorm(par[3], 0, 3*34.625/par[5]^2, log = TRUE)
  d4 = dtnorm(par[4], 0, 3*51.94605/par[6]^2, log = TRUE)
  d5 = sum(dunif(par[5:6], 2, 10, log = TRUE))
  d6 = dtnorm(par[7], mean = 0, sd = 5, log = TRUE)
  d7 = sum(dlnorm(par[8:14],log(1),log(syst),log = TRUE))
  d8 = sum(dtnorm(par[15:16], mean = 0, sd = 1E4, log = TRUE))
  d9 = sum(dunif(par[17:(N + 16)],obsy - 2*erry,obsy + 2*erry,log = TRUE))
  d10 = dunif(par[231],1,100,log=TRUE)
  return(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9 + d10)
}

sampler = function(){
   c(runif(1, 0.1, 0.35),
   runif(1, 0, 1),
    exp(runif(2, log(1e-3), log(2))),
    runif(2, 2.5, 7),
    runif(1, 0, 5),
    rlnorm(7, log(1), log(syst)), #ynorm
    runif(2, 0, 300),
    runif(N, obsy - 2*erry,obsy + 2*erry),
    runif(1,1,100)
)
}

out <- createPrior(density = density, sampler = sampler, lower = lower, upper = upper, best = best)
return(out)
}

prior <- createHedPrior(lower = low, upper = up)


#prior <- createUniformPrior(lower = low,
#                            upper = up)



setup <- createBayesianSetup(likelihood = likelihood,prior = prior,
names = c("e0","er","gd2","gp2","ad","ap","sigma",to("scale", 7),to("ue", 2),to("y", N),"nu"))

#setup <- createBayesianSetup(likelihood = likelihood,lower = low,upper = up,
#names = c("e0","er","gd2","gp2","ad","ap","sigma",to("scale", 7),to("ue", 2),to("y", N)))



settings <- list(iterations = 2E6,thin=5,
                 burnin = 1E4, message = T,nrChains = 1)




system.time(
res <- runMCMC(bayesianSetup = setup,  settings = settings,sampler = "DEzs")
)


# Quick diagnostic
#summary(res)
#tracePlot(sampler = res, thin = 1, start = 2E4, whichParameters = c(1,2,3,4,5,6,15,16,231))


codaObject = getSample(res, start = 2E4, coda = TRUE)

getmcmc_var <- function(outjags=outjags,vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
#write.matrix(ssDat,"He3dp_DREAM.dat")


sDat <- getmcmc_var(codaObject,vars = c("e0","er","gd2","gp2","ad","ap","sigma",to("scale", 7),to("ue", 2),to("y", N)))
index <- sample(seq(1:nrow(sDat)),1E4,replace=FALSE)
ssDat <- sDat[index,]


# Plot fit
xx <- exp(seq(log(min(obsx)),log(max(obsx)),length.out = 250))
plot_Sfactor_DREAM(ssDat)
#





dream_dat <- read.table("He3dp_DREAM.dat",header=T)
ssDat <- dream_dat[,c("e0","er","gd2","gp2","ad","ap","ue1","ue2")]




xx <- 1e-2
S0 <- sfactor3Hedp_5p(xx,ssDat[,1],ssDat[,2],ssDat[,3],ssDat[,4],ssDat[,5],ssDat[,6],0)
quantile(S0,prob=c(0.16, 0.5, 0.84))


# Read chains
dream_dat <- read.table("data/He3dp_DREAM.dat",header=T)


# Case II
SpII <- ggs(as.mcmc(ssDat[,c("e0","er","gd2","gp2","ad","ap","ue1","ue2")]))
Sp0II <- SpII %>% as_tibble()
Sp0II$Parameter <- ordered(Sp0II$Parameter, levels =c("e0","er","gd2","gp2","ad","ap","ue1","ue2"))
levels(Sp0II$Parameter) <- as.factor(c("E[0]~(MeV)","E[r]~(MeV)","gamma[d]^2~(MeV)", "gamma[p]^2~(MeV)","a[d]~(fm)",
                                       "a[p]~(fm)","U[e1]~(eV)", "U[e2]~(eV)"))
#

# Corner plot
Corr_chain <- ggs(as.mcmc(ssDat[,c("e0","er","gd2","gp2","ad","ap","ue1","ue2")])) %>%
as_tibble() %>%
mutate(Parameter = factor(Parameter, levels = c("e0","er","gd2","gp2","ad","ap","ue1","ue2"))) %>%
mutate(Parameter  = factor(Parameter, labels = c("E[0]~(MeV)","E[r]~(MeV)","gamma[d]^2~(MeV)",
                                                 "gamma[p]^2~(MeV)","a[d]~(fm)",
                                                 "a[p]~(fm)","U[e1]~(eV)", "U[e2]~(eV)")))

pdf("plot/He3dp_corr.pdf",height = 9,width = 9)
pair_wise_plot(Corr_chain)
dev.off()



# Normalization Factors

pdf("plot/He3dp_scale_syst.pdf",height = 5,width = 4.5)
plot_normfactors_DREAM(ssDat)
dev.off()



# Estimate S_0
quantile(sfactor3Hedp_5p(1e-4,ssDat[,1],ssDat[,2],ssDat[,3],ssDat[,4],
                         ssDat[,5],ssDat[,6],0),probs = c(0.16,0.5,0.84))
quantile2 <- function(x) quantile(x, probs = c(0.16,0.5,0.84))
su <- apply(sDat, 2, quantile2)
#





# Estimate reaction rates

require(nuclear)

Tgrid = c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,
          0.013,0.014,0.015,0.016,0.018,0.020,0.025,0.030,0.040,0.050,0.060,0.070,
          0.080,0.090,0.100,0.110,0.120,0.130,0.140,0.150,0.160,0.180,0.200,0.250,0.300,
          0.350,0.400,0.450,0.500,0.600,0.700,0.800,0.900,1.000)


NAI <-  table_reaction_He3dp(ssDat, vars=c("e0","er","gd2","gp2","ad","ap"),N=5000,T9=Tgrid )


#write.csv(NAI ,"NA_I.csv",row.names = F)
#NAI <- read.csv("NA_I.csv")

Norm <- NAI$mean


NAI_new <- NAI[,c("T9","mean","lower","upper")]  %>%
  set_colnames(c("T9","Adopted","Lower","Upper")) %>%
  mutate(Adopted = Adopted/Norm) %>%
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)  %>%
  mutate(data="presentI")


old <- read.csv("tabula-tab_he3dp.csv",header = TRUE)
old <- old[,c("T9","Adopted","Lower","Upper")]  %>%
  mutate(data="previous") %>%
  mutate(Adopted = Adopted/Norm) %>%
  mutate(Lower = Lower/Norm) %>%
  mutate(Upper = Upper/Norm)


joint <- rbind(old,NAI_new)
joint$data <- as.factor(joint$data)
joint$data <- factor(joint$data, levels = c("previous","presentI"))



#write.csv(joint ,"joint_rate.csv",row.names = F)

#joint <- read.csv("joint_rate.csv",header = T)

jointf <- joint
jointf$data <- factor(jointf$data, levels = c("presentI","previous"))

pdf("rate_ratio_he3dp.pdf",height = 0.75*7,width = 0.75*10)
ggplot(jointf,aes(x=T9,y=Adopted, group=data,fill=data,linetype=data)) +
  #  geom_rect(aes(xmin=0.045, xmax=0.356, ymin=-1, ymax=22), fill="#F0F8FF",alpha=0.4) +
  geom_ribbon(aes(x=T9,ymin=Lower, ymax=Upper,alpha=data),show.legend=FALSE) +
  scale_alpha_manual(values=c(1,0.375))+
  geom_line(size=0.5) +
  coord_cartesian(ylim=c(0.9,1.06),xlim=c(0.00125,1)) +
  theme_bw() + xlab("Temperature (GK)") + ylab("Reaction rate ratio") +
  scale_fill_manual(values=c("gray50","#984ea3"),name="") +
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



