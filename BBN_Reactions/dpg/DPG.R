# Dpg analysis
#
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
require(ggridges);require(plyr);require(MASS);
require(R2jags);require(scales)
source("..//auxiliar_functions/jagsresults.R")
source("..//auxiliar_functions/pair_wise_plot.R")

######################################################################
# DATA INPUT
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the 
# latter is the individual statistical error of each datum [i];
# energy is in MeV, S-factor in MeVb

ensamble <- read.csv("ensamble_DPG.csv",header = T) 

re <- as.numeric(ensamble$lab)
Nre <- length(unique(ensamble$lab))

N <- nrow(ensamble)
obsy <- ensamble$S*1e6    # Response variable
obsx <-  ensamble$E   # Predictors
erry <- ensamble$Stat*1e6
set <- ensamble$lab

syst = c(0.08617770,0.07696104,0.08617770,0.04401689)
M <- 500
xx <- seq(0.9*min(obsx),1.1*max(obsx),length.out = M)

theory <- read.table("Marcucci2005.dat", header = FALSE)

interp.x <- theory[,1]
interp.y <- theory[,2]

model.data <- list(obsy = obsy,    # Response variable
                   obsx =  obsx,   # Predictors
                   erry = erry,
                   N = nrow(ensamble), # Sample size
                   syst = syst,
                   Nre = Nre,
                   re = re,
                   M = M,
                   xx = xx,
                   interp.x = interp.x,
                   interp.y = interp.y
)

######################################################################
Model <- "model{
# LIKELIHOOD
for (i in 1:N) {
obsy[i] ~ dnorm(ya[i], pow(erry[i], -2))
ya[i] ~ dnorm(y.norm[re[i]]*mut[i],pow(tau[re[i]], -2))
mut[i] <- a.scale*interp.lin(obsx[i], interp.x, interp.y)
}

# Predicted values

for (j in 1:M){

# Bare...

mux[j] <- a.scale*interp.lin(xx[j], interp.x, interp.y)
}


### scaling factor of theory 
a.scale ~ dt(0, pow(5,-2), 1)T(0,)    

for (k in 1:Nre){
y.norm[k] ~ dlnorm(log(1.0),pow(syst[k],-2))
tau[k] ~  dt(mt, pow(5,-2), 1)T(0,)
}

mt ~  dt(0, pow(5,-2), 1)T(0,)

}"


inits <- function(){list(a.scale = runif(1,0.5,1.5),y.norm=runif(4,0.5,1.5),tau=runif(4,0.01,10)) }

# JAGS model with R2Jags;
Normfit <- jags(data = model.data,
                inits = inits,
                parameters.to.save  = c("a.scale", "y.norm","tau","mux"),
                model.file  = textConnection(Model),
                n.thin = 5,
                n.chains = 3,
                n.burnin = 10000,
                n.iter = 20000)


jagsresults(x=Normfit , params=c('tau',"a.scale", "y.norm"),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))


y <- jagsresults(x=Normfit , params=c('mux'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))


gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],lwr3=y[,"0.5%"],upr1=y[,"75%"],
                    upr2=y[,"97.5%"],upr3=y[,"99.5%"])
gobs <- data.frame(obsx,obsy,erry,set)



ggplot(gobs,aes(x=obsx,y=obsy))+
  
   #  Bare
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),fill=c("gray80"),show.legend=FALSE)+
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),  fill = c("gray60"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),fill=c("gray40"),show.legend=FALSE) +
  #
  #  
  
  geom_point(data=gobs,aes(x=obsx,y=obsy,group=set,color=set,shape=set),size=2.75)+
  geom_errorbar(show.legend=FALSE,data=gobs,aes(x=obsx,y=obsy,ymin=obsy-erry,ymax=obsy+erry,group=set,color=set),
                width=0.01,alpha=0.4)+

  scale_colour_stata(name="")+
  scale_shape_manual(values=c(0,19,8,10,4,17,3),name="") +
  theme_bw() + xlab("Energy (MeV)") + ylab("S-Factor (eV b)") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  +
  scale_y_log10() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.1,0.675),
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size=14,colour = "black"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(size=18.5),
        axis.text  = element_text(size=13),
        axis.ticks = element_line(size = 0.75),
        axis.line = element_line(size = 0.5, linetype = "solid")) 





Sp <- ggs(as.mcmc(Normfit)[,c("a.scale", "y.norm[1]","y.norm[2]","y.norm[3]","y.norm[4]",
                              "tau[1]","tau[2]","tau[3]","tau[4]")]) %>% as_tibble()

pair_wise_plot(Sp)