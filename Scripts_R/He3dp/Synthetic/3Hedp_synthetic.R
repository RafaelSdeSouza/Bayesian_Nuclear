# He3dp analysis
#
# purpose: ARTIFICIAL DATA
#
# - 3 parameters are assumed: Er, gamma_d^2, gamma_n^2 [e1, gin, gout]
#
# - uses the function sfactor3Hedp(obsx1[i], e1, gin, gout), which
#   is a C++ version of a Fortran code that includes Coulomb wave
#   function calculations; JAGS has been recompiled with this C++ function
#
#
######################################################################
# preparation: remove all variables from the work space
#rm(list=ls())
set.seed(123)
######################################################################
# data input
# format: obsx, obsy, errobsy; the latter are the individual statistical
# errors of each datum [i]
#
# energy is in units of MeV, and the S-factor in MeVb;
######################################################################
######################################################################
# import jags package
library(rjags);library(R2jags);library(mcmcplots)
require(RcppGSL);require(ggplot2);require(ggthemes)
require(nuclear);require(ggmcmc);require(plyr)
source("https://raw.githubusercontent.com/johnbaums/jagstools/master/R/jagsresults.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")
load.module("nuclear")

######################################################################
## ARTIFICIAL DATA GENERATION

N <- 100

#obsx1 <- runif(N,0,0.7)
tau <- 1.5
obsx1 <- exp(runif(N,log(1e-3),log(1)))
errobsy1 <- runif(N,0.25,1.25)
y1 <- rnorm(N, sfactor3Hedp(obsx1 ,0.35779,1.0085,0.025425),errobsy1)
obsy1 <- rnorm(N,y1,tau)

M <- 150
xx <- seq(min(obsx1),max(obsx1),length.out = M)
model.data <- list(obsy = obsy1,    # Response variable
                   obsx =  obsx1,   # Predictors
                   erry = errobsy1,
                   N = N, # Sample size
                   M = M,
                   xx = xx
)
#
######################################################################
cat('model {

    # LIKELIHOOD
    for (i in 1:N) {
    obsy[i] ~ dnorm(y1[i], pow(erry[i], -2))
    y1[i] ~ dnorm(sfactor3Hedp(obsx[i], e1, gin, gout,6,5,0),pow(tau,-2))
    }

    # Predicted values

    for (j in 1:M){
    mux[j] <- sfactor3Hedp(xx[j], e1, gin, gout,6,5,0)
    yx[j] ~ dnorm(mux[j],pow(tau,-2))
    }

    # PRIORS
    # e1, gin, gout are defined as in tdn.f (by Alain Coc):
    # resonance energy, initial reduced width, final reduced
    # width;

    tau ~ dunif(0,100)
    e1 ~ dunif(0,20)
    gin ~ dunif(1e-4,10)
    gout ~ dunif(1e-4,10)
    }', file={f <- tempfile()})
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


n.burnin  <- 5000
n.iter   <- 10000

n.chains <- 3
n.thin   <- 5
inits <- function () { list(e1 = runif(1,0.1,0.75),gin=runif(1,2,10),gout=runif(1,0.01,1)) }
# "f": is the model specification from above;

# JAGS model with R2Jags;
out <- jags(data = model.data,
            inits = inits,
            parameters = c("e1", "gin", "gout","tau","mux","yx"),
            model.file = f,
            n.thin = n.thin,
            n.chains = n.chains,
            n.burnin = n.burnin,
            n.iter = n.iter)
jagsresults(x=out, params=c("e1", "gin", "gout","tau"),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))


traplot(out ,c("e1", "gin", "gout"),style="plain")
denplot(out ,c("e1", "gin", "gout"),style="plain")


# Plot
y <- jagsresults(x=out, params=c('mux'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))
x <- xx
gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],lwr3=y[,"0.5%"],upr1=y[,"75%"],
                    upr2=y[,"97.5%"],upr3=y[,"99.5%"])
gobs <- data.frame(obsx1,obsy1,errobsy1)


# Import the kitten
#img <- image_read('http://thecatapi.com/api/images/get?size=med')
pdf("plot/He3dp_synthetic.pdf",height = 7,width = 8)
ggplot(gobs,aes(x=obsx1,y=obsy1))+
#  annotation_custom(rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc")), 
#                    -Inf, Inf, -Inf, Inf) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr3, ymax=upr3,y=NULL), fill=c("#deebf7"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), fill = c("#9ecae1"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL), fill=c("#3182bd"),show.legend=FALSE) +
  geom_point(colour="#969696")+
  geom_errorbar(data=gobs,mapping=aes(x=obsx1,y=obsy1,ymin=obsy1-errobsy1,ymax=obsy1+errobsy1),alpha=0.85,
                colour="#969696",width=0.05)+
  geom_line(data=gdata,aes(x=xx,y=mean),colour="#ffffff",linetype="dashed",size=1,show.legend=FALSE)+
  theme_wsj() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + scale_x_log10()  +
  theme(legend.position = "none",
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=10),
        strip.background = element_rect("gray85")) +
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He"))) 
dev.off()




S <- ggs(as.mcmc(out)[,c("e1", "gin", "gout","tau")])

S$Parameter <- revalue(S$Parameter, c("e1"="E[r]", "gin"="Gamma['in']",
                                      "gout"="Gamma['out']","tau"="tau"))
vline.dat <- data.frame(Parameter=factor(c("E[r]","Gamma['in']","Gamma['out']","tau")), vl=c(0.35779,1.0085,0.025425,1.5))

pdf("plot/He3dp_synthetic_posterior.pdf",height = 7,width = 8)
ggplot(data=S,aes(x=value,group=Parameter,fill=Parameter)) +
  geom_density(alpha=0.75) + 
  theme_wsj() +
  scale_fill_tableau()+
  theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="#666666", face="bold", size=15),
        axis.text  = element_text(size=12),
        strip.text = element_text(size=15),
        strip.background = element_rect("white")) + 
  geom_vline(aes(xintercept=vl), data=vline.dat,linetype="dashed",color="gray75",size=1) +
  facet_wrap(~Parameter,scales="free",ncol=2,nrow=2,labeller=label_parsed) +
  ylab("Posterior probability") + xlab("Parameter value")+
  ggtitle(expression(paste(NULL^"3","He(d,p)",NULL^"4","He"))) 
dev.off()


pdf("plot/He3dp_synthetic_trace.pdf",height = 7,width = 8)
ggplot(data=S,aes(x= Iteration,y=value,group=Parameter,color=factor(Chain))) +
  geom_line(alpha=0.5,size=0.25) + 
  theme_wsj() +
  scale_color_tableau()+
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