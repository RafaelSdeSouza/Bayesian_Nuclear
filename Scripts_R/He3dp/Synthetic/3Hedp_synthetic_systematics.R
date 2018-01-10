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
source("/Users/Rafael/Documents/GitHub/JAGS_UNC/Scripts_R/auxiliar_functions/theme_rafa.R")
## for block updating [we do not need to center predictor variables]
load.module("glm")
load.module("nuclear")

######################################################################
## ARTIFICIAL DATA GENERATION

N <- 90
Nre <- 3
#re <- rep(1:Nre, each = 200/4)
re <- sample( c(1,2,3), N, replace=TRUE )


set.seed(123)
a <- rlnorm(Nre, meanlog = 0, sdlog = 0.05)

print(a,2)
#a <- c(0.97,1.04,0.95,1.06)


#obsx1 <- runif(N,0,0.7)
sigma_t <- 1
obsx1 <- exp(runif(N,log(1e-3),log(1)))
errobsy1 <- runif(N,0,0.25*sfactor3Hedp(obsx1 ,0.36,1,0.03))
y1 <- rnorm(N, a[re]*sfactor3Hedp(obsx1 ,0.36,1,0.03),errobsy1)
obsy1 <- rnorm(N,y1,sigma_t)
plot(obsx1,obsy1,col = re,log="x")


M <- 150
xx <- seq(min(obsx1),max(obsx1),length.out = M)
model.data <- list(obsy = obsy1,    # Response variable
                   obsx =  obsx1,   # Predictors
                   erry = errobsy1,
                   N = N, # Sample size
                   M = M,
                   xx = xx,
                   Nre = Nre,
                   re = re, 
                   a = a
)


######################################################################
cat('model {
    
    # LIKELIHOOD
    for (i in 1:N) {
    obsy[i] ~ dnorm(y1[i], pow(erry[i], -2))
    y1[i] ~ dnorm(scale[re[i]]*sfactor3Hedp(obsx[i], e1, gin, gout,6,4,0),pow(tau,-2))
    }
    
    # Predicted values
    
    for (j in 1:M){
    mux[j] <- sfactor3Hedp(xx[j], e1, gin, gout,6,4,0)
    yx[j] ~ dnorm(mux[j],pow(tau,-2))
    }
    
    # PRIORS
    # e1, gin, gout are defined as in tdn.f (by Alain Coc):
    # resonance energy, initial reduced width, final reduced
    # width;
    
    for (k in 1:Nre){
    scale[k] ~ dlnorm(0,pow(log(a[k]),-2))
    }
    
    tau ~ dgamma(0.01,0.01)
    e1 ~ dgamma(0.01,0.01)
    gin ~ dgamma(0.01,0.01)
    gout ~ dgamma(0.01,0.01)
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


n.burnin  <- 10000
n.iter   <- 20000

n.chains <- 3
n.thin   <- 10
inits <- function () { list(e1 = runif(1,0.1,0.75),gin=runif(1,2,10),gout=runif(1,0.01,0.03)) }
# "f": is the model specification from above;

# JAGS model with R2Jags;
out <- jags(data = model.data,
            inits = inits,
            parameters = c("e1", "gin", "gout","tau","mux","yx","scale"),
            model.file = f,
            n.thin = n.thin,
            n.chains = n.chains,
            n.burnin = n.burnin,
            n.iter = n.iter)
jagsresults(x=out, params=c("e1", "gin", "gout","tau","scale"),probs = c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))


jagsresults(x=out, params=c("e1", "gin", "gout","tau","scale"),probs = c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))



traplot(out ,c("e1", "gin", "gout"),style="plain")
denplot(out ,c("e1", "gin", "gout"),style="plain")


# Plot
y <- jagsresults(x=out, params=c('mux'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))
x <- xx
gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],lwr3=y[,"0.5%"],upr1=y[,"75%"],
                    upr2=y[,"97.5%"],upr3=y[,"99.5%"])
gobs <- data.frame(obsx1,obsy1,errobsy1,re)
gobs$re <- as.factor(gobs$re)

# Import the kitten
#img <- image_read('http://thecatapi.com/api/images/get?size=med')
pdf("plot/He3dp_synthetic_syst.pdf",height = 6,width = 8)
ggplot(gobs,aes(x=obsx1,y=obsy1),colour=re)+

  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr3, ymax=upr3,y=NULL), fill=c("#f0f0f0"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), fill = c("#bdbdbd"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL), fill=c("#737373"),show.legend=FALSE) +
  geom_errorbar(show.legend=FALSE,data=gobs,mapping=aes(colour=re,x=obsx1,y=obsy1,ymin=obsy1-errobsy1,ymax=obsy1+errobsy1),alpha=0.6,
                width=0.025)+
  geom_point(aes(shape=re,colour=re),size=4)+
  scale_color_manual(name="",values=c("#cd4f39","#7ebca1","#4e2f19"))+
  scale_shape_manual(values=c(8,17,11),name="")+
#  geom_line(data=gdata,aes(x=xx,y=mean),colour="#ffffff",linetype="dashed",size=1,show.legend=FALSE)+
  theme_rafa() + xlab("Energy (MeV)") + ylab("S-Factor (MeV b)") + 
  scale_x_log10()  +
  theme(legend.position = c(0.9,0.675))
dev.off()




#S <- ggs(as.mcmc(out)[,c("e1", "gin", "gout","tau","scale[1]","scale[2]","scale[3]","scale[4]")])

#S$Parameter <- revalue(S$Parameter, c("e1"="E[r]", "gin"="Gamma['in']",
#                                      "gout"="Gamma['out']","tau"="tau",
#                                      "scale[1]" = "a[1]","scale[2]" = "a[2]","scale[3]" = "a[3]","scale[4]" = "a[4]"))
#vline.dat <- data.frame(Parameter=factor(c("E[r]","Gamma['in']","Gamma['out']","tau","a[1]","a[2]","a[3]","a[4]")), vl=c(0.35779,1.0085,0.025425,1.25,0.95, 0.98, 1.17, 1.01))


S <- ggs(as.mcmc(out)[,c("e1", "gin", "gout","tau")])

S$Parameter <- revalue(S$Parameter, c("e1"="E[r]", "gin"="Gamma['in']",
"gout"="Gamma['out']","tau"="tau"))
vline.dat <- data.frame(Parameter=factor(c("E[r]","Gamma['in']","Gamma['out']","tau")), vl=c(0.36,1,0.03,1))



pdf("plot/He3dp_synthetic_posterior_syst.pdf",height = 7,width = 8)
ggplot(data=S,aes(x=value,group=Parameter,fill=Parameter)) +
  geom_histogram(alpha=0.75) + 
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


my_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(bins = 10,fill="#4271AE",colour="#1F3552",...) +
    theme_void() + theme( panel.grid.minor=element_blank(),
                          panel.grid.major=element_blank())
}

my_bin <- function(data, mapping, ..., low = "#3698BF", high = "#D97C2B") {
  ggplot(data = data, mapping = mapping) +
    geom_bin2d(...) +
    scale_fill_gradient(low = low, high = high) +
    theme_bw()
}

Sp <- ggs(as.mcmc(out)[,c("e1", "gin", "gout")])
levels(Sp$Parameter) <- as.factor(c("E[r]","Gamma[d]", "Gamma[p]"))

pdf("plot/He3dp_synthetic_corr_syst.pdf",height = 6,width = 6)
ggs_pairs(Sp, 
          labeller = "label_parsed",
          diag=list(continuous = my_hist),
          upper = "blank",
          lower = list(continuous = my_bin),
          switch="y",
          showStrips=FALSE
          ) 
dev.off()









Sa <- ggs(as.mcmc(out),family="scale")

Sa$Parameter <- revalue(Sa$Parameter, c("scale[1]" = "a[1]","scale[2]" = "a[2]","scale[3]" = "a[3]"))
vlinea.dat <- data.frame(Parameter=factor(c("a[1]","a[2]","a[3]")), vl=a)


pdf("plot/He3dp_synthetic_scale_syst.pdf",height = 7,width = 6)
ggs_caterpillar(Sa) + aes(color=Parameter) +
  theme_wsj() +
  scale_color_manual(name="",values=c("gray60","gray60","gray60"))+
  geom_point(data=vlinea.dat ,aes(y=Parameter, x = vl),size=4,shape=c(8,17,11),color=c("#cd4f39","#7ebca1","#4e2f19")) +
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
  xlab("Highest Probability Interval") 

dev.off()



# Reaction rates

Nsamp <- 500
mcdat <- as.data.frame(rbind(as.mcmc(out)[,2:4][[1]],
                       as.mcmc(out)[,2:4][[2]],as.mcmc(out)[,2:4][[3]]))

index <- sample(1:nrow(mcdat),size=Nsamp,replace=FALSE)
mcdat <- mcdat[index,]
Tgrid <- exp(seq(log(1e-3),log(10),length.out =  100)) 


gdat <- list()
for(i in 1:500){
x <- Tgrid
y <- sapply(Tgrid,numrate3Hedp,ER = mcdat[i,1],gi = mcdat[i,2],gf = mcdat[i,3])
dd <- data.frame(y)
gdat[[i]] <- dd
}

gg <-  as.data.frame(gdat)
gg$x <- Tgrid

gg2<-apply(gg, 1, quantile, probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995), na.rm=TRUE) 

gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],lwr3=y[,"0.5%"],upr1=y[,"75%"],
                    upr2=y[,"97.5%"],upr3=y[,"99.5%"])

gg2data <- data.frame(x =Tgrid, mean = gg2["50%",],lwr1=gg2["25%",],
                      lwr2 = gg2["2.5%",],lwr3=gg2["0.5%",],upr1=gg2["75%",],
                      upr2=gg2["97.5%",],upr3=gg2["99.5%",])


g1 <- ggplot(gg2data,aes(x=x,y=mean))+
  theme_bw()  +
  
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr3, ymax=upr3,y=NULL), fill=c("#deebf7"),show.legend=FALSE) +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr2, ymax=upr2,y=NULL), fill=c("#9ecae1"),show.legend=FALSE) +
  geom_ribbon(data=gg2data,aes(x=Tgrid,ymin=lwr1, ymax=upr1,y=NULL), fill=c("#3182bd"),show.legend=FALSE) +
  geom_line(size=1,colour="#ffffff",linetype="dashed",size=1,show.legend=FALSE) +
  theme_wsj() + xlab("Temperature (GK)") + ylab(TeX('$N_A\\sigma v$')) + 
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
