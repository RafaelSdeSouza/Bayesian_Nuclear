
require(stats4)
require(bbmle)


sfactorTdn_fast<-function(x=x, e1=e1, gin=gin, gout=gout){
res <- vector()
y <- vector()
tdn_AD.in <- vector()
res[2] <- e1   # resonance energy
res[3] <-  gin    # reduced width incoming
res[4] <- gout   # reduced width outgoing

for (i in 1:length(x)){
  res[1] <- x[i]
  write.table(res, file="tdn_AD.in", quote=TRUE, 
              row.names=FALSE, col.names=FALSE)

  # Load the fortran code needed to calculate S-factor curve
  if(!is.loaded("tdn_AD_Sub"))
   dyn.load("tdn_AD.so") 
  .Fortran("tdn_AD_Sub")
  
  tab1 <- read.table("tdn_AD.out", header=FALSE)
  y[i] <- tab1[1,2]
}
return(y)
}



N <- 100
obsx1 <- log(runif(N, exp(0.0), exp(0.2)))
obsy1 <- sfactorTdn_fast(obsx1,e1=0.0912,gin=2.93,gout=0.0794)
plot(obsx1,obsy1)  
  
  
chi2 <- function(x, e1, gin, gout){
 sum((obsy1 -sfactorTdn_fast(x=obsx1,e1, gin, gout))^2)
}

fit<-mle2(chi2, start = list(e1 = 01,gin= 1,gout = 0.1))  




confint(fit)


