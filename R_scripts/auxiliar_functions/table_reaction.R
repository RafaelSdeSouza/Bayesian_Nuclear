table_reaction <- function(mat, vars=vars, N = 1000){
  mcdat_I <- as.data.frame(do.call(rbind, as.mcmc(mat)[,vars]))
  index <- sample(1:nrow(mcdat_I ),size=N,replace=FALSE)
  mcdat_I  <- mcdat_I [index,]
  
  gdat <- vector('list',N)
  system.time(for(i in 1:N){
    y <- sapply(Tgrid,nuclear_rate3Hedp_5p,ER = mcdat_I[i,1],gi = mcdat_I[i,2],gf = mcdat_I[i,3],r_i=mcdat_I[i,4],r_f=mcdat_I[i,5] )
    dd <- data.frame(y)
    gdat[[i]] <- dd
  }
  )
  gg <-  as.data.frame(gdat)
  
  gg2 <- apply(gg, 1, quantile, probs=c(0.16, 0.5, 0.84), na.rm=TRUE)
  
  fu <- function(x){exp(sqrt(log(1+var(x)/mean(x)^2)))}       
  
  fu_I<-apply(gg, 1, fu)
  
  gg2data <- data.frame(x =Tgrid, mean = gg2["50%",],lower = gg2["16%",], upper = gg2["84%",] )
  gg2data$fu <- fu_I
  return(gg2data)
}