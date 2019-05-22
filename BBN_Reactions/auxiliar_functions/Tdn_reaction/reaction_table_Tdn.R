reaction_table_Tdn <- function(mat, vars=vars, N = 1000, T9){
#  mcdat_I <- as.data.frame(as.mcmc(mcmcChain)[,vars])
  mcdat_I <- mat[,vars]
  index <- sample(1:nrow(mcdat_I ),size=N,replace=FALSE)
  mcdat_I  <- mcdat_I[index,]
  
 gdat <- vector('list',N)
  system.time(for(i in 1:N){
    y <- sapply(T9,nuclear_rate_Tdn,e0 = mcdat_I[i,1], er = mcdat_I[i,2],gi = mcdat_I[i,3],gf = mcdat_I[i,4],ri=mcdat_I[i,5],
            rf=mcdat_I[i,6], ue = mcdat_I[i,7] )
    dd <- data.frame(y)
    gdat[[i]] <- dd
  }
  )
  gg <-  as.data.frame(gdat)
  
  gg2 <- apply(gg, 1, quantile, probs=c(0.16, 0.5, 0.84), na.rm=TRUE)
  
  fu <- function(x){exp(sqrt(log(1+var(x)/mean(x)^2)))}       
  
  fu_I<-apply(gg, 1, fu)
  
  gg2data <- data.frame(T9 =T9, lower = gg2["16%",],mean = gg2["50%",], upper = gg2["84%",] )
  gg2data$fu <- fu_I
  return(gg2data)
}