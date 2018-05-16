
  mcdat_I <- as.data.frame(do.call(rbind, as.mcmc(Normfit )[,c("e1", "gin", "gout","ri","rf")]))
  mcdat_II <- as.data.frame(do.call(rbind, as.mcmc(Normfit )[,c("e1_2", "gin_2", "gout_2","ri_2","rf_2")]))
  
  index <- sample(1:nrow(mcdat_I ),size=Nsamp,replace=FALSE)
  mcdat_I  <- mcdat_I [index,]
  mcdat_II  <- mcdat_II [index,]
  
  Tgrid <- 10^(seq(log(1e-3,10),log(10,10),length.out =  500))
  
   gdat <- vector('list',10)
   for(i in 1:10){
    y_I <- sapply(Tgrid,nuclear_rate3Hedp_5p,ER = mcdat_I[i,1],gi = mcdat_I[i,2],gf = mcdat_I[i,3],r_i=mcdat_I[i,4],r_f=mcdat_I[i,5] )
    y_II <- sapply(Tgrid,nuclear_rate3Hedp_5p,ER = mcdat_II [i,1],gi = mcdat_II [i,2],gf = mcdat_II[i,3],r_i=mcdat_II[i,4],r_f=mcdat_II[i,5] )
    
    dd <- data.frame(y_I/y_II)
    gdat[[i]] <- dd
  }

  gg <-  as.data.frame(gdat)
  
  gg2 <- apply(gg, 1, quantile, probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995), na.rm=TRUE)
  
 
  Drate <- data.frame(T =Tgrid, mean = gg2["50%",],lwr1 = gg2["25%",],lwr2 = gg2["2.5%",],
                      lwr3 =  gg2["0.5%",],upr1=gg2["75%",],
                      upr2=gg2["97.5%",],upr3=gg2["99.5%",])
