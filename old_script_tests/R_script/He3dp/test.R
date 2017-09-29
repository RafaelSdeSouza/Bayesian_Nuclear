## ARTIFICIAL DATA GENERATION 

N <- 300

obsx1 <- exp(seq(log(1e-3), log(1),length.out=N))

res <- vector()
obsy1 <- vector()
errobsy1 <- vector()

# Barker values:
# Er  = 0.0912 MeV
# g^2_in = 2.93 MeV         ! reduced width of deuteron
# g^2_out = 0.0794 MeV      ! reduced width of neutron

#res[2] <-  0.35779  # resonance energy
res[2] <-  0.35  # resonance energy
res[3] <-  1.0085    # reduced width incoming
res[4] <- 0.025425   # reduced width outgoing

for (i in 1:length(obsx1)){
  res[1] <- obsx1[i]
  write.table(res, file="He3dp_AD.in", quote=TRUE, 
              row.names=FALSE, col.names=FALSE)
  
  # Load the fortran code needed to calculate S-factor curve
  if(!is.loaded("He3dp_AD_Sub"))
    dyn.load("He3dp_AD.so") 
  .Fortran("He3dp_AD_Sub")
  
  tab1 <- read.table("He3dp_AD.out", header=FALSE)
  
  errobsy1[i] <- 0.1
  obsy1[i] <- rnorm( 1, tab1[1,2], errobsy1[i] )
}


# Test


plot(obsx1,obsy1,col="green",log="x",ylim=c(0,25),xlim=c(1e-3,1),cex=1.5)
#par(new=TRUE)
lines(obsx1,Sfactor_He3dp(obsx1,0.35,1.0085,0.025425),ylim=c(0,25),
      xlim=c(1e-3,1),col="red",cex=1.25)

