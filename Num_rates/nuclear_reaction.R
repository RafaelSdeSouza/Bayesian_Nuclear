require(cubature)

fd <- read.table("fort_dat.dat",header = FALSE)


numrate3Hedp <- function(ER,gi,gf,T9){
# Constants  
M0 = 3.01493216; M1 <- 2.01355332  # mass projectile in u
Z0 = 2; Z1 = 1  
  
#   DEFINITIONS  
mue <- (M0*M1)/(M0 + M1)
dpieta <- function(E){0.98951013*Z0*Z1*sqrt(mue/E)}

# Integrand
integrand <- function(E,T9) {exp(-dpieta(E))*sfactor3Hedp(E,ER,gi,gf)*exp(-E/(0.086173324*T9))}

Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 0, upper = Inf,
                      abs.tol=1e-10,subdivisions = 1000,T9 = Temp)$value}

out <- Nasv(T9)
return(out)
}


T_h <- sapply(fd$V1,numrate3Hedp,ER=0.48000449,gi=1.16116126,gf= 0.03976223)



Rd <- data.frame(fd$V1,T_h)

ggplot(fd,aes(x=V1,y=V2)) +
  geom_point() +
  geom_line(data=Rd,aes(x=fd.V1,y=T_h)) + xlab("Temperature") +
  ylab("NA") 

# 0.48000449 1.16116126 0.03976223
# 0.35,1.0085,0.025425
numrate3Hedp(0.35,1.0085,0.025425,0.8175)


test <- read.table("numRates.dat",header=TRUE)
test$T9 <- 2

sapply(test[1:10,],numrate3Hedp)
numrate3Hedp(test[1:10,1],test[1:10,2],test[1:10,3],test[1:10,4])




numrate3Hedp(0.35,1.0085,0.025425,Tgrid)

Tgrid <- exp(seq(log(1e-3),log(10),length.out =  300)) 



int <- function(Temp){3.7318e10/Temp^{3/2}*sqrt(1/mue)*adaptIntegrate(integrand, lower = 0, upper = Inf,T9=Temp)$integral}
T_h<-sapply(Tgrid,int)
plot(Tgrid,T_h)


Na <- function(T){
3.7318e10/T9^{3/2}*sqrt(mr)*integrate(sfactor3Hedp(E,0.35,1.0085,0.025425)*exp(-11.605*E/T9))  
}





if(!is.loaded("numRates_Sub"))
  dyn.load("numRates.so") 
.Fortran("numRates_Sub")



dyn.load("numRates.so") 
.Fortran("numRates_Sub")

testFn0 <- function(x) {
  prod(cos(x))
}

adaptIntegrate(testFn0, 0, 1, tol=1e-4)