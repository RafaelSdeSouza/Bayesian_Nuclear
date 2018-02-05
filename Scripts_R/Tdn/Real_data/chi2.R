library(stats4)
## DATA SETS
tab1 <- read.table("data_jarmie.dat", header=FALSE)
obsx1 <- c(tab1[,1])
obsy1 <- c(tab1[,2])
errobsy1 <- c(tab1[,3])

tab2 <- read.table("data_brown.dat", header=FALSE)
obsx2 <- c(tab2[,1])
obsy2 <- c(tab2[,2])
errobsy2 <- c(tab2[,3])

tab3 <- read.table("data_kobzev.dat", header=FALSE)
obsx3 <- c(tab3[,1])
obsy3 <- c(tab3[,2])
errobsy3 <- c(tab3[,3])

tab4 <- read.table("data_arnold.dat", header=FALSE)
obsx4 <- c(tab4[,1])
obsy4 <- c(tab4[,2])
errobsy4 <- c(tab4[,3])


tab_all <- rbind(tab1,tab2,tab3,tab4)

x <- tab_all$V1
y <- tab_all$V2
erry <- tab_all$V3


likelihood <- function(param){
  "Likelihood function"  
  
  er = param[1]
  gi = param[2]
  gf = param[3]
  LL = sum(((y - Sfactor(x,er,gi,gf,5,3,0))^2)/erry^2) 
  return(-log(LL))
}

# Fit via Maximum Likelihood
LL <- function(er,gi,gf){
  likelihood(c(er,gi,gf))
}

# Fit
fit <- mle2(LL , start = list(er = 0.5, gi = 0.5, gf = 0.05), method = "L-BFGS-B", lower = c(0, 0,0))
summary(fit)