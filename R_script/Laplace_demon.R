library(LaplacesDemon)

N <- 100
J <- 3 #Number of predictors, including the intercept
X <- matrix(1,N,J)
for (j in 2:J) {X[,j] <- rnorm(N,runif(1,-3,3),runif(1,0.1,1))}
beta.orig <- runif(J,-3,3)
e <- rnorm(N,0,0.1)
y <- as.vector(tcrossprod(beta.orig, X) + e)

mon.names <- c("LP","sigma")
parm.names <- c("ER","gi","gf") 
PGF <- function(Data) return(c(rnormv(Data$J,0,0.01),
                               log(rhalfcauchy(1,1))))
MyData <- list(J=J, PGF=PGF, N=N, mon.names=mon.names,
               parm.names=parm.names) #Notice that X and y are not included here
write.table(cbind(y,X), "X.csv", sep=",", row.names=FALSE,
            col.names=FALSE)
Model <- function(parm, Data)
{
  ### Parameters
  ER <- parm[1]
  gi <- parm[2]
  gf <- parm[3]
  sigma <- exp(parm[4])
  ### Log(Prior Densities)
  ER.prior <- sum(dnorm(beta, 0, 1000, log=TRUE))
  ER.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
  ER.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
  
  sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
  ### Log-Likelihood
  LL <- BigData(file="X.csv", nrow=Data$N, ncol=Data$J+1, size=1000,
                Method="add", CPUs=1, Type="PSOCK",
                FUN=function(x, beta, sigma) sum(dnorm(x[,1], tcrossprod(x[,-1],
                                                                         t(beta)), sigma, log=TRUE)), beta, sigma)
  ### Log-Posterior
  LP <- LL + beta.prior + sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma),
                   yhat=0,#rnorm(length(mu), mu, sigma),
                   parm=parm)
  return(Modelout)
}


Initial.Values <- GIV(Model, MyData, PGF=TRUE)
########################  Laplace Approximation  ##########################
Fit1 <- LaplaceApproximation(Model, Initial.Values, Data=MyData,
                             Iterations=10000)
Fit1

