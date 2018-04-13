require(GA)
SF<- function(x1, x2){
sfactorTdn(ER=0.0912,gi=x1,gf=x2)
}

x <- runif(100, 0.01, 0.3)
erry <- 1

y <- rnorm( 100, sfactorTdn(x,  0.0912, 2.93 , 0.0794), erry^2) 

likelihood <- function(param){
  "Likelihood function"  
  
  er = param[1]
  gi = param[2]
  gf = param[3]
  LL = -sum(dnorm(y, mean=sfactorTdn(x,er,gi,gf),sd=sigma,log=T)) 
  return(-LL)
}

LK2 <- function(x1,x2) {return(
  likelihood(c(0.0912,x1,x2)))
  }



x1 <-  seq(0.001, 3.5, by = 0.01)
x2 <-  seq(0.001, 0.2, by = 0.01)

slopevalues <- function(x1,x2){return(likelihood(c(0.0912,x1,x2)))}
slopevalues <- Vectorize(slopevalues)

z <- outer(x1, x2, slopevalues)

persp3D(x1, x2, z, theta = 60, phi = 15)

filled.contour(x1, x2, z, color.palette=colorRampPalette(c('white','blue','yellow','red','darkred')),nlevels = 100)


GA <- ga(type = "real-valued", fitness =  function(x) -slopevalues(x[1],x[2]),
         min = c(0.0001, 0.0001), max = c(10, 10),
         popSize = 100, maxiter = 1000)
summary(GA)
plot(GA)