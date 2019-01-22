require(GA)

x <- runif(100, 0.01, 0.3)
sd <- 0.2
y <- rnorm( 100, sfactor3Hedp_5p(x,  0.35779,1.0085,0.025425,5,6), sd = sd) 

likelihood <- function(param){
  "Likelihood function"  
  
  er = param[1]
  gi = param[2]
  gf = param[3]
  ri = param[4]
  rf = param[5]
  xb = sfactor3Hedp_5p(x,er,gi,gf,ri,rf)
  LL = dnorm(y,mean=xb, sd = sd,log=T) 
  return(sum(LL))
}

x1 <-  seq(0.01, 3, by = 0.01)
x2 <-  seq(0.01, 1, by = 0.01)
x3 <-  seq(0.01,1,by = 0.01)

slopevalues <- function(x1,x2){return(likelihood(c(0.357792,x1,x2,5,6)))}
slopevalues <- Vectorize(slopevalues)

z <- outer(x1, x2, slopevalues)

persp3D(x1, x2, z, theta = 50, phi = 20,
        xlab="gi",ylab="gf",zlab="Likelihood")

filled.contour(x1, x2, z, color.palette=colorRampPalette(c('white','blue','yellow','red','darkred')),nlevels = 100)


# Now for Er
slopevalues <- function(x3,x1){return(likelihood(c(x3,x1,0.025425,5,6)))}
slopevalues <- Vectorize(slopevalues)

z <- outer(x3, x1, slopevalues)

persp3D(x3, x1, z, theta = 50, phi = 20,
        xlab="gi",ylab="gf",zlab="Likelihood")
filled.contour(x3, x1, z, color.palette=colorRampPalette(c('white','blue','yellow','red','darkred')),nlevels = 100)



GA <- ga(type = "real-valued", fitness =  function(x) -slopevalues(x[1],x[2]),
         min = c(0.0001, 0.0001), max = c(10, 10),
         popSize = 100, maxiter = 1000)
summary(GA)
plot(GA)