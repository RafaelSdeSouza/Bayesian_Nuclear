# Numrates benckmark
require(nuclear)
require(ggplot2)
require(parallel)
Cdat <- read.table("numRates1.out",header=F)

Tgrid <- exp(seq(log(1e-3),log(10),length.out =  100))


ptm0 <- proc.time()
Nasv <- sapply(Tgrid ,nuclear_rate3Hedp,ER=0.48000449,gi=1.16116126,gf= 0.03976223)
proc.time() - ptm0

ptm <- proc.time()
pna <-  simplify2array(mclapply(Tgrid, nuclear_rate3Hedp,ER=0.48000449,gi=1.16116126,gf= 0.03976223, mc.preschedule = FALSE,
                        mc.set.seed = FALSE))
proc.time() - ptm

Rd <- data.frame(x=Tgrid ,y=Nasv)


ggplot(Cdat,aes(x=V1,y=V2)) +
  geom_point() +
  geom_line(data=Rd,aes(x=x,y=y)) + xlab("Temperature") +
  ylab("NA") + scale_y_log10() +
  scale_x_log10()
