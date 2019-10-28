# Summary Statistics

samp <- read.csv("MCMC_ApJ_ultimaterun.csv",header = T) 
en <- samp[,c('e0_1','e0_2',
              'e0_3', 'e0_4',
              'e0_5', 'e0_6', 
              'e0_7')]

gan <- samp[,c('ga_1','ga_2', 'ga_3', 'ga_4', 'ga_5', 'ga_6', 'ga_7')]
gbn <- samp[,c('gb_1','gb_2', 'gb_3', 'gb_4', 'gb_5', 'gb_6', 'gb_7')]
norm <- samp[,c("y.norm.1.","y.norm.2.",
                 "y.norm.3.","y.norm.4.",
                 "y.norm.5.","y.norm.6.",
                 "y.norm.7.",
                 "y.norm.8.",
                 "y.norm.9.","y.norm.10.")]
sigscat <- samp[,c("y.scat.1.","y.scat.2.",
                   "y.scat.3.","y.scat.4.",
                   "y.scat.5." ,"y.scat.6.",
                   "y.scat.7.",
                   "y.scat.8.",
                   "y.scat.9.","y.scat.10.")]


 probBe7 <- function(x){
 quantile(x,probs=c(0.16, 0.5, 0.84))
}


tabEr <- apply(en,2,probBe7)
 tabga <- apply(gan,2,probBe7)
 tabgb <- apply(gbn,2,probBe7)
 
 formtab <- function(tab){
 low <-  tab[2,] - tab[1,]
 up <- tab[3,] - tab[2,]
 mean <- tab[2,]
 out <- data.frame(low,mean,up)
 return(out)
 }
 
 
 round(formtab(tabEr),5)
 round(formtab( tabga ),5)
 round(formtab( tabgb ),5)
 round(formtab(apply(norm,2,probBe7)),2)
 round(formtab(apply(sigscat,2,probBe7)),2)
 
tab2 <- as.data.frame(tab2)
tab2$low <- tab2[,4] - tab2[,3]
tab2$hi <-  tab2[,5] - tab2[,4]