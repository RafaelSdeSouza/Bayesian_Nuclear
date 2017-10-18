
colpal <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")

source("plot_dist.R")

png("Normal.png",width = 1.25*180, height = 1.25*60)
plot_dist(dists$normal, labels = expression(S[ij]^{lat}, sigma[i]^2),color=colpal[1],plot_dist_name=F)
dev.off()

png("Normal2.png",width = 1.25*180, height = 1.25*60)
plot_dist(dists$normal, labels = expression(kappa[j]*S[ij]^{th}, tau^{2}),color=colpal[1],plot_dist_name=F)
dev.off()

png("logNormal.png",width = 1.25*180, height = 1.25*60)
plot_dist(dists$log_normal, labels = expression(mu,epsilon[j]^2),color=colpal[3],plot_dist_name=F)
dev.off()

png("Gamma.png",width = 1.25*180, height = 1.25*60)
plot_dist(dists$gamma, labels = expression(alpha[1]~beta[1]), color=colpal[4],plot_dist_name=F)
dev.off()

for(i in 1:3){
png(paste("Uniform",i,".png",sep =""),width = 1.25*180, height = 1.25*60)
plot_dist(dists$uniform, labels = expression(a,b), color=colpal[5],plot_dist_name=F)
dev.off()
}