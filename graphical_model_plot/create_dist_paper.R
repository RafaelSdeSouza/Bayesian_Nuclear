source("plot_dist.R")
plot_dist(dists$normal, labels=expression(S[ij]^{lat}, sigma[i]^2),color="red",fill="green",plot_dist_name=F)

plot_dist(dists$log_normal, labels=expression(S[ij]^{lat}, sigma[i]^2),color="red",fill="green",plot_dist_name=F)