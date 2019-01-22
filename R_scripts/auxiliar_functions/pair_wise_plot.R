require(viridis)
pair_wise_plot <- function(s){

my_hist <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
   
   stat_density_ridges(aes(y=0,fill=factor(..quantile..),alpha=factor(..quantile..)),geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                        quantiles = c(0.025, 0.16,  0.84, 0.975)) + 
 #   stat_density_ridges(geom = "density_ridges_gradient",aes(y=0, color="red",fill = 0.5-abs(0.5 - ..ecdf..),
 #                       alpha=factor(..quantile..)),
 #                       calc_ecdf = TRUE, quantile_lines = TRUE,
 #                       quantiles = c(0.025, 0.25,  0.75, 0.975)) +
#    scale_fill_viridis() +    
   scale_fill_manual(name = "Probability", values = c("#440154BF",  "#21908CBF","#FDE725BF",
                                                      "#21908CBF","#440154BF")) +
#    scale_fill_viridis(discrete=TRUE)+
#     geom_histogram(bins = 10,fill="#bf812d",colour="#1F3552",...) +
#    geom_density(adjust = 1.5,fill="#bf812d",colour="#1F3552",...) +

      theme_void() + theme( panel.grid.minor=element_blank(),
                            panel.grid.major = element_blank())
  }
  
my_bin <- function(data, mapping, ..., low = "#74dbef", high = "#0074e4") {
  # get the x and y data to use the other code
 
#   data2 <- probs_CI(data=data)
   ggplot(data = data, mapping = mapping) +
   stat_density_2d(aes(fill =..level..,alpha =..level..),n=500,
                   geom="polygon",...) +
 #   scale_alpha_continuous(range=0.4,1) +
#    stat_ellipse(type = "norm",level = 0.95,geom = "polygon",fill="#feb24c",...) +
#    stat_ellipse(type = "norm",level = 0.68,geom = "polygon",fill="#f03b20",...) +
#    geom_point() +
#   stat_ellipse(type = "t",level = 0.997,alpha=0.95,color="#dfc27d",...) +
  stat_ellipse(type = "t",level = 0.95,color="black",linetype="dashed",alpha=0.5,...) +
  stat_ellipse(type = "t",level = 0.68,color="black",linetype="dashed",alpha=0.5,...) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
     scale_fill_viridis() +
    scale_alpha(range=c(0.2,1)) +
#    scale_fill_gradient(high=high,low=low) +
    theme_bw() +
    theme( text = element_text(size=60),axis.text = element_text(size=9.5),
           axis.text.x = element_text(size=8, angle=30),
           axis.text.y = element_text(size=8),
           strip.text.x = element_text(colour = 'red',size = 20),
           strip.text.y = element_text(colour = 'red',size = 20),
           panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
}


 ggs_pairs(s, 
           labeller = "label_parsed",
           diag=list(continuous = my_hist),
           upper = "blank",
           lower = list(continuous = my_bin),
           switch="both",
           showStrips=TRUE
 ) 
 
 }