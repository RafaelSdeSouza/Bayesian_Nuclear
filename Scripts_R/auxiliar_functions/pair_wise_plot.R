pair_wise_plot <- function(s){
 
my_hist <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
#     geom_histogram(bins = 10,fill="#bf812d",colour="#1F3552",...) +
    geom_density(adjust = 1.5,fill="#bf812d",colour="#1F3552",...) +

      theme_void() + theme( panel.grid.minor=element_blank(),
                            panel.grid.major = element_blank())
  }
  
my_bin <- function(data, mapping, ..., low = "gray80", high = "gray25") {
  # get the x and y data to use the other code
  ggplot(data = data, mapping = mapping) +
    stat_density_2d(aes(fill =..level..), geom="polygon",...) +
#    stat_ellipse(type = "norm",level = 0.95,geom = "polygon",fill="#feb24c",...) +
#    stat_ellipse(type = "norm",level = 0.68,geom = "polygon",fill="#f03b20",...) +
#    geom_bin2d(bins = 50,...) +
#   stat_ellipse(type = "t",level = 0.997,alpha=0.95,color="#dfc27d",...) +
    stat_ellipse(type = "t",level = 0.95,alpha=0.95,color="#bf812d",...) +
   stat_ellipse(type = "t",level = 0.68,alpha=0.95,color="#8c510a",...) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    
    scale_fill_gradient(low = low, high = high) +
    theme_bw() +
    theme( text = element_text(size=60),axis.text = element_text(size=7.5),
           strip.text.x = element_text(colour = 'red',size = 20),
           panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
}


 ggs_pairs(s, 
           labeller = "label_parsed",
           diag=list(continuous = my_hist),
           upper = "blank",
           lower = list(continuous = my_bin),
           switch="both",
           showStrips=FALSE
 ) }