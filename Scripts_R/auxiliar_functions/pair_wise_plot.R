pair_wise_plot <- function(s){
  
my_hist <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_histogram(bins = 10,fill="#b2182b",colour="#1F3552",...) +
      theme_void() + theme( panel.grid.minor=element_blank(),
                            panel.grid.major=element_blank())
  }
  
my_bin <- function(data, mapping, ..., low = "#d1e5f0", high = "#b2182b") {
  ggplot(data = data, mapping = mapping) +
    #    geom_bin2d(...) +
    stat_density2d(aes(fill=..level..,alpha=..level..), geom="polygon",...) +
    scale_fill_gradient(low = low, high = high) +
    theme_bw()
}


 ggs_pairs(s, 
           labeller = "label_parsed",
           diag=list(continuous = my_hist),
           upper = "blank",
           lower = list(continuous = my_bin),
           switch="both",
           showStrips=FALSE
 ) }