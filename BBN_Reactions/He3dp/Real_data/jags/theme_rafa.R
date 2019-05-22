theme_rafa <- function () { 
  theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(colour = "white", fill = "white"),
          plot.background = element_rect(colour = "white", fill = "white"),
          panel.background = element_rect(size=0.75, fill = "white"),
          legend.key = element_rect(colour = "white", fill = "white"),
          axis.title = element_text(size=18.5),
          axis.text  = element_text(size=13),
          axis.ticks = element_line(size = 0.75),
          axis.line = element_line(size = 0.5, linetype = "solid"))
}