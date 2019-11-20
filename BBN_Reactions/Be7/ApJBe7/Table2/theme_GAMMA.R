theme_GAMMA  <- function() {
  theme_economist_white()  +
    theme(legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        plot.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        axis.title = element_text(color="black", size=16),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x  = element_text(vjust = -1),
        strip.text = element_text(size=16),
        strip.background = element_rect("gray85"),
        axis.ticks = element_line(size = 0.85),
        axis.line = element_line(size = 0.5, linetype = "solid")) 
}
