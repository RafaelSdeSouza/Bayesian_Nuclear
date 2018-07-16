plotAreaInROPE <- 
  function(paramSampleVec, credMass = 0.95, compVal = 0, maxROPEradius,
           n = 250,...) {
    # Plots the probability mass included in the ROPE as a function of
    #   the half-width of the ROPE.
    
    # Sanity checks:
    if(missing(maxROPEradius))
      stop("maxROPEradius is missing with no default.")
    if(!isTRUE(maxROPEradius > 0))
      stop("maxROPEradius must be > 0.")
    if(!isTRUE(is.finite(compVal)))
      stop("A finite value for compVal is needed.")
    
    ropeRadVec = seq( 0 , maxROPEradius , length = n ) # arbitrary comb
    areaInRope = rep( NA , n )
    for ( rIdx in 1:n ) {
      areaInRope[rIdx] <- mean( paramSampleVec > (compVal-ropeRadVec[rIdx])
                                & paramSampleVec < (compVal+ropeRadVec[rIdx]) )
    }
    
  
      defaultArgs <- list(xlab=bquote("Radius of ROPE around "*.(compVal)),
                          ylab = expression(paste(E[0]," Posterior in ROPE")))
      
      HDIlim = hdi( paramSampleVec , credMass=credMass )
      farHDIlim = HDIlim[which.max(abs(HDIlim-compVal))]
      ropeRadHDI = abs(compVal-farHDIlim)
      areaInFarHDIlim <- mean(paramSampleVec > (compVal-ropeRadHDI)
                               & paramSampleVec < (compVal+ropeRadHDI) )
      
    dd <- data.frame(x=ropeRadVec, y = areaInRope)
    dd2 <- filter(dd,areaInRope <= 0.95)
    y2 <- 0.95
    x2 <- max(dd2$x)
    
    gg <- ggplot(dd, aes(x=x,y=y)) + geom_line() +
      geom_segment(x=ropeRadHDI,xend= ropeRadHDI,y=0,yend=areaInFarHDIlim,
                   linetype="dashed") +

        annotate(geom="text", x = 1.2*ropeRadHDI , y = 0.1, size=5,
                 label = paste(100*credMass,"% HDI limit \n",
                               "farthest from ",compVal))  +
      
      geom_segment(x=-0.1,xend = ropeRadHDI,y=areaInFarHDIlim,yend=areaInFarHDIlim,
                   linetype="dashed") +
      
      annotate(geom="text", x = 0.02 , y = 0.9*areaInFarHDIlim, size=5,
               label = paste(round(100*areaInFarHDIlim,1),"% \n","probability"))  +
      
        geom_ribbon(data=subset(dd , x < ropeRadHDI),aes(ymax=y),ymin=0,
                    fill="#756bb1",colour=NA,alpha=0.6) +
      geom_ribbon(data=subset(dd , x < 0.01),aes(ymax=y),ymin=0,
                  fill="skyblue",colour=NA) +
      annotate(geom="text", x = 0.01, y = 0.075, size=5,
               label = paste(4.5,"%"))  +
      geom_segment(x=-0.1,xend = 0.01,y=0.045,yend=0.045,
                   linetype="dashed") +
        theme_bw() +
        coord_cartesian(ylim=c(0,1)) +
        xlab(defaultArgs$xlab) + ylab(defaultArgs$ylab) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              legend.background = element_rect(colour = "white", fill = "white"),
              #       plot.background = element_rect(size = 3, linetype='dashed',colour = "white", fill = "white"),
              panel.background = element_rect(colour = "white", fill = "white"),
              legend.key = element_rect(colour = "white", fill = "white"),
              axis.title = element_text(size=22),
              axis.text  = element_text(size=18),
              axis.ticks = element_line(size = 0.45),
              #        axis.line = element_line(size = 0.45, linetype = "solid"),
              axis.text.y = element_text(size = 20, margin = unit(c(t = 0, r = 5, b = 0, l = 0), "mm")),
              axis.text.x = element_text(size = 20, margin = unit(c(t = 5, r = 0, b = 0, l = 0), "mm")),
              axis.ticks.length = unit(-3, "mm"))

           

    return(gg)
    }
      
        