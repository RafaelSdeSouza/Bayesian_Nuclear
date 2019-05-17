######################################################################
# plotpdf.R
#
# input: numRates_samp.out are the rate samples created by numRates
#        for a given temperature
#
# the samples are plotted as a histogram, together with the lognormal
# approximation of the rate; the lognormal parameters were directly
# obtained from the rate samples [as the mean and standard deviation
# of ln(x)]
#
######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
require(pracma)
######################################################################
# input
data <- read.table("numRates_samp.out", header=FALSE)
temp <- unname(unlist(data[1,]))
logmu <- unname(unlist(data[2,]))
logsigma <- unname(unlist(data[3,]))

data1 <- data$V1[-c(1:3)]
data2 <- data$V2[-c(1:3)]
data3 <- data$V3[-c(1:3)]
data4 <- data$V4[-c(1:3)]
data5 <- data$V5[-c(1:3)]
data6 <- data$V6[-c(1:3)]

# plot histogram of rate samples
pdf("plotpdf.pdf")
par(mfcol=c(3,2), mar=c(3.5,3.5,0.5,0.5), oma=c(0.2,0.2,0.2,0.2), tck=0.05, 
    las=1, mgp=c(1.8,0.2,0))

# for linear x-axis scale, take out log="x"

# first
histRes <- hist(data1, plot=FALSE, breaks=200)
xvals <- histRes$breaks
yvals <- histRes$counts
xvals <- c(xvals,xvals[length(xvals)])
yvals <- c(0,yvals,0)
yvals <- yvals/trapz(xvals,yvals)       # normalizes by area
plot(xvals, yvals, type="S", col="red", xlab="Reaction rate", ylab="", log="x")
# don't touch the 'x'
curve(dlnorm(x, logmu[1], logsigma[1]), add=TRUE, col="blue", lwd=2)
legend("topright", c(paste("T9=", as.character(temp[1]), sep=""), 
                     paste("mu=", as.character(logmu[1]), sep=""), 
                     paste("sig=", as.character(logsigma[1]), sep=""))) 

# third
histRes <- hist(data3, plot=FALSE, breaks=490)
xvals <- histRes$breaks
yvals <- histRes$counts
xvals <- c(xvals,xvals[length(xvals)])
yvals <- c(0,yvals,0)
yvals <- yvals/trapz(xvals,yvals)       # normalizes by area
plot(xvals, yvals, type="S", col="red", xlab="Reaction rate", ylab="", log="x")
# don't touch the 'x'
curve(dlnorm(x, logmu[3], logsigma[3]), add=TRUE, col="blue", lwd=2)
legend("topright", c(paste("T9=", as.character(temp[3]), sep=""), 
                     paste("mu=", as.character(logmu[3]), sep=""), 
                     paste("sig=", as.character(logsigma[3]), sep=""))) 

# fifth
histRes <- hist(data5, plot=FALSE, breaks=490)
xvals <- histRes$breaks
yvals <- histRes$counts
xvals <- c(xvals,xvals[length(xvals)])
yvals <- c(0,yvals,0)
yvals <- yvals/trapz(xvals,yvals)       # normalizes by area
plot(xvals, yvals, type="S", col="red", xlab="Reaction rate", ylab="", log="x")
# don't touch the 'x'
curve(dlnorm(x, logmu[5], logsigma[5]), add=TRUE, col="blue", lwd=2)
legend("topright", c(paste("T9=", as.character(temp[5]), sep=""), 
                     paste("mu=", as.character(logmu[5]), sep=""), 
                     paste("sig=", as.character(logsigma[5]), sep=""))) 

# second
histRes <- hist(data2, plot=FALSE, breaks=490)
xvals <- histRes$breaks
yvals <- histRes$counts
xvals <- c(xvals,xvals[length(xvals)])
yvals <- c(0,yvals,0)
yvals <- yvals/trapz(xvals,yvals)       # normalizes by area
plot(xvals, yvals, type="S", col="red", xlab="Reaction rate", ylab="", log="x")
# don't touch the 'x'
curve(dlnorm(x, logmu[2], logsigma[2]), add=TRUE, col="blue", lwd=2)
legend("topright", c(paste("T9=", as.character(temp[2]), sep=""), 
                     paste("mu=", as.character(logmu[2]), sep=""), 
                     paste("sig=", as.character(logsigma[2]), sep=""))) 

# fourth
histRes <- hist(data4, plot=FALSE, breaks=490)
xvals <- histRes$breaks
yvals <- histRes$counts
xvals <- c(xvals,xvals[length(xvals)])
yvals <- c(0,yvals,0)
yvals <- yvals/trapz(xvals,yvals)       # normalizes by area
plot(xvals, yvals, type="S", col="red", xlab="Reaction rate", ylab="", log="x")
# don't touch the 'x'
curve(dlnorm(x, logmu[4], logsigma[4]), add=TRUE, col="blue", lwd=2)
legend("topright", c(paste("T9=", as.character(temp[4]), sep=""), 
                     paste("mu=", as.character(logmu[4]), sep=""), 
                     paste("sig=", as.character(logsigma[4]), sep=""))) 

# sixth
histRes <- hist(data6, plot=FALSE, breaks=490)
xvals <- histRes$breaks
yvals <- histRes$counts
xvals <- c(xvals,xvals[length(xvals)])
yvals <- c(0,yvals,0)
yvals <- yvals/trapz(xvals,yvals)       # normalizes by area
plot(xvals, yvals, type="S", col="red", xlab="Reaction rate", ylab="", log="x")
# don't touch the 'x'
curve(dlnorm(x, logmu[6], logsigma[6]), add=TRUE, col="blue", lwd=2)
legend("topright", c(paste("T9=", as.character(temp[6]), sep=""), 
                     paste("mu=", as.character(logmu[6]), sep=""), 
                     paste("sig=", as.character(logsigma[6]), sep=""))) 

### single plot
#histRes <- hist(data, plot=FALSE, breaks=1000)
#xvals <- histRes$breaks
#yvals <- histRes$counts
#xvals <- c(xvals,xvals[length(xvals)])
#yvals <- c(0,yvals,0)
#yvals <- yvals/trapz(xvals,yvals)       # normalizes by area
#plot(xvals, yvals, type="S", col="red", xlab="Reaction rate", ylab="Probability")
# don't touch the 'x'
#curve(dlnorm(x, param[1], param[2]), add=TRUE, col="blue", lwd=2)
            
dev.off()