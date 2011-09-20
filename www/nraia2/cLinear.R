### R code from vignette source '/home/bates/Documents/books/nraia2/cLinear.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width=69,show.signif.stars=FALSE,str=strOptions(strict.width="cut"))
library(lattice)
library(ggplot2)
library(NRAIA)
#lattice.options(default.theme = function() standard.theme())
lattice.options(default.theme = function() standard.theme(color=FALSE))


###################################################
### code chunk number 2: PCBdata
###################################################
#print(xyplot(conc ~ age, PCB, type=c("g","p"),
#             xlab="Age (years)", ylab="PCB concentration (ppm)"))
p <- qplot(age, conc, data=PCB, xlab="Age (years)", ylab="PCB concentration (ppm)") + theme_bw()
print(p)


###################################################
### code chunk number 3: PCBlogs
###################################################
print(p + scale_y_log2(), TRUE, vp=viewport(0.25, 0.5, 0.5, 1))
print(p + scale_y_log2() + scale_x_sqrt(), FALSE, vp=viewport(0.75, 0.5, 0.5, 1))


