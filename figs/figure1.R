# generate figure 1....

# layout:
#  3 panels
#     1 - detection function humpback from Williams and Thomas w. hist
#     2 - average detection function pilot whale from Dan Pike w. hist
#     3 - levels of BSS covar for pilot whale from Dan Pike w.o. hist

postscript(file="figure1.eps",width=7.5,height=3,paper="special",horizontal=FALSE)

par(mfrow=c(1,3),mar=c(4,4.2,2,2))

# 1
# read in the plot data. 3 columns: histogram x coords, hist y coords, detfct
hb.plot<-read.table("humpback.txt")

# make the histogram
plot(hb.plot$V1,hb.plot$V2,type="l",axes=F,xlab="Distance (m)",
     ylab="Probability of detection",main="Humpback")
axis(2,at=c(0,0.5,1))
axis(1,at=c(0,1000,2000))
# add in the baseline and the first bar line
lines(x=c(2000,0,0),y=c(0,0,0.39636))

lines(x=hb.plot$V1,y=hb.plot$V3)

# 2
dp.df<-read.table("danpike-detfct.txt")

plot(dp.df$V1,dp.df$V2,type="l",axes=F,xlab="Distance (m)",
     ylab="Probability of detection",main="Long-finned pilot whale")
axis(2,at=c(0,0.5,1))
axis(1,at=c(0,1000,2000,3000))
# add in the baseline and the first bar line
lines(x=c(3000,0,0),y=c(0,0,1.1197))

lines(x=dp.df$V1,y=dp.df$V3)

# 3
# quantiles at v2, v3, v4  are BSS=1.5,2.0,3.0 
dp.df<-read.table("danpike-covar.txt")

plot(dp.df$V1,dp.df$V2,type="l",axes=F,xlab="Distance (m)",ylim=c(0,1),
     ylab="Probability of detection",main="Long-finned pilot whale\nQuantiles of Beaufort sea state",col=grey(0.25))
axis(2,at=c(0,0.5,1))
axis(1,at=c(0,1000,2000,3000))

lines(x=dp.df$V1,y=dp.df$V3,col=grey(0.5))
lines(x=dp.df$V1,y=dp.df$V4,col=grey(0.75))

dev.off()

