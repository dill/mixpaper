# amakihi analysis

library(mmds)

set.seed(12499)

amakihi<-read.csv(file="amakihi.csv")

amakihi$mas2<-amakihi$mas
mas.sd <- sd(amakihi$mas)
amakihi$mas<-(amakihi$mas)/mas.sd


aics<-c()

# no covariates
amakihi.fit1<-fitmix(amakihi,mix.terms=1,width=82.5,pt=TRUE)
amakihi.best<-step.ds.mixture(amakihi.fit1)
aics<-c(aics,amakihi.best$aic)

# covariates
amakihi.fit1.o<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~as.factor(obs)",pt=TRUE)
amakihi.o.best<-step.ds.mixture(amakihi.fit1.o)
aics<-c(aics,amakihi.o.best$aic)


amakihi.fit1.h<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~as.factor(has)",pt=TRUE)
amakihi.h.best<-step.ds.mixture(amakihi.fit1.h)
aics<-c(aics,amakihi.h.best$aic)


amakihi.fit1.m<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~mas",pt=TRUE)
amakihi.m.best<-step.ds.mixture(amakihi.fit1.m)
aics<-c(aics,amakihi.m.best$aic)


amakihi.fit1.oh<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~as.factor(obs)+as.factor(has)",pt=TRUE)
amakihi.oh.best<-step.ds.mixture(amakihi.fit1.oh)
aics<-c(aics,amakihi.oh.best$aic)


amakihi.fit1.om<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~as.factor(obs)+mas",pt=TRUE)
amakihi.om.best<-step.ds.mixture(amakihi.fit1.om)
aics<-c(aics,amakihi.om.best$aic)


amakihi.fit1.mh<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~mas+as.factor(has)",pt=TRUE)
amakihi.mh.best<-step.ds.mixture(amakihi.fit1.mh)
aics<-c(aics,amakihi.mh.best$aic)


amakihi.fit1.omh<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~mas+as.factor(has)+as.factor(obs)",pt=TRUE)
amakihi.omh.best<-step.ds.mixture(amakihi.fit1.omh)
aics<-c(aics,amakihi.omh.best$aic)

#############
# observer + mas is best model


# plot pdf
pdf(file="amakihi-om-pdf.pdf",width=4.2,height=3.9)
par(cex.lab=1.2)
plot(amakihi.om.best,breaks=c(seq(0,82.5,len=20)),main="PDF of distances",style="comp",pdf=TRUE)
dev.off()

# load a custom plot.ds.mixutre
switchpars <- mmds:::switchpars
getpars <- mmds:::getpars
integrate.hn <- mmds:::integrate.hn
integrate.hn.pt <- mmds:::integrate.hn.pt
keyfct.hn <- mmds:::keyfct.hn
source("plot.ds.mixture.R")

# need to monkey around to get things looking like Marques et al (2007
# figure 7. They set:
#  OBS -- fix mas=0900HST
#  mas -- fix obs="TJS"

# using http://www.timeanddate.com/worldclock/astronomy.html?n=103&month=7&year=1994&obj=sun&afl=-11&day=1
# sunrise was about 0715 => 0900HST is about 105 minutes

## 3 plots to make
pdf(file="amakihi-om-hh.pdf",width=11,height=4.7)
par(cex.lab=1.6, cex.axis=1.4,cex.main=1.4,mfrow=c(1,3))

# 1. mixtures plot

plot(amakihi.om.best,breaks=c(seq(0,82.5,len=20)),main=c("Detection function","Levels of observer","Quartiles of minutes after sunrise"),style="comp",hide.hist=TRUE,pseq=1,nomf=TRUE,lty=1,lcol="black")

# 2. obs
am2 <- amakihi.om.best
dat<-data.frame(object=1:3,
                distance=c(20,41,61),
                obs=factor(c("SGF","TJS","TKP"),levels=levels(amakihi.om.best$data$obs)),
                mas=rep(105,3)/sd(amakihi$mas2),
                observed=rep(1,3))
amakihi.om.best$data<-dat
z <- amakihi.om.best$z[[1]]
z <- z[1:3,]
z[,4] <- (rep(105,3))/sd(amakihi$mas2)
z[,2] <- c(0,1,0)
z[,3] <- c(0,0,1)
amakihi.om.best$z[[1]]<-z

plot(amakihi.om.best,breaks=c(seq(0,82.5,len=20)),main=c("Detection function","Levels of observer","Quartiles of minutes after sunrise"),style="comp",hide.hist=TRUE,pseq=2,nomf=TRUE,lty=1:3,lcol=rep("black",3))

legend(0,0.2,legend=c("SGF","TJS","TKP"),lty=1:3)


# 3. mas

amakihi.om.best <- am2
dat<-data.frame(object=1:6,
                distance=c(1,21,30,60,70,80),
                obs=factor(rep("TJS",6),levels=levels(amakihi.om.best$data$obs)),
                mas=(seq(0,300,by=60))/sd(amakihi$mas2),
                observed=rep(1,6))
amakihi.om.best$data<-dat
z <- amakihi.om.best$z[[1]]
z <- z[1:6,]
z[,4] <- (seq(0,300,by=60))/sd(amakihi$mas2)
z[,2] <- rep(1,6)
z[,3] <- rep(0,6)
amakihi.om.best$z[[1]]<-z


plot(amakihi.om.best,breaks=c(seq(0,82.5,len=20)),main=c("Detection function","Levels of observer","Minutes after sunrise 0-300"),style="comp",hide.hist=TRUE,pseq=3,nomf=TRUE,lty=rep(1,6),lcol=rep("black",6))

text(17,0.56,"mas=300",cex=1.5)
text(57,0.5,"mas=0",cex=1.5)

dev.off()


# generate results table output
# needs to be processed for delta AIC
source("grabresults.R")

grab_results(amakihi.best)
grab_results(amakihi.o.best)
grab_results(amakihi.h.best)
grab_results(amakihi.m.best)
grab_results(amakihi.oh.best)
grab_results(amakihi.om.best)
grab_results(amakihi.mh.best)
grab_results(amakihi.omh.best)

