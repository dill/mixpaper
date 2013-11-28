# amakihi analysis

library(mmds)

set.seed(12499)

amakihi<-read.csv(file="amakihi.csv")

amakihi$mas2<-amakihi$mas
amakihi$mas<-(amakihi$mas-mean(amakihi$mas))/sd(amakihi$mas)

aics<-c()

### no covariates
##amakihi.fit1<-fitmix(amakihi,mix.terms=1,width=82.5,pt=TRUE)
##amakihi.best<-step.ds.mixture(amakihi.fit1)
##aics<-c(aics,amakihi.best$aic)
##
### covariates
##amakihi.fit1.o<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~as.factor(obs)",pt=TRUE)
##amakihi.o.best<-step.ds.mixture(amakihi.fit1.o)
##aics<-c(aics,amakihi.o.best$aic)
##
##
##amakihi.fit1.h<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~as.factor(has)",pt=TRUE)
##amakihi.h.best<-step.ds.mixture(amakihi.fit1.h)
##aics<-c(aics,amakihi.h.best$aic)
##
##
##amakihi.fit1.m<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~mas",pt=TRUE)
##amakihi.m.best<-step.ds.mixture(amakihi.fit1.m)
##aics<-c(aics,amakihi.m.best$aic)
##
##
##amakihi.fit1.oh<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~as.factor(obs)+as.factor(has)",pt=TRUE)
##amakihi.oh.best<-step.ds.mixture(amakihi.fit1.oh)
##aics<-c(aics,amakihi.oh.best$aic)
##
##
amakihi.fit1.om<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~as.factor(obs)+mas",pt=TRUE)
amakihi.om.best<-step.ds.mixture(amakihi.fit1.om)
aics<-c(aics,amakihi.om.best$aic)
##
##
##amakihi.fit1.mh<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~mas+as.factor(has)",pt=TRUE)
##amakihi.mh.best<-step.ds.mixture(amakihi.fit1.mh)
##aics<-c(aics,amakihi.mh.best$aic)
##
##
##amakihi.fit1.omh<-fitmix(amakihi,mix.terms=1,width=82.5,model.formula="~mas+as.factor(has)+as.factor(obs)",pt=TRUE)
##amakihi.omh.best<-step.ds.mixture(amakihi.fit1.omh)
##aics<-c(aics,amakihi.omh.best$aic)

#############

##pdf(file="amakihi-om.pdf",width=10.7,height=4.7)
##par(cex.lab=1.6, cex.axis=1.4,cex.main=1.4)
##plot(amakihi.om.best,breaks=c(seq(0,82.5,len=20)),main=c("Detection function","Levels of observer","Quantiles of minutes after sunrise"),style="comp")
##dev.off()


##pdf(file="amakihi-om-hh.pdf",width=10.7,height=4.7)
switchpars <- mmds:::switchpars
getpars <- mmds:::getpars
integrate.hn <- mmds:::integrate.hn
integrate.hn.pt <- mmds:::integrate.hn.pt
keyfct.hn <- mmds:::keyfct.hn
source("~/current/mmds/R/plot.ds.mixture.R")
par(cex.lab=1.6, cex.axis=1.4,cex.main=1.4)
plot(amakihi.om.best,breaks=c(seq(0,82.5,len=20)),main=c("Detection function","Levels of observer","Quantiles of minutes after sunrise"),style="comp",hide.hist=TRUE)
##dev.off()


##pdf(file="amakihi-om-pdf.pdf",width=4.2,height=3.9)
##par(cex.lab=1.2)
##plot(amakihi.om.best,breaks=c(seq(0,82.5,len=20)),main="PDF of distances",style="comp",pdf=TRUE)
##dev.off()


## to plot the last bit need to edit the code and then run this
## change line 300 plot.ds.mixture.R to read
## this.col.levels<-unique(this.col)
#
## replace the data and z matrix
#dat<-data.frame(object=1:6,
#                distance=c(1,21,30,60,70,80)
#                obs="TJS",
#                mas=(seq(0,300,len=6)-mean(amakihi$mas2))/sd(amakihi$mas2),
#                observed=rep(1,6))
#amakihi.om.best$data<-dat
#z<-amakihi.om.best$z[[1]][1:6,]
#z[,4]<-(seq(0,300,len=6)-mean(amakihi$mas2))/sd(amakihi$mas2)
#amakihi.om.best$z[[1]]<-z
#source("~/current/mmds/R/plot.ds.mixture.R")
#plot(amakihi.om.best,main=c("Detection function","","Minutes after sunrise"),style="comp",hide.hist=TRUE)
#text(0.92,0.56,"mas=0")
#text(0.84,0.48,"mas=300")
#dev.copy2eps(file="amakihi-om-pdf-again.eps",width=4.2,height=3.9)

###source("grabresults.R")
###
###grab_results(amakihi.best)
###grab_results(amakihi.o.best)
###grab_results(amakihi.h.best)
###grab_results(amakihi.m.best)
###grab_results(amakihi.oh.best)
###grab_results(amakihi.om.best)
###grab_results(amakihi.mh.best)
###grab_results(amakihi.omh.best)

