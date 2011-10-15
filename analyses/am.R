# amakihi analysis

library(mmds)

set.seed(12499)

amakihi<-read.csv(file="amakihi.csv")

amakihi$mas<-(amakihi$mas-mean(amakihi$mas))/sd(amakihi$mas)

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

#plot(amakihi.om.best,breaks=c(seq(0,82.5,len=20)),main=c("Detection function","Levels of observer","Quantiles of minutes after sunrise"),style="comp")
#dev.copy2eps(file="amakihi-om.eps",width=10.7,height=4.7)
#plot(amakihi.om.best,breaks=c(seq(0,82.5,len=20)),main="PDF of distances",style="comp",pdf=TRUE)
#dev.copy2eps(file="amakihi-om-pdf.eps",width=4.2,height=3.9)




source("grabresults.R")

grab_results(amakihi.best)
grab_results(amakihi.o.best)
grab_results(amakihi.h.best)
grab_results(amakihi.m.best)
grab_results(amakihi.oh.best)
grab_results(amakihi.om.best)
grab_results(amakihi.mh.best)
grab_results(amakihi.omh.best)

