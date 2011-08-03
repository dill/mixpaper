# amakihi analysis

library(mmds)

set.seed(1241)

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


cat("%CV pa\n")
cat(round(100*summary(amakihi.best)$average.p.cv,2),"\n")
cat(round(100*summary(amakihi.o.best)$average.p.c,2),"\n")
cat(round(100*summary(amakihi.h.best)$average.p.c,2),"\n")
cat(round(100*summary(amakihi.m.best)$average.p.c,2),"\n")
cat(round(100*summary(amakihi.oh.best)$average.p.c,2),"\n")
cat(round(100*summary(amakihi.om.best)$average.p.c,2),"\n")
cat(round(100*summary(amakihi.mh.best)$average.p.c,2),"\n")
cat(round(100*summary(amakihi.omh.best)$average.p.c,2),"\n")

cat("KS p\n")
cat(round(amakihi.best$ks$p,2),"\n")
cat(round(amakihi.o.best$ks$p,2),"\n")
cat(round(amakihi.h.best$ks$p,2),"\n")
cat(round(amakihi.m.best$ks$p,2),"\n")
cat(round(amakihi.oh.best$ks$p,2),"\n")
cat(round(amakihi.om.best$ks$p,2),"\n")
cat(round(amakihi.mh.best$ks$p,2),"\n")
cat(round(amakihi.omh.best$ks$p,2),"\n")
