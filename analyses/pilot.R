# pilot whale analysis

library(mmds)

fin<-read.csv(file="danpike.csv")

aics<-c()

# start with non-covariates
fin.fit1<-fitmix(fin,mix.terms=1,width=3000)
fin.best<-step.ds.mixture(fin.fit1)
aics<-c(aics,fin.best$aic)

# covariate models

# BSS
fin.fit1.BSS<-fitmix(fin,mix.terms=1,width=3000,model.formula="~as.factor(BSS)",opt.method="BFGS+SANN")
fin.BSS.best<-step.ds.mixture(fin.fit1.BSS)
aics<-c(aics,fin.BSS.best$aic)


# BSS2
fin.fit1.BSS2<-fitmix(fin,mix.terms=1,width=3000,model.formula="~as.factor(BSS2)",opt.method="BFGS+SANN")
fin.BSS2.best<-step.ds.mixture(fin.fit1.BSS2)
aics<-c(aics,fin.BSS2.best$aic)


# BSS3
fin.fit1.BSS3<-fitmix(fin,mix.terms=1,width=3000,model.formula="~as.factor(BSS3)",opt.method="BFGS+SANN")
fin.BSS3.best<-step.ds.mixture(fin.fit1.BSS3)
aics<-c(aics,fin.BSS3.best$aic)


# BSS - continuous
fin.fit1.BSSc<-fitmix(fin,mix.terms=1,width=3000,model.formula="~BSS",opt.method="BFGS+SANN")
fin.BSSc.best<-step.ds.mixture(fin.fit1.BSSc)
aics<-c(aics,fin.BSSc.best$aic)


# best model is ???

cat("Min at",which.min(aics),"\n")

plot(fin.BSS2.best,breaks=c(seq(0,1000,len=10),seq(1250,3000,250)),main=c("Average detection function","Levels of Beaufort sea state"),style="comp") 

