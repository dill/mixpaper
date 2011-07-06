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


plot(fin.BSS2.best,breaks=c(0,250,500,750,1000,1250,1500,1750,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2000,2250,2500,2750,3000))




