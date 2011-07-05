# pilot whale analysis

library(mmds)


fin<-read.csv(file="danpike.csv")

# start with non-covariates
fin.fit1<-fitmix(fin,mix.terms=1,width=3000)
fin.best<-step.ds.mixture(fin.fit1)

aics<-c()

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



# best model is BSS2


