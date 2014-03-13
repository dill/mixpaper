# pilot whale analysis

library(mmds)

set.seed(125)

fin<-read.csv(file="danpike.csv")

aics<-c()

# start with non-covariates
fin.fit1<-fitmix(fin,mix.terms=1,width=3000)
fin.best<-step.ds.mixture(fin.fit1)
aics<-c(aics,fin.best$aic)

# covariate models

# recode the one observation of BSS=3.5 as 4
fin$BSS5 <- as.character(fin$BSS)
fin$BSS5[fin$BSS5=="3,5"] <- "4"
# BSS
fin.fit1.BSS<-fitmix(fin,mix.terms=1,width=3000,model.formula="~as.factor(BSS5)",opt.method="BFGS+SANN")
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
levels(fin$BSS)[5]<-"3.5"
fin$BSS<-as.numeric(as.character(fin$BSS))
fin.fit1.BSSc<-fitmix(fin,mix.terms=1,width=3000,model.formula="~BSS",opt.method="BFGS+SANN")
fin.BSSc.best<-step.ds.mixture(fin.fit1.BSSc)
aics<-c(aics,fin.BSSc.best$aic)


# best model is ???

cat("Min at",which.min(aics),"\n")

#plot(fin.BSS2.best,breaks=c(seq(0,1000,len=10),seq(1250,3000,250)),main=c("Average detection function","Levels of Beaufort sea state"),style="comp") 

# best model was continuous BSS
pdf(file="danpike-bssc.pdf",width=7.6, height=4.1)
plot(fin.BSSc.best,breaks=c(seq(0,1000,len=10),seq(1250,3000,250)),main=c("Average detection function","Quantiles of Beaufort sea state"),style="comp") 
#dev.copy2eps(file="danpike-bssc.eps",width=7.6, height=4.1)
dev.off()

# plot without histograms too, stitch together in tex...
pdf(file="danpike-bssc-hh.pdf",width=7.6, height=4.1)
plot(fin.BSSc.best,breaks=c(seq(0,1000,len=10),seq(1250,3000,250)),main=c("Average detection function","Quantiles of Beaufort sea state"),style="comp",hide.hist=TRUE) 
#dev.copy2eps(file="danpike-bssc-hh.eps",width=7.6, height=4.1)
dev.off()
#####
source("grabresults.R")

grab_results(fin.best)
grab_results(fin.BSS.best)
grab_results(fin.BSS2.best)
grab_results(fin.BSS3.best)
grab_results(fin.BSSc.best)


## % CV of Pa
#cat("%CV\n")
#cat(round(100*summary(fin.best)$average.p.cv,2),"\n")
#cat(round(100*summary(fin.BSS.best)$average.p.cv,2),"\n")
#cat(round(100*summary(fin.BSS2.best)$average.p.cv,2),"\n")
#cat(round(100*summary(fin.BSS3.best)$average.p.cv,2),"\n")
#cat(round(100*summary(fin.BSSc.best)$average.p.cv,2),"\n")
#
## KS
#cat("\nKS\n")
#cat(round(fin.best$ks$p,2),"\n")
#cat(round(fin.BSS.best$ks$p,2),"\n")
#cat(round(fin.BSS2.best$ks$p,2),"\n")
#cat(round(fin.BSS3.best$ks$p,2),"\n")
#cat(round(fin.BSSc.best$ks$p,2),"\n")




