# analyse the ants data

library(mmds)
set.seed(1451)

# load some data
ants<-read.csv(file="woodants.csv")

# non-covariate fit
ants.fit1<-fitmix(ants,mix.terms=1,width=25)
ants.best<-step.ds.mixture(ants.fit1)
#plot(ants.best,style="comp")

# standardize nest size first
ants$ns<-(ants$nest.size-mean(ants$nest.size))/sd(ants$nest.size)

### covariates

aics<-c()

# 1 - habitat 
ants.fit1.habitat<-fitmix(ants,mix.terms=1,width=25,model.formula="~as.factor(habitat)")
ants.habitat.best<-step.ds.mixture(ants.fit1.habitat)
aics<-c(aics,ants.habitat.best$aic)

# 2 - species
ants.fit1.species<-fitmix(ants,mix.terms=1,width=25,model.formula="~as.factor(species)")
ants.species.best<-step.ds.mixture(ants.fit1.species)
aics<-c(aics,ants.species.best$aic)

# 3 - nest size
ants.fit1.nest<-fitmix(ants,mix.terms=1,width=25,model.formula="~ns")
ants.nest.best<-step.ds.mixture(ants.fit1.nest)
aics<-c(aics,ants.nest.best$aic)

# 4 - nest + habitat
ants.fit1.nesthab <- fitmix(ants, mix.terms = 1, width = 25, model.formula = "~ns + as.factor(habitat)")
ants.nesthab.best<-step.ds.mixture(ants.fit1.nesthab)
aics<-c(aics,ants.nesthab.best$aic)

# 5 - habitat + species
ants.fit1.habspecies<-fitmix(ants, mix.terms = 1, width = 25, model.formula = "~as.factor(habitat) + as.factor(species)")
ants.habspecies.best<-step.ds.mixture(ants.fit1.habspecies)
aics<-c(aics,ants.habspecies.best$aic)

# 6 - nest + species
ants.fit1.nestspecies<-fitmix(ants, mix.terms = 1, width = 25, model.formula = "~ns + as.factor(species)")
ants.nestspecies.best<-step.ds.mixture(ants.fit1.nestspecies)
aics<-c(aics,ants.nestspecies.best$aic)

# 7 - nest + habitat + species
ants.fit1.nesthabspecies<-fitmix(ants, mix.terms = 1, width = 25, model.formula = "~ns + as.factor(habitat) + as.factor(species)")
ants.nesthabspecies.best<-step.ds.mixture(ants.fit1.nesthabspecies)
aics<-c(aics,ants.nesthabspecies.best$aic)

#cat("min AIC=",min(aics),"at",which.min(aics),"\n")


# plot nest+hab model
#plot(ants.nesthab.best,breaks=c(seq(0,5,len=10),seq(6,25,1)),main=c("Detection function","Quantiles of nest size","Levels of habitat"))
# again without the histograms
#plot(ants.nesthab.best,breaks=c(seq(0,5,len=10),seq(6,25,1)),main=c("Detection function","Quantiles of nest size","Levels of habitat"),hide.hist=TRUE)
#dev.copy2eps(file="ants-nesthab.eps",width=8.6,height=3.4)
### can stitch these together later

#save.image("ants-run.RData")

## Model & Covariates & AIC & $\hat{P_a}$ & $\% CV \hat{P_a}$ & K-S $p$\\
#grab_results<-function(model){
#
#   mod<-paste("Hn ",model$mix.terms,"-point",sep="")
#   covars<-model$formula   
#
#   cv.p<-100*summary(model)$average.p.cv
#
#cat(mod," & ",covars," & ",round(model$aic,2)," & ",round(model$pa,3),
#    " & ",round(cv.p,2)," & ",round(model$ks$p,2),"\\\\n")
#
#}

source("grabresults.R")


grab_results(ants.best)
grab_results(ants.habitat.best)
grab_results(ants.species.best)
grab_results(ants.nest.best)
grab_results(ants.habspecies.best)
grab_results(ants.nesthab.best)
grab_results(ants.nestspecies.best)
grab_results(ants.nesthabspecies.best)


#cat("%CV p\n")
#cat(round(100*summary(ants.best)$average.p.cv,2),"\n")
#cat(round(100*summary(ants.habitat.best)$average.p.cv,2),"\n")
#cat(round(100*summary(ants.species.best)$average.p.cv,2),"\n")
#cat(round(100*summary(ants.nest.best)$average.p.cv,2),"\n")
#cat(round(100*summary(ants.habspecies.best)$average.p.cv,2),"\n")
#cat(round(100*summary(ants.nesthab.best)$average.p.cv,2),"\n")
#cat(round(100*summary(ants.nestspecies.best)$average.p.cv,2),"\n")
#cat(round(100*summary(ants.nesthabspecies.best)$average.p.cv,2),"\n")
#
#cat("\nKS p\n")
#cat(round(ants.best$ks$p,2),"\n")
#cat(round(ants.habitat.best$ks$p,2),"\n")
#cat(round(ants.species.best$ks$p,2),"\n")
#cat(round(ants.nest.best$ks$p,2),"\n")
#cat(round(ants.nesthab.best$ks$p,2),"\n")
#cat(round(ants.habspecies.best$ks$p,2),"\n")
#cat(round(ants.nestspecies.best$ks$p,2),"\n")
#cat(round(ants.nesthabspecies.best$ks$p,2),"\n")
