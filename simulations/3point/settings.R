# 3-point mixtures 
library(mmds)

# fit 2 point mixtures to these too?

parmat<-matrix(NA,2,5)
parmat[1,]<-c(log(0.8),log(0.5),log(0.1),inv.reparam.pi(rep(1/3,3))[1:2])
parmat[2,]<-c(log(15),log(.25),log(0.05),inv.reparam.pi(c(0.1,0.4,0.5))[1:2])

# set the width
width<-1
# number of realisations
n.sims<-200
# sample sizes
n.samps<-c(30,60,120,480,960)
# number of mixture components
mix.terms<-3
# no covariates
model.formula<-"~1"
# starting values (NULL to get fitmix to calculate them)
starting.vals<-NULL

