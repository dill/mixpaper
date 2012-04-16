# settings for the no-covar sims - point transect
library(mmds)
library(Distance)

# so, the real parameter matrix would be...
#parmat[1,]<-c(0.8,0.2, 0.6400597)
#parmat[2,]<-c(0.8,0.01, 0.7050089)

#parmat<-matrix(NA,2,3)
#parmat[1,]<-c(log(0.8),log(0.2), 0.5)
#parmat[2,]<-c(log(0.8),log(0.01),0.6)

# as with lt nocov....
parmat<-matrix(NA,4,3)
parmat[1,]<-c(log(0.8),log(0.15), inv.reparam.pi(0.3))
parmat[2,]<-c(log(0.6),log(0.1), inv.reparam.pi(1-0.3))
parmat[3,]<-c(log(10),log(0.2), inv.reparam.pi(0.15))
parmat[4,]<-c(log(0.7),log(0.05),inv.reparam.pi(0.6))


# set the width
width<-1
# number of realisations
n.sims<-200
# sample sizes
n.samps<-c(30,60,120,480,960)
# number of mixture components
mix.terms<-2
# no covariates
model.formula<-"~1"
# starting values (NULL to get fitmix to calculate them)
starting.vals<-NULL
pt<-TRUE # point transect time!

