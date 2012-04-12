# settings for the no-covar sims - line transect
library(mmds)
library(Distance)

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

############# SETTINGS DONE ##############

# plot the detection functions
#par(mfrow=c(2,3))
#for(i in 1:dim(parmat)[1]){
#   plot.ds.mixture(list(pars=parmat[i,],width=width,mix.terms=2),style="comp",xlim=c(0,1))
#   text(80,0.75,paste(expression(sigma),"=(",round(exp(parmat[i,1]),3),",",
#                                              round(exp(parmat[i,2]),3),")\n",
#                       expression(pi),"=(",parmat[i,3],")",sep=""))
#}
