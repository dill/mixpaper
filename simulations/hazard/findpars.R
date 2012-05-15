# this routine is to "optimize" the parameters of the hazard rate so 
# that we have the same average p as for the half-normal case


library(mmds)

# parameters
width<-1
# set this pretty big to begin with
n.samps<-10000
mix.terms<-2
model.formula<-"~1"
par(mfrow=c(3,2))

hr.par<-7

keyfct.hr<-function(x,keysc,hr.shape){
  return(1-exp(-(x/exp(keysc))^-hr.shape))
}

mix.hr<-function(x,pars){
  pis<-rep(0.5,2)#reparam.pi(pars[1])
  keysc1<-pars[1]
  keysc2<-pars[2]
  #hr.shape<-pars[4]
  hr.shape<-hr.par
  pis[1]*keyfct.hr(x,keysc1,hr.shape)+pis[2]*keyfct.hr(x,keysc2,hr.shape)
}


f.int<-function(pars){
#cat("pars=",pars,"\n")
  abs(integrate(mix.hr,lower=0,upper=1,pars=pars)$value-true.p)
}


xvec<-seq(0,1,len=100)
true.p<-0.5

# do the optimisation
ll<-optim(log(c(0.1,0.1)),f.int)

# what did we get?

#cat("pi=",reparam.pi(ll$par[1]),"\n")
cat("scales=",ll$par[1:2],"\n")
#cat("shapes=",ll$par[4],"\n")

sim.data <- sim.mix(c(ll$par[1:2],inv.reparam.pi(0.5)),mix.terms=2,n=500,
                    width=1,key="hr",hr.shape=hr.par)
g.eval<-mix.hr(xvec,ll$par)
g.eval1<-keyfct.hr(xvec,ll$par[1],hr.par)
g.eval2<-keyfct.hr(xvec,ll$par[2],hr.par)

a<-hist(sim.data$distance,plot=FALSE)#,breaks=seq(0,1,len=50))
hist.sum<-sum(a$density*diff(a$breaks))
a$density<-a$density*(f.int(ll$par)+0.5)/hist.sum


plot(a,xlab="Distance",freq=FALSE,main="shape=7")
lines(xvec,g.eval)
lines(xvec,g.eval1,col="blue")
lines(xvec,g.eval2,lty=2,col="red")

fit<-try(fitmix(sim.data,mix.terms=2,ftype="hn",width=1,
                model.formula="~1",usegrad=TRUE))
plot(fit)






## for a hazard parameter of 1,holding the other pars as-is
#hr.par<-1
#
#sim.data <- sim.mix(c(ll$par[1:2],inv.reparam.pi(0.5)),mix.terms=2,n=500,
#                    width=1,key="hr",hr.shape=hr.par)
#g.eval<-mix.hr(xvec,ll$par)
#g.eval1<-keyfct.hr(xvec,ll$par[1],hr.par)
#g.eval2<-keyfct.hr(xvec,ll$par[2],hr.par)
#
#a<-hist(sim.data$distance,plot=FALSE)#,breaks=seq(0,1,len=50))
#hist.sum<-sum(a$density*diff(a$breaks))
#a$density<-a$density*(f.int(ll$par)+0.5)/hist.sum
#
#plot(a,xlab="Distance",freq=FALSE,main="shape=1")
#lines(xvec,g.eval)
#lines(xvec,g.eval1,col="blue")
#lines(xvec,g.eval2,lty=2,col="red")
#
#fit<-try(fitmix(sim.data,mix.terms=2,ftype="hn",width=1,
#                model.formula="~1",usegrad=TRUE))
#plot(fit)




# with scale=1 but selecting pars again...


ll<-optim(log(c(0.1,0.1)),f.int)
cat("scales=",ll$par[1:2],"\n")
sim.data <- sim.mix(c(ll$par[1:2],inv.reparam.pi(0.5)),mix.terms=2,n=500,
                    width=1,key="hr",hr.shape=hr.par)
g.eval<-mix.hr(xvec,ll$par)
g.eval1<-keyfct.hr(xvec,ll$par[1],hr.par)
g.eval2<-keyfct.hr(xvec,ll$par[2],hr.par)

a<-hist(sim.data$distance,plot=FALSE)#,breaks=seq(0,1,len=50))
hist.sum<-sum(a$density*diff(a$breaks))
a$density<-a$density*(f.int(ll$par)+0.5)/hist.sum

plot(a,xlab="Distance",freq=FALSE,main="shape=1, new scale")
lines(xvec,g.eval)
lines(xvec,g.eval1,col="blue")
lines(xvec,g.eval2,lty=2,col="red")

fit<-try(fitmix(sim.data,mix.terms=2,ftype="hn",width=1,
                model.formula="~1",usegrad=TRUE))
plot(fit)

