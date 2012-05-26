# this routine is to "optimize" the parameters of the hazard rate so 
# that we have the same average p as for the half-normal case
library(mmds)

width<-1
n.samps<-10000
mix.terms<-2
par(mfrow=c(1,2))

keyfct.hr<-function(x,keysc,hr.shape){
  return(1-exp(-(x/keysc)^-hr.shape))
}

# with scale=1 but selecting pars again...
mix.hr<-function(x,pars){
  keysc1<-pars[1]
  keysc2<-pars[2]
  pis<-reparam.pi(pars[3])
  hr.shape1<-pars[4]
  hr.shape2<-pars[5]
  pis[1]*keyfct.hr(x,keysc1,hr.shape1)+pis[2]*keyfct.hr(x,keysc2,hr.shape2)
}

# want to reproduce what's in nocov (4)
pars4<-c(log(0.7),log(0.05),inv.reparam.pi(0.6))
pa4<-mmds:::mu.calc(pars4,2,1)
cat("p (4)=",pa4,"\n")

f.int.c<-function(hr.shape,keysc,true.p){
  abs(integrate(keyfct.hr,lower=0,upper=1,keysc=keysc,
                hr.shape=hr.shape)$value-true.p)
}

# component 1
keysc1<-0.7
pa1<-sqrt(2*pi)*keysc1*(pnorm(1,sd=keysc1)-pnorm(0,sd=keysc1))
ll1<-optimize(f.int.c,interval=c(0,20),true.p=pa1,keysc=keysc1)
# component 2
keysc2<-0.05
pa2<-sqrt(2*pi)*keysc2*(pnorm(1,sd=keysc2)-pnorm(0,sd=keysc2))
ll2<-optimize(f.int.c,interval=c(0,20),true.p=pa2,keysc=keysc2)


pars<-c(0.7,0.05,inv.reparam.pi(0.6),ll1$minimum,ll2$minimum)

cat("p =",integrate(mix.hr,lower=0,upper=1,pars=pars)$value,"\n")

f.int.p<-function(p,hr.shape,keysc,true.p){
  pars<-c(keysc,inv.reparam.pi(p),hr.shape)
  abs(integrate(mix.hr,lower=0,upper=1,pars=pars)$value-true.p)
}
llp<-optimize(f.int.p,interval=c(0,1),true.p=pa4,
              keysc=c(0.7,0.05),hr.shape=c(ll1$minimum,ll2$minimum))

pars[3]<-inv.reparam.pi(llp$minimum)

xvec<-seq(0,1,len=100)
g.eval<-mix.hr(xvec,pars)
g.eval1<-keyfct.hr(xvec,keysc1,pars[4])
g.eval2<-keyfct.hr(xvec,keysc2,pars[5])


sim.data <- sim.mix(c(log(pars[1:2]),inv.reparam.pi(0.6)),mix.terms=2,n=960,
                    width=1,key="hr",hr.shape=pars[4:5])

a<-hist(sim.data$distance,plot=FALSE,breaks=seq(0,1,len=20))
hist.sum<-sum(a$density*diff(a$breaks))
a$density<-a$density*(
                integrate(mix.hr,lower=0,upper=1,pars=pars)$value)/hist.sum

plot(a,xlab="Distance",freq=FALSE,main="shape=1, new scale",ylim=c(0,1))
lines(xvec,g.eval)
#plot(xvec,g.eval,type="l")
lines(xvec,g.eval1,col="blue")
lines(xvec,g.eval2,lty=2,col="red")

fit<-try(fitmix(sim.data,mix.terms=2,ftype="hn",width=1,
                model.formula="~1",usegrad=TRUE))
plot(fit)


cat("pars = ", pars,"\n")
