# function for exponential power series
# lambda,p >0

# eps detection function
g.eps <- function(y,p,lambda){
  exp(-(y/lambda)^p)
}

# eps pdf
f.eps <- function(y,p,lambda){
  g.eps(y,p,lambda)/(lambda*gamma(1+1/p))
}

# simulate from eps
sim.eps <- function(n,pars){
  samples <- c()
  while(length(samples)<n){
    proposal <- runif(2*n)
    p.in <- g.eps(proposal,pars[1],pars[2])
    samples <- c(samples,proposal[runif(2*n)<=p.in])
  }
  return(samples[1:n])
}

# test
#hh <- hist(sim.eps(100000,c(1.5,0.586)),plot=FALSE)
#hh$density <- hh$density*(integrate(g.eps,lower=0,upper=1,p=1.5,lambda=0.586)$value/sum(hh$density*diff(hh$breaks)))
#plot(hh,freq=FALSE)
#lines(xx,g.eps(xx,1.5,0.586),col="red")


# conclusion from the below is we use:
# p =1.5
# lambda = 0.586
### 3.95435,0.6519915 <- p~~0.75
#plot(xx,g.eps(xx,1.5,0.586),col="red",asp=1,ylim=c(0,1),type="l")
#integrate(g.eps,lower=0,upper=1,p=1.5,lambda=0.586)$value
# [1] 0.5000914

#ff <- function(p,l){
#plot(xx,g.eps(xx,p,l),col="red",asp=1,ylim=c(0,1),type="l")
#print(integrate(g.eps,lower=0,upper=1,p=p,lambda=l)$value)
#}

## try some plots -- pars from Roccio
#xx <- seq(0,30,len=1000)
#plot(xx,g.eps(xx,20.125,1.5),type="l",col="red")
#lines(xx,g.eps(xx,19.453,20),col="red",lty=2)
#lines(xx,g.eps(xx,30.975,1.5),col="blue")
#lines(xx,g.eps(xx,23.825,20),col="blue",lty=2)
#lines(xx,g.eps(xx,55.18,1.5),col="green")
#lines(xx,g.eps(xx,27.512,20),col="green",lty=2)
#
#
## manually look at what happens when we vary the parameters
#l.start <- 1
#xx <- seq(0,1,len=1000)
#plot(xx,g.eps(xx,1,l.start),type="l",col="red")
#for(pval in seq(0.3,2,len=20)){
#  lines(xx,g.eps(xx,pval,l.start),col="green")
#}
#
#ob.p <- function(p){
#  abs(integrate(g.eps,lower=0,upper=1,p=p,lambda=l.start)$value-0.5)
#}
#
#oo<-optimize(ob.p,interval=c(0.01,20))
#
#lines(xx,g.eps(xx,oo$minimum,l.start))
#
#ob.l <- function(l,p){
#  abs(integrate(g.eps,lower=0,upper=1,p=p,lambda=l)$value-0.5)
#}
#
#oo.l <- optimize(ob.l,interval=c(0.01,20),p=oo$minimum)
#
#
#lines(xx,g.eps(xx,oo$minimum,oo.l$minimum))
#
#
#par(mfrow=c(1,3))
## set l and re-optimize p
#l.start <- oo.l$minimum
#oo<-optimize(ob.p,interval=c(0.01,5))
#plot(xx,g.eps(xx,oo$minimum,oo.l$minimum),type="l",col="red",ylim=c(0,1))
#
#ob.w <- function(p,l){
#  abs(g.eps(1,p,l)-0.2)
#}
#
#oo.w <- optimize(ob.w,interval=c(0.01,20),p=oo$minimum)
#
#oo.wl <- optimize(ob.w,interval=c(0.01,20),l=oo.l$minimum)
#
#
#plot(xx,g.eps(xx,oo.w$minimum,oo.l$minimum),type="l",col="red",ylim=c(0,1))
#plot(xx,g.eps(xx,oo$minimum,oo.wl$minimum),type="l",col="red",ylim=c(0,1))
#
#
#oo.wl <- optimize(ob.w,interval=c(0.01,20),p=oo.w$minimum)
#plot(xx,g.eps(xx,oo$minimum,oo.wl$minimum),type="l",col="red",ylim=c(0,1))
#
#
#
#ob.big <- function(pars){
#  p <- pars[1]
#  l <- pars[2]
#  ret <- abs(integrate(g.eps,lower=0,upper=1,p=p,lambda=l)$value-0.5)
#
#  if(g.eps(1,p,l)<0.2 | g.eps(1,p,l)> 0.3){
#    ret <- ret+g.eps(1,p,l)*1000
#  }
#}
#
#library(optimx)
##oo <- optimx(c(oo$minimum,oo.wl$minimum),ob.big,method="nlminb",lower=c(0.01,0.01))
#oo2 <- optimx(c(oo$minimum,oo.wl$minimum),ob.big,method="nlminb",lower=c(2.5,0.1))
#
#
#plot(xx,g.eps(xx,oo2$p1,oo2$p2),type="l",col="red",ylim=c(0,1))





