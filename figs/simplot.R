# new sampling function

library(mmds)

plotcomp<-function(pars,mix.terms,pt=FALSE,asp=1,this.id){
   gp<-mmds:::getpars(pars,mix.terms)
   sigmas<-gp$key.scale
   pis<-gp$mix.prop

   if(pt){
      intfcn<-mmds:::integrate.hn.pt
   }else{
      intfcn<-mmds:::integrate.hn
   }

   x<-seq(0,1,len=1000)
   pv<-matrix(0,mix.terms,length(x))
   mu<-c()

   for(i in 1:mix.terms){
      pv[i,]<-mmds:::keyfct.hn(x,sigmas[i])
      mu<-c(mu,intfcn(sigmas[i],width))
   }
   ptotal<-mmds:::detfct(x,pars,mix.terms)
   mutotal<-mmds:::mu.calc(pars,mix.terms,width,pt=pt)

   if(pt){
      ptotal<-2*pi*x*ptotal/mutotal
      for(i in 1:mix.terms){
         pv[i,]<-2*pi*x*pv[i,]/mu[i]
      }
   }

   if(!pt){
      yl<-c(0,max(c(1,ptotal,pv)))
   }else{
      yl<-c(0,max(c(ptotal,pv)))
      asp <- 1/yl[2]
   }
   plot(x,ptotal,type="l",ylim=yl,xlim=c(0,1),
        asp=asp,axes=F)
   axis(1,c(0,0.5,1))
   if(!pt){
      axis(2,c(0,0.5,1))
   }else{
      axis(2)
   }
   box()
   for(i in 1:mix.terms){
      lines(x,pv[i,],lty=2)
   }

  # plot the id
  points(x=1,y=yl[2]*0.85,cex=3.75)
  text(x=1,y=yl[2]*0.85,label=this.id)
}


#postscript(file="sim-detfct.eps",width=9,height=9,
#            paper="special",horizontal=FALSE)
pdf(file="sim-detfct.pdf",width=7,height=7.5)
par(mfrow=c(5,4),mar=c(2,2.2,1.8,1.5),las=1,oma=c(2,2,0,0))

### plot for lt and pt
width<-1
mix.terms<-2
parmat<-matrix(NA,4,3)
parmat[1,]<-c(log(0.8),log(0.15), inv.reparam.pi(0.3))
parmat[2,]<-c(log(0.6),log(0.1), inv.reparam.pi(1-0.3))
parmat[3,]<-c(log(10),log(0.2), inv.reparam.pi(0.15))
parmat[4,]<-c(log(0.7),log(0.05),inv.reparam.pi(0.6))

ids <- c(paste0("A",1:4),paste0("B",1:4))
k<-1
for(pt in c(FALSE,TRUE)){
   if(pt) asp<-NA else asp<-1
   for(par.ind in 1:4){
      pars<-parmat[par.ind,]
      plotcomp(pars,mix.terms,pt,asp,ids[k])
      k<-k+1
   }
}

### plot for 3-point
parmat<-matrix(NA,2,5)
parmat[1,]<-c(log(0.8),log(0.5),log(0.1),inv.reparam.pi(rep(1/3,3))[1:2])
parmat[2,]<-c(log(15),log(.25),log(0.05),inv.reparam.pi(c(0.1,0.4,0.5))[1:2])
mix.terms<-3

ids <- c(paste0("C",1:2))
for(par.ind in 1:2){
   pars<-parmat[par.ind,]
   plotcomp(pars,mix.terms,FALSE,asp=1,ids[par.ind])
}

# two blank plots
plot(0,0,axes=F,type="n")
plot(0,0,axes=F,type="n")


# now plot the covariate models

# covsim1
model.formula<-"~as.factor(cov1)"
pars<-c(log(c(0.1,0.75,0.6)),inv.reparam.pi(0.4))
mix.terms<-2
n.samples<-1000
z<-list(matrix(c(rep(1,n.samples),rep(c(0,1),n.samples/2)),n.samples,2))
zdim<-2
testdata<-sim.mix(pars,mix.terms,n.samples,width,zdim,z)
names(testdata)[5]<-"cov1"
fit<-try(fitmix(testdata,mix.terms=2,ftype="hn",
            width=width,model.formula=model.formula))


#source("~/current/mmds/R/plot.ds.mixture.R")
source("plot.ds.mixture.R")
switchpars <- mmds:::switchpars
getpars <- mmds:::getpars
integrate.hn <- mmds:::integrate.hn
keyfct.hn <- mmds:::keyfct.hn
plot(fit,nomf=TRUE,hide.hist=TRUE,style="comp",main=c(" "," "),x.axis=c(0,0.5,1),this.id="D1")

# covsim2
model.formula<-"~cov1"
pars<-c(log(0.2),log(0.8),log(0.4), inv.reparam.pi(0.4))
z<-list(matrix(c(rep(1,n.samples),pnorm(seq(-4,4,len=n.samples))),n.samples,2))
zdim<-c(2)
testdata<-sim.mix(pars,mix.terms,n.samples,width,zdim,z)
names(testdata)[5]<-"cov1"
fit<-try(fitmix(testdata,mix.terms=2,ftype="hn",width=width,
                model.formula=model.formula))


plot(fit,nomf=TRUE,hide.hist=TRUE,style="comp",,main=c(" "," "),x.axis=c(0,0.5,1),this.id="D2")


# eps
pars<-c(1.5,0.586)
g.eps <- function(pars){
  exp(-(seq(0,1,len=1000)/pars[2])^pars[1])
}
plot(seq(0,1,len=1000),g.eps(pars),type="l", ylim=c(0,1), xlim=c(0,1),axes=F,asp=1)
axis(1,c(0,0.5,1))
axis(2,c(0,0.5,1))
points(x=1,y=0.85,cex=3.75)
text(x=1,y=0.85,label="E1")
box()


# hazard
pars<-matrix(NA,2,5)

pars[1,] <- c(-1.6940326, -0.3037095, 0.5,7,7)
pars[2,] <- c(-0.3566749, -2.995732, 0.4, 10.47956, 3.699274)

haz<-function(pars){
  hr.par <- pars[4:5]
  these.pars <- pars[1:3]

  keyfct.hr<-function(x,keysc,hr.shape){
    return(1-exp(-(x/exp(keysc))^-hr.shape))
  }

  mix.hr<-function(x,pars,hr.par){
    pis <- pars[3]
    keysc1<-pars[1]
    keysc2<-pars[2]
    pis*keyfct.hr(x,keysc1,hr.par[1])+(1-pis)*keyfct.hr(x,keysc2,hr.par[2])
  }

  p1 <- mix.hr(seq(0,1,len=1000),pars,hr.par)
#  int1 <- integrate(mix.hr,lower=0,upper=1,pars=pars,hr.par=hr.par)

#  return(p1/int1)
  return(p1)#/int1)
}

plot(seq(0,1,len=1000),haz(pars[1,]),type="l", ylim=c(0,1), xlim=c(0,1),axes=F,asp=1)
axis(1,c(0,0.5,1))
axis(2,c(0,0.5,1))
pars[1,3]<-0
lines(seq(0,1,len=1000),haz(pars[1,]),type="l",lty=2)
pars[1,3]<-1
lines(seq(0,1,len=1000),haz(pars[1,]),type="l",lty=2)
points(x=1,y=0.85,cex=3.75)
text(x=1,y=0.85,label="E2")
box()

#plot(seq(0,1,len=1000),haz(pars[2,]),type="l", ylim=c(0,1), xlim=c(0,1),axes=F,asp=1)
#axis(1,c(0,0.5,1))
#axis(2,c(0,0.5,1))
#pars[2,3]<-0
#lines(seq(0,1,len=1000),haz(pars[2,]),type="l",lty=2)
#pars[2,3]<-1
#lines(seq(0,1,len=1000),haz(pars[2,]),type="l",lty=2)
#points(x=1,y=0.85,cex=3.75)
#text(x=1,y=0.85,label="E2")
#box()


# y-axis labels
#text.y <- c(0.2,0.5,0.8)
#mtext("Probability of detection",side=2,outer=TRUE,at=text.y[1],las=3,cex=0.7)
##mtext("Probability of detection",side=2,outer=TRUE,at=text.y[1],las=3,cex=0.7)
#mtext("Probability density",side=2,outer=TRUE,at=text.y[2],las=3,cex=0.7)
#mtext("Probability of detection",side=2,outer=TRUE,at=text.y[3],las=3,cex=0.7)

text.y <- rev(c(0.1,0.3,0.5,0.7,0.9))
mtext("Probability of detection",side=2,outer=TRUE,at=text.y[1],las=3,cex=0.7)
#mtext("Probability of detection",side=2,outer=TRUE,at=text.y[1],las=3,cex=0.7)
mtext("Probability density",side=2,outer=TRUE,at=text.y[2],las=3,cex=0.7)
for(i in 3:5)
mtext("Probability of detection",side=2,outer=TRUE,at=text.y[i],las=3,cex=0.7)

# x-axis label
mtext("Distance",side=1,outer=TRUE,at=c(0.5),cex=0.7)


dev.off()
