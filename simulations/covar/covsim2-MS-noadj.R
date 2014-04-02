# covariate sim 2
# continuous
library(mmds)
library(Distance)
# let's do this in parallel
library(foreach)
library(doMC)
options(cores=6)
registerDoMC()

### SETUP
model.formula<-"~cov1"

### parameters
true.pars<-c(log(0.2),log(0.8),log(0.4), inv.reparam.pi(0.4))
mix.terms<-2
width<-1

opt.method<-"BFGS+SANN"

# probably don't touch below here...

# sample sizes
n.samps<-c(30,60,120,480,960)
# number of realisations
n.sims<-200
# set the seed
set.seed(1028)
showit<-0
starting.vals<-NULL # unknown starting values

big.res<-c()
# loop over sample sizes
for(n.samples in n.samps){
  # loop over sims
  #for(sim in 1:n.sims){
  results<-foreach(sim = 1:n.sims, .combine=rbind,
                   .inorder=FALSE, .init=c()) %dopar% {

     res<-c()

     ### # This changes per simulation parameter setting   
     # construct the covariate matrix
     z<-list(matrix(c(rep(1,n.samples),
                      pnorm(seq(-4,4,len=n.samples))),n.samples,2))
     zdim<-c(2)
     testdata<-sim.mix(true.pars,mix.terms,n.samples,width,zdim,z)
     names(testdata)[5]<-"cov1"

     # calculate the true N
     gp<-mmds:::getpars(true.pars,mix.terms,zdim,z)
     sigma<-gp$key.scale
     pis<-gp$mix.prop
     mu<-apply(sigma,1,mmds:::integrate.hn,width)
     mus<-mu%*%matrix(pis,length(pis),1)
     pas<-mus/width
     true.N<-sum(1/pas)


     ######################################################## 
     # MCDS - hn + as.factor(cov1)

     fit<-try(ds(testdata,width,formula=as.formula(model.formula),
                 adjustment=NULL))
     if(all(class(fit)!="try-error")){
       res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
                        n.samples/sum(1/fitted(fit$ddf)),
                        fit$ddf$Nhat,true.N,NA,"hn+cov1"))
     }else{
       res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"hn+cov1"))
     }

     ######################################################## 
     # MCDS - hr + as.factor(cov1)

     fit<-try(ds(testdata,width,formula=as.formula(model.formula),key="hr",
                 adjustment=NULL))
     if(all(class(fit)!="try-error")){
       res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
                        n.samples/sum(1/fitted(fit$ddf)),
                        fit$ddf$Nhat,true.N,NA,"hr+cov1"))
     }else{
       res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"hr+cov1"))
     }


#write.table(res,file=paste("covsim2-noadj",n.samples,"-",opt.method,".csv",sep=""),append=TRUE,col.names=FALSE)
     return(res)
  }

  big.res<-rbind(big.res,results)
}
write.csv(big.res,file=paste("covsim2-noadj",opt.method,".csv",sep=""))
