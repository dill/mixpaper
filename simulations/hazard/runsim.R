# runsim
library(mmds)
library(Distance)
library(foreach)
library(doMC)
options(cores=6)
registerDoMC()


pars<-matrix(NA,2,5)

pars[1,] <- c(-1.6940326, -0.3037095, inv.reparam.pi(0.5),7,NA)
#pars[2,] <- c(-1.814902, -0.8734648, inv.reparam.pi(0.5),1)
pars[2,] <- c(-0.3566749, -2.995732, 0.439975, 10.47956, 3.699274)

# set the width
width<-1
# number of realisations
n.sims<-200
# sample sizes
n.samps<-c(30,60,120,480,960)
n.samps<-c(30,60,480,960)
#n.samps<-c(120)
# number of mixture components
mix.terms<-2
# no covariates
model.formula<-"~1"
# starting values (NULL to get fitmix to calculate them)
starting.vals<-NULL
showit<-0
opt.method<-"BFGS+SANN"

calc.true.N<-function(pars,hr.par,n.samples){
  keyfct.hr<-function(x,keysc,hr.shape){
    return(1-exp(-(x/exp(keysc))^-hr.shape))
  }
  
  mix.hr<-function(x,pars,hr.par){
#    pis<-rep(0.5,2) # for sim 1
    pis<-c(0.6,0.4) # for sim 2
    keysc1<-pars[1]
    keysc2<-pars[2]
    pis[1]*keyfct.hr(x,keysc1,hr.par[1])+pis[2]*keyfct.hr(x,keysc2,hr.par[2])
#    pis[1]*keyfct.hr(x,keysc1,hr.par)+pis[2]*keyfct.hr(x,keysc2,hr.par)
  }
  
  return(n.samples*width/integrate(mix.hr,lower=0,upper=1,pars=pars,hr.par=hr.par)$value)
}

nsim<-200

#for(pari in 1:nrow(pars)){
pari <- 2

  these.pars<-pars[pari,]
#  hr.par<-these.pars[4]
  hr.par<-these.pars[4:5]
  these.pars<-these.pars[1:3]

  for(n.samples in n.samps){
#n.samples<-960

    true.N<-calc.true.N(these.pars,hr.par,n.samples)

    results<-foreach(sim = 1:nsim, .combine=rbind,
                  .inorder=FALSE, .init=c()) %dopar% {
#sim<-1
      res<-c()
    
      # generate some data
      sim.data <- sim.mix(these.pars,2,n.samples,1,key="hr",hr.shape=hr.par)

      #### mmds fitting
      seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed
      fit<-try(fitmix(sim.data,mix.terms=1,ftype="hn",width=1,
                      model.formula="~1",usegrad=TRUE))
      fit<-step.ds.mixture(fit)

      # restore the seed   
      assign(".Random.seed",seed,envir=.GlobalEnv)

      # result vector
      if(class(fit)!="try-error"){
         res<-rbind(res,c("mmds-MS",pari,n.samples,sim,fit$likelihood,
                fit$aic,fit$pa,fit$N, true.N,fit$mix.terms))
      }else{
         res<-rbind(res,c("mmds-MS",pari,n.samples,sim,rep(NA,6)))
      }

      # CDS
      ######################################################## 
      # CDS - hn+cos
      fit<-try(ds(sim.data,width,monotonicity="strict"))
      if(all(class(fit$ddf)!="try-error")){
        res<-rbind(res,c("cds-hnc",pari,n.samples,sim,"ll",fit$ddf$criterion,
                          fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
      }else{
        res<-rbind(res,c("cds-hnc",pari,n.samples,sim,rep(NA,6)))
      }

      ######################################################## 
      # CDS - hr+poly

      fit<-try(ds(sim.data,width,monotonicity="strict",key="hr",
                  adjustment="poly"))
      if(all(class(fit$ddf)!="try-error")){
        res<-rbind(res,c("cds-hnc",pari,n.samples,sim,"ll",fit$ddf$criterion,
                          fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
      }else{
        res<-rbind(res,c("cds-hrp",pari,n.samples,sim,rep(NA,7)))
      }
      ######################################################## 
      # CDS - hn+cos (width)
      fit<-try(ds(sim.data,width,monotonicity="strict",scale="width"))
      if(all(class(fit$ddf)!="try-error")){
        res<-rbind(res,c("cds-hnc-w",pari,n.samples,sim,"ll",fit$ddf$criterion,
                          fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
      }else{
        res<-rbind(res,c("cds-hnc-w",pari,n.samples,sim,rep(NA,6)))
      }

      ######################################################## 
      # CDS - hr+poly (width)

      fit<-try(ds(sim.data,width,monotonicity="strict",key="hr",
                  adjustment="poly",scale="width"))
      if(all(class(fit$ddf)!="try-error")){
        res<-rbind(res,c("cds-hnc-w",pari,n.samples,sim,"ll",fit$ddf$criterion,
                          fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
      }else{
        res<-rbind(res,c("cds-hrp-w",pari,n.samples,sim,rep(NA,7)))
      }

      write.table(res,file=paste("hr-",n.samples,"-",pari,"-results.csv",
                                 sep=""),append=TRUE,col.names=FALSE,sep=",")

#      write.table(res,file=paste("width-hr-",n.samples,"-",pari,"-results.csv",
#                                 sep=""),append=TRUE,col.names=FALSE)
    }
  }
#}

