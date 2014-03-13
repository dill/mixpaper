# run simulations for exponential power series

source("eps.R")

library(mmds)
library(Distance)
library(foreach)
library(doMC)
options(cores=3)
registerDoMC()

set.seed(1233)

# setup parameters
pars<-matrix(NA,1,2)
pars[1,] <- c(1.5,0.586)

# set the width
width<-1
# number of realisations
n.sims<-200
# sample sizes
n.samps<-c(30,60,120,480,960)
# no covariates
model.formula<-"~1"
# starting values (NULL to get fitmix to calculate them)
starting.vals<-NULL
showit<-0
opt.method<-"BFGS+SANN"

calc.true.N<-function(pars,n.samples){
  return(n.samples*width/integrate(g.eps,lower=0,upper=1,p=pars[1],l=pars[2])$value)
}

pari <- 1

these.pars<-pars[pari,]
these.pars<-these.pars[1:2]

for(n.samples in n.samps){

  true.N<-calc.true.N(these.pars,n.samples)

  these.sims <- 1:n.sims

  results<-foreach(sim = these.sims, .combine=rbind,
                .inorder=FALSE, .init=c()) %dopar% {
    res<-c()

    # generate some data
    sim.data <- sim.eps(n.samples,these.pars)
    sim.data <- data.frame(observed = rep(1, n.samples),
                           object = 1:n.samples,
                           distance = sim.data)

    #### mmds fitting
    seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed
    fit<-try(fitmix(sim.data,mix.terms=1,ftype="hn",width=1,
                    model.formula="~1",usegrad=TRUE))
    if(class(fit)!="try-error"){
      fit<-try(step.ds.mixture(fit))
    }

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
    fit<-try(ds(sim.data,width,monotonicity="strict",scale="scale"))
    if(all(class(fit)!="try-error")){
      res<-rbind(res,c("cds-hnc",pari,n.samples,sim,"ll",fit$ddf$criterion,
                        fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
    }else{
      res<-rbind(res,c("cds-hnc",pari,n.samples,sim,rep(NA,6)))
    }

    ######################################################## 
    # CDS - hr+poly

    fit<-try(ds(sim.data,width,monotonicity="strict",key="hr",
                adjustment="poly",scale="scale"))
    if(all(class(fit)!="try-error")){
      res<-rbind(res,c("cds-hrp",pari,n.samples,sim,"ll",fit$ddf$criterion,
                        fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
    }else{
      res<-rbind(res,c("cds-hrp",pari,n.samples,sim,rep(NA,6)))
    }
    ######################################################## 
    # CDS - hn+cos (width)
    fit<-try(ds(sim.data,width,monotonicity="strict",scale="width"))
    if(all(class(fit)!="try-error")){
      res<-rbind(res,c("cds-hnc-w",pari,n.samples,sim,"ll",fit$ddf$criterion,
                        fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
    }else{
      res<-rbind(res,c("cds-hnc-w",pari,n.samples,sim,rep(NA,6)))
    }

    ######################################################## 
    # CDS - hr+poly (width)

    fit<-try(ds(sim.data,width,monotonicity="strict",key="hr",
                adjustment="poly",scale="width"))
    if(all(class(fit)!="try-error")){
      res<-rbind(res,c("cds-hrp-w",pari,n.samples,sim,"ll",fit$ddf$criterion,
                        fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
    }else{
      res<-rbind(res,c("cds-hrp-w",pari,n.samples,sim,rep(NA,6)))
    }
    write.table(res,file=paste("eps-",n.samples,"-",pari,"-results.csv",
                               sep=""),append=TRUE,col.names=FALSE,sep=",")

  }
}
