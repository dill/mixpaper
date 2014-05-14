# this file actually runs the simulations - point transect models

# source in the settings
source("settings.R") 

# let's do this in parallel
library(foreach)
library(doMC)
options(cores=6)
registerDoMC()

######################################
# optimisation options 

# optimisation method
opt.method<-"BFGS+SANN"
# for debug
showit<-0

######################################
# actually do the simulation

## loop over possible parameters
results<-foreach(par.ind = 1:dim(parmat)[1], .combine=rbind,
                   .inorder=FALSE, .init=c()) %dopar% {
  res<-c()

  # make sure that the seed is the same on all cores
  set.seed(123)

  ## loop over number of samples to take
  for(n.samp in n.samps){
    ## loop over number of replicates
   for(sim in 1:n.sims){

      # simulate the data
      sim.data<-sim.mix(parmat[par.ind,],mix.terms,n.samp,width,pt=TRUE)
      # calculate the true abundance
      true.N<-(n.samp*pi*width^2)/mmds:::mu.calc(parmat[par.ind,],
                  mix.terms,width,pt=TRUE)

      # because SANN changes the seed, save it first
      seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed

      ##############################################
      # fit the TRUE model -- ie the 2-point mixture
      fit<-try(fitmix(sim.data,initialvalues=starting.vals,
                      mix.terms=mix.terms,
                  ftype="hn",width=width,model.formula=model.formula,
                  usegrad=TRUE,opt.method=opt.method,showit=showit,pt=pt))

      # restore the seed
      assign(".Random.seed",seed,envir=.GlobalEnv)

      # result vector -
      #  parameter set, # samples, sim #, par est, aic
      if((class(fit)!="try-error")){
         # par index, no. samples, sim no., par est(3),likelihood, aic, trueN, Nhat
         res<-rbind(res,c("mmds-2",par.ind,n.samp,sim,fit$pars,
                          fit$likelihood,fit$aic,true.N,fit$N))
      }else{
         res<-rbind(res,c("mmds-2",par.ind,n.samp,sim,rep(NA,7)))
      }

      ##############################################
      # fit the WRONG model -- ie the 1-point mixture

      # because SANN changes the seed, save it first
      seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed
      fit<-try(fitmix(sim.data,initialvalues=starting.vals,mix.terms=1,
                  ftype="hn",width=width,model.formula=model.formula,
                  usegrad=TRUE,opt.method=opt.method,showit=showit,pt=pt))

      # restore the seed
      assign(".Random.seed",seed,envir=.GlobalEnv)

      # result vector -
      #  parameter set, # samples, sim #, par est, aic
      if(class(fit)!="try-error"){
         # par index, no. samples, sim no., par est(1),likelihood, aic, trueN, Nhat
         res<-rbind(res,c("mmds-1",par.ind,n.samp,sim,fit$pars,NA,NA,
                          fit$likelihood,fit$aic,true.N,fit$N))
      }else{
         res<-rbind(res,c("mmds-1",par.ind,n.samp,sim,rep(NA,7)))
      }

      ##############################################
      # CDS - hn+cos (scale scaling)

      fit<-try(ds(sim.data,width,monotonicity="strict",transect="point"))
      if(all(class(fit)!="try-error")){
        res<-rbind(res,c("cds-hnc",par.ind,n.samp,sim,rep(NA,3),
                         fitted(fit$ddf)[1],fit$ddf$criterion,
                         true.N,fit$ddf$Nhat))
      }else{
        res<-rbind(res,c("cds-hnc",par.ind,n.samp,sim,rep(NA,7)))
      }

      ######################################################## 
      # CDS - hr+poly (scale scaling)

      fit<-try(ds(sim.data,width,monotonicity="strict",key="hr",
                  adjustment="poly",transect="point"))
      if(all(class(fit)!="try-error")){
        res<-rbind(res,c("cds-hrp",par.ind,n.samp,sim,rep(NA,3),
                         fitted(fit$ddf)[1],fit$ddf$criterion,
                         true.N,fit$ddf$Nhat))
      }else{
        res<-rbind(res,c("cds-hrp",par.ind,n.samp,sim,rep(NA,7)))
      }

      ##############################################
      # CDS - hn+cos (width scaling)

      fit<-try(ds(sim.data,width,monotonicity="strict",transect="point",
                  scale="width"))
      if(all(class(fit)!="try-error")){
        res<-rbind(res,c("cds-hnc-w",par.ind,n.samp,sim,rep(NA,3),
                         fitted(fit$ddf)[1],fit$ddf$criterion,
                         true.N,fit$ddf$Nhat))
      }else{
        res<-rbind(res,c("cds-hnc-w",par.ind,n.samp,sim,rep(NA,7)))
      }

      ######################################################## 
      # CDS - hr+poly (width scaling)

      fit<-try(ds(sim.data,width,monotonicity="strict",key="hr",
                  adjustment="poly",transect="point",scale="width"))
      if(all(class(fit)!="try-error")){
        res<-rbind(res,c("cds-hrp-w",par.ind,n.samp,sim,rep(NA,3),
                         fitted(fit$ddf)[1],fit$ddf$criterion,
                         true.N,fit$ddf$Nhat))
      }else{
        res<-rbind(res,c("cds-hrp-w",par.ind,n.samp,sim,rep(NA,7)))
      }

      ######################################################## 
      # CDS - unif+cos

      fit<-try(ds(sim.data,width,monotonicity="strict",key="unif",
                  adjustment="cos",transect="point"))
      if(all(class(fit)!="try-error")){
        res<-rbind(res,c("cds-uc",par.ind,n.samp,sim,rep(NA,3),
                         fitted(fit$ddf)[1],fit$ddf$criterion,
                         true.N,fit$ddf$Nhat))
      }else{
        res<-rbind(res,c("cds-uc",par.ind,n.samp,sim,rep(NA,7)))
      }

    }


  }
  # write them out
  write.csv(res,file=paste("width-",opt.method,"-",par.ind,"-pt-results.csv",sep=""))
  return(1)
}
