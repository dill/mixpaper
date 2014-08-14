# this file actually runs the simulation

# source in the settings
source("settings.R") 

# let's do this in parallel
library(foreach)
library(doMC)
options(cores=6)
registerDoMC()

######################################
# optimisation options 

# optimisation methods
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
      sim.data<-sim.mix(parmat[par.ind,],mix.terms,n.samp,width)

      # because SANN changes the seed, save it first
      seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed

      # fit the model for 1-point mixture
      fit<-try(fitmix(sim.data,initialvalues=starting.vals,mix.terms=1,
                  ftype="hn",width=width,model.formula=model.formula,
                  usegrad=TRUE,opt.method=opt.method,showit=showit))

      fit<-step.ds.mixture(fit)

      # restore the seed
      assign(".Random.seed",seed,envir=.GlobalEnv)

      true.N<-n.samp*width/mmds:::mu.calc(parmat[par.ind,],mix.terms,width)

      # result vector -
      #  parameter set, # samples, sim #, par est, aic 
      if(class(fit)!="try-error"){
         # par index, no. samples, sim no., par est(3),likelihood, aic, pa, Nhat, N
         res<-rbind(res,c("mmds-MS",par.ind,n.samp,sim,fit$likelihood,
                          fit$aic,fit$pa,fit$N, true.N,fit$mix.terms))
      }else{
         res<-rbind(res,c("mmds-MS",par.ind,n.samp,sim,rep(NA,6)))
      }

      #### CDS below here!
      ######################################################## 
      # CDS - hn+cos (scale scaling)

      fit<-try(ds(sim.data,width,monotonicity="strict"))
      if(all(class(fit)!="try-error")){
        res<-rbind(res,c("cds-hnc",par.ind,n.samp,sim,"ll",fit$ddf$criterion,
                          fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
      }else{
        res<-rbind(res,c("cds-hnc",par.ind,n.samp,sim,rep(NA,6)))
      }

      ######################################################## 
      # CDS - hr+poly (scale scaling)

      fit<-try(ds(sim.data,width,monotonicity="strict",key="hr",
                  adjustment="poly"))
      if(all(class(fit)!="try-error")){
        res<-rbind(res,c("cds-hnc",par.ind,n.samp,sim,"ll",fit$ddf$criterion,
                          fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
      }else{
        res<-rbind(res,c("cds-hrp",par.ind,n.samp,sim,rep(NA,7)))
      }
      ######################################################## 
      # CDS - hn+cos (width scaling)

      fit<-try(ds(sim.data,width,monotonicity="strict",scale="width"))
      if(all(class(fit)!="try-error")){
        res<-rbind(res,c("cds-hnc-w",par.ind,n.samp,sim,"ll",fit$ddf$criterion,
                          fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
      }else{
        res<-rbind(res,c("cds-hnc-w",par.ind,n.samp,sim,rep(NA,6)))
      }

      ######################################################## 
      # CDS - hr+poly (width scaling)

      fit<-try(ds(sim.data,width,monotonicity="strict",key="hr",
                  adjustment="poly",scale="width"))
      if(all(class(fit)!="try-error")){
        res<-rbind(res,c("cds-hnc-w",par.ind,n.samp,sim,"ll",fit$ddf$criterion,
                          fitted(fit$ddf)[1],fit$ddf$Nhat,true.N,"mt"))
      }else{
        res<-rbind(res,c("cds-hrp-w",par.ind,n.samp,sim,rep(NA,7)))
      }

      ######################################################### 
      ## CDS - unif+cos

      #fit<-try(ds(sim.data,width,monotonicity="strict",key="unif",
      #            adjustment="cos"))
      #if(all(class(fit)!="try-error")){
      #  res<-rbind(res,c("cds-uc",par.ind,n.samp,sim,rep(NA,3),
      #                   fitted(fit$ddf)[1],fit$ddf$criterion,
      #                   true.N,fit$ddf$Nhat))
      #}else{
      #  res<-rbind(res,c("cds-uc",par.ind,n.samp,sim,rep(NA,7)))
      #}
    }
  }
  # write them out
  write.csv(res,file=paste(opt.method,"-",par.ind,"-",n.samp,"-3pt-results.csv",sep=""))
  return(1)
}
