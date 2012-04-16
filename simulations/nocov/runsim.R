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
# for debug
showit<-0
######################################
# actually do the simulation

opt.method<-"BFGS+SANN"
## loop over possible parameters
results<-foreach(par.ind = 1:dim(parmat)[1], .combine=rbind,
                  .inorder=FALSE, .init=c()) %dopar% {
  res<-c()

  # make sure that the seed is the same on all cores
  set.seed(123)

  ## loop over number of samples to take
  for(n.samp in n.samps){
    prefer<-c()
    ## loop over number of replicates
    for(sim in 1:n.sims){
  
      # simulate the data
      sim.data<-sim.mix(parmat[par.ind,],mix.terms,n.samp,width)
      true.N<-n.samp*width/mmds:::mu.calc(parmat[par.ind,],mix.terms,width)

      ######################################################## 
      # fit the "correct" model 
      # because SANN changes the seed, save it first
      seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed

      fit<-try(fitmix(sim.data,initialvalues=starting.vals,mix.terms=mix.terms,
                 ftype="hn",width=width,model.formula=model.formula,
                 usegrad=TRUE,opt.method=opt.method,showit=showit))
    
      # restore the seed   
      assign(".Random.seed",seed,envir=.GlobalEnv) 

      if(class(fit)!="try-error"){

        # result vector -
        res<-rbind(res,c("mmds-c",par.ind,n.samp,sim,fit$pars,fit$pa,
                          fit$aic,fit$N,true.N))
      }else{
        res<-rbind(res,c("mmds-c",par.ind,n.samp,sim,rep(NA,7)))
      }

      ######################################################## 
      # 1-point fit
      seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed

      fit<-try(fitmix(sim.data,initialvalues=starting.vals,mix.terms=1,
                 ftype="hn",width=width,model.formula=model.formula,
                 usegrad=TRUE,opt.method=opt.method,showit=showit))
    
      # restore the seed   
      assign(".Random.seed",seed,envir=.GlobalEnv) 

      # result vector - as above but with 1 par
      if(class(fit)!="try-error"){
        res<-rbind(res,c("mmds-1",par.ind,n.samp,sim,fit$pars,NA,NA,fit$pa,
                          fit$aic,fit$N,true.N))
      }else{
        res<-rbind(res,c("mmds-1",par.ind,n.samp,sim,rep(NA,8)))
      }

      ######################################################## 
      # CDS - hn+cos

      fit<-try(ds(sim.data,width,monotonicity="strict"))
      if(all(class(fit$ddf)!="try-error")){
        res<-rbind(res,c("cds-hnc",par.ind,n.samp,sim,rep(NA,3),
                         fitted(fit$ddf)[1],fit$ddf$criterion,
                         true.N,fit$ddf$Nhat))
      }else{
        res<-rbind(res,c("cds-hnc",par.ind,n.samp,sim,rep(NA,7)))
      }

      ######################################################## 
      # CDS - hr+poly

      fit<-try(ds(sim.data,width,monotonicity="strict",key="hr",
                  adjustment="poly"))
      if(all(class(fit$ddf)!="try-error")){
        res<-rbind(res,c("cds-hrp",par.ind,n.samp,sim,rep(NA,3),
                         fitted(fit$ddf)[1],fit$ddf$criterion,
                         true.N,fit$ddf$Nhat))
      }else{
        res<-rbind(res,c("cds-hrp",par.ind,n.samp,sim,rep(NA,7)))
      }

      ######################################################### 
      ## CDS - unif+cos

      #fit<-try(ds(sim.data,width,monotonicity="strict",key="unif",
      #            adjustment="cos"))
      #if(all(class(fit$ddf)!="try-error")){
      #  res<-rbind(res,c("cds-uc",par.ind,n.samp,sim,rep(NA,3),
      #                   fitted(fit$ddf)[1],fit$ddf$criterion,
      #                   true.N,fit$ddf$Nhat))
      #}else{
      #  res<-rbind(res,c("cds-uc",par.ind,n.samp,sim,rep(NA,7)))
      #}

    }# sim block end
  }# foreach block end

  # write them out
  write.csv(res,file=paste(opt.method,"-",par.ind,"-",n.samp,"-results.csv",sep=""))
  return(1)
}
