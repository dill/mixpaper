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

# optimisation methods
opt.methods<-c("BFGS+SANN","EM")
opt.methods<-"BFGS+SANN"
# for debug
showit<-0

######################################
# actually do the simulation

## loop over optimisation methods
for(opt.method in opt.methods){
   ## loop over possible parameters
   results<-foreach(par.ind = 1:dim(parmat)[1], .combine=rbind,
                      .inorder=FALSE, .init=c()) %dopar% {
      res<-c()
      res1<-c()
   
      # make sure that the seed is the same on all cores
      set.seed(123)
   
      ## loop over number of samples to take
      for(n.samp in n.samps){
         ## loop over number of replicates
         for(sim in 1:n.sims){
      
            # simulate the data
            sim.data<-sim.mix(parmat[par.ind,],mix.terms,n.samp,width,pt=TRUE)
            # calculate the true abundance
            true.N<-(n.samp*pi*width^2)/mu.calc(parmat[par.ind,],
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
               res<-rbind(res,c(par.ind,n.samp,sim,fit$pars,
                                fit$likelihood,fit$aic,true.N,fit$N))
            }else{
               res<-rbind(res,c(par.ind,n.samp,sim,rep(NA,7)))
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
               res1<-rbind(res1,c(par.ind,n.samp,sim,fit$pars,
                                fit$likelihood,fit$aic,true.N,fit$N))
            }else{
               res1<-rbind(res1,c(par.ind,n.samp,sim,rep(NA,5)))
            }

         }
      }
      # write them out
      write.csv(res,file=paste(opt.method,"-",par.ind,"-pt-results.csv",sep=""))
      write.csv(res1,file=paste(opt.method,"-",par.ind,"-1-pt-results.csv",sep=""))
      return(1)
   }
}
