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
opt.methods<-c("BFGS+SANN","EM")
opt.methods<-c("BFGS+SANN")
# for debug
showit<-0

######################################
# actually do the simulation

## loop over optimisation methods
#for(opt.method in opt.methods){
opt.method<-"BFGS+SANN"
   ## loop over possible parameters
   results<-foreach(par.ind = 1:dim(parmat)[1], .combine=rbind,
                      .inorder=FALSE, .init=c()) %dopar% {
      res<-c()
      res1<-c()

      # make sure that the seed is the same on all cores
      set.seed(123)
   
      ## loop over number of samples to take
      for(n.samp in n.samps){
         prefer<-c()
         ## loop over number of replicates
         for(sim in 1:n.sims){
      
            aic1<-aic2<-NULL

            # simulate the data
            sim.data<-sim.mix(parmat[par.ind,],mix.terms,n.samp,width)
            true.N<-n.samp*width/mu.calc(parmat[par.ind,],mix.terms,width)

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
               res<-rbind(res,c(par.ind,n.samp,sim,fit$pars,fit$pa,
                                 fit$aic,fit$N,true.N))
               aic2<-fit$aic
            }else{
               res<-rbind(res,c(par.ind,n.samp,sim,rep(NA,7)))
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
               res1<-rbind(res1,c(par.ind,n.samp,sim,fit$pars,fit$pa,
                                 fit$aic,fit$N,true.N))
               aic1<-fit$aic
            }else{
               res1<-rbind(res1,c(par.ind,n.samp,sim,rep(NA,7)))
            }

            # which was better, 1 or 2 point?
            if(!is.null(aic1) & !is.null(aic2)){
               prefer<-c(prefer,which.min(c(aic1,aic2)))
            }
         }

         cat("TKTKTKTK",n.samp,sum(prefer==1),sum(prefer==2),"\n")
      }

      # write them out
      write.csv(res,file=paste(opt.method,"-",par.ind,"-",n.samp,"-results.csv",sep=""))
      write.csv(res1,file=paste(opt.method,"-",par.ind,"-",n.samp,"-1-results.csv",sep=""))
      return(1)
   }
#}
