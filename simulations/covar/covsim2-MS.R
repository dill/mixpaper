# covariate sim 2
# continuous
library(mmds)
library(msm)
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

     # because SANN changes the seed, save it first
     seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed
   
     fit<-try(fitmix(testdata,initialvalues=starting.vals,mix.terms=1,ftype="hn",
                 showit=showit,width=width,model.formula=model.formula,
                 usegrad=TRUE,opt.method=opt.method))

     fit<-step.ds.mixture(fit)

     # restore the seed   
     assign(".Random.seed",seed,envir=.GlobalEnv)

     if(class(fit)!="try-error"){
       # calculate the true N
       gp<-mmds:::getpars(true.pars,mix.terms,zdim,z)
       sigma<-gp$key.scale
       pis<-gp$mix.prop
       mu<-apply(sigma,1,mmds:::integrate.hn,width)
       mus<-mu%*%matrix(pis,length(pis),1)
       pas<-mus/width
       true.N<-sum(1/pas)

       # results...
       res<-c(n.samples,sim,fit$aic,fit$pa,fit$N,
                    true.N,fit$mix.terms,"cov")
     }else{
       # if it failed...
       res<-c(n.samples,sim,NA,NA,NA,NA,NA,"cov")
     }

     #######################################################
     # with no covariates
     # because SANN changes the seed, save it first
     seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed
   
     fit<-try(fitmix(testdata,initialvalues=starting.vals,mix.terms=1,ftype="hn",
                 showit=showit,width=width,model.formula="~1",
                 usegrad=TRUE,opt.method=opt.method))

     fit<-step.ds.mixture(fit)

     # restore the seed   
     assign(".Random.seed",seed,envir=.GlobalEnv)

     if(class(fit)!="try-error"){
       # calculate the true N
       gp<-mmds:::getpars(true.pars,mix.terms,zdim,z)
       sigma<-gp$key.scale
       pis<-gp$mix.prop
       mu<-apply(sigma,1,mmds:::integrate.hn,width)
       mus<-mu%*%matrix(pis,length(pis),1)
       pas<-mus/width
       true.N<-sum(1/pas)

       # results...
       res<-rbind(res,c(n.samples,sim,fit$aic,fit$pa,fit$N,
                    true.N,fit$mix.terms,"nocov"))
     }else{
       # if it failed...
       res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"nocov"))
     }

     # CDS analyses

     ######################################################## 
     # CDS - hn+cos

     fit<-try(ds(testdata,width,monotonicity="strict"))
     if(all(class(fit$ddf)!="try-error")){
       res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
                        n.samples/sum(1/fitted(fit$ddf)),
                        fit$ddf$Nhat,true.N,NA,"hn+cos"))
     }else{
       res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"hn+cos"))
     }

     ######################################################## 
     # CDS - hr+poly

     fit<-try(ds(testdata,width,monotonicity="strict",key="hr",
                 adjustment="poly"))
     if(all(class(fit$ddf)!="try-error")){
       res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
                        n.samples/sum(1/fitted(fit$ddf)),
                        fit$ddf$Nhat,true.N,NA,"hr+poly"))
     }else{
       res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"hr+poly"))
     }

     ######################################################### 
     ## CDS - unif+cos

     #fit<-try(ds(testdata,width,monotonicity="strict",key="unif",
     #            adjustment="cos"))
     #if(all(class(fit$ddf)!="try-error")){
     #  res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
     #                   n.samples/sum(1/fitted(fit$ddf)),
     #                   fit$ddf$Nhat,true.N,NA,"unif+cos"))
     #}else{
     #  res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"unif+cos"))
     #}

     # MCDS analyses
     
     ######################################################## 
     # MCDS - hn+cos + as.factor(cov1)

     fit<-try(ds(testdata,width,formula=as.formula(model.formula)))
     if(all(class(fit$ddf)!="try-error")){
       res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
                        n.samples/sum(1/fitted(fit$ddf)),
                        fit$ddf$Nhat,true.N,NA,"hn+cos+cov1"))
     }else{
       res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"hn+cos+cov1"))
     }

     ######################################################## 
     # MCDS - hr+poly + as.factor(cov1)

     fit<-try(ds(testdata,width,formula=as.formula(model.formula),key="hr",
                 adjustment="poly"))
     if(all(class(fit$ddf)!="try-error")){
       res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
                        n.samples/sum(1/fitted(fit$ddf)),
                        fit$ddf$Nhat,true.N,NA,"hr+poly+cov1"))
     }else{
       res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"hr+poly+cov1"))
     }

     ######################################################### 
     ## MCDS - unif+cos + as.factor(cov1)

     #fit<-try(ds(testdata,width,formula=as.formula(model.formula),key="unif",
     #            adjustment="cos",scale="width"))
     #if(all(class(fit$ddf)!="try-error")){
     #  res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
     #                   n.samples/sum(1/fitted(fit$ddf)),
     #                   fit$ddf$Nhat,true.N,NA,"unif+cos+cov1-width"))
     #}else{
     #  res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"unif+cos+cov1-width"))
     #}


     ######################################################## 
     # MCDS - hn+cos + as.factor(cov1) width scaling

     fit<-try(ds(testdata,width,formula=as.formula(model.formula),scale="width"))
     if(all(class(fit$ddf)!="try-error")){
       res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
                        n.samples/sum(1/fitted(fit$ddf)),
                        fit$ddf$Nhat,true.N,NA,"hn+cos+cov1-width"))
     }else{
       res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"hn+cos+cov1-width"))
     }

     ######################################################## 
     # MCDS - hr+poly + as.factor(cov1) width scaling

     fit<-try(ds(testdata,width,formula=as.formula(model.formula),key="hr",
                 adjustment="poly",scale="width"))
     if(all(class(fit$ddf)!="try-error")){
       res<-rbind(res,c(n.samples,sim,fit$ddf$criterion,
                        n.samples/sum(1/fitted(fit$ddf)),
                        fit$ddf$Nhat,true.N,NA,"hr+poly+cov1-width"))
     }else{
       res<-rbind(res,c(n.samples,sim,NA,NA,NA,NA,NA,"hr+poly+cov1-width"))
     }

write.table(res,file=paste("covsim2-",n.samples,"-",opt.method,".csv",sep=""),append=TRUE,col.names=FALSE)
     return(res)
  }

  big.res<-rbind(big.res,results)
}
#write.csv(big.res,file=paste("covsim2-",opt.method,".csv",sep=""))
