# covariate sim 2 - 1 continuous covariate
library(mmds)
library(msm)
# let's do this in parallel
library(foreach)
library(doMC)
options(cores=2)
registerDoMC()

### SETUP
# mix 1, sigma=exp(beta_00+beta_01*z_2)
# mix 2, sigma=exp(beta_10)
#model.formula<-list()
#model.formula[[1]]<-"~cov1"
#model.formula[[2]]<-"~1"
model.formula<-"~cov1"

### parameters
true.pars<-c(log(0.2),log(0.8),log(0.4), inv.reparam.pi(0.4))
mix.terms<-2
width<-1

# sample sizes
n.samps<-c(30,60,120,480,960)
# number of realisations
n.sims<-200
# set the seed
set.seed(102837)
showit<-0
starting.vals<-NULL # unknown starting values

opt.methods<-c("BFGS+SANN","EM")
opt.methods<-c("BFGS+SANN")

for(opt.method in opt.methods){

   res<-c()

   # loop over sample sizes
   for(n.samples in n.samps){
      # loop over sims
   #   for(sim in 1:n.sims){
      results<-foreach(sim = 1:n.sims, .combine=rbind,
                       .inorder=FALSE, .init=c()) %dopar% {

         this.line<-c()
   
         z<-list(matrix(c(rep(1,n.samples),pnorm(seq(-4,4,len=n.samples))),n.samples,2))
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

         # calculate the true N
         gp<-getpars(true.pars,mix.terms,zdim,z)
         sigma<-gp$key.scale
         pis<-gp$mix.prop
         mu<-apply(sigma,1,integrate.hn,width)
         mus<-mu%*%matrix(pis,length(pis),1)
         pas<-mus/width
         true.N<-sum(1/pas)
   
         if(class(fit)!="try-error"){
            # results...
            this.line<-c(n.samples,sim,fit$aic,fit$pa,fit$N,true.N,fit$mix.terms,"cov")
         }else{
            # if something went wrong
            return(rbind(c(n.samples,sim,NA,NA,NA,NA,NA,NA),
                         c(n.samples,sim,NA,NA,NA,NA,NA,NA)))
         }
         
#### ignoring the covariates
         # because SANN changes the seed, save it first
         seed <- get(".Random.seed",envir=.GlobalEnv) ## store RNG seed
   
         fit<-try(fitmix(testdata,initialvalues=starting.vals,mix.terms=1,ftype="hn",
                            showit=showit,width=width,model.formula="~1",
                            usegrad=TRUE,opt.method=opt.method))
   
         fit<-step.ds.mixture(fit)

         # restore the seed   
         assign(".Random.seed",seed,envir=.GlobalEnv)

         # calculate the true N
         gp<-getpars(true.pars,mix.terms,zdim,z)
         sigma<-gp$key.scale
         pis<-gp$mix.prop
         mu<-apply(sigma,1,integrate.hn,width)
         mus<-mu%*%matrix(pis,length(pis),1)
         pas<-mus/width
         true.N<-sum(1/pas)
   
         if(class(fit)!="try-error"){
            # results...
            this.line<-rbind(this.line,
                             c(n.samples,sim,fit$aic,fit$pa,fit$N,
                               true.N,fit$mix.terms,"nocov"))
         }else{
            # if something went wrong
            return(rbind(c(n.samples,sim,NA,NA,NA,NA,NA,NA),
                         c(n.samples,sim,NA,NA,NA,NA,NA,NA)))
         }

         return(this.line)
      }
      res<-rbind(res,results)
   }
   write.csv(res,file=paste("covsim2-",opt.method,".csv",sep=""))
}


