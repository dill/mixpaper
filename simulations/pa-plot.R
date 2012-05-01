# Pa boxplots
# makes a table of boxplots for all of the simulations with the 
# proportion of time that the true model was selected by AIC printed
# below each box.

# also of interest are table(aic.winners) and table(cov.aic.winners)
# which give the breakdown of the AIC best models according to model
# type.

library(ggplot2)

# pull all of the results together for the paper graphic...

baf<-data.frame(pa=NA,n=NA,id=NA,model=NA)
aic.winners<-data.frame(mix.terms=0,n=0,model=NA,id=0)
true.ps<-c()

samp.sizes<-c(30,60,120,480,960)


### nocov
# no covariate models
#for(selection in c("none")){#,"AIC")){
for(selection in c("AIC")){

  cat("Model selection:",selection,"\n")

  for(par.ind in 1:4){
    for(n.samps in samp.sizes){
    
      fit.method<-"BFGS+SANN"
      
      dat<-read.csv(file=paste("nocov/",fit.method,"-",par.ind,
                                "-960-results.csv",sep=""))
      dat<-dat[,-1]

      names(dat)<-c("model","par.ind","n","sim","par1","par2",
                    "par3","pa","AIC","Nhat","N")

      dat$Nhat<-as.double(as.character(dat$Nhat))
      dat$N<-as.double(as.character(dat$N))

      dat1<-dat[dat$model=="mmds-1",]
      dat<-dat[dat$model=="mmds-c",]

      dat<-dat[dat$n==n.samps,]
      dat1<-dat1[dat1$n==n.samps,]

      #dat1<-dat1[dat$simno,] #????
      ind<- is.na(dat$Nhat) & is.na(dat1$Nhat) & is.na(dat$N) & is.na(dat1$N)
      dat<-dat[!ind,]
      dat1<-dat1[!ind,]

      if(selection=="AIC"){

         aics<-cbind(dat1$AIC,dat$AIC)
         aic.pick<-apply(aics,1,which.min)

         p<-n.samps/cbind(dat1$Nhat,dat$Nhat)

         this.p<-c(p[,1][aic.pick==1],p[,2][aic.pick==2])
    
         aic.winners<-rbind(aic.winners,
                        data.frame(mix.terms=aic.pick,
                                   n=rep(n.samps,length(aic.pick)),
                                   model=rep("nocov",length(aic.pick)),
                                   id=rep(par.ind,length(aic.pick))))
      }else{

         this.p<-n.samps/dat$Nhat
      }

      baf<-rbind(baf,data.frame(pa=this.p,
                                n=rep(n.samps,length(this.p)),
                                id=rep(par.ind,length(this.p)),
                                model=rep("nocov",length(this.p))
                 ))
      
      true.ps<-c(true.ps,n.samps/as.double(as.character(dat$N[!is.na(dat$N)])))
    }
  }
}
# truth lines -- nocov
true.p<-data.frame(t=rep(unique(round(true.ps,6)),4),
                   id=rep(1:4,4),
                   model=rep("nocov",16))




true.ps<-c()
# point transects
for(selection in c("AIC")){

  cat("Model selection:",selection,"\n")

  for(par.ind in 1:4){
    for(n.samps in samp.sizes){

      n.samps<-as.integer(n.samps)
    
      fit.mthod<-"BFGS+SANN"
      dat<-read.csv(file=paste("pt/",fit.method,"-",par.ind,"-pt-results.csv",sep=""))
    
      dat<-dat[,-1]
      names(dat)<-c("model","par.ind","n.samps","sim","par1","par2",
                    "par3","pall","AIC","N","Nhat")

      dat1<-dat[dat$model=="mmds-1",]
      dat<-dat[dat$model=="mmds-2",]

      dat<-dat[dat$n==n.samps,]
      dat1<-dat1[dat1$n==n.samps,]

      #dat1<-dat1[dat$simno,] # ???
      ind<- is.na(dat$Nhat) & is.na(dat1$Nhat)
      dat<-dat[!ind,]
      dat1<-dat1[!ind,]

      if(selection=="AIC"){

        aics<-cbind(dat1$AIC,dat$AIC)
        aic.pick<-apply(aics,1,which.min)

        p<-n.samps/cbind(dat1$Nhat,dat$Nhat)

        this.p<-c(p[,1][aic.pick==1],p[,2][aic.pick==2])
    
        aic.winners<-rbind(aic.winners,
                       data.frame(mix.terms=aic.pick,
                                  n=rep(n.samps,length(aic.pick)),
                                  model=rep("pt",length(aic.pick)),
                                  id=rep(par.ind,length(aic.pick))))
      }else{
        this.p<-n.samps/dat$Nhat
      }

      baf<-rbind(baf,data.frame(pa=this.p,
                                n=rep(n.samps,length(this.p)),
                                id=rep(par.ind,length(this.p)),
                                model=rep("pt",length(this.p))
                 ))
      
      true.ps<-c(true.ps,n.samps/as.double(as.character(dat$N[!is.na(dat$N)])))
    }
  }
}
# truth lines -- nocov
true.p<-rbind(true.p,data.frame(t=rep(unique(round(true.ps,6)),4),
                   id=rep(1:4,4),
                   model=rep("pt",16)))


### covar
# covariates with model selection
# here models were fit with and without covariates
# both covariate and noncovariate fits were stored, covariate as "covar"
# and non-covariate as "nocov"

true.ps<-c()
true.pp<-c()

# separate cov table stuff to see what's going on
cov.aic.winners<-data.frame(mix.terms=0,n=0,model=NA,id=0)


for(par.ind in 1:2){
  for(n.samps in samp.sizes){
  
    fit.method<-"BFGS+SANN"
    # read in dat
    # format of the file:
    # par index, no. samples, sim no., par est(3),likelihood, aic, pa
    dat<-read.csv(file=paste("covar/covsim",par.ind,"-",fit.method,".csv",sep=""))
  
    dat<-dat[,-1]
    names(dat)<-c("n.samps","sim","AIC","pa","Nhat","N","mix.terms","mod")
  
    dat<-dat[dat$n==n.samps,]

    ind<- !is.na(dat$Nhat) & !is.na(dat$N)
    dat<-dat[ind,]


    # do the selection
    sel<-TRUE
    if(sel){
      dat.cov<-dat[dat$mod=="cov",]
      dat.nocov<-dat[dat$mod=="nocov",]

      pa.res<-n.samps/cbind(dat.nocov$Nhat,dat.cov$Nhat)

      aic.res<-cbind(dat.nocov$AIC,dat.cov$AIC)
      aic.pick<-apply(aic.res,1,which.min)

      mixs<-cbind(dat.nocov$mix.terms,dat.cov$mix.terms)
      mixs2<-rep(0,nrow(mixs))
      mixs2[aic.pick==1]<-dat.nocov$mix.terms[aic.pick==1]
      mixs2[aic.pick==2]<-dat.cov$mix.terms[aic.pick==2]
      aic.pick2<-aic.pick
      aic.pick2[aic.pick2==1]<-paste("nocov",mixs2[aic.pick2==1])
      aic.pick2[aic.pick2==2]<-paste("covar",mixs2[aic.pick2==2])

      cov.aic.winners<-rbind(cov.aic.winners,
                             data.frame(mix.terms=aic.pick2,
                                        n=rep(n.samps,length(aic.pick)),
                                        model=rep(1,length(mixs2)),
                                        id=rep(par.ind,length(aic.pick))))

      # make aic.winners be right vs wrong model...
      # recode cov model when the number of mix terms != 2 
      # to be "1"
      ind.aic<-aic.pick==2
      mt<-dat.cov$mix.terms[ind.aic]
      mt[mt!=2]<-1
      aic.pick[ind.aic]<-mt

      aic.winners<-rbind(aic.winners,
                     data.frame(mix.terms=aic.pick,
                                n=rep(n.samps,length(aic.pick)),
                                model=rep("covar",length(aic.pick)),
                                id=rep(par.ind,length(aic.pick))))

      this.p<-c(pa.res[,1][aic.pick==1],pa.res[,2][aic.pick==2])

    }else{
      this.p<-n.samps/dat$Nhat
    }

    baf<-rbind(baf,data.frame(pa=this.p,
                              n=rep(n.samps,length(this.p)),
                              id=rep(par.ind,length(this.p)),
                              model=rep("covar",length(this.p))
               ))
    
    true.pp<-c(true.pp,n.samps/dat$N)

  }
  true.ps<-c(true.ps,median(true.pp))
  true.pp<-c()
}
# truth lines -- covar
true.p<-rbind(true.p,data.frame(t=true.ps,
                                id=1:2,
                                model=rep("covar",2)))


### 3 point
true.ps<-c()

for(par.ind in 1:2){
  for(n.samps in samp.sizes){

    n.samps<-as.integer(n.samps)
  
    fit.method<-"BFGS+SANN"
    # read in dat
    dat<-read.csv(file=paste("3point/",fit.method,"-",par.ind,"-960-3pt-results.csv",sep=""))
  
    dat<-dat[,-1]
    names(dat)<-c("model","par.ind","n.samps","sim","ll","aic",
                  "pa","Nhat","N","mix.terms")

    dat<-dat[dat$model=="mmds-MS",]
    dat<-dat[dat$n==n.samps,]
    dat$Nhat<-as.double(as.character(dat$Nhat))
    ind<- !is.na(dat$Nhat)
    dat<-dat[ind,]

    this.p<-n.samps/dat$Nhat

    aic.winners<-rbind(aic.winners,
                   data.frame(mix.terms=dat$mix.terms,
                              n=rep(n.samps,length(dat$mix.terms)),
                              model=rep("3pt",length(dat$mix.terms)),
                              id=rep(par.ind,length(dat$mix.terms))))
    

    baf<-rbind(baf,data.frame(pa=this.p,
                              n=rep(n.samps,length(this.p)),
                              id=rep(par.ind,length(this.p)),
                              model=rep("3pt",length(this.p))
               ))
    
#    true.ps<-c(true.ps,n.samps/dat$N)
    true.ps<-c(true.ps,n.samps/as.double(as.character(dat$N[!is.na(dat$N)])))
  }
}
# truth lines -- covar
true.p<-rbind(true.p,data.frame(t=unique(round(true.ps,6)),
                   id=1:2,
                   model=rep("3pt",2)))



baf<-baf[-1,]
#baf<-as.data.frame(baf)
#names(baf)<-c("pa","n","id","model")
#baf$pa<-as.numeric(as.character(baf$pa))
baf$model<-as.factor(baf$model)


# make the annotations of the proportion of results that had
# the true number of mixture terms
aic.winners<-aic.winners[-1,]
cov.aic.winners<-cov.aic.winners[-1,]

#(in true model)
itm<-data.frame(n=0,prop=0,model=NA,id=0)

# nocov
for(i in 1:4){
   # calculate the proportions
   this.prop<-table(aic.winners)[2,,"nocov",i]/
                (table(aic.winners)[1,,"nocov",i]+table(aic.winners)[2,,"nocov",i])
   itm<-rbind(itm,cbind(n=samp.sizes,
                        prop=round(this.prop,2),
                        model=rep("nocov",5),
                        id=rep(i,5)))
}
# pt
for(i in 1:4){
   # calculate the proportions
   this.prop<-table(aic.winners)[2,,"pt",i]/
                (table(aic.winners)[1,,"pt",i]+table(aic.winners)[2,,"pt",i])
   itm<-rbind(itm,cbind(n=samp.sizes,
                        prop=round(this.prop,2),
                        model=rep("pt",5),
                        id=rep(i,5)))
}
# covar
for(i in 1:2){
   # calculate the proportions
   this.prop<-table(aic.winners)[2,,"covar",i]/
                (table(aic.winners)[1,,"covar",i]+table(aic.winners)[2,,"covar",i])
   itm<-rbind(itm,cbind(n=samp.sizes,
                        prop=round(this.prop,2),
                        model=rep("covar",5),
                        id=rep(i,5)))
}
# 3pt
for(i in 1:2){
   # calculate the proportions
   this.prop<-table(aic.winners)[3,,"3pt",i]/
                (table(aic.winners)[1,,"3pt",i]+
                 table(aic.winners)[2,,"3pt",i]+
                 table(aic.winners)[3,,"3pt",i])
   itm<-rbind(itm,cbind(n=samp.sizes,
                        prop=round(this.prop,2),
                        model=rep("3pt",5),
                        id=rep(i,5)))
}



itm<-itm[-1,]

# reorder the models
baf$model<-factor(baf$model,levels=c("nocov","pt","3pt","covar"))
itm$model<-factor(itm$model,levels=c("nocov","pt","3pt","covar"))
true.p$model<-factor(true.p$model,levels=c("nocov","pt","3pt","covar"))

# rename the models
levels(baf$model)<-c("No covariates","Point transect","3-point","Covariate")
levels(itm$model)<-c("No covariates","Point transect","3-point","Covariate")
levels(true.p$model)<-c("No covariates","Point transect","3-point","Covariate")


# actually do the plotting here
p<-ggplot(baf,aes(x=factor(n),y=pa))
p<-p+geom_boxplot(outlier.size=1)
p<-p+facet_grid(model~id)
#p<-p+facet_grid(.~.)
#p<-p+geom_text(aes(x=factor(n),y=-0.1,label=prop),size=3,data=itm)


p<-p+labs(x="Sample size",y="Probability of detection")#,fill="Fitting algorithm")
p<-p+opts(panel.grid.major=theme_blank(),
          panel.grid.minor=theme_blank(),
          legend.background=theme_blank(),
          legend.key=theme_blank(),
          panel.background=theme_rect())

# transparency not supported by eps!
#p<-p+geom_hline(aes(yintercept=t),true.p,col="red",alpha=0.4)
p<-p+geom_hline(aes(yintercept=t),true.p,col="grey")
p


#dev.copy2eps(file="pa-plot.eps")



# second plot of AIC winners
#pp<-ggplot(aic.winners)
##pp<-pp+geom_histogram(aes(y=mix.terms,x=n),binwidth=1)
#pp<-pp+geom_histogram(aes(factor(mix.terms)),binwidth=1)
#pp<-pp+facet_grid(n~id)
#pp<-pp+opts(panel.grid.major=theme_blank(),
#          panel.grid.minor=theme_blank(),
#          legend.background=theme_blank(),
#          legend.key=theme_blank(),
#          panel.background=theme_rect())
#pp



