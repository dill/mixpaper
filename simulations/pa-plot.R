# Pa boxplots
# makes a table of boxplots for all of the simulations with the 
# proportion of time that the true model was selected by AIC printed
# below each box.

# also of interest is table(aic.winners)
# which give the breakdown of the AIC best models according to model
# type.

library(ggplot2)

model.names <- c("A. No covariates","B. Point transect",
                 "C. 3-point","D. Covariate","E. Hazard mixture")

samp.sizes<-c(30,60,120,480,960)

#for(set in c("mmds")){#,"cds","combined")){
#for(set in c("combined")){
for(set in c("mmds","cds","combined")){
  baf<-data.frame(pa=NA,n=NA,id=NA,model=NA)
  aic.winners<-data.frame(mix.terms=0,n=0,model=NA,id=0)
  true.ps<-c()

  ### nocov
  for(par.ind in 1:4){
    for(n.samps in samp.sizes){

      dat<-read.csv(file=paste("nocov/BFGS+SANN-",par.ind,
                                "-960-results.csv",sep=""))
      dat<-dat[,-1]

      names(dat)<-c("model","par.ind","n","sim","par1","par2",
                    "par3","pa","AIC","Nhat","N")

      dat<-dat[dat$n==n.samps,]

      dat$Nhat<-as.double(as.character(dat$Nhat))
      dat$N<-as.double(as.character(dat$N))

      if(set=="mmds"){
        models<-c("mmds-1","mmds-c")
      }else if(set=="cds"){
        models<-c("cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w")
      }else{
        models<-c("mmds-1","mmds-c","cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w")
      }

      pa.res<-matrix(NA,200,length(models))
      aic.res<-matrix(NA,200,length(models))

      for(modi in seq_along(models)){
        if(sum(dat$mod==models[modi])>0){
          pa.res[,modi][dat$sim[dat$mod==models[modi]]]<-n.samps/
                                              dat$Nhat[dat$mod==models[modi]]
          aic.res[,modi][dat$sim[dat$mod==models[modi]]]<-
                                              dat$AIC[dat$mod==models[modi]]
        }
      }

      aic.res[is.na(aic.res)]<-Inf

      aic.pick<-apply(aic.res,1,which.min)
      pa.cols<-ncol(pa.res)
      pa.aic<-cbind(pa.res,aic.pick)
      this.p<-apply(pa.aic,1,function(x){x[x[pa.cols+1]]})

      ind <- !is.na(this.p)
      this.p <- this.p[ind]
      aic.pick <- aic.pick[ind]

      if(set=="mmds" | set == "combined"){
        aic.winners<-rbind(aic.winners,
                       data.frame(mix.terms=aic.pick,
                                  n=rep(n.samps,length(aic.pick)),
                                  model=rep("nocov",length(aic.pick)),
                                  id=rep(par.ind,length(aic.pick))))
      }
      baf<-rbind(baf,data.frame(pa=this.p,
                                n=rep(n.samps,length(this.p)),
                                id=rep(par.ind,length(this.p)),
                                model=rep("nocov",length(this.p))
                 ))

      true.ps<-c(true.ps,n.samps/dat$N[!is.na(dat$Nhat)])
    }
  }

  # truth lines -- nocov
  true.p<-data.frame(t=rep(unique(round(true.ps,6)),4),
                     id=rep(1:4,4),
                     model=rep("nocov",16))

  true.ps<-c()

  ################################################################
  # point transects

  for(par.ind in 1:4){
    for(n.samps in samp.sizes){

      n.samps<-as.integer(n.samps)

      fit.mthod<-"BFGS+SANN"
      dat<-read.csv(file=paste("pt/BFGS+SANN-",par.ind,
                               "-pt-results.csv",sep=""))

      dat<-dat[,-1]
      names(dat)<-c("model","par.ind","n.samps","sim","par1","par2",
                    "par3","pall","AIC","N","Nhat")

      dat<-dat[dat$n==n.samps,]
      dat$Nhat<-as.double(as.character(dat$Nhat))
      dat$N<-as.double(as.character(dat$N))

      if(set=="mmds"){
        models<-c("mmds-1","mmds-2")
      }else if(set=="cds"){
        models<-c("cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w")
      }else{
        models<-c("mmds-1","mmds-2","cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w")
      }

      pa.res<-matrix(NA,200,length(models))
      aic.res<-matrix(NA,200,length(models))

      for(modi in seq_along(models)){
        if(sum(dat$mod==models[modi])>0){
          pa.res[,modi][dat$sim[dat$mod==models[modi]]]<-n.samps/
                                              dat$Nhat[dat$mod==models[modi]]
          aic.res[,modi][dat$sim[dat$mod==models[modi]]]<-
                                              dat$AIC[dat$mod==models[modi]]
        }
      }

      # one of the results is way out for CDS, just remove that p
      if(set=="cds" & n.samps==30){
        ind<-FALSE
        for(this.col in 1:ncol(pa.res)){
          ind<-ind|pa.res[,this.col]<6
        }
        pa.res<-pa.res[ind,]
        aic.res<-aic.res[ind,]
      }
      aic.res[is.na(aic.res)]<-Inf

      aic.pick<-apply(aic.res,1,which.min)
      pa.cols<-ncol(pa.res)
      pa.aic<-cbind(pa.res,aic.pick)
      this.p<-apply(pa.aic,1,function(x){x[x[pa.cols+1]]})

      if(set=="mmds" | set=="combined"){
        aic.winners<-rbind(aic.winners,
                       data.frame(mix.terms=aic.pick,
                                  n=rep(n.samps,length(aic.pick)),
                                  model=rep("pt",length(aic.pick)),
                                  id=rep(par.ind,length(aic.pick))))
      }
      baf<-rbind(baf,data.frame(pa=this.p,
                                n=rep(n.samps,length(this.p)),
                                id=rep(par.ind,length(this.p)),
                                model=rep("pt",length(this.p))
                 ))

      true.ps<-c(true.ps,n.samps/dat$N[!is.na(dat$N)])
    }
  }
  # truth lines
  true.p<-rbind(true.p,data.frame(t=rep(unique(round(true.ps,6)),4),
                     id=rep(1:4,4),
                     model=rep("pt",16)))

  ##################################################################
  ### covar
  true.ps<-c()

  for(par.ind in 1:2){
    true.pp<-c()
    for(n.samps in samp.sizes){

      dat<-read.csv(file=paste("covar/covsim",par.ind,"-BFGS+SANN.csv",sep=""))
      dat<-dat[,-1]
      names(dat)<-c("n","sim","AIC","pa","Nhat","N","mix.terms","mod")

      # load the covar data with no adjustments
      dat2<-read.csv(file=paste("covar/covsim",par.ind,
                                "-noadj-BFGS+SANN.csv",sep=""))
      dat2<-dat2[,-1]
      names(dat2)<-c("n","sim","AIC","pa","Nhat","N","mix.terms","mod")

      # bind that on the end and delete the data
      dat<-rbind(dat,dat2)
      rm(dat2)

      dat<-dat[dat$n==n.samps,]

      if(set == "mmds"){
        models<-c("nocov","cov")
      }else if(set=="cds"){
        #models<-c("hn+cos","hr+poly","hn+cos+cov1","hr+poly+cov1",
        #          "hn+cos+cov1-width","hr+poly+cov1-width")
        models<-c("hn+cos","hr+poly","hn+cov1","hr+cov1")
      }else{
        #models<-c("nocov","cov","hn+cos","hr+poly","hn+cos+cov1","hr+poly+cov1",
        #          "hn+cos+cov1-width","hr+poly+cov1-width")
        models<-c("nocov","cov","hn+cos","hr+poly","hn+cov1","hr+cov1")
      }

      pa.res<-matrix(NA,200,length(models))
      aic.res<-matrix(NA,200,length(models))

      for(modi in seq_along(models)){
        if(sum(dat$mod==models[modi])>0){
          pa.res[,modi][dat$sim[dat$mod==models[modi]]]<-n.samps/
                                              dat$Nhat[dat$mod==models[modi]]
          aic.res[,modi][dat$sim[dat$mod==models[modi]]]<-
                                              dat$AIC[dat$mod==models[modi]]
        }
      }

      aic.pick<-apply(aic.res,1,which.min)

      if(set=="mmds" | set=="combined"){


        if(set=="combined"){
          # actually do a recode here...
          # models where # mixtures >1 cov and nocov
          ind.aic<-aic.pick==2
          dat.cov<-dat[dat$mod=="cov",]
          mt <- as.numeric(dat.cov$mix.terms[ind.aic])
          mt[mt==1]<-4
          aic.pick[ind.aic]<-mt

          ind.aic<-aic.pick==1
          dat.cov<-dat[dat$mod=="nocov",]
          mt<- as.numeric(dat.cov$mix.terms[ind.aic])
          mt[mt==1]<-3
          aic.pick[ind.aic]<-mt
        }

        aic.winners<-rbind(aic.winners,
                       data.frame(mix.terms=aic.pick,
                                  n=rep(n.samps,length(aic.pick)),
                                  model=rep("covar-norecode",length(aic.pick)),
                                  id=rep(par.ind,length(aic.pick))))

        # make aic.winners be right vs wrong model...
        # recode cov model when the number of mix terms != 2 
        # to be "1"
        ind.aic<-aic.pick==2
        dat.cov<-dat[dat$mod=="cov",]
        mt<-dat.cov$mix.terms[ind.aic]
        mt[mt!=2]<-1
        aic.pick[ind.aic]<-mt

        aic.winners<-rbind(aic.winners,
                       data.frame(mix.terms=aic.pick,
                                  n=rep(n.samps,length(aic.pick)),
                                  model=rep("covar",length(aic.pick)),
                                  id=rep(par.ind,length(aic.pick))))
      }


      # pick the best -- I really like this idiom
      pa.cols<-ncol(pa.res)
      pa.aic<-cbind(pa.res,aic.pick)
      this.p<-apply(pa.aic,1,function(x){x[x[pa.cols+1]]})

      if(set!="mmds"){
        baf<-rbind(baf,data.frame(pa=this.p,
                                  n=rep(n.samps,length(this.p)),
                                  id=rep(par.ind,length(this.p)),
                                  model=rep("covar",length(this.p))
                   ))
      }else{
        # covariates observed
        baf<-rbind(baf,data.frame(pa=pa.res[,1],
                                  n=rep(n.samps,length(this.p)),
                                  id=rep(par.ind,length(this.p)),
                                  model=rep("covar",length(this.p))
                   ))
        # covariates not observed
        baf<-rbind(baf,data.frame(pa=pa.res[,2],
                                  n=rep(n.samps,length(this.p)),
                                  id=rep(par.ind+2,length(this.p)),
                                  model=rep("covar",length(this.p))
                   ))

      }
    }
    true.pp<-c(true.pp,n.samps/dat$N[!is.na(dat$Nhat)])
    true.ps<-c(true.ps,median(true.pp))
  }
  # truth lines -- covar
  if(set!="mmds"){
    true.p<-rbind(true.p,data.frame(t=true.ps,
                                    id=1:2,
                                    model=rep("covar",2)))
  }else{
    true.p<-rbind(true.p,data.frame(t=rep(true.ps,c(2,2)),
                                    id=1:4,
                                    model=rep("covar",4)))
    baf$id[baf$model=="covar" & baf$id==2] <- 5
    baf$id[baf$model=="covar" & baf$id==3] <- 2
    baf$id[baf$model=="covar" & baf$id==5] <- 3
  }

  ##################################################################
  ### 3 point
  true.ps<-c()

  for(par.ind in 1:2){
    for(n.samps in samp.sizes){

      n.samps<-as.integer(n.samps)

      dat<-read.csv(file=paste("3point/BFGS+SANN-",par.ind,
                               "-960-3pt-results.csv",sep=""))

      dat<-dat[,-1]
      names(dat)<-c("model","par.ind","n.samps","sim","ll","AIC",
                    "pa","Nhat","N","mix.terms")

      dat<-dat[dat$n==n.samps,]
      dat$Nhat<-as.double(as.character(dat$Nhat))
      dat$N<-as.double(as.character(dat$N))

      if(set=="mmds"){
        dat1<-dat[dat$model=="mmds-MS",]

        this.p<-n.samps/dat1$Nhat
        aic.pick<-rep(1,200)

      }else if(set=="cds"){
        models<-c("cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w")

        pa.res<-matrix(NA,200,length(models))
        aic.res<-matrix(NA,200,length(models))

        for(modi in seq_along(models)){
          if(sum(dat$mod==models[modi])>0){
            pa.res[,modi][dat$sim[dat$mod==models[modi]]]<-n.samps/
                                                dat$Nhat[dat$mod==models[modi]]
            aic.res[,modi][dat$sim[dat$mod==models[modi]]]<-
                                                dat$AIC[dat$mod==models[modi]]
          }
        }
      }else{

        models<-c("mmds-MS","cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w")

        pa.res<-matrix(NA,200,length(models))
        aic.res<-matrix(NA,200,length(models))

        for(modi in seq_along(models)){
          if(sum(dat$mod==models[modi])>0){
            pa.res[,modi][dat$sim[dat$mod==models[modi]]]<-n.samps/
                                                dat$Nhat[dat$mod==models[modi]]
            aic.res[,modi][dat$sim[dat$mod==models[modi]]]<-
                                                dat$AIC[dat$mod==models[modi]]
          }
        }
      }

      if(set=="combined" | set=="cds"){
        aic.res[is.na(aic.res)]<-Inf
        aic.pick<-apply(aic.res,1,which.min)
        pa.cols<-ncol(pa.res)
        pa.aic<-cbind(pa.res,aic.pick)
        this.p<-apply(pa.aic,1,function(x){x[x[pa.cols+1]]})
      }

      if(set=="mmds" | set=="combined"){

        # without recoding
        aic.winners<-rbind(aic.winners,
                       data.frame(mix.terms=aic.pick,
                                  n=rep(n.samps,length(aic.pick)),
                                  model=rep("3pt-norecode",length(aic.pick)),
                                  id=rep(par.ind,length(aic.pick))))

        # since we used model selection, recode mix.terms!=3
        # as the wrong model
        mt<-rep(NA,200)
        mt[dat$sim[dat$mod=="mmds-MS"]] <- dat$mix.terms[dat$mod=="mmds-MS"]
        aic.pick2<-aic.pick
        aic.pick2[mt!=3 & aic.pick==1]<-4

        # calculate the aic winners...
        aic.winners<-rbind(aic.winners,
                       data.frame(mix.terms=aic.pick2,
                                  n=rep(n.samps,length(aic.pick2)),
                                  model=rep("3pt",length(aic.pick2)),
                                  id=rep(par.ind,length(aic.pick2))))
      }

      baf<-rbind(baf,data.frame(pa=this.p,
                                n=rep(n.samps,length(this.p)),
                                id=rep(par.ind,length(this.p)),
                                model=rep("3pt",length(this.p))
                 ))

      true.ps<-c(true.ps,
                 n.samps/dat$N[!is.na(dat$Nhat)])
    }
  }

  true.p<-rbind(true.p,data.frame(t=unique(round(true.ps,6)),
                     id=1:2,
                     model=rep("3pt",2)))

  ##################################
  # hazard rate mixtures
  for(par.ind in 1:2){
    for(n.samps in samp.sizes){

      n.samps<-as.integer(n.samps)

      dat<-read.csv(file=paste("hazard/hr-",par.ind,
                               "-results.csv",sep=""))

      dat<-dat[,-1]
      names(dat)<-c("mod","par.ind","n.samp","sim","ll","AIC",
                    "pa","Nhat","N","mixterms")

      dat<-dat[dat$n==n.samps,]
      dat$Nhat<-as.double(as.character(dat$Nhat))
      dat$pa<-as.double(as.character(dat$pa))
      dat$N<-as.double(as.character(dat$N))

      if(set=="mmds"){
        models<-c("mmds-MS")
      }else if(set=="cds"){
        models<-c("cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w")
      }else{
        models<-c("mmds-MS","cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w")
      }
      pa.res<-matrix(NA,200,length(models))
      aic.res<-matrix(NA,200,length(models))

      for(modi in seq_along(models)){
        if(sum(dat$mod==models[modi])>0){
          pa.res[,modi][dat$sim[dat$mod==models[modi]]]<-dat$pa[dat$mod==models[modi]]
          aic.res[,modi][dat$sim[dat$mod==models[modi]]]<-
                                              dat$AIC[dat$mod==models[modi]]
        }
      }

      aic.res[is.na(aic.res)]<-Inf

      # CDS has two outliers
      if(set=="cds"){
        rem.ind <- which(pa.res>2,arr.ind=TRUE)[,1]
        if(length(rem.ind)>0){
          pa.res<-pa.res[-rem.ind,]
          aic.res<-aic.res[-rem.ind,]
        }
      }

      aic.pick<-apply(aic.res,1,which.min)
      pa.cols<-ncol(pa.res)
      pa.aic<-cbind(pa.res,aic.pick)
      this.p<-apply(pa.aic,1,function(x){x[x[pa.cols+1]]})

      if(set=="mmds" | set=="combined"){
        aic.winners<-rbind(aic.winners,
                       data.frame(mix.terms=aic.pick,
                                  n=rep(n.samps,length(aic.pick)),
                                  model=rep("haz",length(aic.pick)),
                                  id=rep(par.ind,length(aic.pick))))
      }

      baf<-rbind(baf,data.frame(pa=this.p,
                                n=rep(n.samps,length(this.p)),
                                id=rep(par.ind,length(this.p)),
                                model=rep("haz",length(this.p))
                 ))

      true.ps<-c(true.ps,n.samps/dat$N[!is.na(dat$N)])
    }
  }
  # truth lines 
  true.p<-rbind(true.p,
#                data.frame(t=unique(round(true.ps,2)),
#                data.frame(t=rep(0.5,2),
                data.frame(t=c(0.5,0.4708531),
                           id=1:2,
                           model=rep("haz",2)))


  ##################################
  ##################################
  ##################################

  baf<-baf[-1,]
  baf$model<-as.factor(baf$model)

  ##################################

  # make the annotations of the proportion of results that had
  # the true number of mixture terms
  if(set=="mmds" | set=="combined"){  
    aic.winners<-aic.winners[-1,]

    # (time) in true model
    itm<-data.frame(n=0,prop=0,model=NA,id=0)

    # nocov
    for(i in 1:4){
       # calculate the proportions
       this.prop<-table(aic.winners)[2,,"nocov",i]/200

       itm<-rbind(itm,cbind(n=samp.sizes,
                            prop=round(this.prop,2),
                            model=rep("nocov",5),
                            id=rep(i,5)))
    }
    # pt
    for(i in 1:4){
       # calculate the proportions
       this.prop<-table(aic.winners)[2,,"pt",i]/200
       itm<-rbind(itm,cbind(n=samp.sizes,
                            prop=round(this.prop,2),
                            model=rep("pt",5),
                            id=rep(i,5)))
    }
    # covar
    for(i in 1:2){
       # calculate the proportions
       this.prop<-table(aic.winners)[2,,"covar",i]/200
       itm<-rbind(itm,cbind(n=samp.sizes,
                            prop=round(this.prop,2),
                            model=rep("covar",5),
                            id=rep(i,5)))
    }
    # 3pt
    for(i in 1:2){
       # calculate the proportions
       this.prop<-table(aic.winners)[1,,"3pt",i]/200
       itm<-rbind(itm,cbind(n=samp.sizes,
                            prop=round(this.prop,2),
                            model=rep("3pt",5),
                            id=rep(i,5)))
    }

    itm <- itm[-1,]
    itm$model <- factor(itm$model,levels=c("nocov","pt","3pt","covar"))
    levels(itm$model) <- model.names[1:4]

    ######################
    # time in mmds (>=2 pt mix)
    aic.winners<-aic.winners[-1,]

    immds<-data.frame(n=0,prop=0,model=NA,id=0)

    # nocov
    for(i in 1:4){
       # calculate the proportions
#       this.prop<-colSums(table(aic.winners)[1:2,,"nocov",i]/200)
         this.prop<-table(aic.winners)[2,,"nocov",i]/200

       immds<-rbind(immds,cbind(n=samp.sizes,
                            prop=round(this.prop,2),
                            model=rep("nocov",5),
                            id=rep(i,5)))
    }
    # pt
    for(i in 1:4){
       # calculate the proportions
#       this.prop<-colSums(table(aic.winners)[1:2,,"pt",i]/200)
       this.prop<-table(aic.winners)[2,,"pt",i]/200
       immds<-rbind(immds,cbind(n=samp.sizes,
                            prop=round(this.prop,2),
                            model=rep("pt",5),
                            id=rep(i,5)))
    }
    # covar
    for(i in 1:2){
       # calculate the proportions
       this.prop<-colSums(table(aic.winners)[1:2,,"covar-norecode",i]/200)
       immds<-rbind(immds,cbind(n=samp.sizes,
                            prop=round(this.prop,2),
                            model=rep("covar",5),
                            id=rep(i,5)))
    }
    # 3pt
    for(i in 1:2){
       # calculate the proportions
       this.prop<-table(aic.winners)[1,,"3pt-norecode",i]/200
       immds<-rbind(immds,cbind(n=samp.sizes,
                            prop=round(this.prop,2),
                            model=rep("3pt",5),
                            id=rep(i,5)))
    }
    # hazard
    for(i in 1:2){
       # calculate the proportions
#       this.prop<-colSums(table(aic.winners)[1,,"haz",i]/200)
       this.prop<-table(aic.winners)[1,,"haz",i]/200
       immds<-rbind(immds,cbind(n=samp.sizes,
                            prop=round(this.prop,2),
                            model=rep("haz",5),
                            id=rep(i,5)))
    }

    immds<-immds[-1,]
    immds$model<-factor(immds$model,levels=c("nocov","pt","3pt","covar","haz"))
    levels(immds$model) <- model.names
  }



  #######################################
  # reorder the models
  baf$model<-factor(baf$model,levels=c("nocov","pt","3pt","covar","haz"))
  true.p$model<-factor(true.p$model,levels=c("nocov","pt","3pt","covar","haz"))

  # rename the models
  levels(baf$model) <- model.names
  levels(true.p$model) <- model.names


  # make little circles with labels in them...
  if(set!="mmds"){
    circ.labs <- data.frame(model=rep(model.names,c(4,4,2,2,2)),
                            id=c(1:4,1:4,1:2,1:2,1:2),
                            circ.label=c(paste0("A",1:4),
                                         paste0("B",1:4),
                                         paste0("C",1:2),
                                         paste0("D",1:2),
                                         paste0("E",1:2)))
  }else{
    circ.labs <- data.frame(model=rep(model.names,c(4,4,2,4,2)),
                            id=c(1:4,1:4,1:2,1:4,1:2),
                            circ.label=c(paste0("A",1:4),
                                         paste0("B",1:4),
                                         paste0("C",1:2),
                                         paste0("D",c(1:2,1:2)),
                                         paste0("E",1:2)))
    nocov.labs <- data.frame(model=model.names[c(4,4)],
                            id=c(3,4),
                            circ.label=rep("No cov.",2))
  }

  #######################################
  # actually do the plotting here
  p<-ggplot(baf,aes(x=factor(n),y=pa))
  p<-p+geom_boxplot(outlier.size=1)
  p<-p+facet_grid(model~id)

  if(set=="mmds" | set=="combined"){
    p<-p+geom_text(aes(x=factor(n),y=-0.1,label=prop),size=3,data=itm)
  }
  if(set=="combined"){
    p<-p+geom_text(aes(x=factor(n),y=1.1,label=prop),size=3,data=immds)
  }
  if(set=="mmds"){
    p <- p + geom_text(aes(x=3.5,y=0.85,label=circ.label),size=3,data=nocov.labs)
  }


  # plot the label and circle
  p <- p + geom_point(aes(x=5,y=0.87),colour="black",
                      size=7,data=circ.labs,shape=1)
  p <- p + geom_text(aes(x=5,y=0.87,label=circ.label),size=3,data=circ.labs)

  p<-p+labs(x="Sample size",
            y=substitute(paste(lab, hat(P[a])),
                         list(lab="Probability of detection, ")))
  p<-p+theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            legend.background=element_blank(),
            legend.key=element_blank(),
            strip.text.x = element_blank(),
            strip.background = element_blank(),
            panel.background=element_blank())

  # transparency not supported by eps!
  p<-p+geom_hline(aes(yintercept=t),true.p,col="grey")

  quartz()
  print(p)

  dev.copy2pdf(file=paste("pa-plot-",set,".pdf",sep=""))

}
