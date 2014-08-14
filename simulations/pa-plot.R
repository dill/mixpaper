# Pa boxplots
# makes a table of boxplots for all of the simulations with the 
# proportion of time that the true model was selected by AIC printed
# below each box.

# also of interest is table(aic.winners)
# which give the breakdown of the AIC best models according to model
# type.

library(ggplot2)
library(plyr)
library(reshape2)
options(stringsAsFactors=FALSE)


model.names <- c("A. No covariates","B. Point transect",
                 "C. 3-point","D. Covariate","E. Other")

samp.sizes<-c(30,60,120,480,960)

# function to extract the boxplot data
get.bp.data <- function(dat,type){
  this.baf <- ddply(dat,.(par.ind,n,sim),
                    function(x){
                       res<- c(x$pa[which.min(x$AIC)], min(x$AIC),
                               x$model[which.min(x$AIC)])
                       names(res) <- c("pa","aic","model")
                       return(res)
                    })
  # set the type and make sure that the pa and aic are numeric
  this.baf$type <- type
  this.baf$pa <- as.numeric(this.baf$pa)
  this.baf$aic <- as.numeric(this.baf$aic)
  return(this.baf)
}

# function to get the aic winners data
get.aicw <- function(dat,true.model,type){
  this.aicw <- ddply(dat,.(par.ind,n,sim),
                     function(x){
                        # pick out the best model
                        bm <- x$model[which.min(x$AIC)]
                        names(bm) <- "model"
                        return(bm)
                     })
  this.aicw <- ddply(this.aicw,.(par.ind,n),
                     function(x){
                       prop <- round(sum(x$model %in% true.model)/nrow(x),2)
                       names(prop)<-"prop"
                       return(prop)
                     })
  this.aicw$type <- type
  return(this.aicw)
}

for(set in c("mmds","cds","combined")){

  # set up containers for results
  itm <- c() # time in true model
  immds <- c() # time in mmds models
  true.p <- c() # true p_a
  true.ps <- c()
  baf <- c() # boxplot data


  ################################################################
  ### Models with no covariates
  for(par.ind in 1:4){

    # read and correctly format the data
    dat<-read.csv(file=paste("nocov/BFGS+SANN-",par.ind,
                              "-960-results.csv",sep=""))
    dat<-dat[,-1]
    names(dat)<-c("model","par.ind","n","sim","par1","par2",
                  "par3","pa","AIC","Nhat","N")
    dat$Nhat<-as.double(dat$Nhat)
    dat$N<-as.double(dat$N)

    models <- switch(set,
                     mmds = c("mmds-1","mmds-c"),
                     cds  = c("cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w"),
                     combined = c("mmds-1","mmds-c","cds-hnc","cds-hrp",
                                  "cds-hnc-w","cds-hrp-w")
                    )

    dat <- dat[dat$model%in%models,]

    ## build the frame which will become the boxplots
    this.baf <- get.bp.data(dat,"nocov")
    baf <- rbind(baf,this.baf)

    ## build the frame of AIC winners
    if(set=="mmds" | set == "combined"){
      # get time in true model
      this.itm <- get.aicw(dat,"mmds-c","nocov")
      itm <- rbind(itm,this.itm)
      # get time in mmds
      this.immds <- get.aicw(dat,"mmds-c","nocov")
      immds <- rbind(immds,this.immds)
    }

    # what is the true p?
    true.ps <- c(true.ps, median(dat$n/dat$N,na.rm=TRUE))
  }

  # truth lines -- nocov
  true.p<-data.frame(t=true.ps,
                     par.ind=1:4,
                     type=rep("nocov",4))

  true.ps<-c()

  ################################################################
  ### point transects
  for(par.ind in 1:4){

    # read and correctly format the data
    dat<-read.csv(file=paste("pt/BFGS+SANN-",par.ind,
                             "-pt-results.csv",sep=""))
    dat<-dat[,-1]
    names(dat)<-c("model","par.ind","n","sim","par1","par2",
                  "par3","ll","AIC","N","Nhat")
    dat$Nhat<-as.numeric(dat$Nhat)
    dat$N<-as.numeric(dat$N)
    dat$pa <- dat$n/dat$Nhat

    models <- switch(set,
                     mmds = c("mmds-1","mmds-2"),
                     cds  = c("cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w"),
                     combined = c("mmds-1","mmds-2","cds-hnc","cds-hrp",
                                  "cds-hnc-w","cds-hrp-w")
                    )
    dat <- dat[dat$model%in%models,]

    this.baf <- get.bp.data(dat,"pt")
    baf <- rbind(baf,this.baf)

    ## build the frame of AIC winners
    if(set=="mmds" | set == "combined"){
      # get time in true model
      itm <- rbind(itm,get.aicw(dat,"mmds-2","pt"))
      # get time in mmds
      immds <- rbind(immds,get.aicw(dat,"mmds-2","pt"))
    }

    # what is the true p?
    true.ps <- c(true.ps, median(dat$n/dat$N,na.rm=TRUE))
  }

  # truth lines -- pt
  true.p<-rbind(true.p,
                data.frame(t=true.ps,
                           par.ind=1:4,
                           type=rep("pt",4)))

  ##################################################################
  ### covar
  true.ps<-c()

  for(par.ind in 1:2){
    true.pp<-c()

    ## read in the data
    dat<-read.csv(file=paste("covar/covsim",par.ind,"-BFGS+SANN.csv",sep=""))
    dat<-dat[,-1]
    names(dat)<-c("n","sim","AIC","pa","Nhat","N","mix.terms","model")

    # select only data with sample size n.samps
    dat$Nhat <- as.numeric(dat$Nhat)
    dat$N <- as.numeric(dat$N)
    dat$par.ind <- par.ind

    # which models are we interested in?
    models <- switch(set,
                     mmds = c("nocov","cov"),
                     cds  = c("hn+cos","hr+poly","hn+cov1","hr+cov1"),
                     combined = c("nocov","cov","hn+cos","hr+poly","hn+cov1",
                                  "hr+cov1")
                    )
    dat <- dat[dat$model%in%models,]

    if(set=="combined"){
      # drop mixture models with 1 term -- just a hn+cos model
      dat <- dat[!(dat$model=="cov" & dat$mix.terms==1),]
      dat <- dat[!(dat$model=="nocov" & dat$mix.terms==1),]
    }

    ## build the frame which will become the boxplots
    baf <- rbind(baf,get.bp.data(dat,"covar"))

    ## build the frame of AIC winners
    if(set=="mmds" | set == "combined"){
      # get time in true model
      itm <- rbind(itm,get.aicw(dat,"cov","covar"))
      # get time in mmds
      immds <- rbind(immds,get.aicw(dat,c("cov","nocov"),"covar"))
    }

    ## separate plot for when covariates are not included
    if(set == "mmds"){
      this.baf <- get.bp.data(dat[dat$model%in%"nocov",],"covar")
    }else if(set=="cds"){
      this.baf <- get.bp.data(dat[dat$model%in%c("hn+cos","hr+poly"),],"covar")
    }else{
      this.baf <- get.bp.data(dat[dat$model%in%c("nocov","hn+cos","hr+poly"),],
                              "covar")
    }

    this.baf$par.ind<-this.baf$par.ind+2
    baf <- rbind(baf,this.baf)

    true.ps<-c(true.ps,median(dat$n/dat$N,na.rm=TRUE))
  }
  # truth lines -- covar
  true.p<-rbind(true.p,data.frame(t=rep(true.ps,c(2,2)),
                                  par.ind=1:4,
                                  type=rep("covar",4)))

  ##################################################################
  ### 3 point
  true.ps<-c()

  for(par.ind in 1:2){
    dat<-read.csv(file=paste("3point/BFGS+SANN-",par.ind,
                             "-960-3pt-results.csv",sep=""))
    dat<-dat[,-1]
    names(dat)<-c("model","par.ind","n","sim","ll","AIC",
                  "pa","Nhat","N","mix.terms")
    dat$Nhat<-as.numeric(dat$Nhat)
    dat$N<-as.numeric(dat$N)

    if(set=="mmds" | set == "combined"){
      dat$model <- paste0(dat$model,ifelse(dat$mix.terms=="mt","",
                                           paste0("-",dat$mix.terms)))
    }

    models <- switch(set,
                     mmds = c("mmds-MS-1","mmds-MS-2","mmds-MS-3"),
                     cds  = c("cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w"),
                     combined = c("mmds-MS-1","mmds-MS-2","mmds-MS-3",
                                  "cds-hnc","cds-hrp",
                                  "cds-hnc-w","cds-hrp-w")
                    )
    dat <- dat[dat$model%in%models,]
    baf <- rbind(baf,get.bp.data(dat,"3pt"))

    ## build the frame of AIC winners
    if(set=="mmds" | set == "combined"){
      itm <- rbind(itm,get.aicw(dat,"mmds-MS-3","3pt"))
      immds <- rbind(immds,get.aicw(dat,c("mmds-MS-1","mmds-MS-2","mmds-MS-3"),
                                    "3pt"))
    }

    true.p<-rbind(true.p,data.frame(t=median(dat$n/dat$N,na.rm=TRUE),
                                    par.ind=par.ind,
                                    type="3pt"))
  }

  ##################################
  # hazard rate mixtures and eps (counted as hazard for simplicity of this code)
  for(par.ind in 2:1){
    dat<-read.csv(file=paste("hazard/hr-",par.ind,
                             "-results.csv",sep=""),header=FALSE)

    # eps is listed as par.ind=1 for file but should be plotted second
    par.ind <- ifelse(par.ind==1,2,1)
    dat<-dat[,-1]
    names(dat)<-c("model","par.ind","n","sim","ll","AIC",
                  "pa","Nhat","N","mixterms")
    dat$Nhat<-as.double(dat$Nhat)
    dat$pa<-as.double(dat$pa)
    dat$N<-as.double(dat$N)
    dat$par.ind <- par.ind

    # since we used model selection recode
    if(set=="combined"){
      dat$model[dat$model=="mmds-MS" & dat$mixterms==1]<-"cds-hnc-w"
    }

    models <- switch(set,
                     mmds = c("mmds-MS"),
                     cds  = c("cds-hnc","cds-hrp","cds-hnc-w","cds-hrp-w"),
                     combined = c("mmds-MS","cds-hnc","cds-hrp",
                                  "cds-hnc-w","cds-hrp-w")
                    )
    dat <- dat[dat$model%in%models,]

    baf <- rbind(baf,get.bp.data(dat,"other"))

    ## build the frame of AIC winners
    if(set=="mmds" | set == "combined"){
      # get time in mmds
      immds <- rbind(immds,get.aicw(dat,"mmds-MS","other"))
    }

    true.p<-rbind(true.p,data.frame(t=median(dat$n/dat$N,na.rm=TRUE),
                                    par.ind=par.ind,
                                    type="other"))
  }


  #######################################
  ## reorder the models
  fixfactors<-function(xx){
    xx <- factor(as.factor(xx),levels=c("nocov","pt","3pt","covar","other"))
    levels(xx) <- model.names
    xx
  }

  baf$type <- fixfactors(baf$type)
  true.p$type <- fixfactors(true.p$type)

  if(set=="mmds" | set == "combined"){
    itm$type <- fixfactors(itm$type)
    immds$type <- fixfactors(immds$type)
  }

  # make little circles with labels in them...
  circ.labs <- data.frame(type=rep(model.names,c(4,4,2,4,2)),
                          par.ind=c(1:4,1:4,1:2,1:4,1:2),
                          circ.label=c(paste0("A",1:4),
                                       paste0("B",1:4),
                                       paste0("C",1:2),
                                       paste0("D",c(1:2,1:2)),
                                       paste0("E",1:2)))
  nocov.labs <- data.frame(type=factor(rep(model.names[4],2),
                                        levels(baf$type)),
                          par.ind=c(3,4),
                          circ.label=rep("No cov.",2))

  circ.labs$type <- factor(circ.labs$type,model.names)

  #######################################
  # actually do the plotting here
  p<-ggplot(baf,aes(x=factor(n),y=pa))
  p<-p+geom_boxplot(outlier.size=1)
  p<-p+facet_grid(type~par.ind)

  # add the annotations for time in true model (itm) and
  # time in >1pt mixture (immds)
  if(set=="mmds" | set=="combined"){
    p <- p+geom_text(aes(x=factor(n),y=-0.1,label=prop),size=3,data=itm)
  }
  if(set=="combined"){
    p <- p+geom_text(aes(x=factor(n),y=1.1,label=prop),size=3,data=immds)
  }

  # label the covariate sims where the covariates were not included
  p <- p + geom_text(aes(x=3.5,y=0.85,label=circ.label),size=3,data=nocov.labs)


  # plot the label and circle
  p <- p + geom_point(aes(x=5,y=0.88),colour="black",
                      size=5,data=circ.labs,shape=1)
  p <- p + geom_text(aes(x=5,y=0.88,label=circ.label),size=2,data=circ.labs)

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

  # add a scale, limits and nice breaks
  p <- p + scale_y_continuous(breaks=c(0,0.5,1),limits=c(-0.1,1.1))

  # transparency not supported by eps!
  p<-p+geom_hline(aes(yintercept=t),true.p,col="grey")

  quartz()
  print(p)

  dev.copy2pdf(file=paste("pa-plot-",set,".pdf",sep=""))

}
