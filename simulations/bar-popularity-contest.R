# what won?
library(ggplot2)

sim.combs<-data.frame(types=c("nocov","pt","3point","covar","hazard"),
                      number=c(4,4,2,2,2),
                      start=c("BFGS+SANN-","BFGS+SANN-","BFGS+SANN-",
                              "covsim","hr-"),
                      end=c("-960-results.csv","-pt-results.csv",
                            "-960-3pt-results.csv","-BFGS+SANN.csv",
                            "-results.csv"),
                      n.y=c(-2.0,-3.8,-2.2,-1.1,-1.1))

headers<-list(
c("model","par.ind","n.samp","sim","par1","par2","par3","pa","aic","Nhat","N"),
c("model","par.ind","n.samp","sim","par1","par2","par3","pall","aic","N","Nhat"),
c("model","par.ind","n.samp","sim","ll","aic","pa","Nhat","N","mixterms"),
c("n.samp","sim","aic","pa","Nhat","N","mixterms","model"),
c("model","par.ind","n.samp","sim","ll","aic","pa","Nhat","N","mixterms")
)

todo<-1:5

big.pop<-c()

for(combi in todo){

  simn<-sim.combs$number[combi]
  type<-sim.combs$types[combi]
  start<-sim.combs$start[combi]
  end<-sim.combs$end[combi]
  n.y<-sim.combs$n.y[combi]

  for(pari in 1:simn){

    dat<-read.csv(paste(type,"/",start,pari,end,sep=""))
    dat<-dat[,-1]

    names(dat)<-headers[[combi]]

    if(combi==4){
      # add in the extra models with no adjustments
      dat2<-read.csv(paste(type,"/",start,pari,"-noadj",end,sep=""))
      dat2<-dat2[,-1]
      names(dat2)<-headers[[combi]]
      dat<-rbind(dat,dat2)
      rm(dat2)

      # some things got messed up, not sure why...
      dat$model[which(dat$pa==Inf)]<-dat$mixterms[which(dat$pa==Inf)]

      # ignore the adjustment models
      models<-c("nocov","cov","hn+cos","hr+poly","hn+cov1",
                 "hr+cov1")#,"hn+cos+cov1-width","hr+poly+cov1-width")

    }else{
      models<-as.character(unique(dat$model))
    }

    dat.aic<-c()
    dat.N<-c()
    dat.Nhat<-c()
    dat.n.samp<-c()
    dat.mixterms<-c()

    n.samples<-c(30,60,120,480,960)

    for(n in n.samples){
      dat.aic.t<-matrix(NA,200,length(models))
      dat.N.t<-matrix(NA,200,length(models))
      dat.Nhat.t<-matrix(NA,200,length(models))
      dat.n.samp.t<-matrix(NA,200,length(models))
      dat.mixterms.t<-rep(NA,200)

      datn<-dat[dat$n.samp==n,]

      # put some data together
      for(modi in seq_along(models)){
        ind <- datn$sim[datn$model==models[modi]]
        ind2 <- datn$model==models[modi] & datn$n.samp==n

        # this throws some NA coercion warnings, those are for the NAs
        # don't panic
        dat.aic.t[,modi][ind]<-as.numeric(as.character(datn$aic[ind2]))
        dat.N.t[,modi][ind]<-as.numeric(as.character(datn$N[ind2]))
        dat.Nhat.t[,modi][ind]<-as.numeric(as.character(datn$Nhat[ind2]))
        dat.n.samp.t[,modi][ind]<-as.numeric(as.character(datn$n.samp[ind2]))
      }

      if(combi==3 | combi==5){
        dat.mixterms.t[datn$sim[datn$model=="mmds-MS"]]<-
               as.numeric(as.character(datn$mixterms[datn$model=="mmds-MS"]))
      }
      if(combi==4){
        dat.mixterms.t[datn$sim[datn$model=="cov"]]<-
               as.numeric(as.character(datn$mixterms[datn$model=="cov"]))
      }

      dat.aic<-rbind(dat.aic,dat.aic.t)
      dat.N<-rbind(dat.N,dat.N.t)
      dat.Nhat<-rbind(dat.Nhat,dat.Nhat.t)
      dat.n.samp<-rbind(dat.n.samp,dat.n.samp.t)
      #dat.mixterms<-rbind(dat.mixterms,dat.mixterms.t)
      dat.mixterms<-c(dat.mixterms,dat.mixterms.t)
    }


    # find the AIC-best models
    dat.aic[is.na(dat.aic)]<-Inf
    winners <- apply(dat.aic,1,which.min)
    these.winners <- models[winners]

    if(combi==1 | combi==2){
      these.winners[these.winners=="mmds-c"]<-"MMDS 2-pt"
      these.winners[these.winners=="mmds-2"]<-"MMDS 2-pt"
      these.winners[these.winners=="mmds-1"]<-"hn"
      these.winners[these.winners=="cds-hrp"]<-"K+A hr+poly"
      these.winners[these.winners=="cds-hrp-w"]<-"K+A hr+poly (w)"
      these.winners[these.winners=="cds-hnc"]<-"K+A hn+cos"
      these.winners[these.winners=="cds-hnc-w"]<-"K+A hn+cos (w)"
    }else if(combi==3){
      these.winners[these.winners=="mmds-MS" & (dat.mixterms==3)]<-"MMDS 3-pt"
      these.winners[these.winners=="mmds-MS" & (dat.mixterms==1)]<-"cds-hnc"
      these.winners[these.winners=="mmds-MS"]<-"MMDS 2-pt"
      these.winners[these.winners=="mmds-1"]<-"hn"
      these.winners[these.winners=="cds-hrp"]<-"K+A hr+poly"
      these.winners[these.winners=="cds-hrp-w"]<-"K+A hr+poly (w)"
      these.winners[these.winners=="cds-hnc"]<-"K+A hn+cos"
      these.winners[these.winners=="cds-hnc-w"]<-"K+A hn+cos (w)"
    }else if(combi==4){
      these.winners[these.winners=="cov" & dat.mixterms!=2]<-"hn+cov1"
      these.winners[these.winners=="nocov" & dat.mixterms!=2]<-"hr+poly"
      these.winners[these.winners=="cov" & (dat.mixterms==3)]<-"MMDS 3-pt (cov)"
      these.winners[these.winners=="cov" & (dat.mixterms==2)]<-"MMDS 2-pt (cov)"
      these.winners[these.winners=="cov" & (dat.mixterms==1)]<-"hn (cov)"
      these.winners[these.winners=="nocov" & (dat.mixterms==3)]<-"MMDS 3-pt"
      these.winners[these.winners=="nocov" & (dat.mixterms==2)]<-"MMDS 2-pt"
      these.winners[these.winners=="nocov" & (dat.mixterms==1)]<-"hn"
      these.winners[these.winners=="hr+cov1"]<-"K+A hr+poly (cov)"
      these.winners[these.winners=="hn+cov1"]<-"K+A hn+cos (cov)"
      these.winners[these.winners=="hn+cos"]<-"K+A hn+cos"
      these.winners[these.winners=="hr+poly"]<-"K+A hr+poly"
      these.winners[these.winners=="hn+cos+cov1-width"] <- "K+A (cov) hn+cos (w)"
      these.winners[these.winners=="hr+poly+cov1-width"]<- "K+A (cov) hr+poly (w)"
    }else if(combi==5){
      these.winners[these.winners=="mmds-MS" & (dat.mixterms==3)]<-"MMDS 3-pt"
      these.winners[these.winners=="mmds-MS" & (dat.mixterms==2)]<-"MMDS 2-pt"
      these.winners[these.winners=="mmds-MS" & (dat.mixterms==1)]<-"hn"
      these.winners[these.winners=="cds-hrp"]<-"K+A hr+poly"
      these.winners[these.winners=="cds-hrp-w"]<-"K+A hr+poly (w)"
      these.winners[these.winners=="cds-hnc"]<-"K+A hn+cos"
      these.winners[these.winners=="cds-hnc-w"]<-"K+A hn+cos (w)"
    }


    # cheatcode to pull the winners out of the other data
    # bind the winners onto the end, then use the element 
    #  in that col to lookup
    this.N<-apply(cbind(dat.N,winners),1,function(x){x[x[length(x)]]})
    this.Nhat<-apply(cbind(dat.Nhat,winners),1,function(x){x[x[length(x)]]})
    this.n.samp<-apply(cbind(dat.n.samp,winners),1,function(x){x[x[length(x)]]})

    # ignore whether width was the scaling or not
    these.winners[these.winners=="K+A hr+poly (w)"] <- "K+A hr+poly"
    these.winners[these.winners=="K+A hn+cos (w)"] <- "K+A hn+cos"

    # make a factor
    these.winners<-as.factor(these.winners)

    pop<-data.frame(winner=these.winners,
                    pari=rep(pari,length(winners)),
                    type=rep(type,length(winners)),
                    n.samp=this.n.samp,
                    N=this.N
                   )

    big.pop <- rbind(big.pop,pop)
  }
}

big.pop$winner <- as.factor(big.pop$winner)
big.pop$n.samp <- as.factor(big.pop$n.samp)


# give things nice names... 
model.names <- c("A. No covariates","B. Point transect",
                 "C. 3-point","D. Covariate","E. Hazard mixture")
big.pop$type <- factor(big.pop$type,
                       levels=c("nocov","pt","3point","covar","hazard"))
levels(big.pop$type) <- model.names




# actually do the plotting
p <- ggplot(big.pop)
p <- p + geom_bar(aes(n.samp,fill=winner))
p <- p + facet_grid(type~pari)
p <- p + theme(panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               legend.background=element_blank(),
               legend.key=element_blank(),
               axis.text.x=element_text(angle=-90),
               strip.text.x = element_blank(),
               strip.background = element_blank(),
               panel.background=element_blank())
p<-p+labs(x="Sample size",
          fill="AIC best model",
          y="Number of simulations")
print(p)
ggsave("bar-winners.pdf")

