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

for(combi in todo){

  big.pop<-c()
  win.labels<-c()

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

        dat.aic.t[,modi][ind]<-as.numeric(as.character(datn$aic[ind2]))
        dat.N.t[,modi][ind]<-as.numeric(as.character(datn$N[ind2]))
        dat.Nhat.t[,modi][ind]<-as.numeric(as.character(datn$Nhat[ind2]))
        dat.n.samp.t[,modi][ind]<-as.numeric(as.character(datn$n.samp[ind2]))
      }

      if(combi==3){
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
      dat.mixterms<-rbind(dat.mixterms,dat.mixterms.t)
    }


    # find the AIC-best models
    dat.aic[is.na(dat.aic)]<-Inf
    winners <- apply(dat.aic,1,which.min)
    these.winners <- models[winners]

    # now re-code the models so that we are just comparing mmds and (m)cds
    if(combi==1 | combi==2){
      these.winners<-as.factor(these.winners)
      levels(these.winners)[levels(these.winners)=="mmds-c"]<-"MMDS"
      levels(these.winners)[levels(these.winners)=="mmds-2"]<-"MMDS"
      levels(these.winners)[levels(these.winners)=="mmds-1"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hrp"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hrp-w"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hnc"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hnc-w"]<-"CDS"
    }else if(combi==3){
      these.winners<-as.factor(these.winners)
      these.winners[these.winners=="mmds-MS" & 
                            (dat.mixterms==1)]<-"cds-hnc"
      levels(these.winners)[levels(these.winners)=="mmds-MS"]<-"MMDS"
      levels(these.winners)[levels(these.winners)=="mmds-1"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hrp"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hrp-w"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hnc"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hnc-w"]<-"CDS"
    }else if(combi==4){
      these.winners<-as.factor(these.winners)
      these.winners[these.winners=="cov" & dat.mixterms!=2]<-"hn+cov1"
      these.winners[these.winners=="nocov" & dat.mixterms!=2]<-"hr+poly"
      levels(these.winners)[levels(these.winners)=="cov"]<-"CMMDS"
      levels(these.winners)[levels(these.winners)=="nocov"]<-"MMDS"
      levels(these.winners)[levels(these.winners)=="hr+cov1"]<-"MCDS"
      levels(these.winners)[levels(these.winners)=="hn+cov1"]<-"MCDS"
      levels(these.winners)[levels(these.winners)=="hn+cos"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="hr+poly"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="hn+cos+cov1-width"]<-"MCDS"
      levels(these.winners)[levels(these.winners)=="hr+poly+cov1-width"]<-"MCDS"
    }else if(combi==5){
      these.winners<-as.factor(these.winners)
      levels(these.winners)[levels(these.winners)=="mmds-MS"]<-"MMDS"
      levels(these.winners)[levels(these.winners)=="cds-hrp"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hrp-w"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hnc"]<-"CDS"
      levels(these.winners)[levels(these.winners)=="cds-hnc-w"]<-"CDS"
    }

    # cheatcode to pull the winners out of the other data
    # bind the winners onto the end, then use the element 
    #  in that col to lookup
    this.N<-apply(cbind(dat.N,winners),1,function(x){x[x[length(x)]]})
    this.Nhat<-apply(cbind(dat.Nhat,winners),1,function(x){x[x[length(x)]]})
    this.n.samp<-apply(cbind(dat.n.samp,winners),1,function(x){x[x[length(x)]]})
  
    pop<-data.frame(winner=these.winners,
                    pari=rep(pari,length(winners)),
                    type=rep(type,length(winners)),
                    n.samp=this.n.samp,
                    N=this.N,
                    Nhat=this.Nhat,
                    prb=(this.N-this.Nhat)/this.N,
                    p.prb=(1/this.Nhat-1/this.N)*this.N
                   )
  
    big.pop <- rbind(big.pop,pop)

    this.win <- c()

    if(combi==1 | combi==2 | combi==3 | combi == 5){
      ermcds<-"CDS"
      for(win.n.samp in unique(pop$n.samp)){

        this.win<-rbind(this.win,
                        data.frame(
                                winner = c(ermcds,"MMDS"),
                                pari = rep(pari,2),
                                type = rep(type,2), 
                                n.samp=rep(win.n.samp,2),
                                n.winner=c(sum(as.character(pop$winner[pop$n.samp==win.n.samp])==ermcds),
                                           sum(as.character(pop$winner[pop$n.samp==win.n.samp])=="MMDS"))
                                  )
                       )
      }
    }else{
      for(win.n.samp in unique(pop$n.samp)){

        this.win<-rbind(this.win,
                        data.frame(
                                winner = c("CDS","CMMDS","MCDS","MMDS"),
                                pari = rep(pari,4),
                                type = rep(type,4),
                                n.samp=rep(win.n.samp,4),
                                n.winner=c(
              sum(as.character(pop$winner[pop$n.samp==win.n.samp])=="CDS"),
              sum(as.character(pop$winner[pop$n.samp==win.n.samp])=="CMMDS"),
              sum(as.character(pop$winner[pop$n.samp==win.n.samp])=="MCDS"),
              sum(as.character(pop$winner[pop$n.samp==win.n.samp])=="MMDS")
                                  ))
                       )
      }
    }

    win.labels<-rbind(win.labels,this.win)
  }

  # actually do the plotting
  p <- ggplot(big.pop)
  p <- p + geom_boxplot(aes(y=p.prb,x=as.factor(winner)))
  p <- p + facet_grid(pari~n.samp)
  p<-p+opts(panel.grid.major=theme_blank(),
            panel.grid.minor=theme_blank(),
            legend.background=theme_blank(),
            legend.key=theme_blank(),
            axis.text.x=theme_text(angle=-90),
            panel.background=theme_rect())
  p<-p+labs(x="Sample size",y="rel. bias in p")
  p<-p+geom_hline(yintercept=0,col="grey")
  p<-p+geom_text(aes(x=factor(winner),y=n.y,label=n.winner),
                      size=3,data=win.labels,col="black")

  quartz()
  print(p)
  ggsave(paste("boxplots-",type,".pdf",sep=""))

}

