# what won?
library(ggplot2)

sim.combs<-data.frame(types=c("nocov","pt","3point","covar"),
                      number=c(4,4,2,2),
                      start=c("BFGS+SANN-","BFGS+SANN-","BFGS+SANN-","covsim"),
                      end=c("-960-results.csv","-pt-results.csv",
                            "-960-3pt-results.csv","-BFGS+SANN.csv"),
                      n.y=c(-2.1,-3.6,-1.75,-0.5))

headers<-list(
c("model","par.ind","n.samp","sim","par1","par2","par3","pa","aic","Nhat","N"),
c("model","par.ind","n.samp","sim","par1","par2","par3","pall","aic","N","Nhat"),
c("model","par.ind","n.samp","sim","ll","aic","pa","Nhat","N","mixterms"),
c("n.samp","sim","aic","pa","Nhat","N","mixpart","model")
)

todo<-1:4 

for(combi in todo){

  big.pop<-c()
  win.labels<-c()

  simn<-sim.combs$number[combi]
  type<-sim.combs$types[combi]
  start<-sim.combs$start[combi]
  end<-sim.combs$end[combi]
  n.y<-sim.combs$n.y[combi]
  
  for(pari in 1:simn){
    # nocov
    dat<-read.csv(paste(type,"/",start,pari,end,sep=""))
    dat<-dat[,-1]

    names(dat)<-headers[[combi]]
    
    # some things got messed up, not sure why...
    if(combi==4){
      dat$model[which(dat$pa==Inf)]<-dat$mixpart[which(dat$pa==Inf)]
    }

    dat.aic<-c()
    dat.N<-c()
    dat.Nhat<-c()
    dat.n.samp<-c()

    # put some data together
    for(model.names in unique(dat$model)){
      dat.aic<-cbind(dat.aic,
                   as.numeric(as.character(dat$aic[dat$model==model.names])))
      dat.N<-cbind(dat.N,
                 as.numeric(as.character(dat$N[dat$model==model.names])))
      dat.Nhat<-cbind(dat.Nhat,
                    as.numeric(as.character(dat$Nhat[dat$model==model.names])))
      dat.n.samp<-cbind(dat.n.samp,
                  as.numeric(as.character(dat$n.samp[dat$model==model.names])))
    }
    # find the AIC-best models
    winners <- apply(dat.aic,1,which.min)
    these.winners <- unique(dat$model)[winners]

    # now re-code the models so that we are just comparing mmds and (m)cds
    if(combi==1 | combi==2 ){
      levels(these.winners)<-c("CDS","CDS","MMDS","MMDS")
    }else if(combi==3){
      levels(these.winners)<-c("CDS","CDS","MMDS")
    }else if(combi==4){
      levels(these.winners)<-c(rep("MCDS",4),"CMMDS","CDS","MCDS","MCDS",
                                "CDS","MCDS","MCDS","MMDS")
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
                    p.prb=(1/this.N-1/this.Nhat)*this.N
                   )
  
    big.pop <- rbind(big.pop,pop)

    this.win <- c()

    if(combi==1 | combi==2 | combi==3){
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
                      size=3,data=win.labels,col="grey")

#  quartz()
  print(p)
  ggsave(paste("boxplots-",type,".pdf",sep=""))

}








