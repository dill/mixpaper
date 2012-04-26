# what won?
library(ggplot2)



sim.combs<-data.frame(types=c("nocov","pt","3point","covar"),
                      number=c(4,4,2,2),
                      start=c("BFGS+SANN-","BFGS+SANN-","BFGS+SANN-","covsim"),
                      end=c("-960-results.csv","-pt-results.csv",
                            "-960-3pt-results.csv","-BFGS+SANN.csv"))

headers<-list(
c("model","par.ind","n.samp","sim","par1","par2","par3","pa","aic","Nhat","N"),
c("model","par.ind","n.samp","sim","par1","par2","par3","pall","aic","N","Nhat"),
c("model","par.ind","n.samp","sim","ll","aic","pa","Nhat","N","mixterms"),
c("n.samp","sim","aic","pa","Nhat","N","mixpart","model")
)


todo<-c(1,2,3,4)
#todo<-1:4 

for(combi in todo){

  big.pop<-c()

  simn<-sim.combs$number[combi]
  type<-sim.combs$types[combi]
  start<-sim.combs$start[combi]
  end<-sim.combs$end[combi]
  
  for(pari in 1:simn){
    # nocov
    dat<-read.csv(paste(type,"//",start,pari,end,sep=""))
    dat<-dat[,-1]

    names(dat)<-headers[[combi]]
    
    dat$Nhat<-as.double(as.character(dat$Nhat))
    
    
    # some things got messed up, not sure why...
    if(combi==4){
      dat$model[which(dat$pa==Inf)]<-dat$mixpart[which(dat$pa==Inf)]
    }

    dat.aic<-c()
    dat.N<-c()
    dat.Nhat<-c()
    dat.n.samp<-c()

    for(model.names in unique(dat$model)){
      dat.aic<-cbind(dat.aic,as.numeric(as.character(dat$aic[dat$model==model.names])))
      dat.N<-cbind(dat.N,as.numeric(as.character(dat$N[dat$model==model.names])))
      dat.Nhat<-cbind(dat.Nhat,as.numeric(as.character(dat$Nhat[dat$model==model.names])))
      dat.n.samp<-cbind(dat.n.samp,as.numeric(as.character(dat$n.samp[dat$model==model.names])))
    }

    winners<-apply(dat.aic,1,which.min)

    # cheatcode to pull the winners out of the other data
    # bind the winners onto the end, then use the element 
    #  in that col to lookup
    this.N<-apply(cbind(dat.N,winners),1,function(x){x[x[length(x)]]})
    this.Nhat<-apply(cbind(dat.Nhat,winners),1,function(x){x[x[length(x)]]})
    this.n.samp<-apply(cbind(dat.n.samp,winners),1,function(x){x[x[length(x)]]})
  
    pop<-data.frame(winner=unique(dat$model)[winners],
                    pari=rep(pari,length(winners)),
                    type=rep(type,length(winners)),
                    n.samp=this.n.samp,
                    N=this.N,
                    Nhat=this.Nhat,
                    prb=(this.N-this.Nhat)/this.N,
                    p.prb=(1/this.N-1/this.Nhat)*this.N
                   )
  
    big.pop <- rbind(big.pop,pop)
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
  print(p)
#  ggsave(paste("boxplots-",type,".pdf",sep=""))
  quartz()

}








