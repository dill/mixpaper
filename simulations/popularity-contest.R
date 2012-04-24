# what won?

big.pop<-c()


sim.combs<-data.frame(types=c("nocov","pt","3point","covar"),
                      number=c(4,4,2,2),
                      start=c("BFGS+SANN-","BFGS+SANN-","BFGS+SANN-","covsim"),
                      end=c("-960-results.csv","-pt-results.csv",
                            "-960-3pt-results.csv","-BFGS+SANN.csv"))

headers<-list(
c("model","par.ind","n.samp","sim","par1","par2","par3","pa","aic","Nhat","N"),
c("model","par.ind","n.samp","sim","par1","par2","par3","pall","aic","N","Nhat"),
c("model","par.ind","n.samp","sim","par1","par2","par3","pa","aic","N","Nhat"),
c("n.samp","sim","aic","pa","Nhat","N","??","model")
)


todo<-c(1,2,4)
#todo<-1:4 

for(combi in todo){

  simn<-1#sim.combs$number[combi]
  type<-sim.combs$types[combi]
  start<-sim.combs$start[combi]
  end<-sim.combs$end[combi]
  
  for(pari in 1:simn){
    # nocov
    dat<-read.csv(paste(type,"//",start,pari,end,sep=""))
    dat<-dat[,-1]

    names(dat)<-headers[[combi]]
    
    dat$Nhat<-as.double(as.character(dat$Nhat))
    
    dat.sel<-c()
    
    # some things got messed up, not sure why...
    if(combi==4){
      dat$model[which(dat$pa==Inf)]<-dat$'??'[which(dat$pa==Inf)]
    }

    for(model.names in unique(dat$model)){
      dat.sel<-cbind(dat.sel,dat$aic[dat$model==model.names])
    }
  
    winners<-apply(dat.sel,1,which.min)
  
    pop<-data.frame(winner=unique(dat$model)[winners],
                    pari=rep(pari,length(winners)),
                    type=rep(type,length(winners)))#,
#                    n.samp=dat$n.samp[unique(dat$sim)])
  
    big.pop <- rbind(big.pop,pop)
  }
}


library(ggplot2)


p <- ggplot(big.pop)

p <- p + geom_histogram(aes(x=winner))#,fill=n.samp))

p <- p + facet_wrap(type~pari)

p


