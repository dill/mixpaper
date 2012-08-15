# how many mixture components do we use for the hazard rate models?

library(ggplot2)


baf<-data.frame(n.samp=0,par.ind=0,mix.terms=0)


for(par.ind in 1:2){

  fit.mthod<-"BFGS+SANN"
  dat<-read.csv(file=paste("hr-",par.ind,"-results.csv",sep=""))

  dat<-dat[,-1]
  names(dat)<-c("mod","par.ind","n.samp","sim","ll","AIC",
                "pa","Nhat","N","mixterms")

  dat<-dat[dat$mod=="mmds-MS",]

  baf<-rbind(baf,cbind(n.samp=dat$n.samp,
                       par.ind=dat$par.ind,
                       mix.terms=dat$mixterms))

}

baf<-baf[-1,]
baf$ind<-1:nrow(baf)
baf$n.samp<-as.factor(baf$n.samp)
baf$mix.terms<-as.factor(baf$mix.terms)
baf$par.ind<-as.factor(baf$par.ind)


p<-ggplot(baf)
p<-p+labs(x="Mixture components",y="Count")
p<-p+opts(panel.grid.major=theme_blank(),
          panel.grid.minor=theme_blank(),
          legend.background=theme_blank(),
          legend.key=theme_blank(),
          panel.background=theme_rect())
 
p<-p+geom_histogram(aes(x=mix.terms,y=..count..))
p<-p+facet_grid(par.ind~n.samp)

# pull out the counts and write them...
count_text<-ggplot_build(p)[[1]][[1]]
count_text$y<-count_text$y+40
count_text<-cbind(count_text,
                  par.ind=c(rep(1,5),rep(2,5))[count_text$PANEL])
count_text<-cbind(count_text,
                  n.samp=as.factor(c(30,60,120,480,960)[c(
                       as.numeric(count_text$PANEL[count_text$par.ind==1]),
                       as.numeric(count_text$PANEL[count_text$par.ind==2])-5)]))

p<-p+geom_text(aes(x=x,y=y,label=count),data=count_text)

print(p)

ggsave(file="hazard-mixcomponents.pdf")



