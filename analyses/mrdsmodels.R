library(mrds)
set.seed(145)

# get the CVs using the methods in mrds, rather than Distance

### Ants
# load some data
# scale log(1.481), 0.4113e-3, -1.475, -1.020, -0.1488
# shape 1.108
ants<-read.csv(file="woodants.csv")
# nest.size+habitat hazard-rate
ants$ns<-(ants$nest.size-mean(ants$nest.size))/sd(ants$nest.size)
antss<-ddf(dsmodel = ~mcds(key = "hr", formula = ~ns+as.factor(habitat)), 
           data=ants, method = "ds", meta.data = list(width = 25),
           control=list(initial=list(scale=c(log(1.481), 0.4113e-3, -1.475, -1.020, -0.1488),shape=1.108),refit=FALSE))
summary(antss)


### Pilot whales
# 4040, -0.5123, 0.3909
fin<-read.csv(file="danpike.csv")
# half-normal + cos(2) BSS (continuous)
levels(fin$BSS)[5]<-"3.5"
fin$BSS<-as.numeric(as.character(fin$BSS))
finn<-ddf(dsmodel = ~mcds(key = "hn", adj.series="cos", adj.order=c(2),
          adj.scale="width",formula = ~BSS), data=fin, method = "ds", 
          control=list(initial=list(scale=c(log(4040),-0.5123),adjustment=0.3909),showit=TRUE,refit=FALSE),
#         control=list(showit=TRUE),
          meta.data = list(width = 3000))

summary(finn,se=TRUE)

# Amakihi
# 25.51, 2.305
# obs -0.1037, 0.4505
# mas -0.8272e-3
amakihi<-read.csv(file="amakihi.csv")
amakihi$mas<-(amakihi$mas-mean(amakihi$mas))/sd(amakihi$mas)
# hr obs+mas
am<-ddf(dsmodel=~mcds(key="hr",formula=~as.factor(obs)+mas), 
        data=amakihi, method = "ds", 
        meta.data = list(width = 82.5,point=TRUE),
        control=list(initial=list(scale=c(log(25.51),-0.1037, 
                                          0.4505,-0.8272e-3),
                                  shape=2.305)))

summary(am)

