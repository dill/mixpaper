# Analysis for the Williams and Thomas data
# See JCRM paper for details

library(mmds)

par(mfrow=c(1,3))

###########################
# Harbour seal
# load the data
hs<-read.csv(file="hs.csv")
# fit a 1-point mixture, truncation at 500m
hs.fit1<-fitmix(hs,mix.terms=1,width=500)
# find best AIC model
hs.best<-step.ds.mixture(hs.fit1)

# plot that
plot(hs.best,style="hist", breaks=c(0,11,22,33,43,65,87,109,130,152,174,196,217,239,261,283,304,326,348,370,391,413,435,457,478,500),main="Harbour seal (in water)")


###########################
# Harbour porpoise
hp<-read.csv(file="hp.csv")
# fit a 1-point mixture, truncation at 500m
hp.fit1<-fitmix(hp,mix.terms=1,width=500)
# find best AIC model
hp.best<-step.ds.mixture(hp.fit1)

# plot that
plot(hp.best,style="hist", breaks=c(0,20,40,61,82,124,166,208,250,291,333,375,416,458,500),main="Harbour porpoise")


###########################
# Humpback
hb<-read.csv(file="hb.csv")
# fit a 1-point mixture, truncation at 2000m
hb.fit1<-fitmix(hb,mix.terms=1,width=2000)
# find best AIC model
hb.best<-step.ds.mixture(hb.fit1)

# plot that
plot(hb.best,style="hist", breaks=c(0, 83.5, 167, 250.5, 333, 500, 667, 833, 1000, 1167, 1333, 1500, 1667, 1833, 2000),main="Humpback")

