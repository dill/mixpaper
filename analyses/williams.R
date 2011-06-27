# Analysis for the Williams and Thomas data
# See JCRM paper for details

library(mmds)

par(mfrow=c(1,3))

###########################
# Dall's porpoise
# load the data
dp<-read.csv(file="dp.csv")
# fit a 1-point mixture, truncation at 700m
dp.fit1<-fitmix(dp,mix.terms=1,width=700)
# find best AIC model
dp.best<-step.ds.mixture(dp.fit1)

# plot that
plot(dp.best,style="hist", breaks=c(0, 23.5, 47, 70.5, 93, 140, 187, 233, 280, 327, 373, 420, 467, 513, 560, 607, 653, 700),main="Dall's porpoise")


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

