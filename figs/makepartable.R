# make the table of parameters


Model & $\beta_1$ & $\beta_2$ & $\beta_3$ & $\pi_1$ & $\pi_2$ & $P_a$ \\
LT &1 & -0.223 & -1.897 & &  0.3 & &  0.369\\
 &2 & -0.511 & -2.303 & &  0.7 & &  0.514\\
 &3 &  2.303 & -1.609 & & 0.15 & &  0.363\\
 &4 & -0.357 & -2.996 & &  0.6 & &  0.471\\
PT &1 & -0.223 & -1.897 & &  0.3 & &  0.24\\
 &2 & -0.511 & -2.303 & &  0.7 & &  0.384\\
 &3 &  2.303 & -1.609 & & 0.15 & &  0.218\\
 &4 & -0.357 & -2.996 & &  0.6 & &  0.378\\
3pt & 1 &  -0.22 &  -0.69 &  -2.3 & 0.3 & 0.3 & 0.505\\  
 & 2 &   2.71 &  -1.39 &  -3.0 & 0.1 & 0.4 & 0.257\\
COV & 1 & -2.303 -0.288 -0.511 & 0.4 & & 0.422
 & 2 & -1.609 -0.223 -0.916 0.4 & & 0.389




#mod<-"nocov"
#parmat[1,]<-c(log(0.8),log(0.15), inv.reparam.pi(0.3))
#parmat[2,]<-c(log(0.6),log(0.1), inv.reparam.pi(1-0.3))
#parmat[3,]<-c(log(10),log(0.2), inv.reparam.pi(0.15))
#parmat[4,]<-c(log(0.7),log(0.05),inv.reparam.pi(0.6))
#mod<-"pt"
#parmat[1,]<-c(log(0.8),log(0.15), inv.reparam.pi(0.3))
#parmat[2,]<-c(log(0.6),log(0.1), inv.reparam.pi(1-0.3))
#parmat[3,]<-c(log(10),log(0.2), inv.reparam.pi(0.15))
#parmat[4,]<-c(log(0.7),log(0.05),inv.reparam.pi(0.6))
#"3pt"
#parmat[1,]<-c(log(0.8),log(0.5),log(0.1),inv.reparam.pi(rep(1/3,3))[1:2])
#parmat[2,]<-c(log(15),log(.25),log(0.05),inv.reparam.pi(c(0.1,0.4,0.5))[1:2])


"covsim 1"
mixp<-inv.reparam.pi(0.4)
true.pars<-c(log(c(0.1,0.75,0.6)),mixp)

"covsim2"
true.pars<-c(log(0.2),log(0.8),log(0.4), inv.reparam.pi(0.4))




#%> true.p
#%          t id          model
#%1  0.368835  1  No covariates
#%2  0.513678  2  No covariates
#%3  0.362814  3  No covariates
#%4  0.470853  4  No covariates
#%17 0.239692  1 Point transect
#%18 0.384326  2 Point transect
#%19 0.217625  3 Point transect
#%20 0.378057  4 Point transect
#%33 0.422207  1      Covariate
#%34 0.389230  2      Covariate
#%35 0.504756  1        3-point
#%36 0.256582  2        3-point



