# Model & Covariates & AIC & $\hat{P_a}$ & $\% CV \hat{P_a}$ & K-S $p$\\
grab_results<-function(model){
   
   mod<-paste("Hn ",model$mix.terms,"-pt",sep="")
   covars<-model$model.formula   

   cv.p<-100*summary(model)$average.p.cv

   cat(mod," & ",covars," & ",round(model$aic,2)," & ",round(model$pa,3),
       " & ",round(cv.p,2)," & ",round(model$ks$p,2),"\\\\\n")

}

