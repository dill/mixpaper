# make the line to go in the table
ml<-function(obj){

   pp<-paste(obj$model.formula,
             round(obj$aic,2),
             round(obj$pa,2),
             round(obj$pa.se/obj$pa*100,2),
             round(obj$ks$p,2),sep=" & ")

   cat(pp,"\n")

}
