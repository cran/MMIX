varSelec<-function(data,family,maxVar=10,trace=0){
#data is a data frame including the response variable (first column) and the
#explanatory variables
#Var is a vector with the explanatory variables names
nom<-names(data)
Var<-nom[2:length(nom)]
Y<-nom[1]
#nobs is the number of data 
nobs=length(data[,1])

#Formula for the full model
Vec_Fin <- vector()
mod <- Var [1]
if (length(Var)>1) {
  for (l in c(2:length(Var))) {
    mod <- paste(mod,Var[l],sep="+")
  }
}
Vec_Fin <- c(Vec_Fin,paste (Y,mod,sep="~"))

#Formula of the constant model 
a<-paste(nom[1],1,sep="~")

modele   <-  glm(formula(a),family=family,data=data)
stepwise <- step(modele,direction=c("forward"),data=data,k=2,steps=maxVar,
scope=list(lower=~1,upper=formula(Vec_Fin)),trace=trace)


VarSelec<-matrix(ncol=(length(stepwise$coef)-1),nrow=nobs)
VarSelec<-data.frame(VarSelec)
for ( i in 2:length(stepwise$coef)){
VarSelec[,i-1]=data[, which(names(stepwise$coef[i])==names(data))]}
names(VarSelec)<-names(stepwise$coef)[2:length(stepwise$coef)]
data2<-cbind(data[,1],VarSelec)
names(data2)<-c(names(data)[1],names(VarSelec))

return(data2)
}



