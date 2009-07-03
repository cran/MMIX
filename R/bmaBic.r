
bmaBic=function(data,family,maxVar=10){
#data is a data frame including the response variable (first column) and the 
#explanatory variables 

#Var0 is a vector with the explanatory variables names
nom0<-names(data)
Var0<-nom0[2:length(nom0)]
#nomcoef0 contains the parameter names 
nomcoef0<-nom0
nomcoef0[1]<-"(Intercept)"

#If the number of explanatory variables is higher than maxVar, the function 
#varSelec selects variables by stepwise selection
if ((dim(data)[2]-1)>maxVar){
data2<-varSelec(data=data,family=family,maxVar=maxVar)}else{data2=data}

#Var is a vector of the explanatory variable names
nom<-names(data2)
Var<-nom[2:length(nom)]
#nomcoef contains the names of parameters 
nomcoef<-nom
nomcoef[1]<-"(Intercept)"
nobs=dim(data2)[1] 
#Creation of all the linear models from all combinations of the explanatory 
#variables
#Sorties[[i]] contains all the combinations of Var with i elements
Sorties <- list()
for (i in c(1:length(Var))) {
     Sorties [[i]] <- combn(Var,i)
}

#Vec_Fin contains the formulas of the models : Y~Xl+...+Xk
Vec_Fin <- vector()
#each element of the list label contains a vector of a model explanatory 
#variables 
label<-list()
for (i in c(1:length(Var))) {
  Tab <- Sorties [[i]]
  for (j in c(1:dim(Tab)[2])) {
    mod <- Tab [1,j]
    label<-c(label,list(Tab[,j]))
    if (i>1) {
      for (n in c(2:i)) {
        mod <- paste(mod,Tab[n,j],sep="+")
      }
    } 
    Vec_Fin <- c(Vec_Fin,paste(nom[1],mod,sep="~"))
  }
}

#The constant model is added
a<-paste(nom[1],1,sep="~")
Vec_Fin<-c(a,Vec_Fin)
label<-c(NA,label)

#Table for BIC
TabBIC <-matrix(nrow=length(Vec_Fin),ncol=1)
#Table for parameter estimates
Tabparam<-matrix(0,ncol=length(nomcoef0),nrow=length(Vec_Fin))
#Table for standard deviation of the parameter estimates
TabSE<-matrix(0,ncol=length(nomcoef0),nrow=length(Vec_Fin))

#Glm estimation, BIC calculation
Model0<- glm(formula(a),family=family,data=data2)

for (n in c(1:(length(Vec_Fin)))) {
  Model <-  glm(formula(Vec_Fin[n]),family=family,data=data2)

  p=length(Model$coef)
  BIC=Model$aic-2*p+p*log(nobs)-(Model0$aic)
  TabBIC[n]<-BIC
     
  for(j in 1:length(Model$coef)){
    Tabparam[n,which(names(Model$coef)[j]==nomcoef0)]=Model$coef[j]
    TabSE[n,which(nomcoef0==names(Model$coef)[j])]=summary(Model)[["coefficients"]][,2][j]
  }
}

#Model posterior weights are estimated by BIC
Weight_Model<-exp(-TabBIC/2)/sum(exp(-TabBIC/2))

#Factor weights : proba(factor(i) != 0)  = sum(weight(model n), factor(i) 
#belongs to model n )
Pne0<-matrix(0,ncol=1,nrow=length(nomcoef0))
for( j in 1:length(nomcoef0)){
  Pne0[j]=sum(Weight_Model[which(Tabparam[,j]!=0)])
}

#Estimator ponderation by model performances
ParamBMA<-matrix(nrow=1,ncol=length(nomcoef0))
ParamBMA=t(Weight_Model)%*%Tabparam
ParamBMA<-as.vector(ParamBMA)
names(ParamBMA)=nomcoef0

#Estimated standard deviation of coefficients
Stdhat<-matrix(ncol=1,nrow=length(nomcoef0))
for(i in 1:length(nomcoef0)){
  Stdhat[i]=sqrt(sum(Weight_Model*(TabSE[,i]^2+Tabparam[,i]^2))-ParamBMA[i]^2)
}

#Predictions
varexplicatives<-as.matrix(cbind(1,data[,2:dim(data)[2]]))
yhat<-varexplicatives%*%ParamBMA
if(family$link=='logit'){
  yhat<-plogis(yhat)
} 
selec_bic<-sort(TabBIC,decreasing=FALSE)[1:3]
BIC_selec<-TabBIC[which(TabBIC <= selec_bic[3])]
Mod_selec<-Vec_Fin[which(TabBIC<=selec_bic[3])]
BestModels<-cbind(Mod_selec,BIC_selec)
BestModels[rank(BIC_selec),]=BestModels

#TabResult is a data frame including the variable names, the estimated 
#coefficients by BMA and the posterior probability of the coefficients to be 
#non zero
TabResult <-data.frame(coef=ParamBMA,pne0=Pne0,sd=Stdhat)

ListSummary<-list(coef=ParamBMA,pne0=as.vector(Pne0),fitted.values=as.vector(yhat),sd=as.vector(Stdhat),BestModels=BestModels)

VecPlot<-t(as.matrix(Pne0)[-1])
rownames(VecPlot)<-"Factor weights"
names(VecPlot)<-nomcoef0[-1]

ListResult<-list(TabResult,ListSummary,VecPlot,coef=ParamBMA,fitted.values=as.vector(yhat),
pne0=as.vector(Pne0),sd=as.vector(Stdhat),BestModels=BestModels,label=label,modweights=as.vector(Weight_Model),allcoef=Tabparam)

#The first element of the list "ListResult" is visible when the output of 
#the function is printed, the second element is printed when calling "summary" 
#function, the others elements are printed by calling specifically the element
# of the list (ouput$element)

class(ListResult)<-"MMIXclass"

return(ListResult)

}