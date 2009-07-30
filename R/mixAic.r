#Mixing model with AIC ponderation

mixAic=function(data,family,maxVar=10){
#data is a data frame including the response variable (first column) and the 
#explanatory variables
#family is a description of the error distribution  (gaussian, binomial)
#maxVar is the maximum number of explanatory variables to include in the model 

##Var0 is a vector with the explanatory variable names
#nom0<-names(data)
#Var0<-nom0[2:length(nom0)]
##nomcoef0 contains the parameter names 
#nomcoef0<-nom0
#nomcoef0[1]<-"(Intercept)"
##nobs is the number of data 
#nobs=length(data[,1])

#If the number of explanatory variables is higher than maxVar, the function
#varSelec selects variables by stepwise selection
if ((dim(data)[2]-1)>maxVar){
data2<-varSelec(data=data,family=family,maxVar=maxVar)}else{data2=data}

#Var is a vector of the explanatory variable names
nom<-names(data2)
Var<-nom[2:length(nom)]
nomcoef<-c("(Intercept)",nom[2:length(nom)])

#Creation of all the linear models from every combinations of the explanatory
# variables
# Sorties[[i]] contains all the combinations of Var with i elements
Sorties <- list()
for (i in c(1:length(Var))){
  Sorties [[i]] <- combn(Var,i)
}

#Vec_Fin contains the formulas of all the models : Y~Xl+...+Xk
Vec_Fin <- vector()
#The list label contains the vectors of explanatory variables of all the models
label<-list()
for (i in c(1:length(Var))) {
    Tab <- Sorties [[i]]

    for (j in c(1:dim(Tab)[2])) {
        mod <- Tab [1,j]
        label<-c(label,list(Tab[,j]))
        if (i>1) {
            for (k in c(2:i)) {
                mod <- paste(mod,Tab[k,j],sep="+")
            }
        }
        Vec_Fin <- c(Vec_Fin,paste(nom[1],mod,sep="~"))
    }
}
#The constant model is added
a<-paste(nom[1],1,sep="~")
Vec_Fin<-c(a,Vec_Fin)
label<-c(NA,label)

#Definition of matrices including the parameter estimators, their standard 
#errors and AIC
Tabcoef<-matrix(0,nrow=length(Vec_Fin),ncol=length(nomcoef))
TabSE<-matrix(0,ncol=length(nomcoef),nrow=length(Vec_Fin))
TabAIC <-matrix(ncol=1,nrow=length(Vec_Fin))

#implementation of function glm for each model
for (i in c(1:length(Vec_Fin))) {
  model<- glm(formula(Vec_Fin[i]),family=family,data=data2)
  TabAIC[i]<- model$aic
  for (j in 1:length(model$coef)){
    Tabcoef[i,which(nomcoef==names(model$coef)[j])]=model$coef[j]
    TabSE[i,which(nomcoef==names(model$coef)[j])]=summary(model)[["coefficients"]][,2][j]
  }
}

#AIC differences
Tabdifaic<- TabAIC-min(TabAIC)
#AIC weights
Tabvraissemblancegx<-exp(-0.5*Tabdifaic)
#Weights sum
sommevraiss=sum(Tabvraissemblancegx)
#Normalised weights
poids<-Tabvraissemblancegx/sommevraiss

#Weighted sum of the estimated parameters 
coefMix<-as.vector(t(Tabcoef)%*%as.matrix(poids))
names(coefMix)<-nomcoef

#Probabilty of selection of variables
Tabpoids<-Tabcoef
Tabpoids[Tabcoef!=0 ]<- 1
VarPoids=t(poids)%*%Tabpoids

#Standard deviation estimation of coefficients
Stdhat<-matrix(ncol=1,nrow=length(nomcoef))
for(i in 1:length(nomcoef)){
  Stdhat[i]=sqrt(sum(poids*(TabSE[,i]^2+(Tabcoef[,i]-coefMix[i])^2)))
}

#Predictions
varexplicatives<-as.matrix(cbind(1,data[,2:dim(data2)[2]]))
yhat<-varexplicatives%*%coefMix
if(family$link=='logit'){
  yhat<-plogis(yhat)
}
#The three models with the highest AIC
selec_aic<-sort(TabAIC,decreasing=FALSE)[1:3]
Mod_selec<-which(TabAIC<=selec_aic[3])
AIC_selec<-TabAIC[Mod_selec]
Mod_selec<-Vec_Fin[Mod_selec]
BestModels<-cbind(Mod_selec,AIC_selec)
BestModels[rank(AIC_selec),]=BestModels

TabResult<-data.frame(coef=coefMix,pne0=t(VarPoids),sd=Stdhat)

ListSummary<-list(coef=coefMix,fitted.values=as.vector(yhat),pne0=VarPoids,sd=as.vector(Stdhat),
BestModels=BestModels)

VecPlot<-t(as.matrix(VarPoids[,-1]))
rownames(VecPlot)<-"Factor weights"
names(VecPlot)<-nomcoef[-1]

ListResult<-list(TabResult,ListSummary,VecPlot,coef=coefMix,fitted.values=as.vector(yhat),
pne0=as.vector(VarPoids),sd=as.vector(Stdhat),BestModels=BestModels,label=label,
modweights=as.vector(poids),allcoef=Tabcoef)

#The first element of the list "ListResult" is visible when the output of 
#the function is printed, the second element is printed when calling "summary" 
#function, the others elements are printed by calling specifically the element
# of the list (ouput$element)
class(ListResult)<-"MMIXclass"

return(ListResult)

}