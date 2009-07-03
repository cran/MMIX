
fullModel<-function(data,family){
#data is a data frame including the response variable (first column) and the 
#explanatory variables 

#Var is a vector with the explanatory variable names
nom<-names(data)
Var<-nom[2:length(nom)]
#nobs is the number of data 
nobs=length(data[,1])

#Creation of the full model formula
#If there are more than one explanatory variable, "+" is put between them, and
#the response variable is separated by "~" from the other variables
mod <- Var[1]
if (length(Var)>1) {
  for (l in c(2:length(Var))) {
  mod <- paste(mod,Var[l],sep="+")
  }
}
Vec_Fin <- paste (nom[1],mod,sep="~")

#Coefficient estimation with "glm" function  
model  <-  glm(formula(Vec_Fin),family=family,data=data)
Tabparam <- data.frame(nom=names(model$coef),coef=model$coef)

aic<-model$aic
p<-length(model$coef)
bic<-aic+(log(nobs)-2)*p
     
#The function returns a data frame Tabcoef containing the coefficient 
#names and values, and the probabilty to be non-zero
Tabcoef <-data.frame(coef=model$coef)

#Predictions
yhat<-model$fitted.values

#The first element of the list "ListResult" is visible when the output of 
#the function is printed, the second element appears when calling the 
#function summary, the others elements are printed by calling specifically 
#the element of the list "ouput$element"
ListSummary<-list(coef=model$coef,aic=aic,bic=bic,
cv=model$converged,fitted.values=yhat)

ListResult<-list(Tabcoef,ListSummary,coef=model$coef,
aic=aic,bic=bic,cv=model$converged,fitted.values=yhat)

class(ListResult)<-"MMIXclass"

return(ListResult)

}


