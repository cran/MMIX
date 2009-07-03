
stepSel=function(data,family,direction="both",criterion,trace=0){
#data is a data frame including the response variable (first column) and the 
#explanatory variables 

#Var is a vector with the explanatory variable names
nom<-names(data)
Var<-nom[2:length(nom)]
Y<-nom[1]
#nobs is the number of data 
nobs=length(data[,1])

#Formula for the full model
Vec_Fin <- vector()
mod <- Var[1]
if (length(Var)>1) {
  for (n in c(2:length(Var))) {
    mod <- paste(mod,Var[n],sep="+")
  }
}
Vec_Fin <- c(Vec_Fin,paste (Y,mod,sep="~"))

#glm and stepwise selection
if (criterion=='aic')k=2
if(criterion=='bic')k=log(nobs)

if (direction=="backward" | direction=="both")
{
modele   <-  glm(formula(Vec_Fin),family=family,data=data)
stepwise <- step(modele,direction=c(direction),data=data,k=k,trace=trace)
}
if (direction=="forward") 
{ 
modele   <-  glm(formula(paste(nom[1],1,sep="~")),family=family,data=data)
stepwise <- step(modele,direction=c(direction),data=data,k=k,scope=list(lower=~1,upper=formula(Vec_Fin)),trace=trace)
}

#Estimated parameter values
Tabparam<-{}
Tabparam <-data.frame(nom=names(stepwise$coef),coef=stepwise$coef)
#Parameters names
nomcoef<-nom
nomcoef[1] <- "(Intercept)"

#Tabfull contains the estimated coefficient for the selected variables, else 0
Tabfull<-matrix(0,nrow=1,ncol=length(nomcoef))
for(i in c(1:dim(Tabparam)[1]))  {
  Tabfull[1,(Tabparam$nom[i] == nomcoef)]=stepwise$coef[i]
}
Tabfull<-as.vector(Tabfull)
names(Tabfull)<-nomcoef

#Probability for the factors to be non zero

aic<-stepwise$aic
p<-length(stepwise$coef)
bic<-aic+(log(nobs)-2)*p

#Predictions
yhat<-stepwise$fitted.values

TabResult<-data.frame(coef=Tabfull)

ListSummary<-list(coef=Tabfull,aic=aic,bic=bic,fitted.values=yhat)

ListResult<-list(TabResult,ListSummary,coef=Tabfull,fitted.values=yhat,aic=aic,bic=bic)
#The first element of the list "ListResult" is visible when the output of 
#the function is printed, the others elements are printed by calling 
#specifically the element of the list "ouput$element"
class(ListResult)<-"MMIXclass"

return(ListResult)
}



