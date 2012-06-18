
arms<-function(data,  family, nsample = 20, nbest = 20, criterion= "both",
     weight = "aic", maxVar=10){ 
#data is a data frame including the response variable ( first column) and the 
#explanatory variables 
#maxVar is the maximum number of explanatory variables to include in the model 

#Initial data set
data0<-data
#Var0 is a vector with the explanatory variables names
nom0<-names(data)
Var0<-nom0[2:length(nom0)]
#nomcoef0 contains the parameters names 
nomcoef0<-nom0
nomcoef0[1]<-"(Intercept)"
#nobs is the number of data 
nobs=length(data[,1])

#If the number of explanatory variables is higher than maxVar, the function 
#varSelec removes variables with a stepwise selection
if ((dim(data)[2]-1)>maxVar){
data<-varSelec(data=data0,family=family,maxVar=maxVar)}else{data=data0}

#Var is a vector with the explanatory variables names
nom<-names(data)
Var<-nom[2:length(nom)]
#nobs is the number of data 
nobs=length(data[,1])
names.coef=c("(Intercept)",Var)

##Creation of the linear models formulas from all the possible combinations of
#the explanatory variables
#Sorties[[i]] includes all the combinations of Var with i elements
Sorties <- list()
#the number of parameters (including the intercept) must be lower than the size 
#of train.sample (=nobs/2)
maxp<- min(round(nobs/2),length(Var)+1)
for (i in 1:(maxp-1)){
  Sorties [[i]] <- combn(Var,i)
}
#Vec_Fin includes the formulas of the models : Y~Xl+...+Xk
Vec_Fin <- vector()
#Each element of the list label contains a vector of a model explanatory 
#variables 
label<-list()
for (i in 1:(maxp-1)){
  Tab <- Sorties [[i]]
  for (j in c(1:dim(Tab)[2])) {
    mod <- Tab [1,j]
    label<-c(label,list(Tab[,j]))
    if (i>1) {
      for (n in c(2:i)) {
        mod <- paste(mod,Tab[n,j],sep="+")
      }
    }
  Vec_Fin <- c(Vec_Fin,paste (nom[1],mod,sep="~"))
  }
}

#The constant model is added
a<-paste(nom[1],1,sep="~")
Vec_Fin<-c(a,Vec_Fin)
label<-c(NA,label)

#The number of models created must be higher than or equal to nbest
if (length(Vec_Fin) < nbest) nbest=length(Vec_Fin)
#Matrix of model weights
Stockpoids=matrix(0,ncol=length(Vec_Fin),nrow=nsample)

#Number of colinearity with the intercept or between factors problems in the 
#train sample
cons.pb=0
col.pb=0
# which occurs successively 
pbcons=0
pbcol=0 
#number of convergence, prediction separation and model weight problems 
cv.error=0
sep.error=0
w.error=0
#number of convergence, model weight problems
# which occurs successively 
pbcv=0
pbw=0

#Loop of the nsample permutations

iteration=1
while (iteration <= nsample) {  

  # data is splitted in two parts to obtain the training and test samples
  individustire <- sample(1:nobs,round(nobs/2),replace=FALSE,prob=NULL)
  train.sample<-data[individustire,]
  test.sample<-data[-individustire,]
  
  #If some columns in the training sample take the same values
  #another sample is drawn
  colelim={} 
  for(j in c(2:dim(train.sample)[2])){
    if (max(train.sample[,j]-min(train.sample[,j])) == 0) colelim=c(colelim,j)   
  }    
  if (length(colelim)!=0) {
    pbcons=pbcons+1
    cons.pb=cons.pb+1
    if(pbcons>100) stop("one factor may be constant")
    next 
  }else pbcons=0
   
  #If there are more than two explanatory variables, and two of them are 
  #proportionnal, another sample is drawn
   if (dim(train.sample)[2] >= 3){
    for(j in c(2:((dim(train.sample)[2])-1))){
      restcol=(j+1):(dim(train.sample)[2])
      if (length(which(cor(train.sample[,j],train.sample[,restcol])==1))!=0){
        colelim=c(colelim,j)
      }
    }
  }
  if (length(colelim)!=0){
    pbcol=pbcol+1
    col.pb=col.pb+1
    if(pbcol>100) stop("problem of colinearity between factors")
    next 
  }else pbcol=0
  
  ##Preliminary convergence test

  if (maxp==length(Var)+1){
    #if the full model converges, all the models will converge
    test.cv<-fullModel(data=train.sample,family=family)$cv
    if (test.cv==FALSE){ cv.error=cv.error+1
      pbcv = pbcv +1
      if (pbcv > 100) stop("100 samples created consecutively without one of them are estimable")
      next} else pbcv=0
  }

  if (maxp != length(Var)+1){
    #else we check the convergence of all the models with maxp parameters
    mod<-choose(length(Var),0:(maxp-2))
    modpre<-sum(mod)
    test.cv<-matrix(ncol=1,nrow=length(Var))
    for(m in 1:length(Var)){
      test.cv[m]<-glm(formula(Vec_Fin[modpre+m]),data=train.sample,family=family)$conv
    }
    if (min(test.cv)==0){  cv.error=cv.error+1
      pbcv = pbcv +1
      if (pbcv > 100) stop("100 samples created consecutively without one of them are estimable")
      next} else pbcv=0
  }
  
  #Peliminary separation test
  if(family$family=="binomial") {
    if (maxp==length(Var)+1){
      #full model predictions
      yhat<-fullModel(data=train.sample,family=family)$fitted.values
      test.sep<-yhat[yhat<= 0.0001 | yhat >=0.9999]
      if(length(test.sep)!=0){
        sep.error=sep.error+1 
        }
    }
    #else we check all the models with maxp parameters
    if (maxp != length(Var)+1){
      test.sep<-matrix(ncol=1,nrow=length(Var))
      for(m in 1:length(Var)){
        yhat<-glm(Vec_Fin[modpre+m],data=train.sample,family=family)$fitted.values
        test.sep[m]<-length(yhat[yhat<= 0.0001 | yhat >=0.9999])
      }
     if(min(test.sep)!=0){
        sep.error=sep.error+1
        }
    }   
  }

  #Model selection
  TabAIC<-matrix(ncol=1,nrow=length(Vec_Fin))
  TabBIC<-matrix(ncol=1,nrow=length(Vec_Fin))
  for (i in c(1:length(Vec_Fin))){
    model<- glm(formula(Vec_Fin[i]),family=family,data=train.sample)
    TabAIC[i,]=model$aic
    TabBIC[i,]=model$aic+length(model$coef)*(log(nobs)-2)
  }
    
  #Models are selected using the chosen criteria
  #'both' : "nbest" models by the AIC and "nbest" models by the BIC
  #if the model has not been estimated because of a problem of convergence or 
  #separation in the prediction, the criterion is replaced by a value higher 
  #than the maximum in order to not select it
  
  if (criterion=='both'){  
    selec_aic<-sort(TabAIC,decreasing=FALSE)
    selec_bic<-sort(TabBIC,decreasing=FALSE) 
    selec_aic=selec_aic[nbest]
    selec_bic=selec_bic[nbest]
    modsel={}
    for (i in c(1:nrow(TabAIC))) {
      if (TabAIC[i]<=selec_aic) modsel=c(modsel,i) 
      if (TabBIC[i]<=selec_bic) modsel=c(modsel,i) 
    }
    modsel<-modsel[duplicated(modsel)==FALSE]
  }
  
  if (criterion=='aic'){
    selec_aic<-sort(TabAIC,decreasing=FALSE)
    selec_aic=selec_aic[nbest]
    modsel={}
    for (i in c(1:nrow(TabAIC))) {
      if (TabAIC[i]<=selec_aic) modsel=c(modsel,i) 
    }
  }
  
  if (criterion=='bic'){
    selec_bic<-sort(TabBIC,decreasing=FALSE)
    selec_bic=selec_bic[nbest]
    modsel={}
    for (i in c(1:nrow(TabBIC))) {
      if (TabBIC[i]<=selec_bic) modsel=c(modsel,i) 
    }
  }

   #Estimation of the parameters of the selected models with the "glm" function,
   #applied on the train.sample
   Tabcoef<-matrix(0,nrow=length(modsel),ncol=length(names.coef))
   Tabvariance<-matrix(ncol=1,nrow=length(modsel))
   
   for (i in 1:length(modsel)){
     model<- glm(formula(Vec_Fin[modsel][i]),family=family,data=train.sample)
     Tabvariance[i]<-model$deviance/model$df.residual
      
     for (j in 1:length(model$coef)){
      Tabcoef[i,which(names(model$coef)[j]==names.coef)]=model$coef[j]
     }  
   }
  
  #Test.sample prediction for all the selected models
  vartest.sample<-cbind(rep(1,dim(test.sample)[1]),as.matrix(test.sample[,2:dim(test.sample)[2]]) )
  ychapeau<-as.matrix(vartest.sample)%*%t(Tabcoef)
  
  #Likelihoods of the test.sample for all the selected models
  if(family$family=="gaussian"){
    delta<-colSums((test.sample[,1]-ychapeau)^2)
    b<-Tabvariance^(-0.5*dim(test.sample)[1])*exp(-delta/(2*Tabvariance))
  }
  
  if(family$family=="binomial"){
    log.b<-matrix(ncol=1,nrow=length(modsel))
    for (n in 1:length(modsel)){
      log.b[n]<-sum(test.sample[,1]*log(plogis(ychapeau[,n]))+(1-test.sample[,1])
      *log((1-plogis(ychapeau[,n])))) 
    }
    b<-exp(log.b)
  }
  
  #Number of parameters in each model
  Tabmodele={}
  Tabmodele<-Tabcoef
  #Tabmodele contains 1 if the coefficient is included in the model, else 0
  Tabmodele[Tabcoef[,1:(dim(Tabcoef)[2])]!=0 ]<- 1
  #The sum of the row k is the number of parameters in the model k
  nbreparam<-rowSums(Tabmodele)
  
  if (weight=="likeli") Weight<- b
  if (weight=="aic") Weight<- b*exp(-nbreparam)
  
  #Control of the weight values
  control<-Weight
  #garde is the vector of selected models whose weight is not "na" or "nan"  
  garde=numeric()
  for (i in c(1:length(modsel))){
    if ( is.nan(control[i]) == FALSE && is.na(control[i]) == FALSE){
      garde=c(garde,i)
    }
  }
  #If none of the weights are defined or different from zero the parameter 
  #estimation is not possible and another sample is drawn
  if (length(garde) == 0 | max(control[garde])==0) {pbw=pbw+1
  w.error=w.error+1}
  if (pbw > 100) stop("100 samples created consecutively without one of them are estimable")
  if (length(garde) == 0 | max(control[garde])==0) next else pbw=0

  #Model normalised weights (0 if it does not belong to modsel) 
  sommedespoids=sum(Weight[garde])
  poids<-rep(0,length(Vec_Fin))
  poids[modsel[garde]]=Weight[garde]/sommedespoids
  
  #Model weights from each permutation are stored in the matrix Stockpoids
  Stockpoids[iteration,]=poids
  
  iteration = iteration + 1
}

#Average model weights 
poidsfinal<-colMeans(Stockpoids)

##Weighting of the estimators

#Estimation by 'glm' from the initial sample for all the models 
Tabcoef<-matrix(0,nrow=length(Vec_Fin),ncol=length(nom))
for (i in c(1:length(Vec_Fin))) {
    modele <-  glm(formula(Vec_Fin[i]),family=family,data=data)
    for (j in 1:length(modele$coef)){
      Tabcoef[i,which(names(modele$coef)[j]==names.coef)]=modele$coef[j]
    }  
}
#Weighted mean of the estimators 
coefYuang<-as.vector(t(Tabcoef)%*%as.matrix(poidsfinal))
names(coefYuang)=names.coef

#Predictions
varexplicatives<-as.matrix(cbind(1,data[,2:dim(data)[2]]))
yhat<-varexplicatives%*%coefYuang
if(family$link=='logit'){
  yhat<-plogis(yhat)
} 

#Factors weights : proba(coef(i) != 0) = sum(weight(k), factor i belonging to 
#model k)
Tabmodelepoids={}
Tabmodelepoids <- Tabcoef
Tabmodelepoids[Tabcoef!=0 ] <- 1

for (i in c(1:dim(Tabcoef)[2])){
Tabmodelepoids[,i] <- Tabmodelepoids[,i] * poidsfinal
}
Tabsommepoids<-colSums(Tabmodelepoids)

#If convergence, data separation or rounding problems occured during the 
#program running, a message appears
if (cv.error >0) cat("Number of convergence problems",":",cv.error,"\n")
if (sep.error >0) cat("Number of separation problems",":",sep.error,"\n")
if (w.error > 0) cat ("Number of rounding problems",":",w.error,"\n")

#Results 
TabResult<-data.frame(coef=coefYuang,pne0=Tabsommepoids)
 
ListSummary<-list(coef=coefYuang,pne0=Tabsommepoids,fitted.values=as.vector(yhat))

VecPlot<-t(as.matrix(Tabsommepoids)[-1])
names(VecPlot)<-names.coef[-1]
rownames(VecPlot)<-"Factor weights"

ListResult<-list(TabResult,ListSummary,VecPlot,coef=coefYuang,pne0=Tabsommepoids,
fitted.values=yhat,label=label,modweights=poidsfinal,allcoef=Tabcoef)
#The first element of the list "ListResult" is visible when the output of 
#the function is printed, the second element is printed when calling "summary"
#function,the others elements are printed by calling specifically the element of
#the list (ouput$element)

class(ListResult)<-"MMIXclass"

return(ListResult)
}






