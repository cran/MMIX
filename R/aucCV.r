aucCV=function(data,method=1,np,random=TRUE,npermu=100,file=NULL,...) {
#data is a data frame with the response variable (first column) and the 
#explanatory variables
#method=1: fullModel, method=2: stepSelec, method=3: bmaBic, method=4: mixAic, 
#method=5: arms
#"..." refers to the specific arguments of the method  

#np must be higher than or equal to 1 and lower than the number of data minus 
#the number of model parameters

tab0<-data[data[,1]==0,]
tab1<-data[data[,1]==1,]
n0=dim(tab0)[1]
n1=dim(tab1)[1]

if (np > min(n0,n1) -1 ){print(paste("np must be lower than ", min(n0,n1) -1))
 stop}
if (np < 1 ){ print("np must equal to or higher than 1" ) 
stop}

#Number of convergence and separation problems
cv.error=0
sep.error=0
#Number of successive convergence problems
pbcv=0

###AUC estimated by leave-one-pair-out cross validation
if(np ==1 & random==FALSE){
  n.sample=n0*n1
  #combinations of test couples 
  ind.test<-{}
  for (i in 1:n1){              
    ind.test<-rbind(ind.test,cbind(which(data[,1]==0),rep(which(data[,1]==1)[i],n0)))
  }  

  #Tabindref includes the response variable observations, and their predictions
  # derived by cross-validation
  Tabindref<-matrix(ncol=2,nrow=n.sample,dimnames=list(NULL,c("ind","ref")))
  
  auc<-numeric()
  
  
  for(i in 1:n.sample) {

    #the model is estimated on the whole sample, expecting a couple of 
    #individuals (one of each label, 1 and 0)
  
      test.sample<-rbind(data[ind.test[i,1],],data[ind.test[i,2],])
      work.sample<-data[-ind.test[i,],]
   
    #Preliminary tests of the convergence and data separation
    #If the full model can not be estimated for the work.sample, the program 
    #switches to the next iteration
  
      fullResult<-fullModel(family=binomial('logit'),data=work.sample)
      if (fullResult$cv == FALSE ){
        pbcv<-pbcv+1
        cv.error=cv.error+1
      }
      if (pbcv > 100) stop("100 samples created consecutively without one of them are estimable")
      if (fullResult$cv==FALSE) next else pbcv=0
  
    
      yhat<-fullResult$fitted.values
      test.sep<-yhat[yhat<= 0.0001 | yhat >=0.9999]
      if (length(test.sep)!=0){ 
        sep.error=sep.error+1
      }
  
    #Model estimation with the chosen method (the default is fullModel)
    #method==1 : fullModel
    if (method==1) coeffestime<-fullModel(family=binomial('logit'),data=work.sample)$coef
    #method==2 : stepSelection
    if (method==2) coeffestime<-stepSel(data=work.sample,family=binomial('logit'),...)$coef
    #method==3 : BMAmodel
    if (method==3) coeffestime<-bmaBic(family=binomial('logit'),data=work.sample,...)$coef
    #method==4 : mixAic 
    if (method==4) coeffestime<-mixAic(family=binomial('logit'),data=work.sample,...)$coef
    #method==5 : ARMS
    if(method==5)coeffestime<-arms(data=work.sample,family=binomial('logit'),...)$coef
    
    
    #Explanatory variables of the data i with 1 for the intercept
    varexplicatives<-as.matrix(cbind(c(1,1),data[ind.test[i,],2:dim(data)[2]]))
    #Estimated values of the coefficients
    coefficient<-as.matrix(coeffestime)
    #Prediction of the data i
    ychapeau<-varexplicatives%*%coefficient
    # proba(Y =1) = inverse of the logistical function
    proba<-plogis(ychapeau)
    #Tabindref includes the prediction and the observed value 
    Tabindref=cbind(proba,data[ind.test[i,],1])
      
    if(length(file) != 0){
    #Predictions are stored in the file "file" , in case the program would stop
    write.table(Tabindref,file,sep="\t",append=TRUE,row.names=TRUE,col.names=FALSE)
    }
    
    auc=c(auc,as.numeric(eval(proba[1]<proba[2]))+0.5*as.numeric(eval(proba[1]==proba[2])))
  }

}

###AUC estimated  by leave-np-pair-out cross validation on npermu permuations
if(np>1 | random==TRUE){
  
  #Tabindref includes the predictions and references derived by cross-validation
  Tabindref<-{}
  #auc is the vector of the npermu auc values
  auc<-matrix(ncol=1,nrow=npermu)
  
   permu=1
   while(permu <= npermu) {        

    #the model is estimated on the train sample and auc calculated on the test 
    #sample (np couples of individuals 0 and 1)
      ind0.test<-sample(1:n0,np,replace=FALSE,prob=NULL)
      ind1.test<-sample(1:n1,np,replace=FALSE,prob=NULL)
      test.sample<-rbind(tab0[ind0.test,],tab1[ind1.test,])
      train.sample<-rbind(tab1[-ind1.test,],tab0[-ind0.test,])
    
    #Preliminary tests of the convergence and prediction separation
    #If the full model can not be estimated for the work.sample, the program
    # switches to the next iteration
  
      fullResult<-fullModel(family=binomial('logit'),data=train.sample)
      if (fullResult$cv == FALSE ){
        pbcv<-pbcv+1
        cv.error=cv.error+1
      }
      if (pbcv > 100) stop("100 samples created consecutively without one of them are estimable")
      if (fullResult$cv==FALSE) next else pbcv=0
  
    
      yhat<-fullResult$fitted.values
      test.sep<-yhat[yhat<= 0.0001 | yhat >=0.9999]
      if (length(test.sep)!=0){ 
        sep.error=sep.error+1
      }
  
    #Model estimation with the chosen method (the default is fullModel)
    #method==1 : fullModel
    if (method==1) coeffestime<-fullModel(family=binomial('logit'),data=train.sample)$coef
    #method==2 : stepSelection
    if (method==2) coeffestime<-stepSel(data=train.sample,family=binomial('logit'),...)$coef
    #method==3 : BMAmodel
    if (method==3) coeffestime<-bmaBic(family=binomial('logit'),data=train.sample,...)$coef
    #method==4 : mixAic 
    if (method==4) coeffestime<-mixAic(family=binomial('logit'),data=train.sample,...)$coef
    #method==5 : ARMS
    if(method==5)coeffestime<-arms(data=train.sample,family=binomial('logit'),...)$coef
    
    
    #Explanatory variables of the test.sample with 1 for the intercept
    varexplicatives<-as.matrix(cbind(rep(1,np*2),test.sample[,2:dim(data)[2]]))
    #Estimated values of the coefficients
    coefficient<-as.matrix(coeffestime)
    #Prediction of the data i
    ychapeau<-varexplicatives%*%coefficient
    # proba(Y =1) = inverse of the logistic function
    proba<-plogis(ychapeau)
    #Tabindref includes the prediction and the observed value 
    Tabindref=cbind(proba,test.sample[,1])
    
    proba0<-proba[test.sample[,1]==0]
    proba1<-proba[test.sample[,1]==1]
    compteur=0
    for( i in 1:length(proba0)){
       for( j in 1:length(proba1)){
           compteur=compteur+as.numeric(eval(proba0[i]<proba1[j]))+0.5*as.numeric(eval(proba0[i]==proba1[j]))
       }
    }
    
    auc[permu]=compteur/(length(proba1)*length(proba0))
    
    if (length(file) !=0){
    #auc are stored in a file, in case the program would stop
    write.table(auc[permu],file,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE)
   }
   
   permu=permu+1
  }

}

AUC.CV=mean(auc)


##AUC calculated on the whole sample

#Paramters are estimated with the chosen method
#method==1 : fullModel;
if (method==1) proba<-fullModel(family=binomial('logit'),data=data)$fitted.values
#method==2 : stepSelection;
if (method==2) proba<-stepSel(family=binomial('logit'),data=data,...)$fitted.values
#method==3 : bmaBic;
if (method==3) proba<-bmaBic(family=binomial('logit'),data=data,...)$fitted.values
#method==4 : mixAic ;
if (method==4) proba<-mixAic(family=binomial('logit'),data=data,...)$fitted.values
#method==5 : arms;
if(method==5) proba<-arms(data=data,,family=binomial('logit'),...)$fitted.values

# Tabindref includes the predictions and the observed values
Tabindref={}
Tabindref<-data.frame(indicateur=proba,ref=data[,1])

#Estimation of the AUC with Mann-Whitney statistic
yhat1<-Tabindref$indicateur[which(Tabindref$ref==0)]
yhat2<-Tabindref$indicateur[which(Tabindref$ref==1)]

compteur=0
for( i in 1:length(yhat1)){
for( j in 1:length(yhat2)){
compteur=compteur+as.numeric(eval(yhat1[i]<yhat2[j]))+0.5*as.numeric(eval(yhat1[i]==yhat2[j]))
}}

AUC=compteur/(length(yhat2)*length(yhat1))

if (cv.error!=0) cat("Number of convergence problems",":",cv.error,"\n")
if (sep.error!=0) cat("Number of separation problems",":",sep.error,"\n")

tab.AUC<-data.frame(AUC=AUC,AUC.CV=AUC.CV)

return(tab.AUC)

}