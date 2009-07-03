pmseCV=function(data, method=1, np, random=TRUE, npermu=100, file=NULL, ...) {
#data is a data frame including the response variable(first column) and the 
#explanatory variables 
#method=1: fullModel, method=2: stepSelec, method=3: bmaBic, method=4: mixAic, 
#method=5: arms
#"..." refers to the specific arguments of the method  

#np must be higher than 1 and lower than the number of data minus the number of 
#model parameters
if (dim(data)[1]-np < dim(data)[2] | np<1 ){print("wrong value for np")}

#Number of errors
n.error=0
#Number of convergence problems
cv.error=0
#Number of successive convergence problems
pbcv=0

if(np==1 & random==FALSE){
  #Tabindref includes the predictions derived by cross-validation
  Tabindref<-matrix(ncol=2,nrow=dim(data)[1],dimnames=list(NULL,c("ind","ref")))
  for(i in c(1:dim(data)[1])) {
    #Data set without the i-th data 
    tabtravail<-data[-i,]
    
    #Preliminary test of the convergence
    #If one of the model can not be estimated for the work.sample, the program 
    #switches to the next iteration
    fullResult<-fullModel(family=gaussian('identity'),data=tabtravail)
      if (fullResult$cv == FALSE ){
        pbcv<-pbcv+1
        cv.error=cv.error+1
      }
    if (pbcv > 100) stop("100 samples created consecutively without one of them are estimable")
    if (fullResult$cv==FALSE) next else pbcv=0
  
      
      #method==1 : fullModel;
      if (method==1) coeffestime<-fullModel(family=gaussian('identity'),data=tabtravail)$coef
      #method==2 : stepSel ;
      if (method==2)coeffestime<-stepSel(family=gaussian('identity'),data=tabtravail,...)$coef
      #method==3 : bmaBic;
      if (method==3) coeffestime<-bmaBic(family=gaussian('identity'),data=tabtravail,...)$coef
      #method==4 : mixAic ;
      if (method==4) coeffestime<-mixAic(family=gaussian('identity'),data=tabtravail,...)$coef
      #method==5 : arms;
      if(method==5)coeffestime<-arms(data=tabtravail,family=gaussian('identity'),...)$coef
       
      # Explanatory variables of the data i with 1 for the intercept
      varexplicatives<-as.matrix(cbind(1,data[i,2:dim(data)[2]]))
      
      #Estimated values of the coefficients
      coefficient<-as.matrix(coeffestime)
      
      #Prediction of the data i
      ychapeau<-varexplicatives%*%coefficient
      
      #Tabindref includes the prediction and the observed value 
      Tabindref[i,]=c(ychapeau,data[i,1]) 
      
      if (length(file) !=0){
      #Predictions are stored in the file "file" , in case the program would 
      #stop
      write.table(Tabindref,file,sep="\t",row.names=FALSE,dec=".")
      }
  
  }
}

#else the PMSE is calculated on npermu samples of size np, for the model 
#estimated from the other part of the samples

#number of colinearity problems in the train sample
cons.pb=0
col.pb=0
# which occurs successively
pbcons=0
pbcol=0
#number of convergence, prediction separation
cv.error=0
sep.error=0
#number of convergence, prediction separation
# which occurs successively
pbcv=0
pbsep=0

if(np>1 | random==TRUE){
  #Tabindref includes the predictions derived by cross-validation
  Tabindref<-{}
     
  permu=1 
  while(permu <= npermu) {
       
      test.ind<-sample(1:dim(data)[1],np,replace=FALSE,prob=NULL)
      
      test.sample<-data[test.ind,]
      train.sample<-data[-test.ind,]
     
       #If some columns in the training sample take the same values another sample 
       #is drawn
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
      
      #Preliminary test of the convergence
      #If one of the model can not be estimated for the work.sample, the program 
      #switches to the next iteration
      fullResult<-fullModel(family=gaussian('identity'),data=train.sample)
        if (fullResult$cv == FALSE ){
          pbcv<-pbcv+1
          cv.error=cv.error+1
        }
      if (pbcv > 100) stop("100 samples created consecutively without one of them are estimable")
      if (fullResult$cv==FALSE) next else pbcv=0
  
      
      #method==1 : fullModel;
      if (method==1) coeffestime<-fullModel(data=train.sample,family=gaussian('identity'),...)$coef
      #method==2 : stepSel ;
      if (method==2)coeffestime<-stepSel(data=train.sample,family=gaussian('identity'),...)$coef
      #method==3 : bmaBic;
      if (method==3) coeffestime<-bmaBic(data=train.sample,family=gaussian('identity'),...)$coef
      #method==4 : mixAic ;
      if (method==4) coeffestime<-mixAic(data=train.sample,family=gaussian('identity'),...)$coef
      #method==5 : arms;
      if(method==5)coeffestime<-arms(data=train.sample,family=gaussian('identity'),...)$coef
       
      #Explanatory variables of the test.sample with 1 for the intercept
      varexplicatives<-as.matrix(cbind(rep(1,np),test.sample[,2:dim(data)[2]]))
      
      #Estimated values of the coefficients
      coefficient<-as.matrix(coeffestime)
      
      #Prediction of the data i
      ychapeau<-varexplicatives%*%coefficient
      
      #Tabindref includes the prediction and the observed value 
      indref=cbind(ychapeau,test.sample[,1])
      Tabindref=rbind(Tabindref,indref)
    
      if (length(file) != 0){
      #Predictions are stored in a file in case the program would stop
      write.table(indref,file,sep="\t",row.names=FALSE,dec=".",append=TRUE,col.names=FALSE)
      }

      permu=permu+1
  }
}

#Missing values in the predictions are deleted in order to calculate the PMSE 
colelim={}
for (i in c(1:nrow(Tabindref))){
  if ((is.na(Tabindref[i,2])|is.nan(Tabindref[i,2])) == TRUE) colelim=c(colelim,i)
}
if (length(colelim) >= 1) Tabindref<-Tabindref[-colelim,]  

#Squares of the residuals 
scr=(Tabindref[,2]-Tabindref[,1])^2

#Predictive mean square error derived by cross validation
pmsecv<-mean(scr)  

##MSE on the whole sample

#method==1 : fullModel;
if (method==1) coeffestime<-fullModel(family=gaussian('identity'),data=data)$coef
#method==2 : stepSel ;
if (method==2) coeffestime<-stepSel(data=data,family=gaussian('identity'),...)$coef
#method==3 : bmaBic;
if (method==3) coeffestime<-bmaBic(family=gaussian('identity'),data=data,...)$coef
#method==4: mixAic ;
if (method==4) coeffestime<-mixAic(family=gaussian('identity'),data=data,...)$coef
#method==5: arms;
if(method==5)  coeffestime<-arms(data=data,family=gaussian('identity'),...)$coef

coefficient<-as.matrix(coeffestime)
varexplicatives<-as.matrix(cbind(1,data[,2:dim(data)[2]]))
ychapeau<-varexplicatives%*%coefficient
scr=(data[,1]-ychapeau)^2
mse<-mean(scr)

if (cv.error!=0) cat("Number of convergence problems",":",cv.error,"\n") 
if (n.error!=0) cat("Number of errors",":",n.error,"\n")

Tabmse<-data.frame(MSE=mse,PMSE_CV=pmsecv)

return(Tabmse)

}
