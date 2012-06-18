
bootFreq=function(data, family, nboot = 100, method = 1, file =NULL, ...){
#data is a data frame including the response variable (first column) and the 
#explanatory variables

#Variables names
nom = names(data)
nomcoef=nom
nomcoef[1]="Intercept"
#nobs is the number of data 
nobs=dim(data)[1]

#Tabcoef includes the coefficient values estimated from each bootstrap sample
Tabcoef=matrix(0,nrow=ncol(data),ncol=nboot)
#TabPne0 includes the probability of each coefficient to be non zero, for each 
#bootstrap sample
TabPne0<-matrix(0,nrow=dim(data)[2],ncol=nboot)

#number of colinearity problems with the training sample
cons.pb=0
col.pb=0
# which occurs successively 
pbcons=0
pbcol=0 
#Number of convergence and prediction problems
cv.error=0
sep.error=0
#Number of successive convergence and separation problems 
pbcv=0
pbsep=0

#Loop for all the bootstrap samples
iteration=1
while(iteration <= nboot) {            
    
    #Random sampling with replacement of nobs individuals in the data set 
    individustire <- sample(1:nobs,nobs,replace=TRUE,prob=NULL)
    tabtravail<-data[individustire,]
    
    #If some columns in the training sample are constant another sample is drawn
    colelim={} 
    for(j in c(2:dim(tabtravail)[2])){
      if (max(tabtravail[,j]-min(tabtravail[,j])) == 0) colelim=c(colelim,j)   
    }    
    if (length(colelim)!=0) {
      pbcons=pbcons+1
      cons.pb=cons.pb+1
      if(pbcons>100) stop("one factor may be constant")
      next 
    }else pbcons=0
     
    #If there are more than two explanatory variables, and two of them are 
    #proportionnal, another sample is drawn
     if (dim(tabtravail)[2] >= 3){
      for(j in c(2:((dim(tabtravail)[2])-1))){
        restcol=(j+1):(dim(tabtravail)[2])
        if (length(which(cor(tabtravail[,j],tabtravail[,restcol])==1))!=0){
          colelim=c(colelim,j)
        }
      }
    
      if (length(colelim)!=0){
        pbcol=pbcol+1
        col.pb=col.pb+1
        if(pbcol>100) stop("problem of colinearity between factors")
        next 
      }else pbcol=0
    }
    
    #Convergence test
    #If the full model is not estimable with the sample, another sample is 
    #drawn
    coeffestime<-fullModel(family=family,data=tabtravail)
    if (coeffestime$cv == FALSE){
      pbcv<-pbcv+1
      cv.error=cv.error+1
    }
    if (cv.error > 100) stop("100 samples created consecutively without one of them are estimable")
    if (coeffestime$cv == FALSE) next else cv.error=0
    
    #In the logistic case, data set may create a separation problem
    if (family$family=="binomial") {
      #full model predictions
      yhat<- fullModel(data=tabtravail,family=family) 
      test.sep<-yhat$fitted.values[yhat$fitted.values<= 0.0001 |yhat$fitted.values >=0.9999]
      if (length(test.sep)!=0){ 
        pbsep=pbsep+1
        sep.error=sep.error+1
      }
      if (pbsep > 100) stop("100 samples created consecutively without one of them are estimable")
      if (length(test.sep)!=0) next else pbsep=0
    }
    
    #method==1 : fullModel;
    if (method==1) Result<-fullModel(data=tabtravail,family=family)   
    #method==2 : stepSel;
    if (method==2) Result<-stepSel(data=tabtravail,family=family,...)
     #method==3 : bmaBic;
    if (method==3) Result<-bmaBic(data=tabtravail,family=family, ...)
    #method==4 : mixAic;
    if (method==4) Result<-mixAic(data=tabtravail,family=family,...)
    #method==5 : arms;
    if(method==5)  Result<-arms(data=tabtravail,,family=family,...)
    
    Tabcoef[,iteration]<-Result$coef

    #Variable frequencies of selection
    Tabmodele<-Tabcoef
    #Tabmodele includes 1 where the Tabcoef is different from zero
    Tabmodele[which(Tabcoef!=0)]<- 1
    TabFreq<- rowSums(Tabmodele)/iteration
 
    #Pne0
    if (method==3|method==4|method==5){
    TabPne0[,iteration]=Result$pne0
    MeanPne0=matrix(0,ncol=1,nrow=length(nomcoef))
    MeanPne0=rowSums(TabPne0)/iteration}
   if(method==1|method==2){
   MeanPne0=TabFreq}
  
  Tabinter={}
  sdT<-apply(as.data.frame(t(Tabcoef[,1:iteration])),2,"sd")
  Tabinter<-data.frame(names=nomcoef,coef=Tabcoef[,1:iteration],
  mean=rowSums(Tabcoef)/iteration,sd=sdT,
  frequency=TabFreq,Pne0=MeanPne0)

  if(length(file) != 0){
  #Results are stored in a file during the run
  write.table(Tabinter,file,sep="\t",row.names=FALSE,dec=".")
  }

  iteration=iteration +1
}

if (cv.error!=0) cat("Convergence problems",":",cv.error,"\n")
if (sep.error!=0) cat("Separation problems",":",sep.error,"\n")
 
#This function returns a data frame with the names of the parameters, the 
#estimations obtained for each sample, the mean of all the estimations; their
#standard deviation, the variable frequencies of selection 
Tabcoef<-as.data.frame(t(Tabcoef))
names(Tabcoef)<-nomcoef
mean<-as.vector(Tabinter$mean)
names(mean)<-nomcoef

TabResult<-data.frame(mean=mean,sd=Tabinter$sd,frequency=Tabinter$frequency,
pne0=Tabinter$Pne0)

ListSummary<-list(coef=Tabcoef,mean=TabResult$mean,sd=TabResult$sd,
frequency=TabResult$frequency,pne0=TabResult$pne0)

VecPlot<-t(as.matrix(TabResult$pne0)[-1])
names(VecPlot)<-nomcoef[-1]
rownames(VecPlot)<-"Factors weights"

ListResult<-list(TabResult,ListSummary,VecPlot,coef=Tabcoef,mean=TabResult$mean,
sd=TabResult$sd,frequency=TabResult$frequency,pne0=TabResult$pne0)

#The first element of the list "ListResult" is visible when the output of 
#the function is printed, the second element is printed when calling the
# function "summary"
class(ListResult)<-"MMIXclass"

return(ListResult)
}








