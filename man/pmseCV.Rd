\name{pmseCV}
\alias{pmseCV}
\alias{aucCV}

\title{ Model performance indicators PMSE and AUC}

\description{
\code{pmseCV} and \code{aucCV} calculate the Predictive Mean Square Error (PMSE)
by "leave-np-out" cross-validation and Area Under Roc Curve (AUC) by
"leave-np-pair-out" cross-validation. They can be applied to models created 
using \code{\link{fullModel}}, \code{\link{stepSel}}, \code{\link{bmaBic}}, 
\code{\link{mixAic}} and \code{\link{arms}}.   
}

\usage{
pmseCV(data, method = 1, np, random = TRUE, npermu = 100, file = NULL, ...)

aucCV(data, method = 1, np, random = TRUE, npermu = 100, file = NULL, ...)
}                                                     
 
\arguments{
  \item{data}{ a data frame including the response variable (first column) and 
  the explanatory variables. All the variables must be numeric and the response 
  variable value must be 0 or 1 for the logistic model. }
  \item{method}{ the statistical method used to estimate the model parameters. 
     \code{method} = 1 for \code{\link{fullModel}}, \code{method} = 2 for 
  \code{\link{stepSel}}, \code{method} = 3 for \code{\link{bmaBic}}, 
  \code{method} = 4 for \code{\link{mixAic}}, \code{method} = 5 for 
  \code{\link{arms}}. }
  \item{np}{ number of observations (\code{pmseCV}) or pairs of observations 
  (\code{aucCV}) left out for computing the PMSE or AUC. }
  \item{random}{ observations are selected at random if TRUE. \code{random} can
  be FALSE only if \code{np} = 1. In this case all the possible sets are 
  selected. }
  \item{npermu}{ number of random samples of \code{np} observations if 
  \code{random} = TRUE. }
  \item{file}{ the path of the file where the results are stored during the run.
                If \code{file} = NULL no file is created. }
  \item{\dots}{ the specific arguments of the called method }
}

  
\details{

These two cross validation procedures are implemented to assess the accuracy of 
model predictions. Linear models should be evaluated using \code{pmseCV} and 
logistic models should be evaluated using \code{aucCV}.
In the aucCV procedure, data are held out by pair (one data from each class, 0 
and 1) at each iteration in order to calculate the corresponding AUC. The PMSE 
(AUC) estimated by cross validation is the mean of the PMSE (AUC) calculated 
for possible sets of \code{np} observations (pairs).
If \code{file} is not NULL, the \code{np} predictions (first column) and the
 \code{np} corresponding observations (second column) are saved at each 
 iteration. 


}

\value{
  pmseCV (aucCV) returns a one-row data frame including the PMSE (AUC) 
  calculated from the whole sample in the first column and by cross validation 
  in the other one. 

}

\references{
Sing, T., Sander, O., Beerenwinkel, N. and Lengauer, T. (2005) ROCR: 
visualizing classifier performance, \emph{Bioinformatics applications note} 
\bold{21}, 3940-3941.

Yuan, Z. and Ghosh, D. (2008) Combining Multiple Biomarker Models in Logistic 
Regression, \emph{Biometrics} \bold{64}, 431-439.

Yuan, Z. and Yang, Y. (2005) Combining Linear Regression Models: When and How?,
 \emph{Journal of the American Statistical Association} \bold{100}, 1202-1214.
}

\author{ Marie Morfin and David Makowski }


\section{Warning }{These functions do not accept or treat the missing values.} 

\seealso{ \code{\link{fullModel}}, \code{\link{stepSel}}, \code{\link{bmaBic}},
 \code{\link{mixAic}},
\code{\link{arms}} }

\examples{
##Data 
#Explanatory variables 
X1<-c(-0.2,-2.4,-0.7,1.2,0.0,-1.1,-2.1,-0.3,2.0,-1.7,1.4,-1.3,-3.4,0.4,-1.3,
-4.8)
X2<- c(-3,  2,  1, -2, -2, -4,  0,  1,  1, -1, -1, -4,  0,  2,  0, -4)
X3<-c(2,1,0,-2,1,-2, 0, -1, -4, 1, -3, -3, -3, -1, 0, 2)

#Linear model
Y1<- c(8.7, 6, 9.1, 10.4, 7.6 ,10.4,  7.9, 11.9, 18, 10.5, 16.5, 8.8, 7.7, 13.5,
 8.2, 0.8)
data1<-data.frame(Y1,X1,X2,X3)
#Logistic model
Y2<-c(1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1)
data2<-data.frame(Y2,X1,X2,X3)

##Linear models evaluated by pmseCV
#for a stepwise selection
pmseStepBic<-pmseCV(data=data1,method=2,np=1,random=FALSE,direction="both",
criterion="bic",trace=0)
pmseStepBic
#for the BMA method
pmseBMA<-pmseCV(data=data1,method=3,np=1,random=FALSE)
pmseBMA

##Logistic models evaluated by aucCV
#for a stepwise selection
aucStepBic<-aucCV(data=data2,method=2,np=1,random=FALSE,direction="both",
criterion="bic",trace=0)
aucStepBic
#for the BMA method
aucBMA<-aucCV(data=data2,method=3,np=1,random=FALSE)
aucBMA


}

\keyword{ regression}
\keyword{ robust }
\keyword{ iteration }
