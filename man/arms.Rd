\name{arms}
\alias{arms}

\title{ Adaptative Regression by Mixing with model Screening}

\description{
 Apply ARMS for linear and logistic models.
}
\usage{
arms(data,  family, nsample = 20, nbest = 20, criterion= "both",
     weight = "aic", maxVar=10)
}

\arguments{
  \item{data}{ a data frame including the response variable (first column) and 
  the explanatory variables. All the variables must be numeric and the response 
  variable value must be 0 or 1 for the logistic model.  }
   \item{family}{  a description of the error distribution (gaussian("identity")
   or binomial("logit")).}
  \item{nsample}{ number of sample permutations}
  \item{nbest}{ number of models selected at the screening step }
  \item{criterion}{ selection criterion used at the screening step. 
  criterion="aic", "bic" or "both". 
  "both" means that \code{nbest} models are selected according to AIC, and
   \code{nbest} according to "bic" (between \code{nbest} and 2*\code{nbest} 
   models may be selected). }
  \item{weight}{ model weight type, \code{weight} = "likeli" or "aic" 
  (the default). See \code{details}.  }
  \item{maxVar}{maximum number of explanatory variables to include in the model. 
  See \code{\link{varSelec}}. }
}

\details{
A set of models including all possible combinations of at most \code{maxVar}
variables is defined. The models are weighted as follows:

Step 1. The sample is splitted in two parts. The first one is used to fit the 
models.
 
Step 2. Each model is fitted by least squares or maximum likelihood. 

Step 3. A set of models is selected according to AIC and BIC.

Step 4. Weights are computed for the selected models using the second 
sample. Two weight types can be used, "likeli" or "aic".

At each iteration, the sample is permuted. The final weight of each model is the 
mean across all the permutations.
The coefficients are calculated on the whole sample using the averaged model 
weights.
}

\value{
  \code{arms} returns an object of class "MMIXclass". A data frame with the main
   results is printed with the function "print" and a graphic with the 
   weight of the explanation variables is obtained with the function "plot". An 
   \code{arms} object is a list including the following components:
  
  \item{coef }{ a named vector of coefficients.}
   \item{pne0 }{ a vector containing the probability of each 
  explanatory variable to be different from zero. }
   \item{fitted.values }{ the fitted values, obtained by transforming the 
  linear predictors by the inverse of the link function.}
  \item{label}{a list of the explanatory variables used by each model.}
  \item{modweights}{the weights of all models.}
  \item{allcoef}{matrix with one row per model and one column per coefficient
   giving the ML estimate of each coefficient for each model.}
}

\references{ 
Yuan, Zheng and Ghosh, Debashis (2008) Combining Multiple Biomarkers Models in 
Logistic Regression, \emph{Biometrics} \bold{64}, 431-439.


Yuan, Zheng and Yang, Yuhong (2005) Combining Linear Regression Models : When 
and How?, \emph{Journal of the American Statistical Association} \bold{100}, 
1202-1214.  

}

\author{ Marie Morfin and David Makowski }

\section{Warning }{This function does not accept or treat missing values.} 

\seealso{  \code{\link{fullModel}}, \code{\link{bmaBic}}, \code{\link{mixAic}},
 \code{\link{stepSel}},  \code{\link{bootFreq}}}

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

#arms method
armsResult1<-arms(data=data1,family=gaussian('identity'),nbest=5,nsample=10,
criterion="both",weight="aic")
armsResult1
summary(armsResult1)
plot(armsResult1)

armsResult2<-arms(data=data2,family=binomial("logit"),nbest=5,nsample=10,
criterion="both",weight="aic")
armsResult2
summary(armsResult2)
plot(armsResult2)

}

\keyword{ regression }
\keyword{ robust }
\keyword{ iteration }