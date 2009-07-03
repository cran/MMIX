\name{fullModel}
\alias{fullModel}

\title{ Linear and logistic regressions}
\description{
Linear or logistic model without any variable selection.
}

\usage{
fullModel(data, family)
}

\arguments{
  \item{data}{ a data frame including the response variable (first column) and 
  the explanatory variables. All the variables must be numeric and the response 
  variable value must be 0 or 1 for the logistic model.}
  \item{family}{ a description of the error distribution (gaussian("identity")
   or binomial("logit")).}
}


\value{
  \code{fullModel} returns an object of class "MMIXclass". A data frame with the
   main results is printed when calling the function "print", and the detailed 
   results are obtained with the function "summary". A \code{fullModel} object 
   is a list including the following components:
  
  \item{cv }{logical. Was the IWLS algorithm judged to have converged?}
  \item{coef }{a named vector of coefficients estimated by least squares or
  maximum likelihood.}
  \item{aic }{Akaike Information Criterion, minus twice the maximized 
  log-likelihood plus twice the number of coefficients.}
  \item{bic }{Bayesian information criterion, minus twice the maximized 
  log-likelihood plus the logarithm of the number of observation multiplied by 
  the number of coefficients}
  \item{fitted.values }{the fitted values, obtained by transforming the linear 
  predictors by the inverse of the link function.}
}

\references{ 
Venables, W.N. and Ripley, B.D. (2002) \emph{Modern Applied Statistics with S.}, 
Springer, New York.  }

\author{ Marie Morfin and David Makowski }


\seealso{  \code{\link{family}}, \code{\link{arms}}, \code{\link{bmaBic}}, 
\code{\link{mixAic}}, \code{\link{stepSel}},  \code{\link{bootFreq}}   }

\examples{

###Data 
#Explanatory variables 
X1<-c(-0.2,-2.4,-0.7,1.2,0.0,-1.1,-2.1,-0.3,2.0,-1.7,1.4,-1.3,-3.4,0.4,-1.3,
-4.8)
X2<- c(-3,  2,  1, -2, -2, -4,  0,  1,  1, -1, -1, -4,  0,  2,  0, -4)
X3<-c(2,1,0,-2,1,-2, 0, -1, -4, 1, -3, -3, -3, -1, 0, 2)

##Linear model
Y1<- c(8.7, 6, 9.1, 10.4, 7.6 ,10.4,  7.9, 11.9, 18, 10.5, 16.5, 8.8, 7.7, 13.5,
 8.2, 0.8)
data1<-data.frame(Y1,X1,X2,X3)
##Logistic model
Y2<-c(1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1)
data2<-data.frame(Y2,X1,X2,X3)

##fullModel method
fullResult1<-fullModel(data=data1,family=gaussian("identity"))
fullResult1
summary(fullResult1)

fullResult2<-fullModel(data=data2,family=binomial("logit"))
fullResult2
summary(fullResult2)

}


\keyword{models}
\keyword{ regression }
