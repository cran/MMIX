\name{bmaBic}
\alias{bmaBic}
\alias{mixAic}

\title{ Model mixing methods for linear or logistic regression}

\description{
 \code{bmaBic} is a Bayesian Model Averaging method based on the bic criterion, 
 and \code{mixAic} is a model mixing method based on the aic criterion.
}

\usage{
bmaBic(data, family, maxVar=10)
mixAic(data, family, maxVar=10)
}

\arguments{
  \item{data}{ a data frame including the response variable (first column) and 
  the explanatory variables. All the variables must be numeric and the response 
  variable value must be 0 or 1 for the logistic model.  }
   \item{family}{  a description of the error distribution (gaussian("identity")
   or binomial("logit")).}
  \item{maxVar}{maximum number of explanatory variables to include in the model. 
  See \code{\link{varSelec}}. }
}

\details{
A set of models including all possible combinations of at most \code{maxVar} 
variables is defined. The models are weighted using aic (mixAic) or bic 
(bmaBic).

Computation time may be very long if \code{maxVar} is set to a high value. 
When the number of explanatory variables in data is higher than \code{maxVar},
this number is reduced to at most \code{maxVar} variables by stepwise selection
("forward", "aic"). 

}

\value{
\code{bmaBic} and \code{mixAic} return an object of class "MMIXclass". A data 
frame with the main results is printed when calling the function "print", the 
detailed results are obtained with the function "summary" and a graphic of the 
weights of the explanatory variables is plotted with the function "plot".
A bmaBic or mixAic object is a list including the following components :
  
  \item{coef }{ a named vector of coefficients.}
  \item{pne0 }{ a vector containing the probability of each explanatory variable
   to be different from zero. }
  \item{fitted.values }{ the fitted values, obtained by transforming the 
  linear predictors by the inverse of the link function.}
  \item{BestModels }{the best three models selected by AIC (mixAic) or BIC 
  (bmaBic).}
  \item{sd }{ standard deviation of the estimated parameter values.}
  \item{label}{a list of the explanatory variables used by each model.}
  \item{modweights}{the weights of all models.}
  \item{allcoef}{matrix with one row per model and one column per coefficient
   giving the ML estimate of each coefficient for each model.}
}

\references{ 

Barbottin A., Makowski D., Le Bail M., Jeuffroy M.-H., Bouchard C., Barrier C. 
(2008) Comparison of models and indicators for categorizing soft wheat fields
according to their grain protein contents, \emph{European Journal of Agronomy}
\bold{29}, 159-183. 

Burnham, K.P. and Anderson, D.R. (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}, Springer-Verlag, New 
York.

Casagrande, M., David, C., Valantin-Morison, M., Makowski, D. and Jeuffroy, 
M.-H. (2009) Factors limiting protein content of organic winter wheat in south-
eastern France: a mixed-model approach, \emph{Agronomy for Sustainable 
Development, In press}. 

Hoeting, J.A., Madigan, D., Raftery, A.E. and Volinsky, C.T. (1999) Bayesian 
model averaging: a tutorial, \emph{Statistical Science} \bold{14}, 382-417.

Prost, L., Makowski, D. and Jeuffroy, M.-H. (2006) Comparison of stepwise
selection and Bayesian model averaging for yield gap analysis, \emph{Ecological
Modelling} \bold{219}, 66-76. 

Raftery, A.E., Madigan D. and Hoeting, J.A. (1997) Bayesian model averaging for
 linear regression models, \emph{Journal of the American Statistical
  Association} \bold{92}, 179-191.

Raftery, A.E., Painter, I.S. and Volinsky, C.T. (2005) BMA : An R package for 
Bayesian Model Averaging, \emph{R News} \bold{5}, 2-8.

Viallefont, V., Raftery, A.E. and Richardson, S. (2001) Variable selection and 
Bayesian model averaging in case-control studies, \emph{Statistics in medicine} 
\bold{20}, 3215-3230.

}

\author{ Marie Morfin and David Makowski }

\section{Warning }{ These functions do not accept or treat missing values.}


\seealso{ \code{\link{family}}, \code{\link{varSelec}} }

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

##bmaBic method
bmaResult1<-bmaBic(data=data1,family=gaussian("identity"))
bmaResult1
summary(bmaResult1)
plot(bmaResult1)

bmaResult2<-bmaBic(data=data2,family=binomial("logit"))
bmaResult2
summary(bmaResult2)
plot(bmaResult2)

##mixAic method
mixResult1<-mixAic(data=data1,family=gaussian("identity"))
mixResult1
summary(mixResult1)
plot(mixResult1)

mixResult2<-mixAic(data=data2,family=binomial("logit"))
mixResult2
summary(mixResult2)
plot(mixResult2)

}

\keyword{ regression }
\keyword{ robust }
