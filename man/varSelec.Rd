\name{varSelec}
\alias{varSelec}


\title{Variable selection}


\description{
Create a data frame from the original data set including a subset of variables.
}


\usage{
varSelec(data,family,maxVar=10,trace=0)
}



\arguments{
 \item{data }{ a data frame including the response variable (first column)
  and the explanatory variables. All the variables must be numeric and the 
  response variable value must be 0 or 1 for the logistic model.}
  \item{family }{ a description of the error distribution (gaussian("identity")
   or binomial("logit")).}
  \item{maxVar}{maximum number of explanatory variables to include.}
  \item{trace}{print information during the run if \code{trace} = 1. Larger
   values may give more information. If \code{trace} = 0 no information is 
   printed. }
}

\details{
This function implements a stepwise regression for linear and logistic models, 
in the direction "forward" and with the criterion "aic". The procedure stops if 
the model includes more than \code{maxVar} factors. The returned data frame 
includes no more than \code{maxVar} explanatory variables.
}

\value{
  varSelec returns a data frame including the response variable and the
  selected explanatory variables.

}


\author{ Marie Morfin and David Makowski }


\seealso{  \code{\link{bmaBic}},  \code{\link{mixAic}}, \code{\link{arms}} }

\examples{
##Data 
#Explanatory variables 
X1<-c(-0.2,-2.4,-0.7,1.2,0.0,-1.1,-2.1,-0.3,2.0,-1.7,1.4,-1.3,-3.4,0.4,-1.3,
-4.8)
X2<- c(-3,  2,  1, -2, -2, -4,  0,  1,  1, -1, -1, -4,  0,  2,  0, -4)
X3<-c(2,1,0,-2,1,-2, 0, -1, -4, 1, -3, -3, -3, -1, 0, 2)

#Linear model
Y1<- c(8.7, 6, 9.1, 10.4, 7.6 ,10.4,  7.9, 11.9, 18, 10.5, 16.5, 8.8, 7.7,
 13.5, 8.2, 0.8)
data1<-data.frame(Y1,X1,X2,X3)
#Logistic model
Y2<-c(1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1)
data2<-data.frame(Y2,X1,X2,X3)

##varSelec
data1bis<-varSelec(data=data1,family=gaussian("identity"),maxVar=2)
data2bis<-varSelec(data=data2,family=binomial("logit"),maxVar=2)

}

\keyword{ regression }

