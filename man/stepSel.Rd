\name{stepSel}
\alias{stepSel}

\title{Stepwise selection}


\description{
 Select and fit a model by stepwise regression, for linear and logistic
  models.
}


\usage{
stepSel(data, family, direction = "both", criterion, trace = 0)
}


\arguments{
 \item{data }{ a data frame including the response variable (first column)
  and the explanatory variables. All the variables must be numeric and the 
  response variable value must be 0 or 1 for the logistic model.}
  \item{family }{ a description of the error distribution (gaussian(\dQuote{identity})
   or binomial(\dQuote{logit})). }
  \item{direction }{the type of stepwise search, can be one of \dQuote{both}
   (the default), \dQuote{backward} or \dQuote{forward}}
   \item{criterion}{ selection criterion at each step of the procedure.
    \code{criterion} = \dQuote{aic} for the Akaike Information Criterion, and 
    \code{criterion} = \dQuote{bic} for the Bayesian Information Criterion.  }
   \item{trace}{print information during the run if \code{trace} = 1. Larger
   values may give more information. If \code{trace} = 0 no information is 
   printed. }
}

\details{
  This function uses the function \code{\link{step}}.
}

\value{
  \code{stepSel} returns an object of class \code{"MMIXclass"}. A data frame 
  with the main results is printed whith the function \code{"print"}, and the 
  detailed results are obtained with the function \code{"summary"}. A 
  \code{stepSel} object is a list including the following components:
  
  \item{coef }{a named vector of coefficients estimated by least squares or
  maximum likelihood.}
  \item{aic }{Akaike Information Criterion, minus twice the maximized 
  log-likelihood plus twice the number of coefficients.}
  \item{bic }{Bayesian information criterion, minus twice the maximized 
  log-likelihood plus the logartithm of the number of observation multiplied by 
  the number of coefficients}
  \item{fitted.values }{the fitted values, obtained by transforming the linear 
  predictors by the inverse of the link function.}
}

\references{ 
Akaike H. (1974) A new look at the statistical model identification, \emph{IEEE
Transactions on Automatic Control} \bold{19}, 716-723.

Miller A. (2002) \emph{Subset selection in regression}, 2nd edition Chapman & 
Hall/CRC, New York.

Schwarz, G. (1978) Estimating the dimension of a model,  \emph{Annals of 
Statistics} \bold{6}, 461-464. 

Whittingham M.J., Stephens P., Bradbury R.B.. Freckleton R.P. (2006) Why do we 
still use stepwise modelling in ecology and behaviour?, \emph{J. Anim. Ecol.}
 \bold{75}, 1182-1189.
}

\author{ Marie Morfin and David Makowski }


\seealso{  \code{\link{family}},  \code{\link{step}}}

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

##stepSel method
stepResult1<-stepSel(data=data1,family=gaussian("identity"),criterion="bic",
direction="both")
stepResult1
summary(stepResult1)

stepResult2<-stepSel(data=data2,family=binomial("logit"),criterion="bic",
direction="both")
stepResult2
summary(stepResult2)

}


\keyword{ regression }
