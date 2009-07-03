\name{bootFreq}
\alias{bootFreq}

\title{ Stability of selection and mixing methods }

\description{
This function analyses the stability of the stepwise selection and mixing 
methods using a bootstrap procedure.

}

\usage{
bootFreq(data, family, nboot = 100, method = 1,  file = NULL, ...)
}

\arguments{
  \item{data}{ a data frame including the response variable (first column) and 
  the explanatory variables. All the variables must be numeric and the response 
  variable value must be 0 or 1 for the logistic model. }
  \item{family}{ a description of the error distribution (gaussian("identity")
   or binomial("logit")). }
  \item{nboot}{ number of bootstrap samples drawn. }
  \item{method}{ the statistical method used to estimate the model parameters.  
  \code{method} = 1 for \code{\link{fullModel}}, \code{method} = 2 for 
  \code{\link{stepSel}}, \code{method} = 3 for \code{\link{bmaBic}}, 
  \code{method} = 4 for \code{\link{mixAic}}, \code{method} = 5 for 
  \code{\link{arms}}. }
  \item{file}{ the path of the file where the results are stored during the run.
   If \code{file} = NULL no file is created. }
  \item{...}{ the specific arguments of the called method.}
}

\details{
\code{nboot} samples are generated with replacement from the original dataset.
For each sample, the parameters are estimated using the chosen method. 
The frequency of selection of a variable is the part of the samples for which 
the estimated value of its coefficient is different from zero. Selection 
frequency is an indicator of the stability of the method. Frequencies close to 
0 or 1 indicates that the corresponding variables was almost never or always 
selected in the bootstrap samples. The bootstrap estimators and their standard
deviations are calculated, as well as the mean and standard deviation of the
estimates across all the bootstrap samples.


}
\value{
bootFreq returns an object of class "classMMIX". A data frame with the main 
results is printed whith the function "print", and a graphic with the weight of 
the explanation variables is plotted with the function "plot". A \code{bootFreq} 
object is a list including the following components:
  
  \item{frequency }{ frequency of selection of each variables across the
   bootstrap samples}
  \item{coef }{ estimated parameter values across the bootstrap samples;
   a matrix with \code{nsample} rows}
  \item{mean }{ mean of the estimated parameter values across the bootstrap 
  samples}
  \item{sd }{ standard deviation of the estimated parameter values across the
   bootstrap samples}
  \item{pne0 }{mean values of the variable weights computed using a model
   mixing method. \code{Pne0} = \code{frequency} if \code{method} = 1 or 2.}
}

\references{ 
Buckland, S.T., Burnham, K.P. and Augustin, N.H. (1997) Model selection: 
an integral part of inference, \emph{Biometrics} \bold{53}, 603-618.

Chatfield, C. (1995) Model uncertainty, data mining and statistical inference, 
\emph{Journal of the Royal Statistical Society /A} \bold{158}, 419-466.

Efron, B. (1979) Bootstrap methods : another look at the jackknife,
  \emph{American Statistician} \bold{7}, 1-26.

Efron, B. and Tibshirani, R.J. (1993) \emph{An introduction to the bootstrap},
 Chapman & Hall.

Hammersley, J.M. and Handscomb, D.C. (1964) \emph{Monte Carlo Methods},
 Chapman & Hall.

Miller A. (2002) \emph{Subset selection in regression}, 2nd edition Chapman & 
Hall/CRC, New York.

Mooney, C.Z. and Duval, R.D. (1993) \emph{Bootstrapping: a nonparametric 
approach to statistical inference}, Sage Publications, London.

Prost, L., Makowski, D. and Jeuffroy, M.-H. (2006) Comparison of stepwise
selection and Bayesian model averaging for yield gap analysis, \emph{Ecological
 Modelling} \bold{219}, 66-76. 
}

\author{ Marie Morfin and David Makowski \email{makowski@grignon.inra.fr}}

\seealso{ \code{\link{fullModel}}, \code{\link{stepSel}}, \code{\link{bmaBic}},
 \code{\link{mixAic}}, \code{\link{arms}}}

\section{Warning }{This function does not accept or treat missing values.} 
 
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

##Analysis of a stepwise selection
bootStep1<-bootFreq(data=data1,family=gaussian("identity"),nboot=50,method=2,
criterion="bic",trace=0)
bootStep1
summary(bootStep1)
plot(bootStep1)

bootStep2<-bootFreq(data=data2,family=binomial("logit"),method=2,
criterion="bic",nboot=20,trace=0)
bootStep2
summary(bootStep2)
plot(bootStep2)

}


\keyword{ regression}
\keyword{ iteration }
