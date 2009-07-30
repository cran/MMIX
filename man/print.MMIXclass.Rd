\name{print.MMIXclass}

\alias{print.MMIXclass}
\alias{summary.MMIXclass}
\alias{plot.MMIXclass}
\alias{MMIXclass}

\title{Summary, print and plot of MMIXclass objects}
\description{
     \code{summary}, \code{print} and \code{plot} methods for MMIXclass objects.
}

\usage{
\method{print}{MMIXclass}(x,...)
\method{summary}{MMIXclass}(object,...)
\method{plot}{MMIXclass}(x,...) 


}

\arguments{
 \item{x}{object of type 'MMIXclass'.}
  \item{object}{object of type 'MMIXclass'.}
  \item{...}{other parameters to be passed to \code{plot.MMIXclass}, 
  \code{print.MMIXclass} or \code{summary.MMIXclass} }
}

\details{
\code{print} displays a data.frame with the main results of the function. 
\code{summary} displays a detailed list of the results. \code{plot} displays 
explanatory variable weights, i.e. the probability of each variable to be 
different from zero.  
}

\author{Marie Morfin and David Makowski \email{makowski@grignon.inra.fr}}


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

##bmaBic method
bmaResult1<-bmaBic(data=data1,family=gaussian("identity"))
bmaResult1
summary(bmaResult1)
plot(bmaResult1)
}

\keyword{print}
