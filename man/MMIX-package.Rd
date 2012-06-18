\name{MMIX-package}
\alias{MMIX-package}
\alias{MMIX}
\docType{package}

\title{Model selection uncertainty and model mixing}

\description{
Implement different types of model mixing and model selection methods for 
linear or logistic models.
}

\details{
\tabular{ll}{
Package: \tab MMIX\cr
Type: \tab Package\cr
Version: \tab 2.15.0.\cr
Date: \tab 2012-06-18\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package implements three model mixing methods, two model selection
procedures and  three functions for model evaluation. 
  
}

\author{
Marie Morfin and David Makowski

Maintainer: Who to complain to <makowski@grignon.inra.fr>
}

\references{
Barbottin A., Makowski D., Le Bail M., Jeuffroy M.-H., Bouchard C., Barrier C. 
(2008) Comparison of models and indicators for categorizing soft wheat fields
according to their grain protein contents, \emph{European Journal of Agronomy}
\bold{29}, 159-183. 

Buckland, S.T., Burnham, K.P. and Augustin, N.H. (1997) Model selection:
 an integral part of inference, \emph{Biometrics} \bold{53}, 603-618.

Burnham, K.P. and Anderson, D.R. (2002) \emph{Model selection and multimodel 
inference: a practical information-theoretic approach}, Springer-Verlag, 
New York.

Casagrande, M., David, C., Valantin-Morison, M., Makowski, D. and Jeuffroy, 
M.-H. (2009) Factors limiting protein content of organic winter wheat in south-
eastern France: a mixed-model approach, \emph{Agronomy for Sustainable 
Development, In press}. 

Chatfield, C. (1995) Model uncertainty, data mining and statistical 
inference, \emph{Journal of the Royal Statistical Society /A} \bold{158}, 
419-466.

Efron, B. (1979) Bootstrap methods : another look at the jackknife,
 \emph{American Statistician} \bold{7}, 1-26.

Efron, B. and Tibshirani, R.J. (1993) \emph{An introduction to the bootstrap},
 Chapman & Hall.

Hammersley, J.M. and Handscomb, D.C. (1964) \emph{Monte Carlo Methods},
 Chapman & Hall, London & New York.

Hoeting, JA., Madigan, D., Raftery, A.E. and Volinsky, C.T. (1999),
 \emph{Bayesian model averaging: a tutorial}, Statistical Science, 14, 382-417.

Miller A. (2002) \emph{Subset selection in regression}, 2nd edition Chapman & 
Hall/CRC, New York. 

Mooney, C.Z. and Duval, R.D. (1993) \emph{Bootstrapping : a nonparametric 
approach to statistical inference}, Sage Publications, London.

Prost, L., Makowski, D. and Jeuffroy, M.-H. (2006) Comparison of stepwise
 selection and Bayesian model averaging for yield gap analysis, \emph{Ecological 
 Modelling} \bold{219}, 66-76. 

Raftery, A.E., Madigan D. and Hoeting J.A. (1997) Bayesian model 
averaging for linear regression models, \emph{Journal of the American 
Statistical Association} \bold{92}, 179-191.

Sing, T., Sander, O., Beerenwinkel, N. and Lengauer, T. (2005) ROCR: 
visualizing classifier performance, \emph{Bioinformatics applications note}
 \bold{21}, 3940-3941.

Viallefont, V., Raftery A.E. and Richardson, S. (2001) Variable selection and 
Bayesian model averaging in case-control studies, \emph{Statistics in medicine} 
\bold{20}, 3215-3230.

Yuan, Z. and Ghosh, D. (2008) Combining Multiple Biomarker Models in Logistic 
Regression, \emph{Biometrics} \bold{64}, 431-439.

Yuan, Z. and Yang, Y. (2005) Combining Linear Regression Models: When and How?, 
\emph{Journal of the American Statistical Association} \bold{100}, 1202-1214.


}


\keyword{ package }


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

##bootFreq
#Stepwise selection method
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

##ARMS method
armsResult1<-arms(data=data1,family=gaussian("identity"),nbest=5,nsample=10,
criterion="both",weight="aic")
armsResult1
summary(armsResult1)
plot(armsResult1)

armsResult2<-arms(data=data2,family=binomial("logit"),nbest=5,nsample=10,
criterion="both",weight="aic")
armsResult2
summary(armsResult2)
plot(armsResult2)


##Root Mean Square Error by cross-validation
#Stepwise selection with BIC
pmseStepBic<-pmseCV(data=data1,method=2,np=1,random=FALSE,direction="both",
criterion="bic",trace=0)
pmseStepBic
#BMA
pmseBMA<-pmseCV(data=data1,method=3,np=1,random=FALSE)
pmseBMA

##Area Under ROC Curve by cross-validation
aucStepBic<-aucCV(data=data2,method=2,np=1,random=FALSE,direction="both",
criterion="bic",trace=0)
aucStepBic

aucBMA<-aucCV(data=data2,method=3,np=1,random=FALSE)
aucBMA


}
