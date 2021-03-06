% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelClass.R
\name{setMeasFunc}
\alias{setMeasFunc}
\alias{setMeasFunc,odeModel,function,missing-method}
\alias{setMeasFunc,odeModel,function,logical-method}
\title{Set the measurement equation for the model}
\usage{
setMeasFunc(odeModel, measFunc, custom)

\S4method{setMeasFunc}{odeModel,`function`,missing}(odeModel, measFunc, custom)

\S4method{setMeasFunc}{odeModel,`function`,logical}(odeModel, measFunc, custom)
}
\arguments{
\item{odeModel}{an object of the class odeModel}

\item{measFunc}{measurement function of the model. Has to be a R functions.}

\item{custom}{custom indexing for the measurement function (used by the baysian method)}
}
\value{
an object of odeModel
}
\description{
For a given model a measurement equation can be set. If no measurement function is set the
states become the output of the system. The function should be defined as in the example below.
}
\examples{

data("uvbModel")

uvbMeasure <- function(x) {

  y1 = 2*x[,5] + x[,4] + x[,8]
  y2 = 2*x[,5] + 2* x[,3] + x[,1]
  y3 = x[,6]
  y4 = x[,11]
  y5 = x[,4]

  return(cbind(y1,y2,y3,y4,y5))
  }

newModel <- setMeasFunc(odeModel = uvbModel, measFunc = uvbMeasure)

}
