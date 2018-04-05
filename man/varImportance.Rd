\name{getVarImp}


\alias{getVarImp}
\alias{getVarImp,sdmModels-method}


\title{variable importance}


\description{
Calculates relative importance of different variables in the models using several approaches.
}

\usage{
getVarImp(x,id,wtest,...)

}

\arguments{
  \item{x}{a numeric vector including the observed values; or a \code{sdmEvaluate} object}
  \item{id}{a numeric vector including the predicted values}
  \item{wtest}{specifies which dataset ('training','test.dep','test.indep') should be used (if exist) to calculate the importance of variables }
  \item{...}{additional arguments (not implemented yet.)}
  }

\details{

\code{getVarImp} function returns an object including different measures of variable importance, and if be put in plot function, a barplot is generated:

}



\references{

Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}



\examples{
\dontrun{
# if m is a sdmModels object (output of sdm function) then:

getVarImp(m,id=1,wtest='training') # variable importance based on training dataset

vi <- getVarImp(m,id=1,wtest='test.dep') 

plot(vi,'auc')

plot(vi,'cor')

}
}

\keyword{spatial,data}