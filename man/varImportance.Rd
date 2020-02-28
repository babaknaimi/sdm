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
  \item{x}{sdmModels object}
  \item{id}{numeric, specify the model (modelID) for which the variable importance values are extracted}
  \item{wtest}{specifies which dataset ('training','test.dep','test.indep') should be used (if exist) to calculate the importance of variables }
  \item{...}{additional arguments as for \code{getModelId} function}
  }

\details{

\code{getVarImp} function returns an object including different measures of variable importance, and if be put in plot function, a barplot is generated. If the ggplot2 package is installed on your machine, the plot is generated using ggplot (unless you turn gg = FALSE), otherwise, the standard barplot is used.

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
# if m is a sdmModels object (output of the sdm function) then:

getVarImp(m,id=1,wtest='training') # variable importance based on training dataset

vi <- getVarImp(m,id=1,wtest='test.dep') 

vi

plot(vi,'auc')

plot(vi,'cor')
#############
# You can get Mean variable importance (and confidence interval) for multiple models:

vi <- getVarImp(m,id=1:10,wtest='test.dep') # specify the modelIDs of the models

vi

plot(vi,'cor')

# you can use the getModelId function to find the id of the specific method, replication, etc.
# or you may put the arguments of the getModelId in the getVarImp function:

vi <- getVarImp(m, method='glm') # Mean variable importance for the method glm

vi

plot(vi)

plot(vi, gg = F) # R standard plot is used instead of ggplot


}
}

\keyword{spatial}