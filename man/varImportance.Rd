\name{getVarImp}


\alias{getVarImp}
\alias{getVarImp,sdmModels-method}


\title{variable importance}


\description{
Calculates relative importance of different variables in the models using several approaches.
}

\usage{
getVarImp(x,id,wtest,setting,...)

}

\arguments{
  \item{x}{sdmModels object}
  \item{id}{numeric, specify the model (modelID) for which the variable importance values are extracted; OR it can be character with "ensemble" specifying that the variable importance should be calculated based on the ensemble of all the model objects}
  \item{wtest}{specifies which dataset ('training','test.dep','test.indep') should be used (if exist) to calculate the importance of variables }
  \item{setting}{an optional list with setting of ensemble function; it is only needed when id = 'ensemble'}
  \item{...}{additional arguments as for \code{getModelId} function including \code{species},\code{method}, \code{replication}, and \code{run}}
  }

\details{

\code{getVarImp} function returns an object including different measures of variable importance, and if be put in plot function, a barplot is generated. If the ggplot2 package is installed on your machine, the plot is generated using ggplot (unless you turn gg = FALSE), otherwise, the standard barplot is used.

if \code{id = "ensemble"} is used in the function, the ensemble function is called to calculate the relative variable importance based on the ensemble prediction of all models. \code{setting} can be specified as an additional argument that will be passed to the ensemble function so check the \code{\link{ensemble}} function to see how can \code{setting} be specified!

}



\references{

Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}



\examples{
\dontrun{
# if m is a sdmModels object (output of the sdm function) then:

getVarImp(m,id=1) # variable importance 

vi <- getVarImp(m,id=1) 

vi

plot(vi,'auc')

plot(vi,'cor')
#############
# You can get Mean variable importance (and confidence interval) for multiple models:

vi <- getVarImp(m,id=1:10) # specify the modelIDs of the models

vi

plot(vi,'cor')
#----
# you can use the getModelId function to find the id of the specific method, replication, etc.
# or you may put the arguments of the getModelId in the getVarImp function:

vi <- getVarImp(m, method='glm') # Mean variable importance for the method glm

vi

plot(vi)
#################

##### Variable Importance based on ENSEMBLE:

# You can get variable importance based on the ensemble of multiple models:

# setting is passed to the ensemble function

vi <- getVarImp(m,id="ensemble",setting=list(method='weighted',stat='auc'))

vi

plot(vi,'cor')

#----------------

# if you want the ensemble based on a subset of models, you can define
# the id within the setting list:

vi <- getVarImp(m,id="ensemble",
        setting=list(method='weighted',stat='auc',id=1:10))

vi

plot(vi,'cor')


plot(vi, gg = F) # R standard plot is used instead of ggplot


}
}

\keyword{spatial}