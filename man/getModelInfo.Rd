\name{get models' outputs}

\alias{getModelInfo}
\alias{getModelInfo,sdmModels-method}
\alias{getModelId}
\alias{getModelId,sdmModels-method}
\alias{getModelObject}
\alias{getModelObject,sdmModels-method}

\title{Get information/modelIDs relevant to fitted models in a sdmModels object}

\description{
When SDMs are fitted using the \code{\link{sdm}} function, a \code{sdmModels} object is generated containing all the information and objects created through fitting and evaluation procedures for all species and methods. To each model, a unique modelID is assigned. \code{getModelInfo} returns a data.frame summarising some information relevant to the fitted models including modelID, method name, whether the model is fitted successfully, whether and what replication procedure is used for data partitioning, etc.  \code{getModelInfo} helps to get the unique model IDs for all or certain models given the parameters that users specify. \code{getModelObject} returns the fitted model object for a single model (specified through \code{id}, or other settings).

}

\usage{
getModelId(x,success,species,method,replication, run)

getModelInfo(x,...)

getModelObject(x,id,species,method,replication, run)


}

\arguments{
  \item{x}{a \code{sdmModel} object}
  \item{success}{logical; specifies whether the info/ids should be returned only for the models that are successfully fitted or not (default is TRUE)}
  \item{species}{optional; a character vector specifies the name of species for which the info should be returned (default is NULL meaning for all species)}
  \item{method}{optional; a character vector specifies the name of methods for which the info should be returned (default is NULL meaning for all methods)}
  \item{replication}{optional; a character vector specifies the name of replication method for which the info should be returned (default is NULL meaning for all species)}
  \item{run}{optional; a numeric vector specifies for which replication runs the info should be returned (default is NULL meaning for all runs)}
  \item{id}{a single numeric value specifying the modelID}
  \item{...}{additional arguments. see details}
  }
\details{
in \code{getModelInfo}, as additional arguments, you can use the arguments in the function \code{getModelId} to specify which records should be returned.

}

\value{
 getModelInfo: data.frame
 getModelId: a numeric vector
 getModelObject: The fitted model object with a class depending on the method
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org}
}

\seealso{# }

\examples{

file <- system.file("external/model.sdm", package="sdm")

m <- read.sdm(file)

getModelInfo(m)

# getModelId(m)

# getModelId(m,method='brt')

obj <- getModelInfo(m, id=3) # obj is the fitted BRT model (through the package of gbm)

class(obj) # The class of the model object

summary(obj)

}

\keyword{model}
\keyword{sdm}
