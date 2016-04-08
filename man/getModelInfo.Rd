\name{get models' outputs}

\alias{getModelInfo}
\alias{getModelInfo,sdmModels-method}

\title{get info out of a sdmModels objects}

\description{
When SDMs are fitted using \code{\link{sdm}} function, a \code{sdmModels} object is generated that contains all the information and objects created through fitting and evaluation procedures for all species and methods. These functions facilitate getting these information or specific objects from the \code{sdmModels}.

}

\usage{
getModelInfo(x,...)


}

\arguments{
  \item{x}{a \code{sdmModel} object}
  \item{...}{additional arguments. see details}
  }
\details{
\code{getModelInfo} generates a data.frame containing the information on the whole information available in the model (i.e., \code{sdmModels} object), including the name of methods, what setting were used (i.e., replication) to evaluate, whether the fitting or different kinds of evaluations were successful. The first column in the data.frame contains the model IDs, the unique code for each model that is used everywhere in the package. That can be used for referring to a specific object in the other functions.

}

\value{
 data.frame
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://biogeoinformatics.org}
}

\seealso{# }

\examples{
\dontrun{
file <- system.file("external/model.sdm", package="sdm")

m <- read.sdm(file)

getModelInfo(m)
}
}

\keyword{spatial,species,data}
