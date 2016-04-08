\name{sdmSetting}

\alias{sdmSetting}
\alias{sdmSetting,ANY,sdmdata,character-method}
\alias{show,sdmSetting-method}

\title{creating sdmSetting object}

\description{
Creates sdmSetting object that holds settings to fit and evaluate the models. It can be used to reproduce a study.
}

\usage{
sdmSetting(formula,data,methods,interaction.depth=1,n=1,replication=NULL,cv.folds=NULL,
     test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
     var.selection=FALSE,ncore=1L,...)


}

\arguments{
  \item{formula}{specify the structure of the model}
  \item{data}{sdm data object or data.frame including species and feature data}
  \item{methods}{character, name of the algorithms}
  \item{interaction.depth}{level of interactions between predictors}
  \item{n}{number of replicates (run)}
  \item{replication}{replication method (e.g., 'subsampling', 'bootstrapping', 'cv')}
  \item{cv.folds}{number of folds if cv (cross-validation) is in the selected replication methods}
  \item{test.percent}{test percentage if subsampling is in the selected replication methods}
  \item{bg}{method to generate background}
  \item{bg.n}{number of background records}
  \item{var.importance}{logical, whether variable importance should be calculated}
  \item{response.curve}{method to calculate variable importance}
  \item{var.selection}{logical, whether variable selection should be considered}
  \item{ncore}{number of cores to parallelize processing}
  \item{...}{additional arguments}
  }
\details{
using sdmSetting, the feature types, interaction.depth and all settings of the model can be defined.
}

\value{
 an object of class \code{sdmSettings}
}

\references{#
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://biogeoinformatics.org}
}

\seealso{ #}

\examples{
#
}

\keyword{spatial,species,data}
