\name{roc}

\alias{roc}
\alias{getRoc}
\alias{roc,vector,vector-method}
\alias{roc,sdmModels-method}
\alias{roc,sdmModels,ANY-method}
\alias{getRoc,vector,vector-method}
\alias{getRoc,sdmModels,ANY-method}

\title{plot ROC curves}

\description{
Plot the Receiver Operating Characteristics (ROC) curve with AUC statistic in the legend.
}

\usage{
roc(x,p=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,
     wtest=NULL,smooth=FALSE,legend=TRUE,...)
     
getRoc(x,p,...)
}

\arguments{
  \item{x}{Either \code{sdmModels}, or \code{sdmEvaluate} object; or a numeric vector including observed binary values of species occurrence}
  \item{p}{if x is sdmModels, p is an optional vector with model ID number(s) that should be plotted (NULL (default means all models)); if x is a numeric vector, p is a vector with the same length including the predicted values}
  \item{species}{the name of species should be specified (required if x is \code{sdmModels} containing models for several species)}
  \item{method}{a character vector with the name of modelling methods that one need to get the roc plot for (if NULL  [default], all methods in the object are considered); only if x is \code{sdmModels}}
  \item{replication}{a character vector with the name of replication methods (i.e., 'sub','cv','boot') that one need to get the roc plot for}
  \item{run}{if x is \code{sdmModels}, and the models are fitted through a replication procedure, run specifies which runs of the partitioning (replications) are required; if NULL, all are considered}
  \item{wtest}{evaluation for which test datasets are required, maximum 2 names from 'training', 'test.dep', 'test.indep' (i.e., evaluation for training data, dependent test dataset, and independent test dataset, respectively)}
  \item{smooth}{logical, specified whether the ROC curves should be smoothed through a spline procedure}
  \item{legend}{logical, specified whether a legend including AUC statistic is required on the plot}
  \item{...}{additional arguments passed to plot function}
  }
\details{
roc generates the plots of roc curves, and getRoc generate the values of ROC
}

\value{
 an object of class \code{matrix}
}

\references{
#
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://biogeoinformatics.org}
}

\seealso{# }

\examples{
#
}

\keyword{spatial,species,data}
