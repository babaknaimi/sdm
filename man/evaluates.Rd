\name{evaluates}

\alias{evaluates}
\alias{getEvaluation}
\alias{evaluates,vector,vector-method}
\alias{getEvaluation,sdmModels-method}

\title{evaluate for accuracy}

\description{
evaluates for accuracy 
}

\usage{
evaluates(x,p,...)

getEvaluation(x,w,wtest,stat,opt,...)
}

\arguments{
  \item{x}{a numeric vector including the observed values; or a \code{sdmModels} object}
  \item{p}{a numeric vector including the predicted values}
  \item{w}{a numeric vector indicates model IDs}
  \item{wtest}{which test, training, dep.test, or indep.test?}
  \item{stat}{statistics that should be extracted from the \code{sdmEvaluate} object}
  \item{opt}{a numeric value indicates which threshold optimisation criteria should be considered if a threshold-based statistic is selected in stat}
  \item{...}{additional arguments (see details)}
  }
\details{
Evaluates the preformance (accuracy) given the obsetved values, and the predicted values. As additional argument, the distribution of data can be specified (through \code{distribution}), that can be either of \code{'binomial'}, \code{'gaussian'}, \code{'laplase'}, or \code{'poisson'}. If not specified, it will be guessed by the function!


}

\value{
 an object of class \code{sdmEvaluate}
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://biogeoinformatics.org}
}

\seealso{# }

\examples{
#
}

\keyword{spatial,species,data}