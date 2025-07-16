\name{boxplot}


\alias{boxplot}
\alias{boxplot,sdmEvaluate-method}


\title{boxplot}


\description{
Make a box plot of model evaluation data, i.e., the model predictions for known presence and absence points.
}


\arguments{
  \item{x} Object of class \code{sdmEvaluate}
  \item{names} Optional, a character vector with the x-axis label(s) for the group of data (e.g., 'Absence', 'Presence')
  \item{...} Additional arguments that can be passed to \link[graphics]{boxplot}
}

\references{

Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}



\examples{



e <- evaluates(x = c(1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0),
          p = c(0.69, 0.04, 0.05, 0.95, 0.04, 0.65, 0.09, 0.61, 0.75, 0.84, 0.15))

boxplot(x = e, names = c("Absence", "Presence"))

}

\keyword{visualise}
\keyword{data}