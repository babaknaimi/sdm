\name{names}

\alias{names}
\alias{names<-}
\alias{names,sdmdata-method}
\alias{names<-,sdmdata-method}


  
\title{Names of species}

\description{
Get or set the names of the species of a \code{sdmdata} object
}

\usage{
\S4method{names}{sdmdata}(x)

\S4method{names}{sdmdata}(x)<-value

}

\arguments{
  \item{x}{A sdm data object (\code{sdmdata}}
  \item{value}{character (vector)}
}

\value{
For \code{names}, a character

For \code{names<-}, the updated object.
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://biogeoinformatics.org}
}


\examples{
#
}

\keyword{spatial}
