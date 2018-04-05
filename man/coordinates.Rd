\name{coordinates}


\alias{coordinates}
\alias{coordinates<-}
\alias{coordinates,sdmdata-method}
\alias{coordinates,sdmModels-method}
\alias{coordinates<-,sdmdata-method}
  
\title{get or set spatial coordinates of species data}

\description{
Get or set spatial coordiates of a \code{sdmdata} object. 
}

\usage{
\S4method{coordinates}{sdmdata}(obj)

\S4method{coordinates}{sdmdata}(object)<-value

}

\arguments{
  \item{obj}{speciesData (either of singleSpecies, multiple Species or SpeciesDataList) object}
  \item{object}{same as obj}
  \item{value}{spatial coordinates either a matrix, or data.frame, or column name as character, or a formula (e.g., ~x+y, or ~coords(x+y))}
  \item{...}{Additional arguments}
}

\value{
matrix, or if the coordinates set, the \code{sdmdata} object is returned.
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}


\examples{
file <- system.file("external/data.sdd", package="sdm")
d <- read.sdm(file)

d # a sdmdata object

coordinates(d)

}

\keyword{spatial}
