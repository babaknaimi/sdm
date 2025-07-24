\name{coords}


\alias{coords}
\alias{coords<-}
\alias{coords,sdmdata-method}
\alias{coords,sdmModels-method}
\alias{coords<-,sdmdata-method}
  
\title{Get or set spatial coordinates of species data}

\description{
Get or set spatial coordinates of a \code{sdmdata} object. 
}

\usage{
\S4method{coords}{sdmdata}(obj,...)

\S4method{coords}{sdmdata}(object)<-value

}

\arguments{
  \item{obj}{speciesData (either of singleSpecies, multiple Species or SpeciesDataList) object}
  \item{object}{same as obj}
  \item{value}{spatial coordinates, either a matrix, or data.frame, or as character to change the names of coordinates}
  \item{...}{Additional arguments}
}

\value{
matrix, or if the coordinates set, the \code{sdmdata} object is returned.
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}


\examples{
file <- system.file("external/data.sdd", package = "sdm")
d <- read.sdm(filename = file)

d # a sdmdata object

coords(d)

}

\keyword{spatial}
\keyword{data}