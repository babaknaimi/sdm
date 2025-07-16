\name{sdmAdapt}

\alias{sdmAdapt}
\alias{sdmAdapt,sdmdata-method}
\alias{sdmAdapt,sdmModels-method}

\title{Adapting sdm* objects in the new version}

\description{
The structure of the \code{sdmdata} and \code{sdmModels} classes were slightly changed in the new version of the package (> 1.2-X). If an \code{sdmdata} or \code{sdmModels} object is created and saved in an old version of the package (e.g., 1.1-8), using the \code{sdmAdapt} function, its structure is modified and adapted to the new version.
}

\usage{
sdmAdapt(x)


}

\arguments{
  \item{x}{an object of \code{sdmdata} or \code{sdmModels}}
  
}


\value{
 an object with the same class as \code{x} 
 
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}

}

\seealso{# }

\examples{
\dontrun{
file <- system.file("external/model.sdm", package = "sdm")

m <- read.sdm(filename = file) # an sdmModels Object (fitted using old version of sdm)

m <- sdmAdaptx = m)

m



}
}

\keyword{spatial}
\keyword{data}
\keyword{accuracy}