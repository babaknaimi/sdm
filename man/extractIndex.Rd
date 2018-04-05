\name{Extract by index}

\docType{methods}


\alias{[,sdmdata,missing,missing-method}
\alias{[,sdmdata,ANY,ANY-method}

\title{Indexing to extract records of a sdmdata object}

\description{

This function extracts records of a \code{sdmdata} object and generates a new object of the same type (if drop=FALSE; otherwise a data.frame).
In \code{sdmdata}, rID is the unique ID for each record.
}

\section{Methods}{
\describe{
\code{x[i]}

Arguments
\tabular{rll}{
\tab \code{x} \tab a Raster* object \cr
\tab \code{i} \tab an index: record id (rID) in sdmdata object \cr
\tab \code{drop} \tab  If \code{TRUE}, a \code{data.frame} is returned, otherwise a \code{sdmdata} object is returned. \cr
}

}}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}


\examples{

file <- system.file("external/data.sdd", package="sdm")

d <- read.sdm(file)

# see the number of records:
d 
 
d2 <- d[1:10]

d2
 
d3 <- d[1:10,drop=TRUE]
 
d3

}

\keyword{methods}
\keyword{spatial}
