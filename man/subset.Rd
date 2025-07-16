\name{subset}

\docType{methods}


\alias{subset}
\alias{subset,sdmModels-method}
\alias{[[,sdmModels,ANY,ANY-method}



\title{Subset models in a sdmModels object}

\description{
This function extracts a subset of models from a \code{sdmModels} object. In generates a new object of the same type as the original object. In \code{sdmModels}, modelID provides the unique IDs.

Instead of using the \code{subset} function, double brackes '[[ ]]' can be used.
}


\section{Methods}{
\describe{

\code{subset(x, subset, drop = TRUE, ...)}

\code{x[[i,...]]}


Arguments:

\code{x} - sdmModels object

\code{i}- integer. Indicates the index/id of the models (modelID) should be extracted from sdmModels object

\code{subset} -  same as \code{i}

\code{drop} - if \code{TRUE}, new modelIDs are generated, otherwise, the original modelIDs are kept in the new object.

\code{...} - additional arguments (Not implemented yet!)

}
}

\details{
#
}

\value{
sdmModels object
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}


\examples{
\dontrun{
file <- system.file("external/model.sdm", package = "sdm")

m <- read.sdm(filename = file)

m

getModelInfo(x = m)


m1 <- m[[3:4]]

m1

getModelInfo(x = m1)

m2 <- m[[3:4, drop = FALSE]]

m2

getModelInfo(x = m2)

#---- the following is the same as previous:

m2 <- subset(x = m, 3:4, drop = FALSE)

m2

getModelInfo(x =  m2)
} 
}

\keyword{ spatial }