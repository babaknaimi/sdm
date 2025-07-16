\name{Arith-methods}

\docType{methods}

\alias{Arith-methods}
\alias{+,sdmModels,sdmModels-method}



\title{Combine (merge) two sdmModels into a single object}

\description{
If two sets of models are fitted in two separate sdmModels objects, they can be merged into a single sdmModels object using `+` operator.
}

\value{
 An object of class \code{sdmModels}.
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}



\examples{
\dontrun{
file <- system.file("external/pa_df.csv", package = "sdm")

df <- read.csv(file)

head(df) 

d <- sdmData(sp ~ b15 + NDVI, train = df)

d
#----


m1 <- sdm(sp ~ b15 + NDVI, data = d, methods = c('glm', 'gbm'))

m1

m2 <- sdm(sp ~ b15 + NDVI, data = d, methods = 'svm')

m2

m <- m1 + m2 # combining two sdmModels objects into one

m

}

}

\keyword{sdm}
\keyword{math}
\keyword{methods}