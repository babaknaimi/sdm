\docType{class}
\name{featuresFrame-class}
\alias{featuresFrame-class}
\title{featureFrame class}
\description{
An S4 class contains the information of features used to fit a model
}
\section{Slots}{

\describe{
\item{\code{vars}}{A character vector containing the name of variabes from the dataset used to generate the features}

\item{\code{feature.types}}{A list containing the definition of features}

\item{\code{response.specific}}{NULL, or a list containing the definition of features that their definitions are according to the response variable (i.e. species)}
}}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}


\keyword{spatial,species,data}