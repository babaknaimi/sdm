\docType{class}
\name{sdmdata-class}
\alias{sdmdata-class}
\alias{show,sdmdata-method}

\title{An S4 class representing sdm dataset}
\description{
An S4 class representing sdm dataset \code{sdmdata}
}
\section{Slots}{

\describe{
\item{\code{species.names}}{the names of species}

\item{\code{species}}{contains the species data}

\item{\code{features.name}}{the names of predictor variables}

\item{\code{features}}{a data.frame containing predictor variables}

\item{\code{factors}}{the names of categorical variables (if any)}

\item{\code{info}}{other information such as coordinates, metadata, etc.}

\item{\code{groups}}{a list including information on groups in the dataset}

\item{\code{sdmFormula}}{an object of class sdmFormula containing the formula and its' terms defined by user}

\item{\code{errorLog}}{reports on errors in the data raised through data cleaning (e.g., NA, duplications, etc.)}
}}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}



\keyword{spatial}
\keyword{methods}
\keyword{data}
\keyword{species}