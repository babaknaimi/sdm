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
\item{\code{species.names}}{The names of species}

\item{\code{species}}{Contains the species data}

\item{\code{features.name}}{The names of predictor variables}

\item{\code{features}}{A data.frame containing predictor variables}

\item{\code{factors}}{The names of categorical variables (if any)}

\item{\code{info}}{Other information such as coordinates, metadata, etc.}

\item{\code{groups}}{A list including information on groups in the dataset}

\item{\code{sdmFormula}}{An object of class sdmFormula containing the formula and its' terms defined by user}

\item{\code{errorLog}}{Reports on errors in the data raised throgh data cleaning (e.g., NA, duplications, etc.)}
}}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://biogeoinformatics.org}
}



\keyword{spatial,species,data}