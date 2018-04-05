\name{density}

\alias{density}
\alias{density,sdmEvaluate-method}


\docType{methods}

\title{density}

\description{
Create a density plots of presence and absence data
}


\section{Methods}{
\describe{
\code{density(x, ...)}
  
\tabular{rll}{  
 \tab  \code{x}  \tab  Object of class 'sdmEvaluate' (or a numeric vector of observed presence/absence)  \cr
 
}
}
}

\value{
A density plot. Presence data are in darkblue, and absence data are in red.
}



\references{

Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}



\examples{



e <- evaluates(x=c(1,1,0,1,0,0,0,1,1,1,0),
          p=c(0.69,0.04,0.05,0.95,0.04,0.65,0.09,0.61,0.75,0.84,0.15))


density(e)


}

\keyword{spatial,data}