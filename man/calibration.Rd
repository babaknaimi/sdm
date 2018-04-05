\name{calibration}

\alias{calibration}
\alias{calibration,vector,vector-method}
\alias{calibration,sdmEvaluate,missing-method}

\title{Calibration}

\description{
evaluates for calibration 
}

\usage{
calibration(x,p,nbin,weight,...)

}

\arguments{
  \item{x}{a numeric vector including the observed values; or a \code{sdmEvaluate} object}
  \item{p}{a numeric vector including the predicted values}
  \item{nbin}{number of bins to discretize the predicted values into the specified bins (default: 10); instead, it can be the keyword of 'seek' to ask for seeking the best number }
  \item{weight}{logical, specifies whether a weight should be calculated based on the number of records at each bin. The weight will be used to summarize the calibration statistic}
  \item{...}{additional arguments (not implemented yet.)}
  }
\details{
The output of this function can be used in the plot function to generate Calibration plot. The calibration statistic is calculated using a method developed by the authors of this package (the journal article is not published yet, but in preparation)

}

\value{
 an object of class \code{.sdmCalibration} 
}

\references{
Naimi, B., Niamir, A., Jimenez-Valverde, A., Araujo, M.B. (In preparation) Measuring calibration capacity of statistical models: a new statistic.

Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}



\examples{



ca <- calibration(x=c(1,1,0,1,0,0,0,1,1,1,0),
          p=c(0.69,0.04,0.05,0.95,0.04,0.65,0.09,0.61,0.75,0.84,0.15))

ca

plot(ca)


}

\keyword{spatial,species,data}