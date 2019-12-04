\name{rcurve}

\alias{rcurve}
\alias{getResponseCurve}
\alias{rcurve,sdmModels-method}
\alias{rcurve,.responseCurve-method}
\alias{plot,.responseCurve-method}
\alias{getResponseCurve,sdmModels-method}

\title{Generate and plot response curves}

\description{
Calculate the response of species to the range of values in each predictor variable based on the fitted models in a sdmModels object.
}

\usage{
rcurve(x,n,id,mean,confidence,gg,...)
     
getResponseCurve(x,id,...)
}

\arguments{
  \item{x}{A \code{sdmModels} object; in the function response, it can be a \code{.responseCurve} object which is the output of the \code{getResponse} function }
  \item{id}{specifies the modelIDs corresponding to the models in the sdmModels object for which the response curves should be generated}
  \item{n}{A vector with the name of variables for which the response curve should be generated}
  \item{mean}{logical, specifies whether a mean should be calculated over responses to a variable when multiple models are specified in ids}
  \item{confidence}{logical, specifies whether a confidence interval should be added to the curve when the mean response curve is calculated based on multiple models}
  \item{gg}{logical, specifies whether the plot should be generated using the ggplot2 package (if the package is installed)}
  \item{...}{additional arguments passed to plot function}
  }
\details{
\code{getResponseCurve} calculates the responses for the models that are specified in \code{id} argument, and put the results in a \code{.responseCurve} object. This object can be used as an input in the \code{plot} function, or \code{rcurve} function.

If you just need the response curve graphs (plots), you can put a \code{sdmModels} object directly in the \code{rcurve} function, and do not need to first use \code{getResponseCurve} function.

In \code{getResponseCurve} function (or in \code{rcurve} when \code{x} is \code{sdmModels}), there are some additional arguments:

- \code{size}: a numeric value; default is 100. Specifies the size of the variable sequence that is used as the x-axis in the response curve plot. Greater number results to a smoother curve.

- \code{includeTest}: a logical value; default is FALSE; when a data object based on which a sdmModels is created containing independent test data; it specifies whether those records should be included into the response curve generation or not.

- \code{...}: additional arguments for the \code{plot} function (e.g., xlab, ylab, main, col, lwd, lty)


}


\value{
 an object of class \code{.responseCurve} or a series of graphs
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}


\examples{
\dontrun{
file <- system.file("external/model.sdm", package="sdm")

m <- read.sdm(file) # a sdmModels Object (fitted using sdm function)



rcurve(m)

rcurve(m,id=1) # for the first model

rcurve(m, id=1:2)

rcurve(m,method = 'glm',smooth = T) # only for models fitted using glm method & with smoothed curve


}
}

\keyword{spatial}
\keyword{sdm}
